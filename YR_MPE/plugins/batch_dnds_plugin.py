# batch_dnds_plugin.py
#
# Copyright (c) 2026 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Batch dN/dS: NG86 counting or IQ-TREE GY/MG MLE."""

from __future__ import annotations

import os
import re
import shutil
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QBrush, QColor
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .components.stdout_summary import SUMMARY_PREFIX, emit_summary
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.methods.selection_utils import (
    GENETIC_CODES,
    PairwiseSelection,
    analyze_file,
    analyze_file_yn00,
    selection_sort_key,
    format_results_csv,
    parse_iqtree_omega,
)

FASTA_SUFFIXES = {".fas", ".fna", ".fa", ".fasta"}


class BatchDnDsThread(BaseProcessThread):
    result_ready = pyqtSignal(object)

    def __init__(self, files, method, code_id=1, threads=1, iqtree_path="",
                 output_dir="", workdir=None):
        super().__init__(iqtree_path or "", list(files), [], workdir=workdir)
        self.files = [Path(path) for path in files]
        self.method = method
        self.code_id = code_id
        self.threads = threads
        self.iqtree_path = Path(iqtree_path) if iqtree_path else None
        self.output_dir = Path(output_dir) if output_dir else None

    def get_tool_name(self) -> str:
        return "Batch dN/dS"

    def execute_commands(self):
        if self.method in ("counting", "yn00"):
            self._run_counting()
        else:
            self._run_mle()

    def _run_counting(self):
        summaries = []
        yn00 = self.method == "yn00"
        for index, path in enumerate(self.files, start=1):
            if not self.is_running:
                break
            self.progress.emit(
                "%s %d/%d: %s" % ("YN00" if yn00 else "Counting",
                                  index, len(self.files), path.name)
            )
            try:
                if yn00:
                    summary = analyze_file_yn00(str(path), code_id=self.code_id)
                else:
                    summary = analyze_file(str(path), code_id=self.code_id)
                summaries.append(summary)
            except Exception as exc:
                self.console_output.emit("%s: %s" % (path.name, exc), "error")
        self.result_ready.emit({"method": self.method, "summaries": summaries})
        self.finished.emit([], [])

    def _run_mle(self):
        if self.iqtree_path is None or not self.iqtree_path.is_file():
            raise FileNotFoundError("IQ-TREE executable is not configured (config.json → IQ-TREE 3).")
        if self.output_dir is None:
            raise ValueError("Output folder is required for MLE mode.")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        model = "GY" if self.method == "gy" else "MG"
        rows = []
        for index, path in enumerate(self.files, start=1):
            if not self.is_running:
                break
            self.progress.emit("IQ-TREE %s %d/%d: %s" % (model, index, len(self.files), path.name))
            stem = path.stem
            prefix = self.output_dir / ("%s_%s" % (stem, model.lower()))
            cmd = [
                str(self.iqtree_path), "-s", str(path),
                "--seqtype", "CODON", "-m", model,
                "-T", str(self.threads), "--prefix", str(prefix), "--redo",
            ]
            row = {
                "alignment": path.name,
                "method": "IQ-TREE " + model,
                "omega": None,
                "logl": None,
                "status": "ok",
                "artifact": str(prefix),
                "seq1": "", "seq2": "",
            }
            try:
                result = self.execute_command(cmd, cwd=str(self.output_dir))
                text = (result.stdout or "") + "\n" + (result.stderr or "")
                report = Path(str(prefix) + ".iqtree")
                if report.exists():
                    text += "\n" + report.read_text(encoding="utf-8", errors="replace")
                omega, logl = parse_iqtree_omega(text)
                row["omega"] = omega
                row["logl"] = logl
                if result.returncode != 0 or omega is None:
                    row["status"] = "failed"
            except Exception as exc:
                row["status"] = "failed: %s" % exc
            rows.append(row)
        self.result_ready.emit({"method": self.method, "rows": rows})
        self.finished.emit([], [])


class BatchDnDsPlugin(BasePlugin):
    """Batch dN/dS with counting or IQ-TREE codon MLE."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self._rows = []
        self._summaries = []
        self.analysis_thread = None
        self._import_temp_dir = None
        super().__init__(import_from, import_data, **kwargs)
        self.run_button.setText("Run batch")

    def init_plugin_info(self):
        self.plugin_name = "Batch dN/dS (NG86 / IQ-TREE)"
        self.tool_name = "IQ-TREE 3"
        self.citation = []
        self.citation_keys = [
            "nei1986", "nei2000molecular", "yang2000",
            "goldman1994", "muse1994", "wong2026iqtree",
        ]
        self.input_types = {"Codon alignments": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"Selection table": [".tsv"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    def config(self):
        # Counting mode does not require IQ-TREE; MLE resolves it lazily.
        return True

    def on_destroyed(self):
        super().on_destroyed()
        if self._import_temp_dir:
            shutil.rmtree(self._import_temp_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # layout
    # ------------------------------------------------------------------
    def setup_input_tab(self):
        outer = QVBoxLayout(self.input_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        outer.addWidget(scroll)
        content = QWidget()
        scroll.setWidget(content)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        title = QLabel("Batch dN/dS")
        title.setStyleSheet("font-size: 18px; font-weight: 600;")
        layout.addWidget(title)
        subtitle = QLabel(
            "One row per pair for NG86 or YN00 counting, or one global omega per "
            "alignment from an IQ-TREE codon model (GY/MG)."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        layout.addWidget(subtitle)

        input_group = QGroupBox("Codon alignments")
        form = QFormLayout(input_group)
        self.folder_edit = QLineEdit()
        self.folder_edit.setPlaceholderText("Folder of in-frame codon FASTA files")
        browse = QPushButton("Browse")
        browse.clicked.connect(self._browse_folder)
        row = QHBoxLayout()
        row.addWidget(self.folder_edit, 1)
        row.addWidget(browse)
        form.addRow("Input folder:", row)
        self.file_info = QLabel("No folder selected.")
        self.file_info.setWordWrap(True)
        self.file_info.setStyleSheet("color: palette(mid);")
        form.addRow("", self.file_info)
        layout.addWidget(input_group)

        params = QGroupBox("Parameters")
        pform = QFormLayout(params)
        self.method_combo = QComboBox()
        self.method_combo.addItem("Counting — NG86 (fast, pairwise)", "counting")
        self.method_combo.addItem("Counting — YN00 (F84 correction)", "yn00")
        self.method_combo.addItem("MLE — IQ-TREE GY", "gy")
        self.method_combo.addItem("MLE — IQ-TREE MG", "mg")
        self.method_combo.currentIndexChanged.connect(self._on_method_changed)
        pform.addRow("Method:", self.method_combo)

        self.code_combo = QComboBox()
        for label, code_id in GENETIC_CODES.items():
            self.code_combo.addItem(label, code_id)
        pform.addRow("Genetic code:", self.code_combo)

        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 16)
        self.threads_spin.setValue(1)
        pform.addRow("MLE threads:", self.threads_spin)

        self.output_edit = QLineEdit()
        self.output_edit.setPlaceholderText("Output folder for MLE runs")
        out_browse = QPushButton("Browse")
        out_browse.clicked.connect(self._browse_output)
        out_row = QHBoxLayout()
        out_row.addWidget(self.output_edit, 1)
        out_row.addWidget(out_browse)
        pform.addRow("Output folder:", out_row)
        self.method_hint = QLabel("")
        self.method_hint.setWordWrap(True)
        self.method_hint.setStyleSheet("color: palette(mid);")
        pform.addRow("", self.method_hint)
        layout.addWidget(params)

        note = QLabel(
            "NG86 and YN00 skip stop/ambiguous codons pairwise; YN00 uses F84 "
            "correction and one kappa per alignment (PAML yn00 equal weighting). "
            "MLE runs IQ-TREE with --seqtype CODON and reports the global dN/dS (omega)."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color: palette(mid);")
        layout.addWidget(note)
        layout.addStretch()
        self._on_method_changed()

    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.summary_label = QLabel("No batch result yet.")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.summary_label)

        self.table = QTableWidget(0, 10)
        self.table.setHorizontalHeaderLabels([
            "Alignment", "Seq 1", "Seq 2", "Method", "dN", "dS",
            "dN/dS / omega", "Z", "p", "Status",
        ])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.table, 1)

        actions = QHBoxLayout()
        self.export_button = QPushButton("Export table (TSV)")
        self.export_button.clicked.connect(self._export_table)
        actions.addWidget(self.export_button)
        actions.addStretch()
        layout.addLayout(actions)

    def setup_control_panel(self):
        super().setup_control_panel()
        self.run_button.setText("Run batch")
        self.control_layout.insertStretch(0, 1)

    def _on_method_changed(self):
        method = self.method_combo.currentData()
        mle = method in ("gy", "mg")
        self.threads_spin.setEnabled(mle)
        self.output_edit.setEnabled(mle)
        if mle:
            text = "MLE runs IQ-TREE for each alignment; one global omega per file."
        elif method == "yn00":
            text = ("YN00 uses Yang-Nielsen counting with F84 correction; "
                    "one kappa per alignment, pairwise rows.")
        else:
            text = "Counting runs in Python; no external binary required. Output rows are pairwise."
        self.method_hint.setText(text)

    # ------------------------------------------------------------------
    # input
    # ------------------------------------------------------------------
    def _browse_folder(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select folder with codon alignments",
            self.folder_edit.text() or os.path.expanduser("~"),
        )
        if path:
            self.folder_edit.setText(path)
            self._scan_folder()

    def _browse_output(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select output folder", self.output_edit.text() or os.path.expanduser("~")
        )
        if path:
            self.output_edit.setText(path)

    def _files(self):
        folder = Path(self.folder_edit.text())
        if not folder.is_dir():
            return []
        return sorted(
            path for path in folder.iterdir()
            if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
        )

    def _scan_folder(self):
        files = self._files()
        if not files:
            self.file_info.setText("⚠ No FASTA alignments found.")
            self.file_info.setStyleSheet("color: #b00020;")
            return False
        self.file_info.setText("✓ %d alignment(s) found." % len(files))
        self.file_info.setStyleSheet("color: #2e7d32;")
        return True

    def handle_import_data(self, import_data):
        items = []
        if isinstance(import_data, dict):
            items = list(import_data.get("dataset_items", []))
        elif isinstance(import_data, list):
            items = [type("Item", (), {"sequences": import_data, "loci_name": "alignment"})()]
        elif hasattr(import_data, "sequences"):
            items = [import_data]
        items = [item for item in items if getattr(item, "sequences", None)]
        if not items:
            return
        folder = Path(tempfile.mkdtemp(prefix="yrmpe_batch_dnds_"))
        self._import_temp_dir = str(folder)
        for index, item in enumerate(items):
            raw_name = getattr(item, "loci_name", "") or "alignment_%d" % (index + 1)
            safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw_name).strip("._")
            path = folder / ("%03d_%s.fas" % (index + 1, safe_name or "alignment"))
            with path.open("w", encoding="utf-8") as handle:
                for record in item.sequences:
                    handle.write(">%s\n%s\n" % (
                        getattr(record, "id", "seq"), getattr(record, "seq", record)
                    ))
        if hasattr(self, "folder_edit"):
            self.folder_edit.setText(str(folder))
            self._scan_folder()

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    def _iqtree_path(self):
        config = self._find_plugin_config()
        if not config:
            return ""
        raw = Path(str(config.get("path", "")))
        if not raw:
            return ""
        if not raw.is_absolute():
            raw = Path(self.plugin_path) / str(config.get("path", "")).lstrip("/\\")
        return str(raw) if raw.is_file() else ""

    def run_analysis(self):
        if not self._scan_folder():
            QMessageBox.warning(self, "Input missing", "Please select a folder of codon FASTA alignments.")
            return
        method = self.method_combo.currentData()
        output_dir = self.output_edit.text().strip()
        if method in ("gy", "mg"):
            iqtree = self._iqtree_path()
            if not iqtree:
                QMessageBox.warning(
                    self, "IQ-TREE missing",
                    "IQ-TREE 3 is not configured. Add its path to config.json or use counting mode.",
                )
                return
            if not output_dir:
                output_dir = str(Path(self.folder_edit.text()).parent / "batch_dnds_results")
                self.output_edit.setText(output_dir)
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running batch dN/dS...")
        self.tab_widget.setCurrentIndex(2)
        self.add_console_message("Starting batch dN/dS (%s)..." % method, "info")

        self.analysis_thread = BatchDnDsThread(
            self._files(), method,
            code_id=int(self.code_combo.currentData() or 1),
            threads=self.threads_spin.value(),
            iqtree_path=self._iqtree_path(),
            output_dir=output_dir,
            workdir=self.workdir,
        )
        self.analysis_thread.progress.connect(self._on_progress)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.start()

    def _on_progress(self, message):
        self.progress_bar.setFormat(message)
        self.add_console_message(message, "info")

    def _on_result(self, payload):
        method = payload.get("method") or self.method_combo.currentData()
        self._summaries = payload.get("summaries", [])
        self._rows = payload.get("rows", [])
        if self._summaries:
            self._rows = []
            for summary in self._summaries:
                for result in sorted(summary.results, key=selection_sort_key):
                    self._rows.append({"pairwise": result, "alignment": summary.alignment})
            results = [result for summary in self._summaries for result in summary.results]
            significant = sum(
                1 for result in results
                if result.p_two_sided is not None and result.p_two_sided < 0.05
            )
            line = emit_summary(
                "batch_dnds",
                method=method,
                alignments=len(self._summaries),
                pairs=len(results),
                significant_p_lt_0_05=significant,
                alignments_detail=[
                    {
                        "alignment": summary.alignment,
                        "sequences": summary.sequences,
                        "pairs": len(summary.results),
                    }
                    for summary in self._summaries
                ],
            )
        else:
            rows = payload.get("rows", [])
            line = emit_summary(
                "batch_dnds",
                method=method,
                alignments=len(rows),
                completed=sum(1 for row in rows if row.get("status") == "ok"),
                failed=sum(1 for row in rows if row.get("status") != "ok"),
            )
        self.add_console_message(SUMMARY_PREFIX + line, "output")

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self._display_results()
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("Batch dN/dS completed.", "info")

    def analysis_error(self, message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Batch dN/dS failed: %s" % message, "error")
        QMessageBox.critical(self, "Batch dN/dS failed", message)

    def stop_analysis(self):
        if self.analysis_thread and self.analysis_thread.isRunning():
            self.analysis_thread.stop()
            self.analysis_thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Stopped by user.", "warning")

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------
    def _display_results(self):
        rows = self._rows
        self.table.setRowCount(len(rows))
        for index, row in enumerate(rows):
            if "pairwise" in row:
                result = row["pairwise"]
                values = [
                    result.alignment, result.seq1, result.seq2, result.method,
                    _fmt(result.dN), _fmt(result.dS), _fmt(result.omega),
                    _fmt(result.z), _fmt(result.p_two_sided), result.conclusion,
                ]
            else:
                values = [
                    row.get("alignment", ""), "", "", row.get("method", ""),
                    "", "", _fmt(row.get("omega")), "", "", row.get("status", ""),
                ]
            for column, value in enumerate(values):
                self.table.setItem(index, column, QTableWidgetItem(value))
            if "pairwise" in row:
                result = row["pairwise"]
                self.table.item(index, 9).setToolTip(
                    "Z = %s · two-sided p = %s" % (_fmt(result.z),
                                                   _fmt(result.p_two_sided))
                )
                if result.p_two_sided is not None and result.p_two_sided < 0.05:
                    item = self.table.item(index, 9)
                    font = item.font()
                    font.setBold(True)
                    item.setFont(font)
                    if result.dN is not None and result.dS is not None:
                        if result.dN > result.dS:
                            item.setBackground(QBrush(QColor(198, 40, 40, 45)))
                        elif result.dN < result.dS:
                            item.setBackground(QBrush(QColor(46, 125, 50, 45)))

        if self._summaries:
            pairs = sum(len(summary.results) for summary in self._summaries)
            self.summary_label.setText(
                "%d alignment(s) · %d pairwise comparison(s)"
                % (len(self._summaries), pairs)
            )
        elif rows:
            ok = sum(1 for row in rows if row.get("status") == "ok")
            self.summary_label.setText(
                "%d alignment(s) · %d/%d MLE runs completed" % (len(rows), ok, len(rows))
            )
        else:
            self.summary_label.setText("No results; see Console.")

    def _export_table(self):
        if not self._rows:
            QMessageBox.warning(self, "No results", "Run the batch first.")
            return
        target, _ = QFileDialog.getSaveFileName(
            self, "Export batch dN/dS", str(Path(self.folder_edit.text()) / "batch_dnds.tsv"),
            "TSV (*.tsv);;CSV (*.csv)",
        )
        if not target:
            return
        if self._summaries:
            text = ""
            for index, summary in enumerate(self._summaries):
                block = format_results_csv(sorted(summary.results,
                                                  key=selection_sort_key))
                text += block if index == 0 else block.split("\n", 1)[1]
        else:
            text = "alignment\tmethod\tomega\tlogL\tstatus\tartifact\n"
            for row in self._rows:
                text += "\t".join(str(row.get(key, "")) for key in
                                  ("alignment", "method", "omega", "logl", "status", "artifact")) + "\n"
        if target.lower().endswith(".csv"):
            text = text.replace("\t", ",")
        Path(target).write_text(text, encoding="utf-8")
        self.add_console_message("Table exported to %s" % target, "info")


def _fmt(value):
    if value is None:
        return "NA"
    try:
        if value != value:
            return "NA"
    except Exception:
        pass
    if isinstance(value, float):
        if value == 0.0:
            value = 0.0  # normalise -0.0
        return "%.6g" % value
    return str(value)


class BatchDnDsPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return BatchDnDsPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
