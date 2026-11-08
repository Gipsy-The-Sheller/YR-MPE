# codon_ztest_plugin.py
#
# Copyright (c) 2026 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Codon-based Z-test (NG86 counting + analytic Z)."""

from __future__ import annotations

import os
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
    QStackedWidget,
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
    AlignmentSummary,
    analyze_alignment,
    analyze_alignment_yn00,
    format_results_csv,
    selection_sort_key,
)


class CodonZTestThread(BaseProcessThread):
    result_ready = pyqtSignal(object)

    def __init__(self, alignment_path: str, code_id: int, method: str = "ng86",
                 workdir: str | None = None):
        super().__init__("", [alignment_path], [], workdir=workdir)
        self.alignment_path = alignment_path
        self.code_id = code_id
        self.method = method

    def get_tool_name(self) -> str:
        return "Codon-based Z-test"

    def execute_commands(self):
        self.progress.emit("Reading codon alignment...")
        records = list(SeqIO.parse(self.alignment_path, "fasta"))
        if len(records) < 2:
            raise ValueError("At least two sequences are required for a pairwise Z-test.")
        self.progress.emit("Counting synonymous/nonsynonymous sites and differences...")
        if self.method == "yn00":
            summary = analyze_alignment_yn00(
                records,
                code_id=self.code_id,
                alignment=Path(self.alignment_path).name,
            )
        else:
            summary = analyze_alignment(
                records,
                code_id=self.code_id,
                alignment=Path(self.alignment_path).name,
            )
        self.result_ready.emit(summary)
        self.finished.emit([], [])


class CodonZTestPlugin(BasePlugin):
    """Pairwise codon-based Z-test for positive/purifying selection."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self._summary = None
        self.analysis_thread = None
        super().__init__(import_from, import_data, **kwargs)
        self.run_button.setText("Run Z-test")

    def init_plugin_info(self):
        self.plugin_name = "Codon-based Z-test (NG86 / YN00)"
        self.tool_name = "Built-in Python Implementation"
        self.citation = []
        self.citation_keys = ["nei1986", "nei2000molecular", "yang2000"]
        self.input_types = {"Codon alignment FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"Selection table": [".tsv"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    def config(self):
        return True

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

        title = QLabel("Codon-based Z-test")
        title.setStyleSheet("font-size: 18px; font-weight: 600;")
        layout.addWidget(title)
        subtitle = QLabel(
            "Pairwise dN/dS with an analytic Z-test, using either Nei-Gojobori "
            "(NG86) counting or Yang-Nielsen (YN00) with F84 correction. "
            "The alignment must be in-frame (length a multiple of 3)."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        layout.addWidget(subtitle)

        input_group = QGroupBox("Input")
        form = QFormLayout(input_group)
        self.alignment_edit = QLineEdit()
        self.alignment_edit.setPlaceholderText("Codon-aligned FASTA")
        browse = QPushButton("Browse")
        browse.clicked.connect(self._browse_alignment)
        row = QHBoxLayout()
        row.addWidget(self.alignment_edit, 1)
        row.addWidget(browse)
        form.addRow("Alignment:", row)
        self.input_info = QLabel("No alignment selected.")
        self.input_info.setWordWrap(True)
        self.input_info.setStyleSheet("color: palette(mid);")
        form.addRow("", self.input_info)
        layout.addWidget(input_group)

        params = QGroupBox("Parameters")
        pform = QFormLayout(params)
        self.method_combo = QComboBox()
        self.method_combo.addItem("NG86 — Nei-Gojobori counting", "ng86")
        self.method_combo.addItem("YN00 — Yang-Nielsen, F84 correction", "yn00")
        self.method_combo.currentIndexChanged.connect(self._on_method_changed)
        pform.addRow("Method:", self.method_combo)
        self.code_combo = QComboBox()
        for label, code_id in GENETIC_CODES.items():
            self.code_combo.addItem(label, code_id)
        pform.addRow("Genetic code:", self.code_combo)
        self.method_hint = QLabel("")
        self.method_hint.setWordWrap(True)
        self.method_hint.setStyleSheet("color: palette(mid);")
        pform.addRow("", self.method_hint)
        self.ztest_hint = QLabel(
            "Z-test: analytic delta method; two-sided and one-sided p-values"
        )
        self.ztest_hint.setWordWrap(True)
        self.ztest_hint.setStyleSheet("color: palette(mid);")
        pform.addRow("", self.ztest_hint)
        layout.addWidget(params)

        note = QLabel(
            "Gaps, ambiguity, and stops are skipped pairwise. YN00 follows PAML "
            "yn00 equal weighting (F84; one kappa). For batch/MLE, use "
            "SELECTION → Batch dN/dS."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color: palette(mid);")
        layout.addWidget(note)
        layout.addStretch()
        self._on_method_changed()

    def _on_method_changed(self):
        yn00 = self.method_combo.currentData() == "yn00"
        if yn00:
            self.method_hint.setText(
                "One kappa per alignment; F84-corrected dS/dN with F84 SEs."
            )
            self.ztest_hint.setText(
                "Z = (dN - dS) / sqrt(SE(dN)^2 + SE(dS)^2)"
            )
        else:
            self.method_hint.setText(
                "NG86 counts synonymous/nonsynonymous sites and differences directly."
            )
            self.ztest_hint.setText(
                "Z-test: analytic delta method; two-sided and one-sided p-values"
            )

    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.summary_label = QLabel("No Z-test yet.")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.summary_label)

        view_row = QHBoxLayout()
        view_row.addWidget(QLabel("View:"))
        self.view_combo = QComboBox()
        self.view_combo.addItem("Table (pairwise rows)", "table")
        self.view_combo.addItem("Matrix (MEGA style)", "matrix")
        self.view_combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.view_combo.currentIndexChanged.connect(self._on_view_changed)
        view_row.addWidget(self.view_combo)
        view_row.addStretch(1)
        self.view_hint = QLabel(
            "Upper triangle: two-sided p-value · lower triangle: dN/dS (ω). "
            "Bold/tinted p < 0.05 (red: dN > dS; green: dN < dS)."
        )
        self.view_hint.setWordWrap(True)
        self.view_hint.setStyleSheet("color: palette(mid);")
        view_row.addWidget(self.view_hint, 2)
        layout.addLayout(view_row)

        self.view_stack = QStackedWidget()

        self.table = QTableWidget(0, 12)
        self.table.setHorizontalHeaderLabels([
            "Seq 1", "Seq 2", "Sites", "S", "N", "Sd", "Nd",
            "dS (SE)", "dN (SE)", "dN/dS", "Z", "Conclusion",
        ])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.view_stack.addWidget(self.table)

        self.matrix = QTableWidget(0, 0)
        self.matrix.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.matrix.setSelectionMode(QAbstractItemView.SingleSelection)
        self.matrix.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.matrix.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.view_stack.addWidget(self.matrix)

        layout.addWidget(self.view_stack, 1)

        actions = QHBoxLayout()
        self.export_button = QPushButton("Export table (TSV)")
        self.export_button.clicked.connect(self._export_table)
        actions.addWidget(self.export_button)
        actions.addStretch()
        layout.addLayout(actions)
        self.view_hint.setVisible(False)

    def setup_control_panel(self):
        super().setup_control_panel()
        self.run_button.setText("Run Z-test")
        self.control_layout.insertStretch(0, 1)

    # ------------------------------------------------------------------
    # input / run
    # ------------------------------------------------------------------
    def _browse_alignment(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select codon alignment", self.alignment_edit.text() or os.path.expanduser("~"),
            "FASTA (*.fas *.fna *.fa *.fasta);;All files (*)",
        )
        if path:
            self.alignment_edit.setText(path)
            self._scan_alignment()

    def _scan_alignment(self):
        path = Path(self.alignment_edit.text())
        if not path.is_file():
            self.input_info.setText("⚠ Alignment file not found.")
            self.input_info.setStyleSheet("color: #b00020;")
            return False
        try:
            records = list(SeqIO.parse(str(path), "fasta"))
        except Exception as exc:
            self.input_info.setText("⚠ Could not read FASTA: %s" % exc)
            self.input_info.setStyleSheet("color: #b00020;")
            return False
        if len(records) < 2:
            self.input_info.setText("⚠ At least two sequences are required.")
            self.input_info.setStyleSheet("color: #b00020;")
            return False
        lengths = {len(record.seq) for record in records}
        if len(lengths) != 1 or lengths.pop() % 3 != 0:
            self.input_info.setText("⚠ Sequences must be aligned and in-frame (length multiple of 3).")
            self.input_info.setStyleSheet("color: #b00020;")
            return False
        self.input_info.setText("✓ %d sequences · %d codons" % (len(records), len(records[0].seq) // 3))
        self.input_info.setStyleSheet("color: #2e7d32;")
        return True

    def handle_import_data(self, import_data):
        sequences = None
        if isinstance(import_data, dict):
            for item in import_data.get("dataset_items", []):
                if getattr(item, "sequences", None):
                    sequences = item.sequences
                    break
        elif isinstance(import_data, list):
            sequences = import_data
        elif hasattr(import_data, "sequences"):
            sequences = import_data.sequences
        if not sequences:
            return
        path = self.create_temp_file(suffix=".fas")
        with open(path, "w", encoding="utf-8") as handle:
            for index, record in enumerate(sequences):
                seq_id = getattr(record, "id", "seq%d" % (index + 1))
                handle.write(">%s\n%s\n" % (seq_id, getattr(record, "seq", record)))
        self.temp_files.append(path)
        if hasattr(self, "alignment_edit"):
            self.alignment_edit.setText(path)
            self._scan_alignment()

    def run_analysis(self):
        if not self._scan_alignment():
            QMessageBox.warning(self, "Input missing", "Please select a valid in-frame codon alignment.")
            return
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running Z-test...")
        self.tab_widget.setCurrentIndex(2)
        self.add_console_message("Starting codon-based Z-test...", "info")

        self.analysis_thread = CodonZTestThread(
            self.alignment_edit.text(),
            int(self.code_combo.currentData() or 1),
            method=self.method_combo.currentData(),
            workdir=self.workdir,
        )
        self.analysis_thread.progress.connect(self._on_progress)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.result_ready.connect(self._on_summary)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.start()

    def _on_progress(self, message):
        self.progress_bar.setFormat(message)
        self.add_console_message(message, "info")

    def _on_summary(self, summary: AlignmentSummary):
        self._summary = summary

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        if self._summary is not None:
            self._display_summary(self._summary)
            summary = self._summary
            significant = sum(
                1 for result in summary.results
                if result.p_two_sided is not None and result.p_two_sided < 0.05
            )
            line = emit_summary(
                "codon_ztest",
                method=summary.method,
                alignment=summary.alignment,
                sequences=summary.sequences,
                pairs=len(summary.results),
                codons=summary.codons,
                usable_codons=summary.usable_codons,
                skipped_stop=summary.skipped_stop,
                skipped_ambiguous=summary.skipped_ambiguous,
                significant_p_lt_0_05=significant,
            )
            self.add_console_message(SUMMARY_PREFIX + line, "output")
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("Z-test completed.", "info")

    def analysis_error(self, message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Z-test failed: %s" % message, "error")
        QMessageBox.critical(self, "Z-test failed", message)

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
    def _display_summary(self, summary: AlignmentSummary):
        self.summary_label.setText(
            "%s — %s · %d sequences, %d pairs · %d/%d usable codons "
            "(skipped stop: %d, ambiguous: %d)"
            % (
                summary.alignment, summary.method, summary.sequences,
                len(summary.results), summary.usable_codons, summary.codons,
                summary.skipped_stop, summary.skipped_ambiguous,
            )
        )
        ordered = sorted(summary.results, key=selection_sort_key)
        self.table.setRowCount(len(ordered))
        for row, result in enumerate(ordered):
            values = [
                result.seq1, result.seq2, str(result.sites),
                _fmt_table(result.synonymous_sites), _fmt_table(result.nonsynonymous_sites),
                _fmt_table(result.synonymous_differences), _fmt_table(result.nonsynonymous_differences),
                _fmt_with_se(result.dS, result.se_dS, compact=True),
                _fmt_with_se(result.dN, result.se_dN, compact=True),
                _fmt_table(result.omega), _fmt_table(result.z), result.conclusion,
            ]
            for column, value in enumerate(values):
                self.table.setItem(row, column, QTableWidgetItem(value))
            self.table.item(row, 11).setToolTip(
                "Method: %s\nZ = %s · two-sided p = %s\n"
                "dS = %s ± %s · dN = %s ± %s"
                % (result.method, _fmt(result.z), _fmt(result.p_two_sided),
                   _fmt(result.dS), _fmt(result.se_dS),
                   _fmt(result.dN), _fmt(result.se_dN))
            )
            if result.p_two_sided is not None and result.p_two_sided < 0.05:
                item = self.table.item(row, 11)
                font = item.font()
                font.setBold(True)
                item.setFont(font)
                if result.dN is not None and result.dS is not None:
                    if result.dN > result.dS:
                        item.setBackground(QBrush(QColor(198, 40, 40, 45)))
                    elif result.dN < result.dS:
                        item.setBackground(QBrush(QColor(46, 125, 50, 45)))

        self._build_matrix(summary)

    def _build_matrix(self, summary: AlignmentSummary):
        names = list(summary.names)
        if not names:
            names = sorted({result.seq1 for result in summary.results}
                           | {result.seq2 for result in summary.results})
        lookup = {}
        for result in summary.results:
            lookup[(result.seq1, result.seq2)] = result
            lookup[(result.seq2, result.seq1)] = result

        self.matrix.clear()
        self.matrix.setRowCount(len(names))
        self.matrix.setColumnCount(len(names))
        self.matrix.setHorizontalHeaderLabels(names)
        self.matrix.setVerticalHeaderLabels(names)
        for row, name1 in enumerate(names):
            for column, name2 in enumerate(names):
                if row == column:
                    item = QTableWidgetItem("—")
                    item.setTextAlignment(Qt.AlignCenter)
                    item.setFlags(item.flags() & ~Qt.ItemIsSelectable)
                    self.matrix.setItem(row, column, item)
                    continue
                result = lookup.get((name1, name2))
                if row < column:
                    value = result.p_two_sided if result else None
                    item = QTableWidgetItem(_fmt(value))
                    if value is not None and value < 0.05:
                        font = item.font()
                        font.setBold(True)
                        item.setFont(font)
                        # A tint reinforces significance; bold + the numeric
                        # p-value keep the signal non-colour-only.
                        if result.dN is not None and result.dS is not None:
                            if result.dN > result.dS:
                                item.setBackground(QBrush(QColor(198, 40, 40, 45)))
                            elif result.dN < result.dS:
                                item.setBackground(QBrush(QColor(46, 125, 50, 45)))
                else:
                    value = result.omega if result else None
                    item = QTableWidgetItem(_fmt(value))
                item.setTextAlignment(Qt.AlignCenter)
                if result is not None:
                    item.setToolTip(
                        "%s vs %s\nZ = %s\np (two-sided) = %s\n"
                        "p (positive) = %s · p (purifying) = %s\n%s"
                        % (
                            result.seq1, result.seq2, _fmt(result.z),
                            _fmt(result.p_two_sided), _fmt(result.p_positive),
                            _fmt(result.p_negative), result.conclusion,
                        )
                    )
                self.matrix.setItem(row, column, item)

    def _on_view_changed(self):
        is_matrix = self.view_combo.currentData() == "matrix"
        self.view_stack.setCurrentWidget(self.matrix if is_matrix else self.table)
        self.view_hint.setVisible(is_matrix)
        self.export_button.setText(
            "Export matrix (TSV)" if is_matrix else "Export table (TSV)"
        )

    def _matrix_tsv(self) -> str:
        names = list(self._summary.names)
        if not names:
            names = sorted({result.seq1 for result in self._summary.results}
                           | {result.seq2 for result in self._summary.results})
        lookup = {}
        for result in self._summary.results:
            lookup[(result.seq1, result.seq2)] = result
            lookup[(result.seq2, result.seq1)] = result
        lines = ["\t".join([""] + names)]
        for row, name1 in enumerate(names):
            values = [name1]
            for column, name2 in enumerate(names):
                if row == column:
                    values.append("")
                    continue
                result = lookup.get((name1, name2))
                if row < column:
                    values.append(_fmt(result.p_two_sided) if result else "NA")
                else:
                    values.append(_fmt(result.omega) if result else "NA")
            lines.append("\t".join(values))
        return "\n".join(lines) + "\n"

    def _export_table(self):
        if self._summary is None:
            QMessageBox.warning(self, "No results", "Run the Z-test first.")
            return
        is_matrix = self.view_combo.currentData() == "matrix"
        default = str(Path(self.alignment_edit.text()).with_suffix(
            ".ztest_matrix.tsv" if is_matrix else ".ztest.tsv"
        ))
        target, _ = QFileDialog.getSaveFileName(
            self, "Export Z-test matrix" if is_matrix else "Export Z-test table",
            default, "TSV (*.tsv);;CSV (*.csv)"
        )
        if not target:
            return
        text = (self._matrix_tsv() if is_matrix
                else format_results_csv(sorted(self._summary.results,
                                               key=selection_sort_key)))
        if target.lower().endswith(".csv"):
            text = text.replace("\t", ",")
        Path(target).write_text(text, encoding="utf-8")
        self.add_console_message("Table exported to %s" % target, "info")


def _fmt_with_se(value, se, compact=False):
    formatter = _fmt_table if compact else _fmt
    text = formatter(value)
    if se is None:
        return text
    return "%s ± %s" % (text, formatter(se))


def _fmt_table(value):
    """Four significant digits keep the wide pairwise table readable."""
    if value is None:
        return "NA"
    if isinstance(value, float):
        if value != value:
            return "NA"
        if value == 0.0:
            value = 0.0
        return "%.4g" % value
    return str(value)


def _fmt(value):
    if value is None:
        return "NA"
    try:
        if value != value:  # NaN
            return "NA"
    except Exception:
        pass
    if isinstance(value, float):
        if value == 0.0:
            value = 0.0  # normalise -0.0
        return "%.6g" % value
    return str(value)


class CodonZTestPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return CodonZTestPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
