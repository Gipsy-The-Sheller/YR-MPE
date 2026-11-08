# modeltest_ng_plugin.py
#
# Copyright (c) 2026 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""ModelTest-NG model selection for aligned nucleotide/protein alignments."""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path
from typing import Dict, List

from Bio import SeqIO
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QMessageBox, QPushButton, QSpinBox, QTableWidget,
    QTableWidgetItem, QVBoxLayout,
)

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread

TABLE_HEADERS = ["Rank", "Model", "K", "lnL", "Score", "Delta", "Weight"]


def detect_data_type(path: str) -> str:
    """Return ModelTest-NG's ``nt`` or ``aa`` token for an alignment."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            record = next(SeqIO.parse(handle, "fasta"))
    except Exception:
        return "nt"
    residues = set(str(record.seq).upper()) - set("-?.~")
    dna_alphabet = set("ACGTUNRYSWKMBDHV")
    return "nt" if residues <= dna_alphabet else "aa"


def parse_modeltest_output(path: str) -> Dict:
    """Parse ModelTest-NG's ``.out`` report into criterion tables."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    criteria: Dict[str, List[Dict]] = {}
    current = None
    for line in text.splitlines():
        match = re.match(
            r"^(BIC|AIC|AICc)\s+model\s+K\s+lnL\s+score\s+delta\s+weight",
            line.strip(),
        )
        if match:
            current = match.group(1)
            criteria[current] = []
            continue
        if current is None:
            continue
        if line.startswith("---"):
            if criteria[current]:
                current = None
            continue
        parts = line.split()
        if len(parts) >= 7 and parts[0].isdigit():
            try:
                criteria[current].append({
                    "Rank": int(parts[0]),
                    "Model": parts[1],
                    "K": int(parts[2]),
                    "lnL": float(parts[3]),
                    "Score": float(parts[4]),
                    "Delta": float(parts[5]),
                    "Weight": float(parts[6]),
                })
            except ValueError:
                continue
        elif criteria[current]:
            current = None

    best: Dict[str, Dict] = {}
    current = None
    field_re = re.compile(
        r"^(Model|lnL|Frequencies|Subst\. Rates|Inv\. sites prop|"
        r"Gamma shape|Score|Weight):\s*(.*)$"
    )
    for line in text.splitlines():
        match = re.match(r"^Best model according to (BIC|AICc|AIC)", line.strip())
        if match:
            current = match.group(1)
            best[current] = {}
            continue
        if current is None:
            continue
        field = field_re.match(line.strip())
        if field:
            key, value = field.group(1), field.group(2).strip()
            if key in ("Model", "Frequencies", "Subst. Rates",
                       "Inv. sites prop", "Gamma shape"):
                best[current][key] = value
            else:
                try:
                    best[current][key] = float(value)
                except ValueError:
                    best[current][key] = value

    for criterion, rows in criteria.items():
        if criterion not in best and rows:
            best[criterion] = {"Model": rows[0]["Model"]}
    return {"criteria": criteria, "best": best}


class ModelTestNGThread(BaseProcessThread):
    """Run ModelTest-NG for each input alignment."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, data_type="auto", topology="ml",
                 schemes="11", threads=2, imported_files=None, workdir=None):
        super().__init__(
            tool_path, input_files, [], imported_files, workdir=workdir
        )
        self.data_type = data_type
        self.topology = topology
        self.schemes = schemes
        self.threads = threads

    def get_tool_name(self):
        return "ModelTest-NG"

    def execute_commands(self):
        output_dir = tempfile.mkdtemp(
            prefix="modeltest-ng-",
            dir=self.workdir if self.workdir and os.path.isdir(self.workdir) else None,
        )
        output_files = []
        parsed = {}
        total = len(self.input_files)
        for index, input_file in enumerate(self.input_files):
            if not self.is_running:
                break
            self.progress.emit("Processing file %d/%d..." % (index + 1, total))
            self.console_output.emit(
                "ModelTest-NG: %s" % os.path.basename(input_file), "info"
            )
            data_type = self.data_type
            if data_type == "auto":
                data_type = detect_data_type(input_file)
            prefix = os.path.join(
                output_dir,
                "%02d_%s.modeltest" % (index + 1, Path(input_file).stem),
            )
            cmd = [
                self.tool_path,
                "-i", input_file,
                "-d", data_type,
                "-t", self.topology,
                "-p", str(self.threads),
                "-o", prefix,
                "--force",
                "--disable-checkpoint",
            ]
            if self.schemes and data_type == "nt":
                cmd.extend(["-s", self.schemes])
            result = self.execute_command(cmd, cwd=output_dir)
            out_path = prefix + ".out"
            if not os.path.exists(out_path):
                detail = (result.stderr or "").strip()
                raise RuntimeError(
                    "ModelTest-NG produced no .out report for %s%s"
                    % (os.path.basename(input_file),
                       ": " + detail if detail else "")
                )
            if result.returncode != 0:
                self.console_output.emit(
                    "ModelTest-NG exited with code %s after writing %s; "
                    "the report is valid (known Cygwin shutdown issue)."
                    % (result.returncode, os.path.basename(out_path)),
                    "warning",
                )
            output_files.append(out_path)
            parsed[out_path] = parse_modeltest_output(out_path)
        self.progress.emit("ModelTest-NG completed")
        self.result_ready.emit({"output_files": output_files, "parsed": parsed})
        self.finished.emit(output_files, [])


class ModelTestNGPlugin(BasePlugin):
    """ModelTest-NG model-selection plugin."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self.output_files = []
        self.parsed_results = {}
        super().__init__(import_from, import_data, **kwargs)
        if import_from == "DATASET_MANAGER" and import_data is not None:
            self.handle_dataset_import(import_data)
        self.run_button.setText("Run ModelTest-NG")

    def init_plugin_info(self):
        self.plugin_name = "ModelTest-NG"
        self.tool_name = "ModelTest-NG"
        self.citation = []
        self.citation_keys = []
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "PHYLIP": ["phy"]}
        self.output_types = {"ModelTest-NG": ".out"}
        self.plugin_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), ".."
        )

    def _set_imported_files(self, paths):
        self.imported_files = [path for path in paths if os.path.exists(path)]
        if len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        elif self.imported_files:
            self.file_path_edit.setText("%d files selected" % len(self.imported_files))
        else:
            self.file_path_edit.setText("")

    def handle_import_data(self, import_data):
        if isinstance(import_data, list):
            super().handle_import_data(import_data)
            if self.imported_files:
                self.file_path_edit.setText(self.imported_files[0])

    def handle_dataset_import(self, import_data):
        items = import_data.get("dataset_items", []) if isinstance(import_data, dict) else []
        paths = []
        for item in items:
            if getattr(item, "item_type", "") not in ("alignment", "sequence"):
                continue
            file_path = getattr(item, "file_path", "") or ""
            if file_path and os.path.exists(file_path):
                paths.append(file_path)
                continue
            sequences = getattr(item, "sequences", None)
            if not sequences:
                continue
            temp_file = self.create_temp_file(suffix=".fas")
            with open(temp_file, "w", encoding="utf-8") as handle:
                for record in sequences:
                    handle.write(">%s\n%s\n" % (record.id, record.seq))
            paths.append(temp_file)
        self._set_imported_files(paths)
        if paths:
            self.add_console_message(
                "Imported %d alignment(s) from the dataset" % len(paths), "info"
            )

    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)

        input_group = QGroupBox("Input alignment")
        input_layout = QFormLayout(input_group)
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select one or more aligned files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)
        layout.addWidget(input_group)

        params_group = QGroupBox("ModelTest-NG parameters")
        params_layout = QFormLayout(params_group)
        self.data_type_combo = QComboBox()
        self.data_type_combo.addItems(["Auto", "DNA", "Protein"])
        params_layout.addRow("Sequence type:", self.data_type_combo)

        self.criterion_combo = QComboBox()
        self.criterion_combo.addItems(["BIC", "AIC", "AICc"])
        self.criterion_combo.currentTextChanged.connect(self.refresh_display)
        params_layout.addRow("Criterion:", self.criterion_combo)

        self.topology_combo = QComboBox()
        self.topology_combo.addItems(["ML starting tree", "MP starting tree"])
        params_layout.addRow("Starting tree:", self.topology_combo)

        self.schemes_combo = QComboBox()
        self.schemes_combo.addItems(["3", "5", "7", "11"])
        self.schemes_combo.setCurrentText("11")
        params_layout.addRow("DNA schemes:", self.schemes_combo)

        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 64)
        self.threads_spinbox.setValue(2)
        params_layout.addRow("Threads:", self.threads_spinbox)
        layout.addWidget(params_group)

        note = QLabel(
            "ModelTest-NG reports BIC, AIC and AICc in one run. The selected "
            "criterion is used for the exported best model and the result table."
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        layout.addStretch(1)

        if getattr(self, "imported_files", None):
            self._set_imported_files(self.imported_files)

    def browse_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select aligned sequence files", "",
            "Alignments (*.fas *.fna *.fa *.fasta *.phy);;All files (*)"
        )
        if file_paths:
            self._set_imported_files(file_paths)

    def prepare_input_files(self):
        if self.imported_files:
            files = [path for path in self.imported_files if os.path.exists(path)]
            if files:
                return files
        path = self.file_path_edit.text().strip()
        if path and os.path.exists(path):
            return [path]
        QMessageBox.warning(self, "Warning", "Please select an alignment file.")
        return []

    def control_bindings(self):
        return [
            {"id": "input", "widget": "file_path_edit",
             "label": "Input alignment", "kind": "path",
             "cli": ["-i", "{value}"],
             "notes": "One command is run per selected file."},
            {"id": "data_type", "widget": "data_type_combo",
             "label": "Sequence type", "kind": "choice",
             "choices": {"Auto": "auto", "DNA": "nt", "Protein": "aa"},
             "cli": ["-d", "{value}"],
             "notes": "Auto is resolved to nt or aa before invocation."},
            {"id": "criterion", "widget": "criterion_combo",
             "label": "Criterion", "kind": "choice",
             "cli": [],
             "notes": "ModelTest-NG computes BIC/AIC/AICc together; this "
                      "selects the exported best model and table."},
            {"id": "topology", "widget": "topology_combo",
             "label": "Starting tree", "kind": "choice",
             "choices": {"ML starting tree": "ml", "MP starting tree": "mp"},
             "cli": ["-t", "{value}"]},
            {"id": "schemes", "widget": "schemes_combo",
             "label": "DNA schemes", "kind": "choice",
             "cli": ["-s", "{value}"],
             "notes": "Only passed for nucleotide data."},
            {"id": "threads", "widget": "threads_spinbox",
             "label": "Threads", "kind": "int",
             "cli": ["-p", "{value}"]},
        ]

    def cli_preview(self):
        parameters = self.get_parameters()
        input_file = (self.imported_files[0] if self.imported_files
                      else self.file_path_edit.text().strip() or "<input>")
        raw_type = parameters["data_type"]
        data_type = "<auto: nt|aa>" if raw_type == "auto" else raw_type
        command = [
            self.tool_path or "modeltest-ng.exe",
            "-i", input_file,
            "-d", data_type,
            "-t", parameters["topology"],
            "-p", str(parameters["threads"]),
            "-o", "<output-prefix>",
            "--force",
            "--disable-checkpoint",
        ]
        if parameters["schemes"] and raw_type in ("nt", "auto"):
            command.extend(["-s", parameters["schemes"]])
        return command

    def get_parameters(self):
        data_type = {"Auto": "auto", "DNA": "nt", "Protein": "aa"}[
            self.data_type_combo.currentText()
        ]
        topology = "ml" if self.topology_combo.currentIndex() == 0 else "mp"
        return {
            "data_type": data_type,
            "topology": topology,
            "schemes": self.schemes_combo.currentText(),
            "threads": self.threads_spinbox.value(),
        }

    def setup_output_tab(self):
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        selector_row = QHBoxLayout()
        selector_row.addWidget(QLabel("Result file:"))
        self.result_combo = QComboBox()
        self.result_combo.currentIndexChanged.connect(self._on_result_file_changed)
        selector_row.addWidget(self.result_combo, 1)
        layout.addLayout(selector_row)

        self.results_table = QTableWidget()
        self.results_table.setColumnCount(len(TABLE_HEADERS))
        self.results_table.setHorizontalHeaderLabels(TABLE_HEADERS)
        self.results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.results_table.setAlternatingRowColors(True)
        layout.addWidget(self.results_table, 1)

        self.output_info = QLabel("Model selection results will appear here.")
        self.output_info.setWordWrap(True)
        layout.addWidget(self.output_info)

    def setup_control_panel(self):
        super().setup_control_panel()
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.export_table_to_mpea)
        self.import_to_platform_btn.setVisible(False)
        self.control_layout.addWidget(self.import_to_platform_btn)

    def _selected_criterion(self):
        return self.criterion_combo.currentText()

    def _on_result_file_changed(self, index):
        self.refresh_display()

    def refresh_display(self, *_args):
        path = self.result_combo.currentData()
        if not path:
            return
        data = self.parsed_results.get(path) or parse_modeltest_output(path)
        criterion = self._selected_criterion()
        rows = data.get("criteria", {}).get(criterion, [])
        self.results_table.setRowCount(len(rows))
        for row_index, row in enumerate(rows):
            for col_index, header in enumerate(TABLE_HEADERS):
                value = row.get(header, "")
                if isinstance(value, float):
                    text = "%.4f" % value
                else:
                    text = str(value)
                item = QTableWidgetItem(text)
                if header != "Model":
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.results_table.setItem(row_index, col_index, item)
        self.results_table.resizeColumnsToContents()
        best = data.get("best", {}).get(criterion, {})
        best_model = best.get("Model", rows[0]["Model"] if rows else "N/A")
        score = best.get("Score", rows[0]["Score"] if rows else "N/A")
        if isinstance(score, float):
            score = "%.4f" % score
        self.output_info.setText(
            "%s · best model (%s): %s · score: %s · %d candidate model(s)"
            % (os.path.basename(path), criterion, best_model, score, len(rows))
        )

    def display_results(self, output_files):
        self.output_files = list(output_files or [])
        self.result_combo.blockSignals(True)
        self.result_combo.clear()
        for path in self.output_files:
            self.result_combo.addItem(os.path.basename(path), path)
        self.result_combo.blockSignals(False)
        if self.output_files:
            self.result_combo.setCurrentIndex(0)
            self.refresh_display()

    def run_analysis(self):
        if self.is_running:
            return
        input_files = self.prepare_input_files()
        if not input_files:
            return
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "ModelTest-NG executable not found.")
            return
        parameters = self.get_parameters()
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.add_console_message(
            "Starting ModelTest-NG (%s, %s schemes)..."
            % (parameters["data_type"], parameters["schemes"]), "info"
        )
        self.analysis_thread = ModelTestNGThread(
            self.tool_path, input_files,
            data_type=parameters["data_type"],
            topology=parameters["topology"],
            schemes=parameters["schemes"],
            threads=parameters["threads"],
            imported_files=self.imported_files,
            workdir=self.workdir,
        )
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()

    def _on_result(self, payload):
        self.parsed_results = dict(payload.get("parsed", {}))
        self.output_files = list(payload.get("output_files", []))

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.display_results(output_files)
        self.tab_widget.setCurrentIndex(1)
        self.import_to_platform_btn.setVisible(
            bool(self.result_has_consumer()) or self.import_from == "YR_MPEA"
        )
        self.add_console_message(
            "ModelTest-NG completed: %d report(s)" % len(output_files), "info"
        )

    def analysis_error(self, error_message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("ModelTest-NG failed: %s" % error_message, "error")
        QMessageBox.critical(self, "Error", "ModelTest-NG failed: %s" % error_message)

    def stop_analysis(self):
        thread = getattr(self, "analysis_thread", None)
        if thread is not None and thread.isRunning():
            for signal in (
                thread.progress, getattr(thread, "result_ready", None),
                thread.finished, thread.error, thread.console_output,
            ):
                if signal is None:
                    continue
                try:
                    signal.disconnect()
                except (TypeError, RuntimeError):
                    pass
            thread.stop()
            thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("ModelTest-NG stopped by user.", "warning")

    def _build_result_spec(self):
        files = list(self.output_files)
        if not files:
            raise ValueError("No ModelTest-NG report is available.")
        criterion = self._selected_criterion()
        items = []
        for path in files:
            data = self.parsed_results.get(path) or parse_modeltest_output(path)
            rows = data.get("criteria", {}).get(criterion, [])
            best = data.get("best", {}).get(criterion, {})
            best_model = best.get("Model") or (rows[0]["Model"] if rows else "")
            items.append(DatasetItemSpec(
                loci_name=Path(path).stem.replace(".modeltest", ""),
                item_type="model",
                file_path=str(path),
                format="modeltest-ng",
                data={
                    "type": "model_table",
                    "criterion": criterion,
                    "headers": TABLE_HEADERS,
                    "data": rows,
                    "best_model": best_model,
                    "best": best,
                    "criteria": data.get("criteria", {}),
                },
                metadata={
                    "source": self.plugin_name,
                    "tool": "ModelTest-NG",
                    "criterion": criterion,
                },
            ))
        dataset = DatasetSpec(
            name="ModelTest-NG models",
            items=items,
            description="%d model report(s)" % len(items),
            source="modeltest_ng",
            summary={"reports": len(items), "criterion": criterion},
            metadata={
                "input": self.imported_files[0] if len(self.imported_files) == 1 else "",
                "parameters": self.get_parameters(),
            },
        )
        return ResultSpec(
            source="modeltest_ng",
            datasets=[dataset],
            artifacts=[str(path) for path in files],
        )

    def export_table_to_mpea(self):
        if not self.result_has_consumer():
            self.add_console_message(
                "No workspace result channel is connected; open this plugin "
                "from YPhylo Workbench.", "warning"
            )
            return
        try:
            result = self._build_result_spec()
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        if self.export_result(result):
            self.import_to_platform_btn.setEnabled(False)
            self.import_to_platform_btn.setText("Sent to workspace")
            QMessageBox.information(
                self, "Success",
                "ModelTest-NG model report sent to the current workspace.",
            )


class ModelTestNGPluginEntry:
    """Plugin entry point."""

    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return ModelTestNGPlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )
