# raxmlng_plugin.py
#
# Copyright (c) 2026 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""RAxML-NG maximum-likelihood tree inference."""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QPushButton, QSpinBox,
    QVBoxLayout,
)

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread

DNA_MODELS = [
    "GTR+G4", "GTR+G", "HKY+G4", "TIM2+G4", "TVM+G4", "TPM2+G4",
    "K80+G4", "JC", "SYM+G4",
]
AA_MODELS = [
    "LG+G4", "WAG+G4", "JTT+G4", "DAYHOFF+G4", "BLOSUM62+G4",
    "cpREV+G4", "mtREV+G4", "VT+G4", "GTR+G4",
]


def normalize_model_name(model: str) -> str:
    """Translate IQ-TREE frequency tokens to RAxML-NG equivalents."""
    parts = []
    for part in str(model).split("+"):
        part = part.strip()
        if not part:
            continue
        if part == "FQ":
            part = "FE"  # equal base frequencies
        parts.append(part)
    return "+".join(parts)


def parse_final_loglh(log_path: str):
    try:
        text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    matches = re.findall(r"Final LogLikelihood:\s*([-+0-9.eE]+)", text)
    if not matches:
        return None
    try:
        return float(matches[-1])
    except ValueError:
        return None


class RAxMLNGThread(BaseProcessThread):
    """Run one RAxML-NG analysis per input alignment."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, model, threads=2,
                 bootstrap=0, auto_mre=False, start_trees="",
                 outgroup="", imported_files=None, workdir=None):
        super().__init__(
            tool_path, input_files, [], imported_files, workdir=workdir
        )
        self.model = model
        self.threads = threads
        self.bootstrap = bootstrap
        self.auto_mre = auto_mre
        self.start_trees = start_trees
        self.outgroup = outgroup

    def get_tool_name(self):
        return "RAxML-NG"

    def execute_commands(self):
        output_dir = tempfile.mkdtemp(
            prefix="raxml-ng-",
            dir=self.workdir if self.workdir and os.path.isdir(self.workdir) else None,
        )
        outputs = []
        total = len(self.input_files)
        for index, input_file in enumerate(self.input_files):
            if not self.is_running:
                break
            self.progress.emit("Processing file %d/%d..." % (index + 1, total))
            self.console_output.emit(
                "RAxML-NG: %s" % os.path.basename(input_file), "info"
            )
            prefix = os.path.join(
                output_dir, "%02d_%s" % (index + 1, Path(input_file).stem)
            )
            bootstrap_mode = self.bootstrap > 0 or self.auto_mre
            cmd = [
                self.tool_path,
                "--all" if bootstrap_mode else "--search",
                "--msa", input_file,
                "--model", self.model,
                "--threads", str(self.threads),
                "--prefix", prefix,
                "--redo",
            ]
            if bootstrap_mode:
                if self.auto_mre:
                    limit = self.bootstrap if self.bootstrap > 0 else 1000
                    cmd.extend(["--bs-trees", "autoMRE{%d}" % limit])
                else:
                    cmd.extend(["--bs-trees", str(self.bootstrap)])
            if self.start_trees:
                start_trees = self.start_trees
                if os.name == "nt":
                    # The MSYS2 build strips unescaped braces before parsing.
                    start_trees = (start_trees.replace("{", r"\{")
                                   .replace("}", r"\}"))
                cmd.extend(["--tree", start_trees])
            if self.outgroup.strip():
                cmd.extend(["--outgroup", self.outgroup.strip()])
            result = self.execute_command(cmd, cwd=output_dir)
            best_tree = prefix + ".raxml.bestTree"
            support_tree = prefix + ".raxml.support"
            if not os.path.exists(best_tree) and not os.path.exists(support_tree):
                detail = (result.stderr or "").strip()
                raise RuntimeError(
                    "RAxML-NG produced no tree for %s%s"
                    % (os.path.basename(input_file),
                       ": " + detail if detail else "")
                )
            log_path = prefix + ".raxml.log"
            bootstraps = prefix + ".raxml.bootstraps"
            outputs.append({
                "input": input_file,
                "prefix": prefix,
                "best_tree": best_tree if os.path.exists(best_tree) else "",
                "support_tree": support_tree if os.path.exists(support_tree) else "",
                "bootstraps": bootstraps if os.path.exists(bootstraps) else "",
                "best_model": prefix + ".raxml.bestModel",
                "log": log_path if os.path.exists(log_path) else "",
                "model": self.model,
                "bootstrap": (
                    "autoMRE{%d}" % (self.bootstrap if self.bootstrap > 0 else 1000)
                    if self.auto_mre else str(self.bootstrap)
                ) if bootstrap_mode else "0",
                "log_likelihood": parse_final_loglh(log_path),
            })
        self.progress.emit("RAxML-NG completed")
        self.result_ready.emit({"outputs": outputs})
        self.finished.emit(
            [path for output in outputs
             for path in (output["support_tree"] or output["best_tree"],)
             if path],
            [],
        )


class RAxMLNGPlugin(BasePlugin):
    """RAxML-NG tree inference plugin."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self.outputs = []
        super().__init__(import_from, import_data, **kwargs)
        if import_from == "DATASET_MANAGER" and import_data is not None:
            self.handle_dataset_import(import_data)
        self.run_button.setText("Run RAxML-NG")

    def init_plugin_info(self):
        self.plugin_name = "RAxML-NG"
        self.tool_name = "RAxML-NG"
        self.citation = []
        self.citation_keys = []
        self.input_types = {
            "FASTA": ["fas", "fna", "fa", "fasta"],
            "PHYLIP": ["phy"],
        }
        self.output_types = {"Newick": ".bestTree", "RAxML-NG": ".raxml.*"}
        self.plugin_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), ".."
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

        params_group = QGroupBox("RAxML-NG parameters")
        params_layout = QFormLayout(params_group)
        self.model_combo = QComboBox()
        self.model_combo.setEditable(True)
        self.model_combo.addItems(DNA_MODELS + AA_MODELS)
        self.model_combo.setCurrentText("GTR+G4")
        params_layout.addRow("Model:", self.model_combo)

        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 64)
        self.threads_spinbox.setValue(2)
        params_layout.addRow("Threads:", self.threads_spinbox)

        self.bootstrap_spinbox = QSpinBox()
        self.bootstrap_spinbox.setRange(0, 10000)
        self.bootstrap_spinbox.setSingleStep(100)
        self.bootstrap_spinbox.setValue(0)
        self.bootstrap_spinbox.setSpecialValueText("None")
        params_layout.addRow("Bootstrap replicates:", self.bootstrap_spinbox)

        self.auto_mre_checkbox = QCheckBox("Stop at autoMRE convergence")
        self.auto_mre_checkbox.setChecked(False)
        self.auto_mre_checkbox.setToolTip(
            "With 0 replicates, autoMRE{1000} is used. With N replicates, "
            "the run stops when convergence is reached or at N."
        )
        params_layout.addRow("", self.auto_mre_checkbox)

        self.tree_combo = QComboBox()
        self.tree_combo.addItems([
            "Default (10 pars + 10 rand)", "pars{10}", "rand{10}", "rand{1}",
        ])
        self.tree_combo.setToolTip(
            "Starting trees. Default lets RAxML-NG use 10 parsimony and 10 "
            "random trees; rand{1} is the quickest but least thorough."
        )
        params_layout.addRow("Starting trees:", self.tree_combo)

        self.outgroup_edit = QLineEdit()
        self.outgroup_edit.setPlaceholderText("comma-separated taxa (optional)")
        params_layout.addRow("Outgroup:", self.outgroup_edit)
        layout.addWidget(params_group)

        note = QLabel(
            "Model syntax follows RAxML-NG, e.g. GTR+G4 or LG+G4. A best-model "
            "object selected in the workspace is applied automatically when "
            "this window opens."
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        layout.addStretch(1)

        if getattr(self, "imported_files", None):
            self._set_imported_files(self.imported_files)

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

    def set_model(self, model):
        if model:
            self.model_combo.setCurrentText(normalize_model_name(model))

    def control_bindings(self):
        return [
            {"id": "msa", "widget": "file_path_edit",
             "label": "Input alignment", "kind": "path",
             "cli": ["--msa", "{value}"],
             "notes": "One RAxML-NG command is run per selected file."},
            {"id": "model", "widget": "model_combo",
             "label": "Model", "kind": "text",
             "cli": ["--model", "{value}"],
             "notes": "RAxML-NG model syntax, e.g. GTR+G4 or LG+G4."},
            {"id": "threads", "widget": "threads_spinbox",
             "label": "Threads", "kind": "int",
             "cli": ["--threads", "{value}"]},
            {"id": "bootstrap", "widget": "bootstrap_spinbox",
             "label": "Bootstrap replicates", "kind": "int",
             "step": 100,
             "cli": ["--bs-trees", "{value}"],
             "notes": "0 means no bootstrap unless autoMRE is enabled; "
                      "--all is used when bootstrapping is active."},
            {"id": "auto_mre", "widget": "auto_mre_checkbox",
             "label": "autoMRE convergence", "kind": "bool",
             "cli": [],
             "notes": "When enabled, --bs-trees autoMRE{N} is used and the "
                      "replicate count becomes the upper limit."},
            {"id": "start_trees", "widget": "tree_combo",
             "label": "Starting trees", "kind": "choice",
             "choices": {"Default (10 pars + 10 rand)": "",
                         "pars{10}": "pars{10}",
                         "rand{10}": "rand{10}",
                         "rand{1}": "rand{1}"},
             "cli": ["--tree", "{value}"],
             "notes": "Empty value omits --tree and uses RAxML-NG defaults."},
            {"id": "outgroup", "widget": "outgroup_edit",
             "label": "Outgroup", "kind": "text",
             "cli": ["--outgroup", "{value}"],
             "notes": "Comma-separated taxa; RAxML-NG treats this as a "
                      "drawing option."},
        ]

    def cli_preview(self):
        parameters = self.get_parameters()
        input_file = (self.imported_files[0] if self.imported_files
                      else self.file_path_edit.text().strip() or "<input>")
        bootstrap_mode = (parameters["bootstrap"] > 0
                          or parameters["auto_mre"])
        command = [
            self.tool_path or "raxml-ng.exe",
            "--all" if bootstrap_mode else "--search",
            "--msa", input_file,
            "--model", parameters["model"] or "<model>",
            "--threads", str(parameters["threads"]),
            "--prefix", "<output-prefix>",
            "--redo",
        ]
        if bootstrap_mode:
            if parameters["auto_mre"]:
                limit = parameters["bootstrap"] or 1000
                command.extend(["--bs-trees", "autoMRE{%d}" % limit])
            else:
                command.extend(["--bs-trees", str(parameters["bootstrap"])])
        if parameters["start_trees"]:
            command.extend(["--tree", parameters["start_trees"]])
        if parameters["outgroup"]:
            command.extend(["--outgroup", parameters["outgroup"]])
        return command

    def get_parameters(self):
        return {
            "model": self.model_combo.currentText().strip(),
            "threads": self.threads_spinbox.value(),
            "bootstrap": self.bootstrap_spinbox.value(),
            "auto_mre": self.auto_mre_checkbox.isChecked(),
            "start_trees": (
                "" if self.tree_combo.currentIndex() == 0
                else self.tree_combo.currentText()
            ),
            "outgroup": self.outgroup_edit.text().strip(),
        }

    def setup_output_tab(self):
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        self.output_text = QPlainTextEdit()
        self.output_text.setReadOnly(True)
        self.output_text.setPlaceholderText(
            "Tree paths and run statistics will appear here."
        )
        layout.addWidget(self.output_text, 1)
        self.output_info = QLabel("No RAxML-NG run yet.")
        self.output_info.setWordWrap(True)
        layout.addWidget(self.output_info)

    def setup_control_panel(self):
        super().setup_control_panel()
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)
        self.control_layout.addWidget(self.import_to_platform_btn)

    def run_analysis(self):
        if self.is_running:
            return
        input_files = self.prepare_input_files()
        if not input_files:
            return
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "RAxML-NG executable not found.")
            return
        parameters = self.get_parameters()
        parameters["model"] = normalize_model_name(parameters["model"])
        if not parameters["model"]:
            QMessageBox.warning(self, "Warning", "Please specify a model.")
            return
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.add_console_message(
            "Starting RAxML-NG with model %s..." % parameters["model"], "info"
        )
        self.analysis_thread = RAxMLNGThread(
            self.tool_path,
            input_files,
            parameters["model"],
            threads=parameters["threads"],
            bootstrap=parameters["bootstrap"],
            auto_mre=parameters["auto_mre"],
            start_trees=parameters["start_trees"],
            outgroup=parameters["outgroup"],
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
        self.outputs = list(payload.get("outputs", []))

    def display_results(self):
        lines = []
        for output in self.outputs:
            tree = output["support_tree"] or output["best_tree"]
            log_lh = output.get("log_likelihood")
            log_lh_text = "%.4f" % log_lh if log_lh is not None else "N/A"
            lines.append(
                "Input: %s\nModel: %s\nFinal logL: %s\nBootstrap: %s\n"
                "Tree: %s\nBootstrap trees: %s\nLog: %s\n"
                % (
                    os.path.basename(output["input"]),
                    output["model"],
                    log_lh_text,
                    output["bootstrap"],
                    tree or "N/A",
                    output["bootstraps"] or "N/A",
                    output["log"] or "N/A",
                )
            )
        self.output_text.setPlainText("\n".join(lines))
        self.output_info.setText(
            "%d RAxML-NG run(s) completed." % len(self.outputs)
        )

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.display_results()
        self.tab_widget.setCurrentIndex(1)
        self.import_to_platform_btn.setVisible(
            bool(self.result_has_consumer()) or self.import_from == "YR_MPEA"
        )
        self.add_console_message(
            "RAxML-NG completed: %d tree(s)" % len(self.outputs), "info"
        )

    def analysis_error(self, error_message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("RAxML-NG failed: %s" % error_message, "error")
        QMessageBox.critical(self, "Error", "RAxML-NG failed: %s" % error_message)

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
        self.add_console_message("RAxML-NG stopped by user.", "warning")

    def _build_result_spec(self):
        if not self.outputs:
            raise ValueError("No RAxML-NG result is available.")
        items = []
        artifacts = []
        for output in self.outputs:
            stem = Path(output["input"]).stem
            tree = output["support_tree"] or output["best_tree"]
            if tree:
                items.append(DatasetItemSpec(
                    loci_name=stem,
                    item_type="phylogeny",
                    file_path=str(tree),
                    format="newick",
                    data={
                        "model": output["model"],
                        "bootstrap": output["bootstrap"],
                        "log_likelihood": output.get("log_likelihood"),
                        "best_tree": output["best_tree"],
                        "support_tree": output["support_tree"],
                    },
                    metadata={"source": self.plugin_name, "tool": "RAxML-NG"},
                ))
                artifacts.append(str(tree))
            if output["bootstraps"]:
                items.append(DatasetItemSpec(
                    loci_name="%s_bootstrap" % stem,
                    item_type="treeset",
                    file_path=str(output["bootstraps"]),
                    format="newick",
                    data={
                        "model": output["model"],
                        "bootstrap": output["bootstrap"],
                        "kind": "bootstrap",
                    },
                    metadata={"source": self.plugin_name, "tool": "RAxML-NG"},
                ))
                artifacts.append(str(output["bootstraps"]))
            for key in ("best_tree", "log", "best_model"):
                if output.get(key):
                    artifacts.append(str(output[key]))
        if not items:
            raise ValueError("RAxML-NG did not produce an importable tree.")
        dataset = DatasetSpec(
            name="RAxML-NG phylogenies",
            items=items,
            description="%d RAxML-NG result(s)" % len(self.outputs),
            source="raxmlng",
            summary={
                "runs": len(self.outputs),
                "bootstrap": self.get_parameters()["bootstrap"],
                "auto_mre": self.get_parameters()["auto_mre"],
            },
            metadata={
                "input": self.imported_files[0] if len(self.imported_files) == 1 else "",
                "parameters": self.get_parameters(),
            },
        )
        return ResultSpec(
            source="raxmlng",
            datasets=[dataset],
            artifacts=artifacts,
        )

    def import_to_platform(self):
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
                "RAxML-NG trees sent to the current workspace.",
            )


class RAxMLNGPluginEntry:
    """Plugin entry point."""

    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return RAxMLNGPlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )
