# clustal_omega_plugin.py
#
# Copyright (c) 2026 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Clustal Omega aligner, standard ResultSpec plugin."""

from __future__ import annotations

import os
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFormLayout, QFrame, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QPushButton,
    QRadioButton, QScrollArea, QSpinBox, QTextEdit, QVBoxLayout, QWidget,
)

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.alignment_common import (
    build_result, prepare_import, records_temp, setup_preview, show_preview,
)


class ClustalOmegaThread(BaseProcessThread):
    """Run Clustal Omega and return validated FASTA alignments."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters,
                 imported_files=None, suffix=None, workdir=None):
        super().__init__(
            tool_path, input_files, parameters, imported_files, workdir=workdir
        )
        self.output_suffix = suffix or ""

    def get_tool_name(self):
        return "Clustal Omega"

    def _output_path(self, input_file):
        if self.output_suffix:
            return Path("%s%s.fas" % (input_file, self.output_suffix))
        return Path(self.create_temp_file(suffix=".fas"))

    def execute_commands(self):
        try:
            output_files = []
            total = len(self.input_files)
            for index, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                self.progress.emit("Processing file %d/%d..." % (index + 1, total))
                self.console_output.emit(
                    "Processing file %d/%d: %s"
                    % (index + 1, total, os.path.basename(input_file)),
                    "info",
                )
                output_file = self._output_path(input_file)
                output_file.parent.mkdir(parents=True, exist_ok=True)
                cmd = [
                    self.tool_path, "-i", input_file, "-o", str(output_file),
                    "--outfmt=fa", *self.parameters,
                ]
                result = self.execute_command(
                    cmd, cwd=str(output_file.parent)
                )
                detail = (result.stderr or "").strip()
                if result.returncode != 0:
                    raise RuntimeError(
                        "Clustal Omega failed for %s%s"
                        % (os.path.basename(input_file),
                           ": " + detail if detail else "")
                    )
                if not output_file.is_file() or output_file.stat().st_size == 0:
                    raise RuntimeError(
                        "Clustal Omega did not write %s" % output_file
                    )
                output_files.append(str(output_file))

            self.progress.emit("All alignments completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])
        except Exception as exc:
            self.error.emit("Alignment exception: %s" % exc)


class ClustalOmegaPlugin(BasePlugin):
    """Clustal Omega alignment panel."""

    # Deprecated compatibility channels; the host uses result_ready.
    import_alignment_signal = pyqtSignal(list)
    batch_import_alignment_signal = pyqtSignal(list)

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self.batch_mode = False
        self._payload = None
        super().__init__(import_from, import_data, **kwargs)
        if import_data is not None and (
            not isinstance(import_data, list) or import_from != "YR_MPEA"
        ):
            self.handle_import_data(import_data)
        self.run_button.setText("Run Clustal Omega")

    def init_plugin_info(self):
        self.plugin_name = "Clustal Omega Aligner"
        self.tool_name = "Clustal Omega"
        self.citation = []
        self.citation_keys = ["sievers2011fast"]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"Alignment FASTA": ["fas"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    # ------------------------------------------------------------------
    # input
    # ------------------------------------------------------------------
    def control_bindings(self):
        """GUI controls mapped to the Clustal Omega command line."""
        return [
            {"id": "input", "widget": "file_path_edit",
             "label": "Input FASTA", "kind": "path",
             "cli": ["-i", "{value}"],
             "notes": "One Clustal Omega run per file in batch mode."},
            {"id": "auto", "widget": "auto_checkbox",
             "label": "Auto options", "kind": "bool", "cli": ["--auto"],
             "notes": "Default on; Clustal chooses guide-tree and HMM settings."},
            {"id": "seq_type", "widget": "seqtype_combo",
             "label": "Sequence type", "kind": "choice",
             "choices": {"Auto": "", "Protein": "Protein",
                         "RNA": "RNA", "DNA": "DNA"},
             "cli": ["--seqtype", "{value}"]},
            {"id": "threads", "widget": "thread_spinbox",
             "label": "Threads", "kind": "int",
             "cli": ["--threads", "{value}"], "omit_values": [1]},
            {"id": "order", "widget": "order_combo",
             "label": "Output order", "kind": "choice",
             "choices": {"input-order": "", "tree-order": "--outorder=tree-order"},
             "cli": ["{value}"]},
        ]

    def setup_input_tab(self):
        outer = QVBoxLayout(self.input_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        outer.addWidget(scroll)
        content = QWidget()
        scroll.setWidget(content)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        title = QLabel("Clustal Omega")
        title.setStyleSheet("font-size: 18px; font-weight: 600;")
        layout.addWidget(title)
        subtitle = QLabel(
            "Align one or more FASTA files. Output is written as aligned FASTA "
            "and can be sent to the workspace as a ResultSpec."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        layout.addWidget(subtitle)

        input_group = QGroupBox("Input")
        form = QFormLayout(input_group)
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select FASTA files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        row = QHBoxLayout()
        row.addWidget(self.file_path_edit, 1)
        row.addWidget(self.file_browse_btn)
        form.addRow("Files:", row)

        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        form.addRow("", self.file_tags_container)

        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste FASTA sequences...")
        self.sequence_text.setMaximumHeight(160)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        form.addRow("Paste:", self.sequence_text)
        layout.addWidget(input_group)

        params = QGroupBox("Parameters")
        pform = QFormLayout(params)
        self.auto_checkbox = QCheckBox("Auto options (--auto)")
        self.auto_checkbox.setChecked(True)
        pform.addRow("", self.auto_checkbox)
        self.seqtype_combo = QComboBox()
        self.seqtype_combo.addItems(["Auto", "Protein", "RNA", "DNA"])
        pform.addRow("Sequence type:", self.seqtype_combo)
        self.thread_spinbox = QSpinBox()
        self.thread_spinbox.setRange(1, 64)
        self.thread_spinbox.setValue(4)
        pform.addRow("Threads:", self.thread_spinbox)
        self.order_combo = QComboBox()
        self.order_combo.addItems(["input-order", "tree-order"])
        pform.addRow("Output order:", self.order_combo)
        layout.addWidget(params)

        output_group = QGroupBox("Output")
        oform = QFormLayout(output_group)
        mode_row = QHBoxLayout()
        self.save_to_cwd = QRadioButton("Beside input")
        self.save_to_cwd.setChecked(True)
        self.save_to_tmp = QRadioButton("Temporary file")
        mode_row.addWidget(self.save_to_cwd)
        mode_row.addWidget(self.save_to_tmp)
        oform.addRow("Mode:", mode_row)
        self.output_suffix_edit = QLineEdit("_clustal")
        self.output_suffix_edit.setPlaceholderText("Output suffix")
        oform.addRow("Suffix:", self.output_suffix_edit)
        self.save_to_cwd.toggled.connect(self.on_output_radio_changed)
        layout.addWidget(output_group)
        layout.addStretch()

        self.file_tags = []

    def on_output_radio_changed(self):
        enabled = self.save_to_cwd.isChecked()
        self.output_suffix_edit.setEnabled(enabled)
        self.output_suffix_edit.setVisible(enabled)

    def setup_output_tab(self):
        setup_preview(self, "The aligned FASTA output will appear here.")

    def setup_control_panel(self):
        super().setup_control_panel()
        self.import_to_platform_btn = QPushButton("Send to workspace")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)
        self.control_layout.addWidget(self.import_to_platform_btn)

    def handle_import_data(self, import_data):
        files, batch = prepare_import(self, import_data)
        self.batch_mode = batch
        self.imported_files = files
        self.import_file = files[0] if len(files) == 1 else None
        if hasattr(self, "file_path_edit"):
            if not files:
                self.file_path_edit.clear()
            elif len(files) == 1:
                self.file_path_edit.setText(files[0])
            else:
                self.file_path_edit.setText("%d files selected" % len(files))
                self.file_path_edit.setEnabled(False)

    # ------------------------------------------------------------------
    # input helpers
    # ------------------------------------------------------------------
    def browse_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA files", "",
            "FASTA files (*.fas *.fasta *.fa *.fna);;All files (*)",
        )
        for path in paths:
            if path not in self.imported_files:
                self.add_file_tag(path)
        if self.imported_files:
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText("%d files selected" % len(self.imported_files))
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)

    def add_file_tag(self, file_path):
        self.imported_files.append(file_path)
        self.import_file = self.imported_files[0] if len(self.imported_files) == 1 else None
        self.batch_mode = len(self.imported_files) > 1
        tag = QFrame()
        tag.setStyleSheet(
            "QFrame { background-color: palette(base); border: 1px solid palette(mid);"
            " border-radius: 12px; margin: 2px; }"
        )
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag.setLayout(tag_layout)
        label = QLabel(self.get_display_name(file_path))
        tag_layout.addWidget(label)
        close = QPushButton("×")
        close.setFlat(True)
        close.setFixedSize(20, 20)
        close.clicked.connect(lambda: self.remove_file_tag(file_path, tag))
        tag_layout.addWidget(close)
        self.file_tags_layout.addWidget(tag)
        self.file_tags.append((file_path, tag))
        self.file_tags_container.setVisible(True)

    def remove_file_tag(self, file_path, tag_widget):
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        tag_widget.setParent(None)
        tag_widget.deleteLater()
        self.batch_mode = len(self.imported_files) > 1
        self.import_file = self.imported_files[0] if len(self.imported_files) == 1 else None
        if self.imported_files:
            self.file_path_edit.setText(
                self.imported_files[0] if len(self.imported_files) == 1
                else "%d files selected" % len(self.imported_files)
            )
        else:
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)

    def clear_all_file_tags(self):
        for _, tag in self.file_tags:
            tag.setParent(None)
            tag.deleteLater()
        self.file_tags = []
        self.imported_files = []
        self.import_file = None
        self.batch_mode = False
        self.file_tags_container.setVisible(False)

    def get_display_name(self, file_path):
        filename = os.path.basename(file_path)
        duplicates = [
            path for path in self.imported_files
            if os.path.basename(path) == filename and path != file_path
        ]
        if duplicates:
            return "%s (%s)" % (filename, os.path.basename(os.path.dirname(file_path)))
        return filename

    def on_text_changed(self):
        text = self.sequence_text.toPlainText().strip()
        if text:
            self.file_path_edit.setVisible(False)
            self.file_browse_btn.setVisible(False)
            self.file_tags_container.setVisible(False)
            self.save_to_tmp.setChecked(True)
            self.save_to_cwd.setEnabled(False)
            self.clear_all_file_tags()
        else:
            self.file_path_edit.setVisible(True)
            self.file_browse_btn.setVisible(True)
            self.save_to_cwd.setEnabled(True)
            if self.imported_files:
                self.file_tags_container.setVisible(True)

    def prepare_input_files(self):
        if self.imported_files:
            return [path for path in self.imported_files if os.path.exists(path)]
        text = self.sequence_text.toPlainText().strip()
        if not text:
            QMessageBox.warning(self, "Input missing", "Select files or paste FASTA sequences.")
            return None
        try:
            records = list(SeqIO.parse(StringIO(text), "fasta"))
        except Exception as exc:
            QMessageBox.warning(self, "Invalid FASTA", str(exc))
            return None
        if len(records) < 2:
            QMessageBox.warning(self, "Need more sequences", "Paste at least two FASTA records.")
            return None
        return [records_temp(self, records, "pasted")]

    def get_parameters(self):
        params = []
        if self.auto_checkbox.isChecked():
            params.append("--auto")
        seqtype = self.seqtype_combo.currentText()
        if seqtype != "Auto":
            params.extend(["--seqtype", seqtype])
        if self.thread_spinbox.value() > 1:
            params.extend(["--threads", str(self.thread_spinbox.value())])
        if self.order_combo.currentText() == "tree-order":
            params.append("--outorder=tree-order")
        return params

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    def run_analysis(self):
        if self.is_running:
            return
        if not self.tool_path or not Path(self.tool_path).is_file():
            QMessageBox.critical(self, "Tool missing", "Clustal Omega executable not found.")
            return
        input_files = self.prepare_input_files()
        if not input_files:
            return

        suffix = None
        if self.save_to_cwd.isChecked():
            suffix = self.output_suffix_edit.text()
            if any(sep in suffix for sep in ("/", "\\")):
                QMessageBox.warning(self, "Invalid suffix", "Suffix must not contain a path separator.")
                return

        self._payload = None
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running Clustal Omega...")
        self.import_to_platform_btn.setEnabled(True)
        self.import_to_platform_btn.setText("Send to workspace")
        self.tab_widget.setCurrentIndex(2)

        self.alignment_thread = ClustalOmegaThread(
            self.tool_path, input_files, self.get_parameters(),
            self.imported_files, suffix=suffix, workdir=self.workdir,
        )
        self.alignment_thread.progress.connect(self.progress_bar.setFormat)
        self.alignment_thread.result_ready.connect(self._on_result)
        self.alignment_thread.finished.connect(self.analysis_finished)
        self.alignment_thread.error.connect(self.analysis_error)
        self.alignment_thread.console_output.connect(self.add_console_message)
        self.alignment_thread.start()

    def _on_result(self, payload):
        self._payload = payload

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.alignment_output_files = output_files
        self._show_alignment_preview(output_files)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("Clustal Omega alignment completed.", "info")
        if self.result_has_consumer() or self.import_from in [
            "YR_MPEA", "seq_viewer", "DATASET_MANAGER"
        ]:
            self.import_to_platform_btn.setVisible(True)
        self.import_to_platform_btn.setEnabled(True)

    def _show_alignment_preview(self, output_files):
        show_preview(self, output_files)

    def analysis_error(self, error_msg):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Alignment failed", error_msg)
        self.tab_widget.setCurrentIndex(2)

    def stop_analysis(self):
        if hasattr(self, "alignment_thread") and self.alignment_thread.isRunning():
            thread = self.alignment_thread
            for signal in (
                thread.progress, thread.result_ready, thread.finished,
                thread.error, thread.console_output,
            ):
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
        self.add_console_message("Alignment stopped by user.", "warning")

    # ------------------------------------------------------------------
    # result
    # ------------------------------------------------------------------
    def _build_result_spec(self):
        if not getattr(self, "alignment_output_files", None):
            raise ValueError("No Clustal Omega output is available.")
        suffix = (self.output_suffix_edit.text()
                  if self.save_to_cwd.isChecked() else "")
        return build_result(
            tool="Clustal Omega",
            source="clustal_omega",
            output_files=self.alignment_output_files,
            input_files=self.imported_files,
            parameters=self.get_parameters(),
            mode=self.seqtype_combo.currentText(),
            batch=self.batch_mode,
            suffix=suffix,
        )

    def import_to_platform(self):
        if not getattr(self, "alignment_output_files", None):
            QMessageBox.warning(self, "No result", "Run Clustal Omega first.")
            return
        if self.result_has_consumer():
            try:
                result = self._build_result_spec()
            except Exception as exc:
                QMessageBox.critical(self, "Result failed", str(exc))
                return
            if self.export_result(result):
                self.import_to_platform_btn.setEnabled(False)
                self.import_to_platform_btn.setText("Sent to workspace")
                QMessageBox.information(
                    self, "Success", "Alignment result sent to the current workspace."
                )
            return

        # Legacy fallback for external callers.
        try:
            from Bio import SeqIO
            if self.batch_mode:
                batch = [list(SeqIO.parse(path, "fasta"))
                         for path in self.alignment_output_files]
                self.batch_import_alignment_signal.emit(batch)
            else:
                sequences = []
                for path in self.alignment_output_files:
                    sequences.extend(list(SeqIO.parse(path, "fasta")))
                if sequences:
                    self.import_alignment_signal.emit(sequences)
        except Exception as exc:
            QMessageBox.critical(self, "Import failed", str(exc))


class ClustalOmegaPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return ClustalOmegaPlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )
