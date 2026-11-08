# mafft_plugin.py
#
# Copyright (c) 2025 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# YR-MPE program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit,
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit, QPlainTextEdit, QDialog,
                             QDoubleSpinBox, QRadioButton)
from PyQt5.QtCore import Qt, pyqtSignal
import os
from pathlib import Path

from Bio import AlignIO, SeqIO

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread


def _alignment_locus_name(output_file, output_suffix=""):
    """Derive a stable locus name from a MAFFT output path.

    MAFFT output keeps the whole input name and appends the GUI suffix and
    ``.fas`` (``OG0000008.cds.fna_mafft.fas``). Dataset items are named with
    compound sequence suffixes stripped (``OG0000008``), so normalize here to
    let alignment results update the item they came from.
    """
    name = os.path.basename(str(output_file))
    name = os.path.splitext(name)[0]
    if output_suffix and name.endswith(output_suffix):
        name = name[: -len(output_suffix)]
    lower = name.lower()
    for suffix in (
        ".cds.fna", ".cds.fasta", ".cds.fa",
        ".fna", ".fasta", ".fas", ".fa",
        ".phy", ".nex", ".nexus",
    ):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _strip_direction_prefix(records):
    """Remove MAFFT's ``_R_`` marker from reversed-direction records.

    ``--adjustdirection`` prefixes the affected sequence names with ``_R_``.
    The marker is useful inside MAFFT but must not become part of the taxon
    label used by downstream dataset/concatenation code.
    """
    for record in records:
        for attribute in ("id", "name", "description"):
            value = str(getattr(record, attribute, "") or "")
            if value.startswith("_R_"):
                setattr(record, attribute, value[3:])
    return records


class MAFFTThread(BaseProcessThread):
    """MAFFT alignment worker."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters, imported_files=None,
                 suffix=None, workdir=None):
        super().__init__(
            tool_path, input_files, parameters, imported_files, workdir=workdir
        )
        self.output_suffix = suffix or ""
        self.clustal_output = "--clustalout" in (parameters or [])

    def get_tool_name(self):
        return "MAFFT"

    def _output_path(self, input_file):
        if self.output_suffix:
            suffix = ".aln" if self.clustal_output else ".fas"
            return Path(f"{input_file}{self.output_suffix}{suffix}")
        suffix = ".aln" if self.clustal_output else ".fas"
        return Path(self.create_temp_file(suffix=suffix))

    def _strip_reversed_headers(self, output_file):
        if not any(str(param).startswith("--adjustdirection")
                   for param in self.parameters):
            return
        try:
            if self.clustal_output:
                alignment = AlignIO.read(str(output_file), "clustal")
                _strip_direction_prefix(alignment)
                AlignIO.write(alignment, str(output_file), "clustal")
            else:
                records = list(SeqIO.parse(str(output_file), "fasta"))
                _strip_direction_prefix(records)
                SeqIO.write(records, str(output_file), "fasta")
        except Exception as exc:
            self.console_output.emit(
                "Could not normalize MAFFT direction prefixes: %s" % exc,
                "warning",
            )

    def execute_commands(self):
        try:
            output_files = []
            total_files = len(self.input_files)
            for index, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break

                self.progress.emit(
                    "Processing file %d/%d..." % (index + 1, total_files)
                )
                self.console_output.emit(
                    "Processing file %d/%d: %s"
                    % (index + 1, total_files, os.path.basename(input_file)),
                    "info",
                )

                output_file = self._output_path(input_file)
                output_file.parent.mkdir(parents=True, exist_ok=True)
                cmd = [self.tool_path, *self.parameters, input_file]

                result = self.execute_command(
                    cmd, cwd=str(output_file.parent)
                )
                detail = (result.stderr or "").strip()
                if result.returncode != 0:
                    raise RuntimeError(
                        "MAFFT failed for %s%s"
                        % (
                            os.path.basename(input_file),
                            ": " + detail if detail else "",
                        )
                    )
                if not (result.stdout or "").strip():
                    raise RuntimeError(
                        "MAFFT produced empty output for %s"
                        % os.path.basename(input_file)
                    )

                output_file.write_text(result.stdout, encoding="utf-8")
                self._strip_reversed_headers(output_file)
                if output_file.stat().st_size == 0:
                    raise RuntimeError(
                        "MAFFT wrote an empty file: %s" % output_file
                    )
                output_files.append(str(output_file))

            self.progress.emit("All alignments completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])

        except Exception as exc:
            self.error.emit("Alignment exception: %s" % exc)


class MAFFTPlugin(BasePlugin):
    """MAFFT plugin: standard ResultSpec output, legacy signals retained."""
    # Deprecated compatibility channels. New callers use BasePlugin.result_ready
    # and the host's _connect_plugin_result(); the legacy path remains only for
    # external callers that still connect these signals.
    import_alignment_signal = pyqtSignal(list)
    batch_import_alignment_signal = pyqtSignal(list)
    def __init__(self, import_from=None, import_data=None, **kwargs):
        """Initialize the MAFFT plugin."""
        self.batch_mode = False
        # When MAFFT is launched from a workspace dataset, write the
        # alignments back into that dataset instead of creating a detached
        # "MAFFT alignment" dataset.
        self._source_dataset_id = None
        super().__init__(import_from, import_data, **kwargs)
        # BasePlugin handles a plain list only when import_from == "YR_MPEA".
        # Structured payloads (dataset items, paths, batched lists) need the
        # MAFFT-specific handler here.
        if import_data is not None and (
            not isinstance(import_data, list) or import_from != "YR_MPEA"
        ):
            self.handle_import_data(import_data)

    def _records_temp(self, records, label="import"):
        temp_file = self.create_temp_file(suffix="_%s.fas" % label)
        with open(temp_file, "w", encoding="utf-8") as handle:
            SeqIO.write(list(records), handle, "fasta")
        self.temp_files.append(temp_file)
        return temp_file

    def _items_to_files(self, items):
        files = []
        for item in items:
            path = getattr(item, "file_path", "") or ""
            if path and os.path.exists(path):
                files.append(path)
                continue
            records = list(getattr(item, "sequences", []) or [])
            if records:
                files.append(self._records_temp(records, "dataset_item"))
        return files

    def handle_import_data(self, import_data):
        """Accept SeqRecords, batched lists, paths, or dataset-item dicts."""
        files = []
        batch_mode = False

        if isinstance(import_data, dict):
            dataset_config = import_data.get("dataset_config", {}) or {}
            self._source_dataset_id = dataset_config.get("dataset_id")
            items = import_data.get("dataset_items")
            if items is not None:
                files = self._items_to_files(items)
                batch_mode = len(files) > 1
            elif import_data.get("file_path"):
                files = [str(import_data["file_path"])]
            elif import_data.get("files"):
                files = [str(path) for path in import_data["files"]]
            elif import_data.get("sequences"):
                files = [self._records_temp(import_data["sequences"], "import")]
        elif isinstance(import_data, str):
            files = [import_data]
        elif isinstance(import_data, (list, tuple)):
            if import_data and isinstance(import_data[0], (list, tuple)):
                batch_mode = True
                for index, records in enumerate(import_data):
                    files.append(self._records_temp(records, "batch_%d" % index))
            elif import_data and isinstance(import_data[0], dict):
                files = self._items_to_files(import_data)
                batch_mode = len(files) > 1
            elif import_data and hasattr(import_data[0], "seq"):
                files = [self._records_temp(import_data, "import")]

        files = [path for path in files if path and os.path.exists(path)]
        self.batch_mode = batch_mode
        self.imported_files = files
        self.import_file = files[0] if len(files) == 1 else None

        if hasattr(self, "file_path_edit") and self.file_path_edit is not None:
            if not files:
                self.file_path_edit.clear()
            elif len(files) == 1:
                self.file_path_edit.setText(files[0])
            else:
                self.file_path_edit.setText("%d files selected" % len(files))
                self.file_path_edit.setEnabled(False)

    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "MAFFT Aligner"
        self.tool_name = "MAFFT"
        self.citation = []  # citation text lives in YR_MPE/citations/*.bib
        self.citation_keys = ["katoh2013mafft"]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"FASTA": ".fas"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # 初始化文件标签列表
        self.file_tags = []
        
        # 输入组
        input_group = QGroupBox("Input / Output")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # 文件输入
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select FASTA files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)
        
        # 文件标签容器
        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # 文本输入
        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste FASTA format sequence...")
        self.sequence_text.setMaximumHeight(200)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        input_layout.addRow("Sequence text:", self.sequence_text)
        
        # 处理导入的数据
        if self.import_file:
            self.file_path_edit.setText(self.import_file)
            input_group.setVisible(False)
        
        # 基本参数组
        basic_params_group = QGroupBox("Basic Parameters")
        basic_params_layout = QFormLayout()
        basic_params_group.setLayout(basic_params_layout)
        layout.addWidget(basic_params_group)
        
        # 比对模式
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Auto mode (--auto)",
            "High speed (default)",
            "High accuracy - localpair (--localpair)",
            "High accuracy - genafpair (--genafpair)", 
            "High accuracy - globalpair (--globalpair)",
            "Fast (--retree 1)"
        ])
        basic_params_layout.addRow("Alignment mode:", self.mode_combo)
        
        # 线程数
        self.thread_spinbox = QSpinBox()
        self.thread_spinbox.setRange(-1, 32)
        self.thread_spinbox.setValue(-1)  # 默认: 自动
        self.thread_spinbox.setToolTip("-1 for automatic thread detection")
        basic_params_layout.addRow("Threads:", self.thread_spinbox)
        
        # 高级参数按钮
        advanced_btn = QPushButton("Advanced Parameters...")
        advanced_btn.clicked.connect(self.show_advanced_dialog)
        basic_params_layout.addRow("", advanced_btn)
        
        # Output settings
        output_settings_group = QWidget()
        output_settings_group.setLayout(QVBoxLayout())
        output_settings_group.layout().setContentsMargins(0, 0, 0, 0)
        output_settings_group.setContentsMargins(0, 0, 0, 0)
        output_radio_group = QWidget()
        output_radio_group.setLayout(QHBoxLayout())
        self.save_to_cwd = QRadioButton("Beside input")
        self.save_to_cwd.setChecked(True)
        self.save_to_tmp = QRadioButton("Temporary File")
        output_radio_group.layout().addWidget(self.save_to_cwd)
        output_radio_group.layout().addWidget(self.save_to_tmp)
        output_settings_group.layout().addWidget(output_radio_group)

        self.save_to_cwd.toggled.connect(self.on_output_radio_changed)
        self.save_to_tmp.toggled.connect(self.on_output_radio_changed)

        self.output_suffix_edit = QLineEdit('_mafft')
        self.output_suffix_edit.setPlaceholderText("Output suffix")
        output_settings_group.layout().addWidget(self.output_suffix_edit)
        input_layout.addRow("Save to:", output_settings_group)
        
        layout.addStretch()

    def setup_output_tab(self):
        """Show a lightweight FASTA preview without embedding a web view."""
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.summary_label = QLabel("No alignment yet.")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.summary_label)
        self.preview_text = QPlainTextEdit()
        self.preview_text.setReadOnly(True)
        self.preview_text.setPlaceholderText(
            "The aligned FASTA output will appear here."
        )
        layout.addWidget(self.preview_text, 1)

    def _show_alignment_preview(self, output_files):
        if not output_files:
            self.summary_label.setText("No alignment output.")
            self.preview_text.clear()
            return
        first = Path(output_files[0])
        self.summary_label.setText(
            "%d alignment file(s) · %s"
            % (len(output_files), first.name)
        )
        try:
            text = first.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            self.preview_text.setPlainText("Cannot read %s: %s" % (first, exc))
            return
        if len(text) > 200_000:
            text = text[:200_000] + "\n... (preview truncated)"
        self.preview_text.setPlainText(text)

    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)
    
    def show_advanced_dialog(self):
        """显示高级参数对话框"""
        dialog = MAFFTAdvancedDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            # 存储高级参数
            self.advanced_params = dialog.get_parameters()

    def on_output_radio_changed(self):
        """处理输出选项的切换"""
        if self.save_to_cwd.isChecked():
            self.output_suffix_edit.setEnabled(True)
            self.output_suffix_edit.setVisible(True)
        elif self.save_to_tmp.isChecked():
            self.output_suffix_edit.setEnabled(False)
            self.output_suffix_edit.setVisible(False)

    def control_bindings(self):
        return [
            {"id": "input", "widget": "file_path_edit",
             "label": "Input sequences", "kind": "path",
             "cli": ["{value}"],
             "notes": "Positional input appended after all options."},
            {"id": "mode", "widget": "mode_combo",
             "label": "Alignment mode", "kind": "choice",
             "choices": {
                 "Auto mode (--auto)": ["--auto"],
                 "High speed (default)": [],
                 "High accuracy - localpair (--localpair)":
                     ["--maxiterate", "1000", "--localpair"],
                 "High accuracy - genafpair (--genafpair)":
                     ["--maxiterate", "1000", "--genafpair"],
                 "High accuracy - globalpair (--globalpair)":
                     ["--maxiterate", "1000", "--globalpair"],
                 "Fast (--retree 1)": ["--retree", "1"],
             },
             "cli": ["{value}"]},
            {"id": "threads", "widget": "thread_spinbox",
             "label": "Threads", "kind": "int",
             "cli": ["--thread", "{value}"],
             "omit_values": [-1],
             "notes": "-1 means automatic and is omitted."},
            {"id": "output_suffix", "widget": "output_suffix_edit",
             "label": "Output suffix", "kind": "text",
             "cli": [],
             "notes": "GUI-only: MAFFT writes to stdout, the plugin writes "
                      "the file beside the input."},
        ]

    def cli_preview(self):
        input_file = (self.imported_files[0] if self.imported_files
                      else self.file_path_edit.text().strip() or "<input>")
        return [self.tool_path or "mafft"] + self.get_parameters() + [input_file]

    def get_parameters(self):
        """获取MAFFT命令参数"""
        params = []
        
        # 基于模式选择的基本参数
        mode_text = self.mode_combo.currentText()
        if "Auto mode" in mode_text:
            params.append("--auto")
        elif "High accuracy - localpair" in mode_text:
            params.extend(["--maxiterate", "1000", "--localpair"])
        elif "High accuracy - genafpair" in mode_text:
            params.extend(["--maxiterate", "1000", "--genafpair"])
        elif "High accuracy - globalpair" in mode_text:
            params.extend(["--maxiterate", "1000", "--globalpair"])
        elif "Fast" in mode_text:
            params.append("--retree")
            params.append("1")
        # 高速(默认) - 无额外参数
        
        # 线程数
        thread_count = self.thread_spinbox.value()
        if thread_count != -1:
            params.extend(["--thread", str(thread_count)])
        
        # 高级参数
        if hasattr(self, 'advanced_params'):
            params.extend(self.advanced_params)
        
        return params
    
    def run_analysis(self):
        """运行MAFFT分析"""
        if self.is_running:
            return
            
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "MAFFT executable file not found!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting MAFFT alignment...", "info")
        
        # 准备输入文件
        input_files = self.prepare_input_files()
        if not input_files:
            return
            
        # 更新UI状态
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 未知进度
        
        # Output suffix: beside the input file. Temporary mode ignores it.
        suffix = None
        if self.save_to_cwd.isChecked():
            suffix = self.output_suffix_edit.text()
            if any(separator in suffix for separator in ("/", "\\")):
                QMessageBox.warning(
                    self, "Invalid suffix",
                    "The output suffix must not contain a path separator.",
                )
                return

        self._mafft_payload = None
        self.import_to_platform_btn.setEnabled(True)
        self.import_to_platform_btn.setText("Import to Current Platform")
        self.progress_bar.setFormat("Running MAFFT...")
        self.alignment_thread = MAFFTThread(
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
        self._mafft_payload = payload

    def analysis_finished(self, output_files, html_files):
        """分析完成处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.alignment_output_files = output_files

        # HTML reports were removed from this plugin; the base report combo
        # no longer exists, so do not touch it here.
        self.reports = []
        self._show_alignment_preview(output_files)
        self.tab_widget.setCurrentIndex(1)

        self.add_console_message("Alignment completed successfully!", "info")
        QMessageBox.information(self, "Success", "Alignment completed successfully!")

        if (self.import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"]
                or self.result_has_consumer()):
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)

    def analysis_error(self, error_msg):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Error", f"Alignment failed: {error_msg}")
        self.tab_widget.setCurrentIndex(2)

    def stop_analysis(self):
        """Stop MAFFT without terminating the QThread object."""
        if hasattr(self, "alignment_thread") and self.alignment_thread.isRunning():
            thread = self.alignment_thread
            # A killed process can surface a non-zero return code; suppress
            # the completion/error callbacks because the user initiated stop.
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

    def browse_files(self):
        """浏览选择文件"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA files", "", "FASTA files (*.fas *.fasta *.fa *.fna);;All files (*)"
        )
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            # 隐藏文本输入框
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)

    def add_file_tag(self, file_path):
        """添加文件标签"""
        self.imported_files.append(file_path)
        self.import_file = (
            self.imported_files[0] if len(self.imported_files) == 1 else None
        )
        self.batch_mode = len(self.imported_files) > 1

        # Create file tag widget
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.NoFrame)
        tag_widget.setStyleSheet(
            "QFrame {"
            "  background-color: palette(base);"
            "  border: 1px solid palette(mid);"
            "  border-radius: 12px;"
            "  margin: 2px;"
            "}"
        )

        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)

        # Get display name (handle duplicate names)
        display_name = self.get_display_name(file_path)

        # File name label
        name_label = QLabel(display_name)
        tag_layout.addWidget(name_label)

        # Close button
        close_btn = QPushButton("×")
        close_btn.setFixedSize(20, 20)
        close_btn.setToolTip("Remove this file")
        close_btn.setFlat(True)
        close_btn.clicked.connect(lambda: self.remove_file_tag(file_path, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append((file_path, tag_widget))
        
        # Show container if it was hidden
        self.file_tags_container.setVisible(True)

    def remove_file_tag(self, file_path, tag_widget):
        """Remove one file from the input list and its tag widget."""
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        self.file_tags = [
            entry for entry in self.file_tags
            if entry[0] != file_path or entry[1] is not tag_widget
        ]
        tag_widget.setParent(None)
        tag_widget.deleteLater()
        self.batch_mode = len(self.imported_files) > 1
        self.import_file = (
            self.imported_files[0] if len(self.imported_files) == 1 else None
        )
        if self.imported_files:
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(
                    "%d files selected" % len(self.imported_files)
                )
        else:
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)

    def clear_all_file_tags(self):
        """Remove every input tag; used when the user switches to text input."""
        for file_path, tag_widget in list(self.file_tags):
            tag_widget.setParent(None)
            tag_widget.deleteLater()
        self.file_tags = []
        self.imported_files = []
        self.import_file = None
        self.batch_mode = False
        self.file_path_edit.clear()
        self.file_tags_container.setVisible(False)

    def on_text_changed(self):
        """Handle text input changes"""
        text = self.sequence_text.toPlainText().strip()
        if text:
            # Hide file input when text is present
            self.file_path_edit.setVisible(False)
            self.file_browse_btn.setVisible(False)
            self.file_tags_container.setVisible(False)
            # Force select Temporary File mode for text input
            self.save_to_tmp.setChecked(True)
            self.save_to_cwd.setEnabled(False)
            # Clear imported files
            self.clear_all_file_tags()
        else:
            # Show file input when text is empty
            self.file_path_edit.setVisible(True)
            self.file_browse_btn.setVisible(True)
            self.save_to_cwd.setEnabled(True)
            if self.imported_files:
                self.file_tags_container.setVisible(True)

    def prepare_input_files(self):
        """Prepare input files from multiple sources"""
        try:
            input_files = []
            
            # If there are imported files, use them individually
            if self.imported_files:
                for file_path in self.imported_files:
                    if os.path.exists(file_path):
                        input_files.append(file_path)
                    else:
                        QMessageBox.warning(self, "Warning", f"File does not exist: {file_path}")
                return input_files
            elif self.import_file:
                return [self.import_file]
            
            # Otherwise, use text input to create a temporary file
            sequence_text = self.sequence_text.toPlainText().strip()
            if not sequence_text and not self.import_file:
                QMessageBox.warning(self, "Warning", "Please input sequence text!")
                return None
            
            if self.save_to_tmp.isChecked():
                # Create temporary file
                temp_file = self.create_temp_file(suffix='.fas')
                with open(temp_file, 'w') as f:
                    f.write(sequence_text)
                return [temp_file]
            
            elif self.save_to_cwd.isChecked():
                if self.import_file:
                    return [self.import_file]
                QMessageBox.warning(
                    self, "Warning",
                    "No input file is available for this output mode.",
                )
                return None

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None

    def get_display_name(self, file_path):
        """Get display name for file, handling duplicates"""
        filename = os.path.basename(file_path)
        
        # Check for duplicates
        duplicate_count = 0
        for existing_path in self.imported_files:
            if os.path.basename(existing_path) == filename and existing_path != file_path:
                duplicate_count += 1
        
        if duplicate_count > 0:
            # Use last directory name for duplicates
            parent_dir = os.path.basename(os.path.dirname(file_path))
            return f"{filename} ({parent_dir})"
        else:
            return filename

    def _read_alignment(self, path):
        records = _strip_direction_prefix(list(SeqIO.parse(str(path), "fasta")))
        if records:
            return records, "fasta"
        try:
            alignment = AlignIO.read(str(path), "clustal")
        except Exception as exc:
            raise ValueError(
                "Cannot parse MAFFT output %s as FASTA or Clustal: %s"
                % (path, exc)
            )
        return _strip_direction_prefix(list(alignment)), "clustal"

    def _build_result_spec(self):
        """Build the standard ResultSpec for one or more alignments."""
        if not getattr(self, "alignment_output_files", None):
            raise ValueError("No MAFFT alignment output is available.")

        items = []
        total_sequences = 0
        max_length = 0
        inputs = list(getattr(self, "imported_files", []) or [])
        suffix = (self.output_suffix_edit.text()
                  if self.save_to_cwd.isChecked() else "")
        for index, output_file in enumerate(self.alignment_output_files):
            records, file_format = self._read_alignment(output_file)
            if len(records) < 2:
                raise ValueError(
                    "%s contains fewer than two aligned sequences." % output_file
                )
            lengths = {len(str(record.seq)) for record in records}
            if len(lengths) > 1:
                raise ValueError(
                    "%s is not a rectangular alignment." % output_file
                )
            total_sequences += len(records)
            max_length = max(max_length, max(lengths))
            path = Path(output_file)
            loci_name = _alignment_locus_name(path, suffix)
            metadata = {"source": self.plugin_name}
            if index < len(inputs):
                metadata["input"] = inputs[index]
            items.append(DatasetItemSpec(
                loci_name=loci_name,
                item_type="alignment",
                sequences=records,
                file_path=str(path),
                format=file_format,
                is_aligned=True,
                length=max(lengths),
                sequence_count=len(records),
                data={
                    "tool": "MAFFT",
                    "parameters": self.get_parameters(),
                },
                metadata=metadata,
            ))

        dataset = DatasetSpec(
            name="MAFFT alignment" if len(items) == 1 else "MAFFT alignments",
            items=items,
            dataset_id=getattr(self, "_source_dataset_id", None),
            description=(
                "%d alignment(s), %d sequences"
                % (len(items), total_sequences)
            ),
            source="mafft",
            summary={
                "alignments": len(items),
                "sequences": total_sequences,
                "length": max_length,
            },
            metadata={
                "input": inputs[0] if len(inputs) == 1 else "",
                "input_files": inputs,
                "parameters": self.get_parameters(),
                "mode": self.mode_combo.currentText(),
                "batch": bool(getattr(self, "batch_mode", False)),
            },
        )
        return ResultSpec(source="mafft", datasets=[dataset])

    def import_to_platform(self):
        """Send results through the standard result channel when connected."""
        if not hasattr(self, 'alignment_output_files') or not self.alignment_output_files:
            QMessageBox.warning(self, "Warning", "No alignment results to import.")
            return

        if self.result_has_consumer():
            try:
                result = self._build_result_spec()
            except Exception as exc:
                QMessageBox.critical(
                    self, "Error",
                    "Failed to build alignment result: %s" % exc,
                )
                return
            if self.export_result(result):
                self.import_to_platform_btn.setEnabled(False)
                self.import_to_platform_btn.setText("Sent to workspace")
                QMessageBox.information(
                    self, "Success",
                    "Alignment result sent to the current workspace.",
                )
            return

        # Legacy fallback for callers that still connect the old signals.
        try:
            is_batch_mode = getattr(self, 'batch_mode', False)
            if is_batch_mode:
                batch_sequences = []
                for output_file in self.alignment_output_files:
                    file_sequences, _ = self._read_alignment(output_file)
                    batch_sequences.append(file_sequences)
                self.batch_import_alignment_signal.emit(batch_sequences)
            else:
                sequences = []
                for output_file in self.alignment_output_files:
                    file_sequences, _ = self._read_alignment(output_file)
                    sequences.extend(file_sequences)

                if not sequences:
                    QMessageBox.warning(self, "Warning", "No sequences found in alignment results.")
                    return
                self.import_alignment_signal.emit(sequences)

            if is_batch_mode:
                QMessageBox.information(self, "Success", f"Successfully imported {len(self.alignment_output_files)} alignment files to the platform.")
            else:
                QMessageBox.information(self, "Success", f"Successfully imported {len(sequences)} sequences to the platform.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import alignment results: {str(e)}")

class MAFFTAdvancedDialog(QDialog):
    """MAFFT高级参数对话框"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MAFFT Advanced Parameters")
        self.setModal(True)
        self.setMinimumSize(500, 400)
        
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # 为参数创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout()
        scroll_widget.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # 间隙惩罚
        gap_group = QGroupBox("Gap Penalties")
        gap_layout = QFormLayout()
        gap_group.setLayout(gap_layout)
        scroll_layout.addWidget(gap_group)
        
        self.op_spinbox = QDoubleSpinBox()
        self.op_spinbox.setRange(0.0, 10.0)
        self.op_spinbox.setSingleStep(0.1)
        self.op_spinbox.setValue(1.53)
        self.op_spinbox.setDecimals(2)
        self.op_spinbox.setToolTip("Gap opening penalty (default: 1.53)")
        gap_layout.addRow("Gap opening penalty (--op):", self.op_spinbox)
        
        self.ep_spinbox = QDoubleSpinBox()
        self.ep_spinbox.setRange(0.0, 5.0)
        self.ep_spinbox.setSingleStep(0.1)
        self.ep_spinbox.setValue(0.0)
        self.ep_spinbox.setDecimals(2)
        self.ep_spinbox.setToolTip("Offset/gap extension penalty (default: 0.0)")
        gap_layout.addRow("Gap extension penalty (--ep):", self.ep_spinbox)
        
        # 迭代设置
        iter_group = QGroupBox("Iteration Settings")
        iter_layout = QFormLayout()
        iter_group.setLayout(iter_layout)
        scroll_layout.addWidget(iter_group)
        
        self.maxiterate_spinbox = QSpinBox()
        self.maxiterate_spinbox.setRange(0, 10000)
        self.maxiterate_spinbox.setValue(0)
        self.maxiterate_spinbox.setToolTip("Maximum number of iterative refinement (default: 0)")
        iter_layout.addRow("Max iterations (--maxiterate):", self.maxiterate_spinbox)
        
        # 输出选项
        output_group = QGroupBox("Output Options")
        output_layout = QFormLayout()
        output_group.setLayout(output_layout)
        scroll_layout.addWidget(output_group)
        
        self.clustalout_checkbox = QCheckBox("Output in Clustal format")
        self.clustalout_checkbox.setToolTip("Output alignment in Clustal format (--clustalout)")
        output_layout.addRow("", self.clustalout_checkbox)
        
        self.reorder_checkbox = QCheckBox("Reorder sequences")
        self.reorder_checkbox.setToolTip("Reorder sequences in aligned order (--reorder)")
        output_layout.addRow("", self.reorder_checkbox)
        
        # 调整方向选项
        self.adjust_direction_combo = QComboBox()
        self.adjust_direction_combo.addItems([
            "No adjustment",
            "Adjust direction (--adjustdirection)",
            "Adjust direction accurately (--adjustdirectionaccurately)"
        ])
        self.adjust_direction_combo.setToolTip("Handle potentially reversed sequencing reads")
        output_layout.addRow("Direction adjustment:", self.adjust_direction_combo)
        
        scroll_layout.addStretch()
        
        # 对话框按钮
        button_layout = QHBoxLayout()
        
        reset_btn = QPushButton("Reset to Defaults")
        reset_btn.clicked.connect(self.reset_to_defaults)
        button_layout.addWidget(reset_btn)
        
        button_layout.addStretch()
        
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)
        
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self.accept)
        button_layout.addWidget(ok_btn)
        
        layout.addLayout(button_layout)
    
    def reset_to_defaults(self):
        """将所有参数重置为默认值"""
        self.op_spinbox.setValue(1.53)
        self.ep_spinbox.setValue(0.0)
        self.maxiterate_spinbox.setValue(0)
        self.clustalout_checkbox.setChecked(False)
        self.reorder_checkbox.setChecked(False)
        self.adjust_direction_combo.setCurrentIndex(0)
    
    def get_parameters(self):
        """将高级参数获取为列表"""
        params = []
        
        # 间隙惩罚
        if self.op_spinbox.value() != 1.53:  # 仅在非默认值时添加
            params.extend(["--op", str(self.op_spinbox.value())])
        
        if self.ep_spinbox.value() != 0.0:  # 仅在非默认值时添加
            params.extend(["--ep", str(self.ep_spinbox.value())])
        
        # 迭代设置
        if self.maxiterate_spinbox.value() != 0:  # 仅在非默认值时添加
            params.extend(["--maxiterate", str(self.maxiterate_spinbox.value())])
        
        # 输出选项
        if self.clustalout_checkbox.isChecked():
            params.append("--clustalout")
        
        if self.reorder_checkbox.isChecked():
            params.append("--reorder")

        # 添加调整方向参数
        adjust_index = self.adjust_direction_combo.currentIndex()
        if adjust_index == 1:  # --adjustdirection
            params.append("--adjustdirection")
        elif adjust_index == 2:  # --adjustdirectionaccurately
            params.append("--adjustdirectionaccurately")

        return params


class MAFFTPluginEntry:
    """MAFFT插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return MAFFTPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )