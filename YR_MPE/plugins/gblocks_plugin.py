# gblocks_plugin.py
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

from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFormLayout, QFrame, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton, QRadioButton,
    QScrollArea, QSpinBox, QTextEdit, QVBoxLayout, QWidget,
)
from PyQt5.QtCore import pyqtSignal
import os
import shutil
from io import StringIO
from pathlib import Path

from Bio import SeqIO

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.alignment_common import (
    build_result, prepare_import, records_temp, setup_preview, show_preview,
)


class GBlocksThread(BaseProcessThread):
    """Run GBlocks and move its -gb.fa output to the requested path."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters,
                 imported_files=None, suffix=None, workdir=None):
        super().__init__(
            tool_path, input_files, parameters, imported_files, workdir=workdir
        )
        self.output_suffix = suffix or ""

    def get_tool_name(self):
        return "GBlocks"

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
                cmd = [
                    self.tool_path, input_file, *self.parameters,
                    "-p=n", "-q=n", "-r=n", "-e=-gb",
                ]
                result = self.execute_command(
                    cmd, cwd=str(output_file.parent)
                )
                detail = (result.stderr or "").strip()
                if result.returncode != 0:
                    raise RuntimeError(
                        "GBlocks failed for %s%s"
                        % (os.path.basename(input_file),
                           ": " + detail if detail else "")
                    )
                generated = Path("%s-gb.fa" % input_file)
                if not generated.is_file():
                    raise RuntimeError(
                        "GBlocks did not write %s" % generated
                    )
                if generated.resolve() != output_file.resolve():
                    output_file.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(generated), str(output_file))
                if output_file.stat().st_size == 0:
                    raise RuntimeError("GBlocks wrote an empty file: %s" % output_file)
                output_files.append(str(output_file))

            self.progress.emit("All trimming completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])
        except Exception as exc:
            self.error.emit("Trimming exception: %s" % exc)


class GBlocksPlugin(BasePlugin):
    """GBlocks插件类"""
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
        self.run_button.setText("Run GBlocks")

    def init_plugin_info(self):
        self.plugin_name = "GBlocks Trimmer"
        self.tool_name = "GBlocks"
        self.citation = []
        self.citation_keys = ["talavera2007improvement"]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"],
                            "NBRF/PIR": ["pir"]}
        self.output_types = {"Trimmed alignment FASTA": ["fas"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

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

    def control_bindings(self):
        """GUI controls mapped to the GBlocks command line."""
        return [
            {"id": "input", "widget": "file_path_edit",
             "label": "Input alignment", "kind": "path", "cli": [],
             "notes": "Positional input; one GBlocks run per file."},
            {"id": "seq_type", "widget": "seq_type_combo",
             "label": "Sequence type", "kind": "choice",
             "choices": {"Auto": "", "Protein": "-t=p",
                         "DNA": "-t=d", "Codons": "-t=c"},
             "cli": ["{value}"]},
            {"id": "b1", "widget": "b1_spinbox",
             "label": "Minimum block length", "kind": "int",
             "cli": ["-b1={value}"], "omit_values": [0]},
            {"id": "b2", "widget": "b2_spinbox",
             "label": "Minimum conserved positions", "kind": "int",
             "cli": ["-b2={value}"], "omit_values": [0]},
            {"id": "b3", "widget": "b3_spinbox",
             "label": "Maximum contiguous non-conserved", "kind": "int",
             "cli": ["-b3={value}"], "omit_values": [8]},
            {"id": "b4", "widget": "b4_spinbox",
             "label": "Minimum flank positions", "kind": "int",
             "cli": ["-b4={value}"], "omit_values": [10]},
            {"id": "b5", "widget": "gap_combo",
             "label": "Gap positions", "kind": "choice",
             "choices": {"None": "", "With Half": "-b5=h",
                         "All": "-b5=a"},
             "cli": ["{value}"]},
            {"id": "b6", "widget": "sim_matrix_checkbox",
             "label": "Similarity matrix", "kind": "bool", "cli": [],
             "notes": "Emitted as -b6=y when checked, -b6=n otherwise."},
        ]

    def setup_input_tab(self):
        """设置输入和参数标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # 输入组
        input_group = QGroupBox("Input")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # 文件输入
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select alignment files...")
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
        self.sequence_text.setPlaceholderText("Or paste alignment in FASTA or NBRF/PIR format...")
        self.sequence_text.setMaximumHeight(200)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        input_layout.addRow("Sequence text:", self.sequence_text)
        
        # 处理导入的数据
        if self.import_file:
            self.file_path_edit.setText(self.import_file)
            input_group.setVisible(False)
        elif hasattr(self, 'imported_files') and self.imported_files:
            # 显示导入的文件
            for file_path in self.imported_files:
                self.add_file_tag(file_path)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            # 隐藏文本输入框
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
        
        # 参数组
        params_group = QGroupBox("GBlocks Parameters")
        params_form_layout = QFormLayout()
        params_group.setLayout(params_form_layout)
        layout.addWidget(params_group)
        
        # 序列类型
        self.seq_type_combo = QComboBox()
        self.seq_type_combo.addItems(["Auto", "Protein", "DNA", "Codons"])
        params_form_layout.addRow("Sequence Type:", self.seq_type_combo)
        
        # Output settings
        output_settings_group = QWidget()
        output_settings_group.setLayout(QVBoxLayout())
        output_settings_group.layout().setContentsMargins(0, 0, 0, 0)
        output_settings_group.setContentsMargins(0, 0, 0, 0)
        output_radio_group = QWidget()
        output_radio_group.setLayout(QHBoxLayout())
        self.save_to_cwd = QRadioButton("Current Directory")
        self.save_to_cwd.setChecked(True)
        self.save_to_tmp = QRadioButton("Temporary File")
        output_radio_group.layout().addWidget(self.save_to_cwd)
        output_radio_group.layout().addWidget(self.save_to_tmp)
        output_settings_group.layout().addWidget(output_radio_group)

        self.save_to_cwd.toggled.connect(self.on_output_radio_changed)
        self.save_to_tmp.toggled.connect(self.on_output_radio_changed)

        self.output_suffix_edit = QLineEdit('_gblocks')
        self.output_suffix_edit.setPlaceholderText("Output suffix")
        output_settings_group.layout().addWidget(self.output_suffix_edit)
        params_form_layout.addRow("Save to:", output_settings_group)
        
        # 保守位置最小序列数
        self.b1_spinbox = QSpinBox()
        self.b1_spinbox.setMinimum(0)
        self.b1_spinbox.setMaximum(1000)
        self.b1_spinbox.setValue(0)  # 0表示使用默认值
        params_form_layout.addRow("Min. Sequences for Conserved Position (-b1):", self.b1_spinbox)
        
        # 侧翼位置最小序列数
        self.b2_spinbox = QSpinBox()
        self.b2_spinbox.setMinimum(0)
        self.b2_spinbox.setMaximum(1000)
        self.b2_spinbox.setValue(0)  # 0表示使用默认值
        params_form_layout.addRow("Min. Sequences for Flank Position (-b2):", self.b2_spinbox)
        
        # 最大连续非保守位置数
        self.b3_spinbox = QSpinBox()
        self.b3_spinbox.setMinimum(0)
        self.b3_spinbox.setMaximum(1000)
        self.b3_spinbox.setValue(8)  # 默认值
        params_form_layout.addRow("Max. Contiguous Nonconserved Positions (-b3):", self.b3_spinbox)
        
        # 区块最小长度
        self.b4_spinbox = QSpinBox()
        self.b4_spinbox.setMinimum(2)
        self.b4_spinbox.setMaximum(1000)
        self.b4_spinbox.setValue(10)  # 默认值
        params_form_layout.addRow("Minimum Length of a Block (-b4):", self.b4_spinbox)
        
        # 允许的缺口位置
        self.gap_combo = QComboBox()
        self.gap_combo.addItems(["None", "With Half", "All"])
        params_form_layout.addRow("Allowed Gap Positions (-b5):", self.gap_combo)
        
        # 使用相似性矩阵
        self.sim_matrix_checkbox = QCheckBox("Use Similarity Matrix (Protein only) (-b6)")
        self.sim_matrix_checkbox.setChecked(True)
        params_form_layout.addRow(self.sim_matrix_checkbox)
        
        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets
    
    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.fas *.fasta *.fa *.fna *.pir);;All files (*)"
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select alignment files", "", file_filter)
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
        
        # Create file tag widget
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.Box)
        tag_widget.setStyleSheet("""
            QFrame {
                background-color: #e9ecef;
                border-radius: 15px;
                margin: 2px;
            }
        """)
        
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)
        
        # Get display name (handle duplicate names)
        display_name = self.get_display_name(file_path)
        
        # File name label
        name_label = QLabel(display_name)
        name_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(name_label)
        
        # Close button
        close_btn = QPushButton("×")
        close_btn.setFixedSize(20, 20)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: transparent;
                border: none;
                color: #6c757d;
                font-weight: bold;
                font-size: 14px;
            }
            QPushButton:hover {
                color: #dc3545;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_file_tag(tag_widget, file_path))
        tag_layout.addWidget(close_btn)
        
        # Add to container
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append(tag_widget)
        
        # Show container
        self.file_tags_container.setVisible(True)
    
    def remove_file_tag(self, tag_widget, file_path):
        """移除文件标签"""
        # Remove from lists
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        if tag_widget in self.file_tags:
            self.file_tags.remove(tag_widget)
        
        # Remove widget
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Update file path display
        if len(self.imported_files) == 0:
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
        elif len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        else:
            self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def get_display_name(self, file_path):
        """获取文件显示名称（处理重名情况）"""
        base_name = os.path.basename(file_path)
        count = sum(1 for f in self.imported_files if os.path.basename(f) == base_name)
        if count > 1:
            name, ext = os.path.splitext(base_name)
            dir_name = os.path.basename(os.path.dirname(file_path))
            return f"{name} ({dir_name}){ext}"
        return base_name
    
    def on_text_changed(self):
        """文本输入变化处理"""
        if self.sequence_text.toPlainText().strip():
            # 如果有文本输入，清空文件选择
            self.imported_files.clear()
            for tag in self.file_tags:
                self.file_tags_layout.removeWidget(tag)
                tag.deleteLater()
            self.file_tags.clear()
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
            # Force select Temporary File mode for text input
            self.save_to_tmp.setChecked(True)
            self.save_to_cwd.setEnabled(False)
        elif not self.imported_files:
            # 如果既没有文本也没有文件，恢复初始状态
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
            self.save_to_cwd.setEnabled(True)
    
    def setup_gblocks_params(self):
        """设置GBlocks参数"""
        # GBlocks参数已在setup_input_tab中设置，此处留空以保持接口一致性
        pass
    
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Send to workspace")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)
    
    def setup_output_tab(self):
        setup_preview(self, "The GBlocks-trimmed FASTA alignment will appear here.")

    def get_parameters(self):
        params = []
        seq_type_index = self.seq_type_combo.currentIndex()
        if seq_type_index == 1:
            params.append("-t=p")
        elif seq_type_index == 2:
            params.append("-t=d")
        elif seq_type_index == 3:
            params.append("-t=c")
        if self.b1_spinbox.value() > 0:
            params.append("-b1=%d" % self.b1_spinbox.value())
        if self.b2_spinbox.value() > 0:
            params.append("-b2=%d" % self.b2_spinbox.value())
        if self.b3_spinbox.value() != 8:
            params.append("-b3=%d" % self.b3_spinbox.value())
        if self.b4_spinbox.value() != 10:
            params.append("-b4=%d" % self.b4_spinbox.value())
        b5 = {"None": "n", "With Half": "h", "All": "a"}[
            self.gap_combo.currentText()
        ]
        if b5 != "n":
            params.append("-b5=%s" % b5)
        params.append("-b6=%s" % ("y" if self.sim_matrix_checkbox.isChecked() else "n"))
        return params

    def on_output_radio_changed(self):
        """处理输出选项的切换"""
        if self.save_to_cwd.isChecked():
            self.output_suffix_edit.setEnabled(True)
            self.output_suffix_edit.setVisible(True)
        elif self.save_to_tmp.isChecked():
            self.output_suffix_edit.setEnabled(False)
            self.output_suffix_edit.setVisible(False)
    
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
                QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
                return None
            
            try:
                records = list(SeqIO.parse(StringIO(sequence_text), "fasta"))
            except Exception as exc:
                QMessageBox.warning(self, "Invalid FASTA", str(exc))
                return None
            if len(records) < 2:
                QMessageBox.warning(self, "Need more sequences", "Provide at least two aligned sequences.")
                return None
            return [records_temp(self, records, "pasted")]

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
    
    def run_analysis(self):
        if self.is_running:
            return
        if not self.tool_path or not Path(self.tool_path).is_file():
            QMessageBox.critical(self, "Tool missing", "GBlocks executable not found.")
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
        self.progress_bar.setFormat("Running GBlocks...")
        self.import_to_platform_btn.setEnabled(True)
        self.import_to_platform_btn.setText("Send to workspace")
        self.tab_widget.setCurrentIndex(2)

        self.analysis_thread = GBlocksThread(
            self.tool_path, input_files, self.get_parameters(),
            self.imported_files, suffix=suffix, workdir=self.workdir,
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()

    def _on_result(self, payload):
        self._payload = payload

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.alignment_output_files = output_files
        show_preview(self, output_files)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("GBlocks trimming completed.", "info")
        if self.result_has_consumer() or self.import_from in [
            "YR_MPEA", "seq_viewer", "DATASET_MANAGER"
        ]:
            self.import_to_platform_btn.setVisible(True)
        self.import_to_platform_btn.setEnabled(True)

    def analysis_error(self, error_message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("GBlocks failed: %s" % error_message, "error")
        QMessageBox.critical(self, "GBlocks failed", error_message)

    def stop_analysis(self):
        if hasattr(self, "analysis_thread") and self.analysis_thread.isRunning():
            thread = self.analysis_thread
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
        self.add_console_message("GBlocks trimming stopped by user.", "warning")

    def _build_result_spec(self):
        if not getattr(self, "alignment_output_files", None):
            raise ValueError("No GBlocks output is available.")
        suffix = (self.output_suffix_edit.text()
                  if self.save_to_cwd.isChecked() else "")
        return build_result(
            tool="GBlocks",
            source="gblocks",
            output_files=self.alignment_output_files,
            input_files=self.imported_files,
            parameters=self.get_parameters(),
            mode=self.seq_type_combo.currentText(),
            batch=self.batch_mode,
            suffix=suffix,
        )

    def import_to_platform(self):
        if not getattr(self, "alignment_output_files", None):
            QMessageBox.warning(self, "No result", "Run GBlocks first.")
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
                    self, "Success", "GBlocks result sent to the current workspace."
                )
            return
        try:
            if self.batch_mode:
                self.batch_import_alignment_signal.emit([
                    list(SeqIO.parse(path, "fasta"))
                    for path in self.alignment_output_files
                ])
            else:
                sequences = []
                for path in self.alignment_output_files:
                    sequences.extend(list(SeqIO.parse(path, "fasta")))
                if sequences:
                    self.import_alignment_signal.emit(sequences)
        except Exception as exc:
            QMessageBox.critical(self, "Import failed", str(exc))

    def reset_to_defaults(self):
        self.seq_type_combo.setCurrentIndex(0)
        self.b1_spinbox.setValue(0)
        self.b2_spinbox.setValue(0)
        self.b3_spinbox.setValue(8)
        self.b4_spinbox.setValue(10)
        self.gap_combo.setCurrentIndex(0)
        self.sim_matrix_checkbox.setChecked(True)


class GBlocksPluginEntry:
    """GBlocks插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return GBlocksPlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )