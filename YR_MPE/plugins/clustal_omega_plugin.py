# cluster_omega_plugin.py
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
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, QRadioButton)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon
import tempfile
import os
from typing import List, Optional
from Bio import SeqIO
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class ClustalOmegaThread(BaseProcessThread):
    """Clustal Omega比对线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None, suffix=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
        if suffix:
            self.output_suffix = suffix
    
    def get_tool_name(self):
        """返回工具名称"""
        return "Clustal Omega"
        
    def execute_commands(self):
        """执行Clustal Omega比对命令"""
        try:
            output_files = []
            html_files = []
            
            # 分别处理每个输入文件
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                    
                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", "info")
                
                # 创建输出文件
                # 如果是从UI传递过来的suffix参数，则使用它来构造输出文件名
                if self.output_suffix:
                    # 构造带后缀的输出文件名
                    output_file = f"{input_file}{self.output_suffix}.fas"
                else:
                    output_file = self.create_temp_file(suffix='.fas')
                
                # 构建命令
                cmd = [
                    self.tool_path,
                    "-i", input_file,
                    "-o", output_file,
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"Clustal Omega execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 从FASTA输出生成HTML报告
                self.console_output.emit(f"Generating HTML report for file {i+1}...", "info")
                
                output_files.append(output_file)
            
            self.progress.emit("All alignments completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Alignment exception: {str(e)}")


class ClustalOmegaPlugin(BasePlugin):
    """Clustal Omega插件类"""
    import_alignment_signal = pyqtSignal(list)
    batch_import_alignment_signal = pyqtSignal(list)
    
    def __init__(self, import_from=None, import_data=None):
        """初始化Clustal Omega插件"""
        super().__init__(import_from, import_data)
        self.batch_mode = False
        if import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"] and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "Clustal Omega Aligner"
        self.tool_name = "Clustal Omega"
        self.citation = ["""Sievers F, Wilm A, Dineen D, et al. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. <i>Mol Syst Biol.</i> 2011;7:539. Published 2011 Oct 11. DOI: <a href="https://doi.org/10.1038/msb.2011.75">10.1038/msb.2011.75</a>"""]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"FASTA": ".fas"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
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
        basic_group = QGroupBox("Basic Parameters")
        basic_layout = QFormLayout()
        basic_group.setLayout(basic_layout)
        layout.addWidget(basic_group)
        
        # 自动选项
        self.auto_checkbox = QCheckBox("Auto options (--auto)")
        self.auto_checkbox.setChecked(True)
        basic_layout.addRow("", self.auto_checkbox)
        
        # 序列类型
        self.seqtype_combo = QComboBox()
        self.seqtype_combo.addItems(["Auto", "Protein", "RNA", "DNA"])
        basic_layout.addRow("Sequence type:", self.seqtype_combo)
        
        # 线程数
        self.thread_spinbox = QSpinBox()
        self.thread_spinbox.setRange(1, 64)
        self.thread_spinbox.setValue(4)
        basic_layout.addRow("Threads:", self.thread_spinbox)
        
        # 输出顺序
        self.order_combo = QComboBox()
        self.order_combo.addItems(["input-order", "tree-order"]) 
        basic_layout.addRow("Output order:", self.order_combo)
        
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

        self.output_suffix_edit = QLineEdit('_clustal')
        self.output_suffix_edit.setPlaceholderText("Output suffix")
        output_settings_group.layout().addWidget(self.output_suffix_edit)
        input_layout.addRow("Save to:", output_settings_group)
        
        layout.addStretch()
    
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)
    
    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 自动选项
        if self.auto_checkbox.isChecked():
            params.append("--auto")
        
        # 序列类型
        seqtype = self.seqtype_combo.currentText()
        if seqtype != "Auto":
            params.extend(["--seqtype", seqtype])
            
        # 线程数
        threads = self.thread_spinbox.value()
        if threads > 1:
            params.extend(["--threads", str(threads)])
            
        # 输出顺序
        order = self.order_combo.currentText()
        if order == "tree-order":
            params.append("--outorder=tree-order")
            
        return params

    def on_output_radio_changed(self):
        """处理输出选项的切换"""
        if self.save_to_cwd.isChecked():
            self.output_suffix_edit.setEnabled(True)
            self.output_suffix_edit.setVisible(True)
        elif self.save_to_tmp.isChecked():
            self.output_suffix_edit.setEnabled(False)
            self.output_suffix_edit.setVisible(False)

    def run_analysis(self):
        """运行Clustal Omega分析"""
        if self.is_running:
            return
            
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "Clustal Omega executable file not found!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting Clustal Omega alignment...", "info")
        
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
        
        # 如果需要添加输出文件后缀，则设置suffix属性
        suffix = None
        if self.save_to_cwd.isChecked():
            suffix = self.output_suffix_edit.text()
        
        # 在单独的线程中运行比对
        self.alignment_thread = ClustalOmegaThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files, suffix=suffix
        )
        self.alignment_thread.progress.connect(self.progress_bar.setFormat)
        self.alignment_thread.finished.connect(self.analysis_finished)
        self.alignment_thread.error.connect(self.analysis_error)
        self.alignment_thread.console_output.connect(self.add_console_message)
        self.alignment_thread.start()

    def analysis_finished(self, output_files, html_files):
        """分析完成处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # Save output files and reports
        self.reports = html_files
        self.update_report_combo()
        
        # Show results
        if html_files:
            self.show_current_report()
            self.tab_widget.setCurrentIndex(1)  # Switch to output preview tab
        
        self.add_console_message("Alignment completed successfully!", "info")
        QMessageBox.information(self, "Success", "Alignment completed successfully!")
        
        # 显示导入按钮（支持多种来源）
        if self.import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"]:
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        
        # 保存输出文件路径供导入使用
        self.alignment_output_files = output_files

    def analysis_error(self, error_msg):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Error", f"Alignment failed: {error_msg}")
        self.tab_widget.setCurrentIndex(2)  # Switch to console tab    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'alignment_thread') and self.alignment_thread.isRunning():
            self.alignment_thread.terminate()
            self.alignment_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "Alignment has been aborted.")

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
                background-color: #dc3545;
                color: white;
                border: none;
                border-radius: 10px;
                font-weight: bold;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_file_tag(file_path, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append((file_path, tag_widget))
        
        # Show container if it was hidden
        self.file_tags_container.setVisible(True)

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
                # directly use input paths
                temp_file = self.import_file
                return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None

    def clear_all_file_tags(self):
        """Clear all file tags"""
        for file_path, tag_widget in self.file_tags:
            self.file_tags_layout.removeWidget(tag_widget)
            tag_widget.deleteLater()
        
        self.imported_files.clear()
        self.file_tags.clear()
        self.file_tags_container.setVisible(False)

    def import_to_platform(self):
        """将比对结果导入到当前平台"""
        if not hasattr(self, 'alignment_output_files') or not self.alignment_output_files:
            QMessageBox.warning(self, "Warning", "No alignment results to import.")
            return
            
        try:
            # 判断是否为批量模式
            is_batch_mode = getattr(self, 'batch_mode', False)
            
            if is_batch_mode:
                # 批量模式：保持分组结构 [[SeqRecord, ...], [SeqRecord, ...]]
                batch_sequences = []
                for output_file in self.alignment_output_files:
                    file_sequences = list(SeqIO.parse(output_file, "fasta"))
                    batch_sequences.append(file_sequences)
                self.batch_import_alignment_signal.emit(batch_sequences)
                
            else:
                # 单文件模式：扁平化列表 [SeqRecord, SeqRecord, ...]
                sequences = []
                for output_file in self.alignment_output_files:
                    # 读取FASTA文件
                    for record in SeqIO.parse(output_file, "fasta"):
                        sequences.append(record)
                
                if not sequences:
                    QMessageBox.warning(self, "Warning", "No sequences found in alignment results.")
                    return
                    
                # 发送信号将数据导入到平台
                self.import_alignment_signal.emit(sequences)
            
            # 显示成功消息
            if is_batch_mode:
                QMessageBox.information(self, "Success", f"Successfully imported {len(self.alignment_output_files)} alignment files to the platform.")
            else:
                QMessageBox.information(self, "Success", f"Successfully imported {len(sequences)} sequences to the platform.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import alignment results: {str(e)}")

    def remove_file_tag(self, file_path, tag_widget):
        """Remove a file tag"""
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        
        # Remove from tags list
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        
        # Remove widget
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Update file path display
        if not self.imported_files:
            self.file_path_edit.setText("")
            self.file_tags_container.setVisible(False)
            # Show text input when no files
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
        elif len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        else:
            self.file_path_edit.setText(f"{len(self.imported_files)} files selected")

    def handle_import_data(self, import_data):
        """处理从YR-MPEA或Dataset Manager导入的数据"""
        if isinstance(import_data, list):
            # 判断是否为批量模式：[[SeqRecord, ...], [SeqRecord, ...]]
            if len(import_data) > 0 and isinstance(import_data[0], list):
                # 批量模式
                self.batch_mode = True
                self.temp_files = []
                self.imported_files = []
                
                for i, seq_list in enumerate(import_data):
                    temp_file = self.create_temp_file(suffix=f'_batch_{i}.fas')
                    with open(temp_file, 'w') as f:
                        for seq in seq_list:
                            f.write(f">{seq.id}\n{seq.seq}\n")
                    self.temp_files.append(temp_file)
                    self.imported_files.append(temp_file)
                    
                # 更新UI显示
                if hasattr(self, 'file_path_edit') and self.file_path_edit:
                    self.file_path_edit.setText(f"{len(import_data)} files selected")
                    self.file_path_edit.setEnabled(False)
                    
            else:
                # 单文件模式
                self.batch_mode = False
                # 创建临时文件来存储导入的序列数据
                temp_file = self.create_temp_file(suffix='.fas')
                with open(temp_file, 'w') as f:
                    for seq in import_data:
                        f.write(f">{seq.id}\n{seq.seq}\n")
                self.temp_files.append(temp_file)
                self.import_file = temp_file
                self.imported_files = [temp_file]
                
                # 更新UI显示导入的文件
                if hasattr(self, 'file_path_edit') and self.file_path_edit:
                    self.file_path_edit.setText(temp_file)
        else:
            self.import_file = None
            self.imported_files = []

class ClustalOmegaPluginEntry:
    """Clustal Omega插件入口点"""
    
    def run(self, import_from=None, import_data=None):
        return ClustalOmegaPlugin(import_from=import_from, import_data=import_data)