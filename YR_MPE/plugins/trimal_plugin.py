# trimal_plugin.py
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
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, QDoubleSpinBox, QRadioButton, QSizePolicy)
from PyQt5.QtCore import Qt, pyqtSignal, QUrl
from PyQt5.QtGui import QIcon
import tempfile
import os
from typing import List, Optional

from Bio import SeqIO

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class TrimAlThread(BaseProcessThread):
    """TrimAl修剪线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None, suffix=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
        if suffix:
            self.output_suffix = suffix
    
    def get_tool_name(self):
        """返回工具名称"""
        return "TrimAl"
        
    def execute_commands(self):
        """执行TrimAl修剪命令"""
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
                
                # 创建HTML输出文件
                html_file = self.create_temp_file(suffix='.html')
                
                # 构建命令
                cmd = [
                    self.tool_path,
                    "-in", input_file,
                    "-out", output_file,
                    "-htmlout", html_file,
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"TrimAl execution failed for file {i+1}: {result.stderr}")
                    return
                
                output_files.append(output_file)
                html_files.append(html_file)
            
            self.progress.emit("All trimming completed")
            self.finished.emit(output_files, html_files)
            
        except Exception as e:
            self.error.emit(f"Trimming exception: {str(e)}")


class TrimAlPlugin(BasePlugin):
    """TrimAl插件类"""
    import_alignment_signal = pyqtSignal(list)
    batch_import_alignment_signal = pyqtSignal(list)
    
    def __init__(self, import_from=None, import_data=None):
        """初始化TrimAl插件"""
        super().__init__(import_from, import_data)
        # 初始化插件信息
        self.init_plugin_info()
        self.batch_mode = False
        
        # 特别处理YR-MPEA、seq_viewer或DATASET_MANAGER导入的数据
        if import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"] and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "TrimAl Trimmer"
        self.tool_name = "TrimAl"
        self.citation = ["""Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T. 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics. 25: 1972-1973. doi: <a href='https://doi.org/10.1093/bioinformatics/btp348'>10.1093/bioinformatics/btp348</a>."""]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], 
                           "CLUSTAL": ["aln"], 
                           "PHYLIP": ["phy"],
                           "NEXUS": ["nex", "nexus"]}
        self.output_types = {"FASTA": ".fas"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    
    def get_citation(self):
        """获取引用信息"""
        return self.citation

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
        self.sequence_text.setPlaceholderText("Or paste alignment in FASTA format...")
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
        
        # 基本参数组
        basic_params_group = QGroupBox("Basic Parameters")
        basic_params_layout = QFormLayout()
        basic_params_group.setLayout(basic_params_layout)
        layout.addWidget(basic_params_group)
        
        # 修剪模式
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Gappyout (-gappyout)",
            "Strict (-strict)",
            "Strict Plus (-strictplus)",
            "Automated 1 (-automated1)"
        ])
        self.mode_combo.setCurrentText("Automated 1 (-automated1)")
        basic_params_layout.addRow("Trimming mode:", self.mode_combo)
        
        # 输出格式
        self.output_format_combo = QComboBox()
        self.output_format_combo.addItems([
            "FASTA (-fasta)",
            "CLUSTAL (-clustal)",
            "NEXUS (-nexus)",
            "PHYLIP (-phylip)",
            "PHYLIP3.2 (-phylip3.2)",
            "MEGA (-mega)",
            "NBRF (-nbrf)"
        ])
        basic_params_layout.addRow("Output format:", self.output_format_combo)
        
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
        self.save_to_cwd = QRadioButton("Current Directory")
        self.save_to_cwd.setChecked(True)
        self.save_to_tmp = QRadioButton("Temporary File")
        output_radio_group.layout().addWidget(self.save_to_cwd)
        output_radio_group.layout().addWidget(self.save_to_tmp)
        output_settings_group.layout().addWidget(output_radio_group)

        self.save_to_cwd.toggled.connect(self.on_output_radio_changed)
        self.save_to_tmp.toggled.connect(self.on_output_radio_changed)

        self.output_suffix_edit = QLineEdit('_trimal')
        self.output_suffix_edit.setPlaceholderText("Output suffix")
        output_settings_group.layout().addWidget(self.output_suffix_edit)
        input_layout.addRow("Save to:", output_settings_group)

        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets
    
    def browse_files(self):
        """浏览选择文件"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select alignment files", "", "Alignment files (*.fas *.fna *.fa *.fasta *.aln *.phy *.nex *.nexus);;All files (*)"
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
    
    def on_text_changed(self):
        """Handle text input changes"""
        text = self.sequence_text.toPlainText().strip()
        if text:
            # Hide file input when text is present
            self.file_path_edit.setVisible(False)
            self.file_browse_btn.setVisible(False)
            self.file_tags_container.setVisible(False)
            # Clear imported files
            self.clear_all_file_tags()
        else:
            # Show file input when text is empty
            self.file_path_edit.setVisible(True)
            self.file_browse_btn.setVisible(True)
            if self.imported_files:
                self.file_tags_container.setVisible(True)
    
    def clear_all_file_tags(self):
        """Clear all file tags"""
        for file_path, tag_widget in self.file_tags:
            self.file_tags_layout.removeWidget(tag_widget)
            tag_widget.deleteLater()
        
        self.imported_files.clear()
        self.file_tags.clear()
        self.file_tags_container.setVisible(False)

    def setup_output_tab(self):
        """设置输出预览标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 说明文字
        info_label = QLabel("Trimmed alignment results will be displayed here after analysis completes.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # 报告选择下拉框（用于批量处理）
        report_layout = QHBoxLayout()
        report_label = QLabel("Select Report:")
        self.report_combo = QComboBox()
        self.report_combo.setEnabled(False)
        self.report_combo.currentIndexChanged.connect(self.on_report_changed)
        report_layout.addWidget(report_label)
        report_layout.addWidget(self.report_combo)
        report_layout.addStretch()
        layout.addLayout(report_layout)
        
        # 结果展示区域 - 使用WebEngineView显示HTML
        from PyQt5.QtWebEngineWidgets import QWebEngineView
        self.result_view = QWebEngineView()
        self.result_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(self.result_view)
        
        # 输出文件信息
        self.output_info = QLabel("Output file information will be displayed here")
        layout.addWidget(self.output_info)
    
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)

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

    def show_advanced_dialog(self):
        """显示高级参数对话框"""
        dialog = TrimAlAdvancedDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            # 存储高级参数
            self.advanced_params = dialog.get_parameters()

    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 修剪模式
        mode_text = self.mode_combo.currentText()
        if "-strict" in mode_text:
            params.append("-strict")
        elif "-strictplus" in mode_text:
            params.append("-strictplus")
        elif "-gappyout" in mode_text:
            params.append("-gappyout")
        elif "-automated1" in mode_text:
            params.append("-automated1")
            
        # 输出格式
        format_text = self.output_format_combo.currentText()
        if "-clustal" in format_text:
            params.append("-clustal")
        elif "-nexus" in format_text:
            params.append("-nexus")
        elif "-phylip" in format_text:
            params.append("-phylip")
        elif "-phylip3.2" in format_text:
            params.append("-phylip3.2")
        elif "-mega" in format_text:
            params.append("-mega")
        elif "-nbrf" in format_text:
            params.append("-nbrf")
        # 默认为fasta格式，无需额外参数
            
        # 添加高级参数
        if hasattr(self, 'advanced_params'):
            params.extend(self.advanced_params)
            
        return params

    def on_output_radio_changed(self):
        """处理输出选项的切换"""
        if self.save_to_cwd.isChecked():
            self.output_suffix_edit.setEnabled(True)
            self.output_suffix_edit.setVisible(True)
        elif self.save_to_tmp.isChecked():
            self.output_suffix_edit.setEnabled(False)
            self.output_suffix_edit.setVisible(False)

    def display_results(self, output_files):
        """在输出标签页中显示结果"""
        if not output_files:
            self.result_view.setHtml("<p>No output files generated.</p>")
            return
            
        # 显示第一个HTML输出文件的内容
        try:
            html_file = getattr(self, 'current_html_files', [None])[0]
            if html_file and os.path.exists(html_file):
                # 加载HTML文件
                from PyQt5.QtCore import QUrl
                self.result_view.load(QUrl.fromLocalFile(html_file))
            else:
                # 如果没有HTML文件，则显示序列文件内容
                with open(output_files[0], 'r') as f:
                    content = f.read()
                    self.result_view.setPlainText(content)
        except Exception as e:
            self.result_view.setHtml(f"<p>Error reading output file: {str(e)}</p>")
        
        # 更新输出信息
        self.output_info.setText(f"Results from: {os.path.basename(output_files[0])}")
    
    def run_analysis(self):
        """运行TrimAl分析"""
        if self.is_running:
            return
            
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "TrimAl executable file not found!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting TrimAl trimming...", "info")
        
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
        
        # 在单独的线程中运行TrimAl
        self.analysis_thread = TrimAlThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files, suffix=suffix
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()
    
    def analysis_finished(self, output_files, html_files):
        """分析完成处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 保存HTML报告文件
        self.reports = html_files
        self.update_report_combo()
        
        # 显示结果
        if html_files:
            self.show_current_report()
            self.tab_widget.setCurrentIndex(1)  # 切换到输出预览标签页
        
        # 添加控制台消息
        self.add_console_message("Trimming completed successfully!", "info")
        QMessageBox.information(self, "Success", "Trimming completed successfully!")
        
        # 显示导入按钮（支持多种来源）
        if self.import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"]:
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        
        # 保存输出文件路径供导入使用
        self.alignment_output_files = output_files
    
    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 添加控制台消息
        self.add_console_message(f"TrimAl failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"TrimAl failed: {error_message}")
    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "TrimAl has been aborted.")
        
    def import_to_platform(self):
        """将修剪结果导入到当前平台"""
        if not hasattr(self, 'alignment_output_files') or not self.alignment_output_files:
            QMessageBox.warning(self, "Warning", "No trimming results to import.")
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
                    QMessageBox.warning(self, "Warning", "No sequences found in trimming results.")
                    return
                    
                # 发送信号将数据导入到平台
                self.import_alignment_signal.emit(sequences)
            
            # 显示成功消息
            if is_batch_mode:
                QMessageBox.information(self, "Success", f"Successfully imported {len(self.alignment_output_files)} trimmed files to the platform.")
            else:
                QMessageBox.information(self, "Success", f"Successfully imported {len(sequences)} sequences to the platform.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import trimming results: {str(e)}")

    def show_current_report(self):
        """
        显示当前选中的报告（使用result_view而不是web_view）
        """
        if not self.reports or self.current_report_index >= len(self.reports):
            self.result_view.setHtml("<h2>No report available</h2>")
            return
        
        report_file = self.reports[self.current_report_index]
        if report_file and os.path.exists(report_file):
            # 转换Windows路径为文件URL格式
            file_url = QUrl.fromLocalFile(report_file)
            self.result_view.load(file_url)
        else:
            self.result_view.setHtml("<h2>Report file not found</h2>")


class TrimAlAdvancedDialog(QDialog):
    """TrimAl高级参数对话框"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.init_ui()
        self.load_previous_params()
        
    def init_ui(self):
        """初始化UI"""
        self.setWindowTitle("TrimAl Advanced Parameters")
        self.setMinimumWidth(500)
        self.setModal(True)
        
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_area.setWidgetResizable(True)
        scroll_layout = QVBoxLayout()
        scroll_widget.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)
        
        # 阈值设置
        threshold_group = QGroupBox("Threshold Settings")
        threshold_layout = QFormLayout()
        threshold_group.setLayout(threshold_layout)
        scroll_layout.addWidget(threshold_group)
        
        # Gap阈值
        self.gap_threshold_spinbox = QDoubleSpinBox()
        self.gap_threshold_spinbox.setRange(0.0, 1.0)
        self.gap_threshold_spinbox.setSingleStep(0.05)
        self.gap_threshold_spinbox.setValue(0.0)  # 默认不设置
        self.gap_threshold_spinbox.setToolTip("1 - (fraction of sequences with a gap allowed)")
        threshold_layout.addRow("Gap threshold (-gt):", self.gap_threshold_spinbox)
        
        # 相似性阈值
        self.sim_threshold_spinbox = QDoubleSpinBox()
        self.sim_threshold_spinbox.setRange(0.0, 1.0)
        self.sim_threshold_spinbox.setSingleStep(0.05)
        self.sim_threshold_spinbox.setValue(0.0)  # 默认不设置
        self.sim_threshold_spinbox.setToolTip("Minimum average similarity allowed")
        threshold_layout.addRow("Similarity threshold (-st):", self.sim_threshold_spinbox)
        
        # 一致性阈值
        self.con_threshold_spinbox = QDoubleSpinBox()
        self.con_threshold_spinbox.setRange(0.0, 1.0)
        self.con_threshold_spinbox.setSingleStep(0.05)
        self.con_threshold_spinbox.setValue(0.0)  # 默认不设置
        self.con_threshold_spinbox.setToolTip("Minimum consistency value allowed")
        threshold_layout.addRow("Consistency threshold (-ct):", self.con_threshold_spinbox)
        
        # 保守比例
        self.cons_percentage_spinbox = QSpinBox()
        self.cons_percentage_spinbox.setRange(0, 100)
        self.cons_percentage_spinbox.setValue(0)  # 默认不设置
        self.cons_percentage_spinbox.setToolTip("Minimum percentage of the positions in the original alignment to conserve")
        threshold_layout.addRow("Conservation percentage (-cons):", self.cons_percentage_spinbox)
        
        # 窗口大小设置
        window_group = QGroupBox("Window Size Settings")
        window_layout = QFormLayout()
        window_group.setLayout(window_layout)
        scroll_layout.addWidget(window_group)
        
        # 通用窗口大小
        self.window_size_spinbox = QSpinBox()
        self.window_size_spinbox.setRange(0, 100)
        self.window_size_spinbox.setValue(0)  # 默认不设置
        self.window_size_spinbox.setToolTip("(half) Window size, score of position i is the average of the window (i - n) to (i + n)")
        window_layout.addRow("Window size (-w):", self.window_size_spinbox)
        
        # Gap窗口大小
        self.gap_window_spinbox = QSpinBox()
        self.gap_window_spinbox.setRange(0, 100)
        self.gap_window_spinbox.setValue(0)  # 默认不设置
        self.gap_window_spinbox.setToolTip("(half) Window size only applies to statistics/methods based on Gaps")
        window_layout.addRow("Gap window size (-gw):", self.gap_window_spinbox)
        
        # 相似性窗口大小
        self.sim_window_spinbox = QSpinBox()
        self.sim_window_spinbox.setRange(0, 100)
        self.sim_window_spinbox.setValue(0)  # 默认不设置
        self.sim_window_spinbox.setToolTip("(half) Window size only applies to statistics/methods based on Similarity")
        window_layout.addRow("Similarity window size (-sw):", self.sim_window_spinbox)
        
        # 一致性窗口大小
        self.con_window_spinbox = QSpinBox()
        self.con_window_spinbox.setRange(0, 100)
        self.con_window_spinbox.setValue(0)  # 默认不设置
        self.con_window_spinbox.setToolTip("(half) Window size only applies to statistics/methods based on Consistency")
        window_layout.addRow("Consistency window size (-cw):", self.con_window_spinbox)
        
        # Overlap设置
        overlap_group = QGroupBox("Overlap Settings")
        overlap_layout = QFormLayout()
        overlap_group.setLayout(overlap_layout)
        scroll_layout.addWidget(overlap_group)
        
        # 残基重叠
        self.res_overlap_spinbox = QDoubleSpinBox()
        self.res_overlap_spinbox.setRange(0.0, 1.0)
        self.res_overlap_spinbox.setSingleStep(0.05)
        self.res_overlap_spinbox.setValue(0.0)  # 默认不设置
        self.res_overlap_spinbox.setToolTip("Minimum overlap of a positions with other positions in the column to be considered a 'good position'")
        overlap_layout.addRow("Residue overlap (-resoverlap):", self.res_overlap_spinbox)
        
        # 序列重叠
        self.seq_overlap_spinbox = QSpinBox()
        self.seq_overlap_spinbox.setRange(0, 100)
        self.seq_overlap_spinbox.setValue(0)  # 默认不设置
        self.seq_overlap_spinbox.setToolTip("Minimum percentage of 'good positions' that a sequence must have in order to be conserved")
        overlap_layout.addRow("Sequence overlap (-seqoverlap):", self.seq_overlap_spinbox)
        
        # 其他选项
        other_group = QGroupBox("Other Options")
        other_layout = QFormLayout()
        other_group.setLayout(other_layout)
        scroll_layout.addWidget(other_group)
        
        # 移除所有含gap的位点
        self.nogaps_checkbox = QCheckBox("Remove all positions with gaps (-nogaps)")
        other_layout.addRow("", self.nogaps_checkbox)
        
        # 移除全部为gap的列
        self.noallgaps_checkbox = QCheckBox("Remove columns composed only by gaps (-noallgaps)")
        other_layout.addRow("", self.noallgaps_checkbox)
        
        # 互补对齐
        self.complementary_checkbox = QCheckBox("Get the complementary alignment (-complementary)")
        other_layout.addRow("", self.complementary_checkbox)
        
        # 列编号
        self.colnumbering_checkbox = QCheckBox("Get the relationship between the columns in the old and new alignment (-colnumbering)")
        other_layout.addRow("", self.colnumbering_checkbox)
        
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
    
    def load_previous_params(self):
        """加载之前的参数设置"""
        # 如果父窗口有保存的高级参数，则加载它们
        if hasattr(self.parent, 'advanced_params') and self.parent.advanced_params:
            params = self.parent.advanced_params
            # 解析参数并设置到UI控件
            self.parse_and_set_params(params)
    
    def parse_and_set_params(self, params):
        """解析参数并设置到UI控件"""
        i = 0
        while i < len(params):
            param = params[i]
            if param == "-gt" and i + 1 < len(params):
                try:
                    self.gap_threshold_spinbox.setValue(float(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-st" and i + 1 < len(params):
                try:
                    self.sim_threshold_spinbox.setValue(float(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-ct" and i + 1 < len(params):
                try:
                    self.con_threshold_spinbox.setValue(float(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-cons" and i + 1 < len(params):
                try:
                    self.cons_percentage_spinbox.setValue(int(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-w" and i + 1 < len(params):
                try:
                    self.window_size_spinbox.setValue(int(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-gw" and i + 1 < len(params):
                try:
                    self.gap_window_spinbox.setValue(int(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-sw" and i + 1 < len(params):
                try:
                    self.sim_window_spinbox.setValue(int(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-cw" and i + 1 < len(params):
                try:
                    self.con_window_spinbox.setValue(int(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-resoverlap" and i + 1 < len(params):
                try:
                    self.res_overlap_spinbox.setValue(float(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-seqoverlap" and i + 1 < len(params):
                try:
                    self.seq_overlap_spinbox.setValue(int(params[i + 1]))
                    i += 2
                    continue
                except ValueError:
                    pass
            elif param == "-nogaps":
                self.nogaps_checkbox.setChecked(True)
            elif param == "-noallgaps":
                self.noallgaps_checkbox.setChecked(True)
            elif param == "-complementary":
                self.complementary_checkbox.setChecked(True)
            elif param == "-colnumbering":
                self.colnumbering_checkbox.setChecked(True)
            
            i += 1
    
    def reset_to_defaults(self):
        """将所有参数重置为默认值"""
        self.gap_threshold_spinbox.setValue(0.0)
        self.sim_threshold_spinbox.setValue(0.0)
        self.con_threshold_spinbox.setValue(0.0)
        self.cons_percentage_spinbox.setValue(0)
        self.window_size_spinbox.setValue(0)
        self.gap_window_spinbox.setValue(0)
        self.sim_window_spinbox.setValue(0)
        self.con_window_spinbox.setValue(0)
        self.res_overlap_spinbox.setValue(0.0)
        self.seq_overlap_spinbox.setValue(0)
        self.nogaps_checkbox.setChecked(False)
        self.noallgaps_checkbox.setChecked(False)
        self.complementary_checkbox.setChecked(False)
        self.colnumbering_checkbox.setChecked(False)
    
    def get_parameters(self):
        """将高级参数获取为列表"""
        params = []
        
        # 阈值设置
        if self.gap_threshold_spinbox.value() > 0.0:
            params.extend(["-gt", str(self.gap_threshold_spinbox.value())])
        
        if self.sim_threshold_spinbox.value() > 0.0:
            params.extend(["-st", str(self.sim_threshold_spinbox.value())])
        
        if self.con_threshold_spinbox.value() > 0.0:
            params.extend(["-ct", str(self.con_threshold_spinbox.value())])
        
        if self.cons_percentage_spinbox.value() > 0:
            params.extend(["-cons", str(self.cons_percentage_spinbox.value())])
        
        # 窗口大小设置
        if self.window_size_spinbox.value() > 0:
            params.extend(["-w", str(self.window_size_spinbox.value())])
        
        if self.gap_window_spinbox.value() > 0:
            params.extend(["-gw", str(self.gap_window_spinbox.value())])
        
        if self.sim_window_spinbox.value() > 0:
            params.extend(["-sw", str(self.sim_window_spinbox.value())])
        
        if self.con_window_spinbox.value() > 0:
            params.extend(["-cw", str(self.con_window_spinbox.value())])
        
        # Overlap设置
        if self.res_overlap_spinbox.value() > 0.0:
            params.extend(["-resoverlap", str(self.res_overlap_spinbox.value())])
        
        if self.seq_overlap_spinbox.value() > 0:
            params.extend(["-seqoverlap", str(self.seq_overlap_spinbox.value())])
        
        # 其他选项
        if self.nogaps_checkbox.isChecked():
            params.append("-nogaps")
        
        if self.noallgaps_checkbox.isChecked():
            params.append("-noallgaps")
        
        if self.complementary_checkbox.isChecked():
            params.append("-complementary")
        
        if self.colnumbering_checkbox.isChecked():
            params.append("-colnumbering")
        
        return params


class TrimAlPluginEntry:
    """TrimAl插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return TrimAlPlugin(import_from, import_data)
