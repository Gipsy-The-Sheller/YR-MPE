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

from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, QDoubleSpinBox)
from PyQt5.QtCore import Qt, pyqtSignal, QUrl
from PyQt5.QtGui import QIcon
from PyQt5.QtWebEngineWidgets import QWebEngineView
import tempfile
import os
from typing import List, Optional

from Bio import SeqIO

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class GBlocksThread(BaseProcessThread):
    """GBlocks修剪线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "GBlocks"
        
    def execute_commands(self):
        """执行GBlocks修剪命令"""
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
                
                # 构建命令
                cmd = [self.tool_path, input_file]
                
                # 添加参数
                cmd.extend(self.parameters)
                
                # 添加HTML输出参数
                cmd.extend(["-p=y", "-e=-gb"])  # 生成HTML报告
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"GBlocks execution failed for file {i+1}: {result.stderr}")
                    return
                
                # GBlocks会自动生成输出文件，通常为原文件名加-gb后缀
                base_name = os.path.splitext(input_file)[0]
                output_file = f"{base_name}-gb"
                
                # HTML文件通常为原文件名加-gb.htm后缀
                gblocks_html_file = f"{base_name}-gb.htm"
                
                # 如果GBlocks生成了HTML文件，则使用它
                if os.path.exists(gblocks_html_file):
                    html_files.append(gblocks_html_file)
                else:
                    # 创建一个简单的HTML报告
                    html_file = self.create_temp_file(suffix='.html')
                    self.create_simple_html_report(html_file, input_file, output_file, cmd, result)
                    html_files.append(html_file)
                
                # 检查输出文件是否存在，如果不存在则使用输入文件名
                if not os.path.exists(output_file):
                    output_file = input_file
                    
                output_files.append(output_file)
            
            self.progress.emit("All trimming completed")
            self.finished.emit(output_files, html_files)
            
        except Exception as e:
            self.error.emit(f"Trimming exception: {str(e)}")
    
    def create_simple_html_report(self, html_file, input_file, output_file, cmd, result):
        """创建简单的HTML报告"""
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>GBlocks Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                h1 {{ color: #2c3e50; }}
                .info {{ background-color: #e8f4f8; padding: 10px; border-radius: 5px; }}
            </style>
        </head>
        <body>
            <h1>GBlocks Trimming Report</h1>
            <div class="info">
                <p><strong>Input file:</strong> {input_file}</p>
                <p><strong>Output file:</strong> {output_file}</p>
                <p><strong>Command:</strong> {' '.join(cmd)}</p>
            </div>
            <h2>Output</h2>
            <pre>{result.stdout}</pre>
            {('<h2>Errors</h2><pre>' + result.stderr + '</pre>') if result.stderr else ''}
        </body>
        </html>
        """
        
        with open(html_file, 'w') as f:
            f.write(html_content)


class GBlocksPlugin(BasePlugin):
    """GBlocks插件类"""
    import_alignment_signal = pyqtSignal(list)
    
    def __init__(self, import_from=None, import_data=None):
        """初始化GBlocks插件"""
        super().__init__(import_from, import_data)
        # 初始化插件信息
        self.init_plugin_info()
        
        # 特别处理YR-MPEA导入的数据
        if import_from in ["YR_MPEA", "seq_viewer"] and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "GBlocks Trimmer"
        self.tool_name = "GBlocks"
        self.citation = ["""Talavera, G., and Castresana, J. 2007. Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments. Systematic Biology 56: 564-577. doi: <a href='https://doi.org/10.1080/10635150701472164'>10.1080/10635150701472164</a>."""]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], 
                           "NBRF/PIR": ["pir"]}
        self.output_types = {"FASTA": ".fas"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    
    def get_citation(self):
        """获取引用信息"""
        return self.citation

    def handle_import_data(self, import_data):
        """处理从YR-MPEA导入的数据"""
        if isinstance(import_data, list):
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
        elif not self.imported_files:
            # 如果既没有文本也没有文件，恢复初始状态
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
    
    def setup_gblocks_params(self):
        """设置GBlocks参数"""
        # GBlocks参数已在setup_input_tab中设置，此处留空以保持接口一致性
        pass
    
    def setup_output_tab(self):
        """设置输出预览标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 说明文字
        info_label = QLabel("Trimmed alignment results will be displayed here after analysis completes.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # 结果展示区域 - 使用WebEngineView显示HTML
        self.result_view = QWebEngineView()
        layout.addWidget(self.result_view)
        
        # 输出文件信息
        self.output_info = QLabel("Output file information will be displayed here")
        layout.addWidget(self.output_info)
        
    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 序列类型
        seq_type_index = self.seq_type_combo.currentIndex()
        if seq_type_index == 1:  # Protein
            params.append("-t=p")
        elif seq_type_index == 2:  # DNA
            params.append("-t=d")
        elif seq_type_index == 3:  # Codons
            params.append("-t=c")
        # 0为Auto，不添加参数
            
        # 保守位置最小序列数
        if self.b1_spinbox.value() > 0:
            params.append(f"-b1={self.b1_spinbox.value()}")
            
        # 侧翼位置最小序列数
        if self.b2_spinbox.value() > 0:
            params.append(f"-b2={self.b2_spinbox.value()}")
            
        # 最大连续非保守位置数
        if self.b3_spinbox.value() != 8:  # 不是默认值才添加
            params.append(f"-b3={self.b3_spinbox.value()}")
            
        # 区块最小长度
        if self.b4_spinbox.value() != 10:  # 不是默认值才添加
            params.append(f"-b4={self.b4_spinbox.value()}")
            
        # 允许的缺口位置
        gap_index = self.gap_combo.currentIndex()
        if gap_index == 1:  # With Half
            params.append("-b5=h")
        elif gap_index == 2:  # All
            params.append("-b5=a")
        # 0为None，是默认值，不添加参数
            
        # 使用相似性矩阵
        if not self.sim_matrix_checkbox.isChecked():
            params.append("-b6=n")
        
        return params
    
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
                self.result_view.load(QUrl.fromLocalFile(html_file))
            else:
                # 如果没有HTML文件，则显示序列文件内容
                with open(output_files[0], 'r') as f:
                    content = f.read()
                    self.result_view.setHtml(f"<pre>{content}</pre>")
        except Exception as e:
            self.result_view.setHtml(f"<p>Error reading output file: {str(e)}</p>")
        
        # 更新输出信息
        self.output_info.setText(f"Results from: {os.path.basename(output_files[0])}")
    
    def run_analysis(self):
        """运行GBlocks分析"""
        if self.is_running:
            return
            
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "GBlocks executable file not found!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting GBlocks trimming...", "info")
        
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
        
        # 在单独的线程中运行GBlocks
        self.analysis_thread = GBlocksThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()

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
        
        # 在单独的线程中运行比对
        self.alignment_thread = GBlocksThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
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
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from in ["YR_MPEA", "seq_viewer"]:
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
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "GBlocks has been aborted.")
    
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
                
            # Create temporary file
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                f.write(sequence_text)
            return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
    
    def analysis_finished(self, output_files, html_files):
        """分析完成处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 保存HTML文件
        self.current_html_files = html_files
        
        # 显示结果
        self.display_results(output_files)
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        # 添加控制台消息
        self.add_console_message(f"GBlocks completed successfully! Found {len(output_files)} result file(s)", "info")
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from in ["YR_MPEA", "seq_viewer"]:
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        
        # 保存输出文件路径供导入使用
        self.alignment_output_files = output_files
        
        QMessageBox.information(self, "Completed", "Alignment trimming completed!")
    
    def reset_to_defaults(self):
        """将所有参数重置为默认值"""
        self.seq_type_combo.setCurrentIndex(0)  # Auto
        self.b1_spinbox.setValue(0)  # Default
        self.b2_spinbox.setValue(0)  # Default
        self.b3_spinbox.setValue(8)  # Default
        self.b4_spinbox.setValue(10)  # Default
        self.gap_combo.setCurrentIndex(0)  # None
        self.sim_matrix_checkbox.setChecked(True)  # Use similarity matrix
    
    def import_to_platform(self):
        """将修剪结果导入到当前平台"""
        if not hasattr(self, 'alignment_output_files') or not self.alignment_output_files:
            QMessageBox.warning(self, "Warning", "No trimming results to import.")
            return
            
        try:
            # 解析修剪结果文件
            sequences = []
            for output_file in self.alignment_output_files:
                # 读取文件
                for record in SeqIO.parse(output_file, "fasta"):
                    sequences.append(record)
            
            if not sequences:
                QMessageBox.warning(self, "Warning", "No sequences found in trimming results.")
                return
                
            # 发送信号将数据导入到平台
            self.import_alignment_to_platform(sequences)
            
            # 显示成功消息
            QMessageBox.information(self, "Success", f"Successfully imported {len(sequences)} sequences to the platform.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import trimming results: {str(e)}")

    def import_alignment_to_platform(self, sequences):
        """将修剪结果导入到平台的工作区"""
        # 发送信号将数据导入到平台
        self.import_alignment_signal.emit(sequences)


class GBlocksPluginEntry:
    """GBlocks插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return GBlocksPlugin(import_from, import_data)