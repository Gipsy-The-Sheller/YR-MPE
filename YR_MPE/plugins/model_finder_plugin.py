# model_finder_plugin.py
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
                             QSpinBox, QCheckBox, QLabel, QComboBox, QTableWidget, QTableWidgetItem,
                             QTextEdit, QToolButton, QApplication, QFrame, QWidget)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon, QFont
import tempfile
import os
import re
from typing import List, Optional
from Bio import SeqIO

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class ModelFinderThread(BaseProcessThread):
    """ModelFinder线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "IQ-TREE ModelFinder"
        
    def execute_commands(self):
        """执行ModelFinder命令"""
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
                cmd = [
                    self.tool_path,
                    "-s", input_file,
                    "-m", "TESTONLY",
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"IQ-TREE ModelFinder execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的.iqtree文件
                iqtree_file = input_file + ".iqtree"
                if os.path.exists(iqtree_file):
                    output_files.append(iqtree_file)
                else:
                    # 尝试其他可能的命名方式
                    alternate_iqtree = input_file + ".phy.iqtree"
                    if os.path.exists(alternate_iqtree):
                        output_files.append(alternate_iqtree)
                    else:
                        self.console_output.emit(f"Warning: Could not find .iqtree file for {input_file}", "warning")
            
            self.progress.emit("Model finding completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Model finding exception: {str(e)}")


class ModelFinderPlugin(BasePlugin):
    """IQ-TREE ModelFinder插件类"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_model_result_signal = pyqtSignal(dict)  # 导出模型结果信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化ModelFinder插件"""
        super().__init__(import_from, import_data)
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "IQ-TREE ModelFinder"
        self.tool_name = "IQ-Tree 3"
        self.citation = [
            """Bui Quang Minh, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., in press. https://doi.org/10.1038/s41586-020-2176-9""",
            """Subha Kalyaanamoorthy, Bui Quang Minh, Thomas KF Wong, Arndt von Haeseler, and Lars S Jermiin (2017) ModelFinder: Fast model selection for accurate phylogenetic estimates. Nature Methods, 14:587–589. https://doi.org/10.1038/nmeth.4285"""
        ]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "PHYLIP": ["phy"]}
        self.output_types = {"IQ-TREE": ".iqtree"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

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
        """设置输入标签页"""
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
        self.file_path_edit.setPlaceholderText("Select sequence files...")
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
        self.sequence_text.setPlaceholderText("Or paste sequence in FASTA format...")
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
        params_group = QGroupBox("ModelFinder Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # 序列类型
        self.seq_type_combo = QComboBox()
        self.seq_type_combo.addItems(["AUTO", "DNA", "PROT"])
        self.seq_type_combo.setCurrentText("AUTO")
        params_layout.addRow("Sequence Type:", self.seq_type_combo)
        
        # 准则
        self.criterion_combo = QComboBox()
        self.criterion_combo.addItems(["BIC", "AICc", "AIC"])
        self.criterion_combo.setCurrentText("BIC")
        params_layout.addRow("Criterion:", self.criterion_combo)
        
        # 线程数
        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 5)
        self.threads_spinbox.setValue(1)
        self.threads_spinbox.setSpecialValueText("AUTO")
        params_layout.addRow("Threads:", self.threads_spinbox)

        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets
    
    def browse_files(self):
        """浏览选择文件"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select sequence files", "", "Sequence files (*.fas *.fna *.fa *.fasta *.phy);;All files (*)"
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
        info_label = QLabel("Model selection results will be displayed here after analysis completes.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # 结果表格
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(8)
        self.results_table.setHorizontalHeaderLabels([
            "Model", "LogL", "AIC", "w-AIC", "AICc", "w-AICc", "BIC", "w-BIC"
        ])
        self.results_table.setSortingEnabled(True)
        layout.addWidget(self.results_table)
        
        # 输出文件信息
        self.output_info = QLabel("Output file information will be displayed here")
        layout.addWidget(self.output_info)
        
        # # 添加导出到YR-MPEA按钮
        # self.export_to_mpea_btn = QPushButton("Export Best Model to YR-MPEA")
        # self.export_to_mpea_btn.clicked.connect(self.export_to_mpea)
        # self.export_to_mpea_btn.setVisible(False)  # 初始隐藏
        # layout.addWidget(self.export_to_mpea_btn)
        
        # # 添加导出完整模型表到YR-MPEA按钮
        # self.export_table_to_mpea_btn = QPushButton("Export Full Model Table to YR-MPEA")
        # self.export_table_to_mpea_btn.clicked.connect(self.export_table_to_mpea)
        # self.export_table_to_mpea_btn.setVisible(False)  # 初始隐藏
        # layout.addWidget(self.export_table_to_mpea_btn)
        
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.export_table_to_mpea)
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
                QMessageBox.warning(self, "Warning", "Please input sequence text!")
                return None
                
            # Create temporary file
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                f.write(sequence_text)
            return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None

    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 序列类型
        seq_type = self.seq_type_combo.currentText()
        seq_type_code = {"dna": "DNA", "prot": "AA"}[seq_type.lower()]
        if seq_type != "AUTO":
            params.extend(["-st", seq_type_code])
        
        # 准则
        criterion = self.criterion_combo.currentText()
        if criterion == "BIC":
            pass
        elif criterion == "AICc":
            params.extend(["-AICc"])
        elif criterion == "AIC":
            params.extend(["-AIC"])
        else:
            pass
        params.extend(["-mset", "ALL"])  # 测试所有模型
        
        # 线程数
        threads = self.threads_spinbox.value()
        if threads > 1:
            params.extend(["-nt", str(threads)])
        
        return params
    
    def parse_iqtree_file(self, iqtree_file_path):
        """解析.iqtree文件并提取模型信息"""
        models_data = []
        
        try:
            with open(iqtree_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 查找模型列表部分
            model_section_pattern = r"List of models sorted by.*?\n\n(.*?)\n\n"
            model_section_match = re.search(model_section_pattern, content, re.DOTALL)
            
            if model_section_match:
                model_lines = model_section_match.group(1).strip().split('\n')
                if len(model_lines) > 1:
                    # 跳过标题行
                    for line in model_lines[1:]:
                        # 使用正则表达式提取模型数据
                        # 格式: Model LogL AIC w-AIC AICc w-AICc BIC w-BIC
                        parts = line.split()
                        if len(parts) >= 10:
                            model_data = {
                                'Model': parts[0],
                                'LogL': parts[1],
                                'AIC': parts[2],
                                'w-AIC': parts[4],
                                'AICc': parts[5],
                                'w-AICc': parts[7],
                                'BIC': parts[8],
                                'w-BIC': parts[10]
                            }
                            models_data.append(model_data)
            
            return models_data
        except Exception as e:
            self.add_console_message(f"Error parsing .iqtree file: {str(e)}", "error")
            return []
    
    def display_results(self, iqtree_files):
        """在输出标签页中显示结果"""
        if not iqtree_files:
            self.results_table.setRowCount(0)
            return
            
        # 解析第一个文件的结果（如果有多个文件，只显示第一个）
        models_data = self.parse_iqtree_file(iqtree_files[0])
        
        # 更新表格
        self.results_table.setRowCount(len(models_data))
        
        for row, model_data in enumerate(models_data):
            for col, (key, value) in enumerate(model_data.items()):
                item = QTableWidgetItem(value)
                if key not in ['Model']:
                    # 数值列右对齐
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.results_table.setItem(row, col, item)
        
        # 调整列宽
        self.results_table.resizeColumnsToContents()
        
        # 更新输出信息
        self.output_info.setText(f"Results from: {os.path.basename(iqtree_files[0])}")
        
        # 显示导出按钮
        # self.export_to_mpea_btn.setVisible(True)
        # self.export_table_to_mpea_btn.setVisible(True)
    
    def run_analysis(self):
        """运行ModelFinder分析"""
        if self.is_running:
            return
            
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "IQ-TREE executable file not found!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting IQ-TREE ModelFinder...", "info")
        
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
        
        # 在单独的线程中运行ModelFinder
        self.analysis_thread = ModelFinderThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
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
        
        # 显示结果
        self.display_results(output_files)
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        # 添加控制台消息
        self.add_console_message(f"ModelFinder completed successfully! Found {len(output_files)} result file(s)", "info")
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from == "YR_MPEA":
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
            
        # 显示导出按钮
        # self.export_to_mpea_btn.setVisible(True)
        
        # 保存输出文件路径供导入使用
        self.alignment_output_files = output_files
        
        QMessageBox.information(self, "Completed", "Model selection completed!")
    
    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 添加控制台消息
        self.add_console_message(f"ModelFinder failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"ModelFinder failed: {error_message}")
    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "ModelFinder has been aborted.")

    def import_alignment_to_platform(self, sequences):
        """将比对结果导入到平台的工作区"""
        # 发送信号将数据导入到平台
        self.import_alignment_signal.emit(sequences)

    def export_to_mpea(self):
        """导出最佳模型到YR-MPEA"""
        if self.results_table.rowCount() == 0:
            QMessageBox.warning(self, "Warning", "No model results to export.")
            return
            
        # 获取第一行（最佳模型）的数据
        best_model_data = {}
        for col in range(self.results_table.columnCount()):
            header = self.results_table.horizontalHeaderItem(col).text()
            value = self.results_table.item(0, col).text()
            best_model_data[header] = value
            
        # 发送信号将模型结果导出到YR-MPEA
        self.export_model_result_signal.emit(best_model_data)
        
        QMessageBox.information(self, "Success", f"Best model '{best_model_data['Model']}' exported to YR-MPEA successfully!")

    def export_table_to_mpea(self):
        """导出完整模型表到YR-MPEA"""
        if self.results_table.rowCount() == 0:
            QMessageBox.warning(self, "Warning", "No model results to export.")
            return
            
        # 获取完整的模型表数据
        model_table_data = []
        headers = []
        
        # 获取表头
        for col in range(self.results_table.columnCount()):
            headers.append(self.results_table.horizontalHeaderItem(col).text())
        
        # 获取所有行数据
        for row in range(self.results_table.rowCount()):
            row_data = {}
            for col in range(self.results_table.columnCount()):
                header = headers[col]
                value = self.results_table.item(row, col).text()
                row_data[header] = value
            model_table_data.append(row_data)
            
        # 发送信号将完整模型表导出到YR-MPEA
        table_export_data = {
            "type": "model_table",
            "headers": headers,
            "data": model_table_data
        }
        self.export_model_result_signal.emit(table_export_data)
        
        QMessageBox.information(self, "Success", "Full model table exported to YR-MPEA successfully!")

    def copy_citation(self):
        """复制引用"""
        clipboard = QApplication.clipboard()
        clipboard.setText("\n".join(self.citation))
        QMessageBox.information(self, "Copied", "Citation copied to clipboard!")


class ModelFinderPluginEntry:
    """ModelFinder插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return ModelFinderPlugin(import_from, import_data)