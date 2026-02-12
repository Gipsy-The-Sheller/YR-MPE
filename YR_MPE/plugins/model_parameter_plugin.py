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
                             QTextEdit, QToolButton, QApplication, QFrame, QWidget, QDialog,
                             QHeaderView, QSplitter, QTabWidget, QProgressBar)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon, QFont, QColor, QBrush
import tempfile
import os
import re
from typing import List, Optional, Dict
from Bio import SeqIO
import json

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.model_parameter_parser import parse_iqtree_file, validate_dna_sequences, get_q_matrix_color_scale, format_q_matrix_for_display
import subprocess


class ModelParameterThread(BaseProcessThread):
    """Model Parameter Estimation线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "IQ-TREE Model Parameter Estimation"
        
    def execute_commands(self):
        """执行IQ-TREE命令进行模型参数估计"""
        try:
            output_files = []
            
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
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"IQ-TREE execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的.iqtree文件
                iqtree_file = input_file + ".iqtree"
                if os.path.exists(iqtree_file):
                    output_files.append(iqtree_file)
                else:
                    self.console_output.emit(f"Warning: Could not find .iqtree file for {input_file}", "warning")
            
            self.progress.emit("Model parameter estimation completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Model parameter estimation exception: {str(e)}")


class ModelParameterPlugin(BasePlugin):
    """IQ-TREE Model Parameter Estimation插件类"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_model_result_signal = pyqtSignal(dict)  # 导出模型结果信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化Model Parameter Estimation插件"""
        super().__init__(import_from, import_data)
        
        # 存储解析的参数数据
        self.parsed_parameters = {}  # {file_path: parameters_dict}
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "Estimate Model Parameters (DNA, ML)"
        self.tool_name = "IQ-TREE 3"
        self.citation = [
            """Bui Quang Minh, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., in press. https://doi.org/10.1093/molbev/msaa015"""
        ]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "PHYLIP": ["phy"]}
        self.output_types = {"IQ-TREE Report": ".iqtree"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

    def handle_import_data(self, import_data):
        """处理从YR-MPEA导入的数据"""
        if isinstance(import_data, list):
            # 验证是否为DNA序列
            is_valid, error_msg = validate_dna_sequences(import_data)
            if not is_valid:
                QMessageBox.warning(self, "Validation Error", error_msg)
                return
            
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
        self.file_path_edit.setPlaceholderText("Select DNA alignment files...")
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
        
        # .iqtree文件直接导入选项
        self.direct_iqtree_checkbox = QCheckBox("Load from existing .iqtree files (skip estimation)")
        self.direct_iqtree_checkbox.stateChanged.connect(self.on_direct_iqtree_toggled)
        input_layout.addRow("", self.direct_iqtree_checkbox)
        
        # 文本输入
        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste DNA alignment in FASTA format...")
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
        params_group = QGroupBox("Model Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # 模型选择
        model_layout = QHBoxLayout()
        self.model_combo = QComboBox()
        self.model_combo.addItems(["AUTO", "JC", "K2P", "HKY", "TVM", "GTR"])
        self.model_combo.setCurrentText("AUTO")
        self.model_combo.currentTextChanged.connect(self.on_model_changed)
        model_layout.addWidget(QLabel("Substitution Model:"))
        model_layout.addWidget(self.model_combo)
        params_layout.addRow(model_layout)
        
        # Gamma分布参数
        gamma_layout = QHBoxLayout()
        self.gamma_checkbox = QCheckBox("+G")
        self.gamma_spinbox = QSpinBox()
        self.gamma_spinbox.setRange(1, 10)
        self.gamma_spinbox.setValue(4)
        self.gamma_spinbox.setEnabled(False)
        self.gamma_checkbox.stateChanged.connect(
            lambda state: self.gamma_spinbox.setEnabled(state == Qt.Checked)
        )
        gamma_layout.addWidget(self.gamma_checkbox)
        gamma_layout.addWidget(self.gamma_spinbox)
        params_layout.addRow("Gamma Distribution:", gamma_layout)
        
        # Invariable Sites参数
        self.invar_checkbox = QCheckBox("+I")
        params_layout.addRow("Invariable Sites:", self.invar_checkbox)
        
        # Empirical Frequencies参数
        self.empirical_checkbox = QCheckBox("+F (Empirical frequencies)")
        params_layout.addRow("Empirical Frequencies:", self.empirical_checkbox)
        
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
    
    def on_direct_iqtree_toggled(self, state):
        """处理直接导入.iqtree文件选项切换"""
        if state == Qt.Checked:
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
            self.params_group.setVisible(False)  # 禁用参数组
            self.add_console_message("Direct .iqtree file import mode enabled. Parameters will be extracted from existing files.", "info")
        else:
            if not self.imported_files:
                self.sequence_text.setVisible(True)
                self.sequence_text.setEnabled(True)
            self.params_group.setVisible(True)  # 启用参数组
            self.add_console_message("Direct .iqtree file import mode disabled. Parameters will be estimated from alignments.", "info")
    
    def browse_files(self):
        """浏览选择文件"""
        if self.direct_iqtree_checkbox.isChecked():
            # 选择.iqtree文件
            file_paths, _ = QFileDialog.getOpenFileNames(
                self, "Select .iqtree files", "", "IQ-TREE Report files (*.iqtree);;All files (*)"
            )
        else:
            # 选择比对文件
            file_paths, _ = QFileDialog.getOpenFileNames(
                self, "Select DNA alignment files", "", "Alignment files (*.fas *.fna *.fa *.fasta *.phy);;All files (*)"
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
    
    def on_model_changed(self):
        """Handle model change"""
        if self.model_combo.currentText() == "AUTO":
            self.gamma_checkbox.setEnabled(True)
            self.invar_checkbox.setEnabled(True)
            self.empirical_checkbox.setEnabled(True)
        else:
            # For specific models, enable options
            self.gamma_checkbox.setEnabled(True)
            self.invar_checkbox.setEnabled(True)
            self.empirical_checkbox.setEnabled(True)

    def setup_output_tab(self):
        """设置输出预览标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 创建分割器
        splitter = QSplitter(Qt.Vertical)
        
        # 上部分：文件选择和参数概览
        top_widget = QWidget()
        top_layout = QVBoxLayout()
        top_widget.setLayout(top_layout)
        
        # 文件选择器
        file_selector_layout = QHBoxLayout()
        file_selector_layout.addWidget(QLabel("Select file:"))
        
        self.output_file_combo = QComboBox()
        self.output_file_combo.setEnabled(False)
        self.output_file_combo.currentIndexChanged.connect(self.on_output_file_changed)
        file_selector_layout.addWidget(self.output_file_combo)
        
        file_selector_layout.addStretch()
        top_layout.addLayout(file_selector_layout)
        
        # 模型名称显示
        self.model_name_label = QLabel("No model loaded")
        self.model_name_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        top_layout.addWidget(self.model_name_label)
        
        # 参数概览表格
        self.params_table = QTableWidget()
        self.params_table.setColumnCount(2)
        self.params_table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self.params_table.verticalHeader().setVisible(False)
        self.params_table.setEditTriggers(QTableWidget.NoEditTriggers)
        top_layout.addWidget(QLabel("Parameter Summary:"))
        top_layout.addWidget(self.params_table)
        
        splitter.addWidget(top_widget)
        
        # 下部分：详细参数显示
        bottom_widget = QWidget()
        bottom_layout = QVBoxLayout()
        bottom_widget.setLayout(bottom_layout)
        
        # 创建标签页用于Q矩阵和碱基频率
        detail_tabs = QTabWidget()
        
        # Q矩阵标签页
        q_matrix_tab = QWidget()
        q_matrix_layout = QVBoxLayout()
        q_matrix_tab.setLayout(q_matrix_layout)
        
        q_matrix_layout.addWidget(QLabel("Q Matrix (Instantaneous Rate Matrix):"))
        
        self.q_matrix_table = QTableWidget()
        self.q_matrix_table.setRowCount(4)
        self.q_matrix_table.setColumnCount(4)
        self.q_matrix_table.setHorizontalHeaderLabels(['A', 'C', 'G', 'T'])
        self.q_matrix_table.setVerticalHeaderLabels(['A', 'C', 'G', 'T'])
        self.q_matrix_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.q_matrix_table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.q_matrix_table.setEditTriggers(QTableWidget.NoEditTriggers)
        q_matrix_layout.addWidget(self.q_matrix_table)
        
        detail_tabs.addTab(q_matrix_tab, "Q Matrix")
        
        # 碱基频率标签页
        freq_tab = QWidget()
        freq_layout = QVBoxLayout()
        freq_tab.setLayout(freq_layout)
        
        freq_layout.addWidget(QLabel("State Frequencies:"))
        
        self.freq_table = QTableWidget()
        self.freq_table.setRowCount(4)
        self.freq_table.setColumnCount(2)
        self.freq_table.setHorizontalHeaderLabels(["Base", "Frequency"])
        self.freq_table.setVerticalHeaderLabels(["A", "C", "G", "T"])
        self.freq_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.freq_table.verticalHeader().setVisible(True)
        self.freq_table.setEditTriggers(QTableWidget.NoEditTriggers)
        freq_layout.addWidget(self.freq_table)
        
        # 添加堆叠条形图（使用简单的进度条表示）
        freq_chart_layout = QHBoxLayout()
        freq_chart_layout.addWidget(QLabel("Frequency Distribution:"))
        freq_layout.addLayout(freq_chart_layout)
        
        freq_progress_layout = QHBoxLayout()
        self.freq_progress_a = QProgressBar()  # A
        self.freq_progress_c = QProgressBar()  # C
        self.freq_progress_g = QProgressBar()  # G
        self.freq_progress_t = QProgressBar()  # T
        
        # 设置样式
        for bar in [self.freq_progress_a, self.freq_progress_c, self.freq_progress_g, self.freq_progress_t]:
            bar.setTextVisible(True)
            bar.setRange(0, 100)
            bar.setStyleSheet("""
                QProgressBar {
                    border: 2px solid grey;
                    border-radius: 5px;
                    text-align: center;
                }
                QProgressBar::chunk {
                    background-color: #3b82f6;
                }
            """)
        
        self.freq_progress_a.setStyleSheet("""
            QProgressBar::chunk {
                background-color: #ef4444;
            }
        """)
        self.freq_progress_c.setStyleSheet("""
            QProgressBar::chunk {
                background-color: #22c55e;
            }
        """)
        self.freq_progress_g.setStyleSheet("""
            QProgressBar::chunk {
                background-color: #eab308;
            }
        """)
        self.freq_progress_t.setStyleSheet("""
            QProgressBar::chunk {
                background-color: #3b82f6;
            }
        """)
        
        freq_progress_layout.addWidget(QLabel("A:"))
        freq_progress_layout.addWidget(self.freq_progress_a)
        freq_progress_layout.addWidget(QLabel("C:"))
        freq_progress_layout.addWidget(self.freq_progress_c)
        freq_progress_layout.addWidget(QLabel("G:"))
        freq_progress_layout.addWidget(self.freq_progress_g)
        freq_progress_layout.addWidget(QLabel("T:"))
        freq_progress_layout.addWidget(self.freq_progress_t)
        
        freq_layout.addLayout(freq_progress_layout)
        
        detail_tabs.addTab(freq_tab, "State Frequencies")
        
        bottom_layout.addWidget(detail_tabs)
        splitter.addWidget(bottom_widget)
        
        layout.addWidget(splitter)
        
        # 导出按钮
        export_layout = QHBoxLayout()
        export_csv_btn = QPushButton("Export to CSV")
        export_csv_btn.clicked.connect(self.export_to_csv)
        export_json_btn = QPushButton("Export to JSON")
        export_json_btn.clicked.connect(self.export_to_json)
        export_layout.addWidget(export_csv_btn)
        export_layout.addWidget(export_json_btn)
        export_layout.addStretch()
        layout.addLayout(export_layout)
    
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
                QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
                return None
                
            # Validate DNA sequences
            sequences = list(SeqIO.parse(tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fas'), 'fasta'))
            # Actually parse from text
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                f.write(sequence_text)
            
            # Validate
            sequences = list(SeqIO.parse(temp_file, 'fasta'))
            is_valid, error_msg = validate_dna_sequences(sequences)
            if not is_valid:
                QMessageBox.warning(self, "Validation Error", error_msg)
                os.remove(temp_file)
                return None
            
            return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None

    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 模型
        model_text = self.model_combo.currentText()
        
        if model_text != "AUTO":
            model = model_text
            
            # 添加模型扩展参数
            if self.gamma_checkbox.isChecked():
                model += f"+G{self.gamma_spinbox.value()}"
            
            if self.invar_checkbox.isChecked():
                model += "+I"
                
            if self.empirical_checkbox.isChecked():
                model += "+F"
            
            params.extend(["-m", model])
        else:
            # 使用ModelFinder自动选择模型
            params.extend(["-m", "MFP"])
        
        # 线程数
        threads = self.threads_spinbox.value()
        if threads > 1:
            params.extend(["-nt", str(threads)])
        
        params.extend(["-redo"])  # 重新运行
        
        return params
    
    def run_analysis(self):
        """运行模型参数估计分析"""
        # 检查是否直接导入.iqtree文件
        if self.direct_iqtree_checkbox.isChecked():
            if not self.imported_files:
                QMessageBox.warning(self, "Warning", "Please select .iqtree files first!")
                return
            
            self.parse_direct_iqtree_files()
            return
        
        # 常规估计流程
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
        self.add_console_message("Starting IQ-TREE model parameter estimation...", "info")
        
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
        
        # 创建并启动线程
        self.analysis_thread = ModelParameterThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
        )
        
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()
    
    def parse_direct_iqtree_files(self):
        """直接解析.iqtree文件"""
        self.parsed_parameters = {}
        
        for iqtree_file in self.imported_files:
            if not os.path.exists(iqtree_file):
                self.add_console_message(f"File not found: {iqtree_file}", "error")
                continue
            
            params = parse_iqtree_file(iqtree_file)
            if params:
                self.parsed_parameters[iqtree_file] = params
                self.add_console_message(f"Successfully parsed: {os.path.basename(iqtree_file)}", "info")
            else:
                self.add_console_message(f"Failed to parse: {os.path.basename(iqtree_file)}", "error")
        
        # 显示结果
        self.display_parsed_results()
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        self.add_console_message(f"Parsed {len(self.parsed_parameters)} file(s)", "info")
    
    def analysis_finished(self, output_files, html_files):
        """分析完成处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 解析.iqtree文件
        self.parsed_parameters = {}
        
        for iqtree_file in output_files:
            if os.path.exists(iqtree_file):
                params = parse_iqtree_file(iqtree_file)
                if params:
                    self.parsed_parameters[iqtree_file] = params
        
        # 显示结果
        self.display_parsed_results()
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        # 添加控制台消息
        self.add_console_message(f"Model parameter estimation completed! Found {len(output_files)} result file(s)", "info")
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from == "YR_MPEA":
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
    
    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 添加控制台消息
        self.add_console_message(f"Model parameter estimation failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"Model parameter estimation failed: {error_message}")
    
    def display_parsed_results(self):
        """显示解析的参数结果"""
        # 更新文件选择器
        self.output_file_combo.clear()
        
        if not self.parsed_parameters:
            self.model_name_label.setText("No model loaded")
            self.output_file_combo.setEnabled(False)
            return
        
        # 添加文件选项
        for iqtree_file in self.parsed_parameters.keys():
            self.output_file_combo.addItem(os.path.basename(iqtree_file), iqtree_file)
        
        self.output_file_combo.setEnabled(True)
        
        # 显示第一个文件的结果
        if self.parsed_parameters:
            self.on_output_file_changed(0)
    
    def on_output_file_changed(self, index):
        """文件选择改变时的处理"""
        if index < 0 or index >= self.output_file_combo.count():
            return
        
        iqtree_file = self.output_file_combo.itemData(index)
        params = self.parsed_parameters.get(iqtree_file)
        
        if not params:
            return
        
        # 显示模型名称
        model_name = params.get('model_name', 'Unknown')
        self.model_name_label.setText(f"Model: {model_name}")
        
        # 显示参数概览
        self.params_table.setRowCount(0)
        
        # 添加速率异质性参数
        rate_hetero = params.get('rate_heterogeneity', {})
        self.add_param_row("Rate Heterogeneity", rate_hetero.get('model', 'N/A'))
        
        if rate_hetero.get('gamma_alpha') is not None:
            self.add_param_row("Gamma Shape (α)", f"{rate_hetero['gamma_alpha']:.4f}")
        
        if rate_hetero.get('invariant_sites') is not None:
            self.add_param_row("Invariant Sites", f"{rate_hetero['invariant_sites']:.4f}")
        
        # 显示Q矩阵
        q_matrix = params.get('q_matrix', [])
        if q_matrix:
            min_val, max_val = get_q_matrix_color_scale(q_matrix)
            
            for i in range(4):
                for j in range(4):
                    item = QTableWidgetItem(f"{q_matrix[i][j]:.4f}")
                    item.setTextAlignment(Qt.AlignCenter | Qt.AlignVCenter)
                    
                    # 应用颜色（对角线为灰色，其他根据数值着色）
                    if i == j:
                        item.setBackground(QColor(240, 240, 240))
                    else:
                        # 计算颜色（从蓝色到红色）
                        value = q_matrix[i][j]
                        normalized = (value - min_val) / (max_val - min_val) if max_val > min_val else 0.5
                        # 蓝色(0,0,255)到红色(255,0,0)
                        red = int(255 * normalized)
                        blue = int(255 * (1 - normalized))
                        green = 0
                        item.setBackground(QColor(red, green, blue))
                        # 根据背景色设置文字颜色
                        if normalized > 0.5:
                            item.setForeground(QColor(255, 255, 255))
                        else:
                            item.setForeground(QColor(0, 0, 0))
                    
                    self.q_matrix_table.setItem(i, j, item)
        
        # 显示碱基频率
        state_freqs = params.get('state_frequencies', {})
        bases = ['A', 'C', 'G', 'T']
        total = sum(state_freqs.values()) if state_freqs else 0
        
        for i, base in enumerate(bases):
            freq_key = f'pi({base})'
            freq = state_freqs.get(freq_key, 0.0)
            
            # 表格显示
            item = QTableWidgetItem(f"{freq:.4f}")
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.freq_table.setItem(i, 0, QTableWidgetItem(base))
            self.freq_table.setItem(i, 1, item)
            
            # 进度条显示
            percentage = (freq / total * 100) if total > 0 else 0
            if base == 'A':
                self.freq_progress_a.setValue(int(percentage))
                self.freq_progress_a.setFormat(f"{percentage:.2f}%")
            elif base == 'C':
                self.freq_progress_c.setValue(int(percentage))
                self.freq_progress_c.setFormat(f"{percentage:.2f}%")
            elif base == 'G':
                self.freq_progress_g.setValue(int(percentage))
                self.freq_progress_g.setFormat(f"{percentage:.2f}%")
            elif base == 'T':
                self.freq_progress_t.setValue(int(percentage))
                self.freq_progress_t.setFormat(f"{percentage:.2f}%")
    
    def add_param_row(self, param_name, param_value):
        """添加参数行到表格"""
        row = self.params_table.rowCount()
        self.params_table.insertRow(row)
        
        name_item = QTableWidgetItem(param_name)
        value_item = QTableWidgetItem(str(param_value))
        
        self.params_table.setItem(row, 0, name_item)
        self.params_table.setItem(row, 1, value_item)
    
    def export_to_csv(self):
        """导出参数到CSV文件"""
        if not self.parsed_parameters:
            QMessageBox.warning(self, "Warning", "No data to export!")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export to CSV", "", "CSV files (*.csv);;All files (*)"
        )
        
        if not file_path:
            return
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                # 导出每个文件的参数
                for iqtree_file, params in self.parsed_parameters.items():
                    f.write(f"File: {os.path.basename(iqtree_file)}\n")
                    f.write(f"Model: {params.get('model_name', 'Unknown')}\n\n")
                    
                    # 速率异质性
                    rate_hetero = params.get('rate_heterogeneity', {})
                    f.write(f"Rate Heterogeneity: {rate_hetero.get('model', 'N/A')}\n")
                    if rate_hetero.get('gamma_alpha') is not None:
                        f.write(f"Gamma Shape (α): {rate_hetero['gamma_alpha']:.4f}\n")
                    if rate_hetero.get('invariant_sites') is not None:
                        f.write(f"Invariant Sites: {rate_hetero['invariant_sites']:.4f}\n")
                    
                    f.write("\n")
                    
                    # Q矩阵
                    q_matrix = params.get('q_matrix', [])
                    if q_matrix:
                        f.write("Q Matrix:\n")
                        f.write("    A        C        G        T\n")
                        for i, row in enumerate(q_matrix):
                            row_label = ['A', 'C', 'G', 'T'][i]
                            f.write(f"{row_label} ")
                            for val in row:
                                f.write(f"{val:8.4f} ")
                            f.write("\n")
                    
                    f.write("\n")
                    
                    # 碱基频率
                    state_freqs = params.get('state_frequencies', {})
                    if state_freqs:
                        f.write("State Frequencies:\n")
                        for base in ['A', 'C', 'G', 'T']:
                            freq_key = f'pi({base})'
                            freq = state_freqs.get(freq_key, 0.0)
                            f.write(f"  {freq_key}: {freq:.4f}\n")
                    
                    f.write("\n" + "="*50 + "\n\n")
            
            self.add_console_message(f"Exported to CSV: {file_path}", "info")
            QMessageBox.information(self, "Success", "Exported to CSV successfully!")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Export failed: {e}")
    
    def export_to_json(self):
        """导出参数到JSON文件"""
        if not self.parsed_parameters:
            QMessageBox.warning(self, "Warning", "No data to export!")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export to JSON", "", "JSON files (*.json);;All files (*)"
        )
        
        if not file_path:
            return
        
        try:
            # 准备导出数据
            export_data = {}
            for iqtree_file, params in self.parsed_parameters.items():
                export_data[os.path.basename(iqtree_file)] = params
            
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(export_data, f, indent=2)
            
            self.add_console_message(f"Exported to JSON: {file_path}", "info")
            QMessageBox.information(self, "Success", "Exported to JSON successfully!")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Export failed: {e}")
    
    def import_to_platform(self):
        """导入到当前平台"""
        # TODO: 实现导入到平台的功能
        self.add_console_message("Import to platform functionality not yet implemented", "warning")


class ModelParameterPluginEntry:
    """插件入口类"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return ModelParameterPlugin(import_from, import_data)