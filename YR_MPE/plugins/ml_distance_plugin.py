# ml_distance_plugin.py
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

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QTextEdit,
                             QTabWidget, QToolButton, QApplication, QFrame, QTableWidget,
                             QHeaderView, QTableWidgetItem)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QTextCursor, QIcon
import tempfile
import os
import re
from typing import List, Optional

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess
from Bio import SeqIO


class MLDistanceThread(BaseProcessThread):
    """ML Distance计算线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "IQ-TREE ML Distance"
        
    def execute_commands(self):
        """执行ML Distance计算命令"""
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
                    # "-m", "DIST",
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"IQ-TREE ML Distance execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的.mldist文件
                mldist_file = input_file + ".mldist"
                if os.path.exists(mldist_file):
                    output_files.append(mldist_file)
                else:
                    # 尝试其他可能的命名方式
                    alternate_mldist = input_file + ".phy.mldist"
                    if os.path.exists(alternate_mldist):
                        output_files.append(alternate_mldist)
                    else:
                        self.console_output.emit(f"Warning: Could not find .mldist file for {input_file}", "warning")
            
            self.progress.emit("ML distance calculation completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"ML distance calculation exception: {str(e)}")


class MLDistancePlugin(BasePlugin):
    """IQ-TREE ML Distance插件类"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_distance_result_signal = pyqtSignal(dict)  # 导出模型结果信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化ML Distance插件"""
        super().__init__(import_from, import_data)
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "IQ-TREE ML Distance"
        self.tool_name = "IQ-TREE 3"
        self.citation = [
            """Bui Quang Minh, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., in press. https://doi.org/10.1093/molbev/msaa015"""
        ]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "PHYLIP": ["phy"]}
        self.output_types = {"ML Distance Matrix": ".mldist"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    
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
        params_group = QGroupBox("IQ-TREE ML Distance Parameters")
        params_form_layout = QFormLayout()
        params_group.setLayout(params_form_layout)
        layout.addWidget(params_group)
        
        # 序列类型
        self.seq_type_combo = QComboBox()
        self.seq_type_combo.addItems(["AUTO", "DNA", "AA"])
        params_form_layout.addRow("Sequence Type:", self.seq_type_combo)
        
        # 模型选择
        self.model_combo = QComboBox()
        # DNA模型
        dna_models = [
            ["JC69", "JC"], ["F81"], ["K2P", "K80"], ["HKY85", "HKY"], ["TNe"], ["TN93", "TN"], ["K3P", "K81"], ["K81u", "K3Pu"], ["TPM2"], ["TPM2u"], 
            ["TPM3"], ["TPM3u"], ["TIM"], ["TIMe"], ["TIM2"], ["TIM2e"], ["TIM3"], ["TIM3e"], 
            ["TVM"], ["TVMe"], ["SYM"], ["GTR"]
        ]
        # 蛋白质模型
        protein_models = [
            "Blosum62", "cpREV", "Dayhoff", "DCMut", "EAL", "ELM", "FLAVI", "FLU", "GTR20", "HIVb", "HIVw", "JTT",
            "JTTDCMut", "LG", "mtART", "mtMAM", "mtREV", "mtZOA", "mtMet", "mtVer", "mtInv", 
            "NQ.bird", "NQ.insect", "NQ.mammal", "NQ.pfam", "NQ.plant", "NQ.yeast",
            "Poisson", "PMB", 
            "Q.bird", "Q.insect", "Q.mammal", "Q.pfam", "Q.plant", "Q.yeast",
            "rtREV", "VT", "WAG"
        ]
        
        # 创建模型映射字典和显示列表
        self.model_map = {}  # 用于存储别名到主名称的映射
        model_display_items = ['auto']
        
        # 处理DNA模型
        for model_entry in dna_models:
            if isinstance(model_entry, list):
                main_model = model_entry[0]
                model_display_items.append(main_model)  # 主名称作为显示项
                
                # 为所有别名建立映射关系
                for alias in model_entry[1:]:
                    self.model_map[alias] = main_model
                    model_display_items.append(f"{main_model} ({alias})")  # 添加带别名的显示项
            else:
                model_display_items.append(model_entry)
        
        # 处理蛋白质模型（暂时没有别名）
        for model in protein_models:
            model_display_items.append(model)
        
        # 添加所有模型到组合框
        self.model_combo.addItems(model_display_items)
        self.model_combo.currentTextChanged.connect(self.on_model_changed)
        params_form_layout.addRow("Substitution Model:", self.model_combo)
        
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
        params_form_layout.addRow("Gamma Distribution:", gamma_layout)
        
        # Invariable Sites参数
        self.invar_checkbox = QCheckBox("+I")
        params_form_layout.addRow("Invariable Sites:", self.invar_checkbox)
        
        # FreeRate参数
        self.freerate_checkbox = QCheckBox("+R")
        params_form_layout.addRow("Free Rate Model:", self.freerate_checkbox)

        # state freq selection
        self.state_freq_combo = QComboBox()
        self.state_freq_combo.addItems(["Estimated", "Empirical (+F)", "ML-optimized (+FO)", "Equal (+FQ)"])
        params_form_layout.addRow("State Freq.:", self.state_freq_combo)
        
        # 线程数
        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 32)
        self.threads_spinbox.setValue(1)
        self.threads_spinbox.setToolTip("Number of threads for calculation")
        params_form_layout.addRow("Threads:", self.threads_spinbox)
        
        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets

    def on_model_changed(self):
        if self.model_combo.currentText() == "auto":
            # disable +I +G +F parameters
            self.gamma_checkbox.setEnabled(False)
            self.invar_checkbox.setEnabled(False)
            self.state_freq_combo.setEnabled(False)
            self.freerate_checkbox.setEnabled(False)
            self.gamma_checkbox.setChecked(False)
            self.invar_checkbox.setChecked(False)
            self.freerate_checkbox.setChecked(False)
            self.state_freq_combo.setCurrentText('Estimated')
        else:
            # enable +I +G +F parameters
            self.gamma_checkbox.setEnabled(True)
            self.invar_checkbox.setEnabled(True)
            self.state_freq_combo.setEnabled(True)
            self.freerate_checkbox.setEnabled(True)
    
    def reset_to_defaults(self):
        """重置参数为默认值"""
        self.seq_type_combo.setCurrentText("AUTO")
        self.model_combo.setCurrentIndex(0)  # First model (JC)
        self.gamma_checkbox.setChecked(False)
        self.gamma_spinbox.setValue(4)
        self.invar_checkbox.setChecked(False)
        self.freerate_checkbox.setChecked(False)
        self.threads_spinbox.setValue(1)
    
    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.fas *.fasta *.fa *.fna *.phy);;All files (*)"
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
    
    def setup_output_tab(self):
        """设置输出标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 创建表格控件用于显示距离矩阵
        self.distance_table = QTableWidget()
        self.distance_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.distance_table.horizontalHeader().setStretchLastSection(True)
        self.distance_table.verticalHeader().setStretchLastSection(True)
        self.distance_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.distance_table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        # self.distance_table.setStyleSheet("""
        #     QTableWidget {
        #         gridline-color: #ddd;
        #         border: 1px solid #ccc;
        #     }
        #     QHeaderView::section {
        #         background-color: #f8f9fa;
        #         padding: 4px;
        #         border: 1px solid #ddd;
        #     }
        #     QTableWidget::item {
        #         padding: 4px;
        #         border: 1px solid #eee;
        #     }
        # """)
        
        # 添加到布局
        layout.addWidget(QLabel("ML Distance Matrix Preview:"))
        layout.addWidget(self.distance_table)
    
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)
    
    def reset_to_defaults(self):
        """重置参数为默认值"""
        self.seq_type_combo.setCurrentText("AUTO")
        self.model_combo.setCurrentText("AUTO (automatic model selection)")
        self.custom_model_edit.clear()
        self.custom_model_edit.setVisible(False)
        self.threads_spinbox.setValue(1)
    
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
    
    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 序列类型
        seq_type = self.seq_type_combo.currentText()
        if seq_type != "AUTO":
            params.extend(["-st", seq_type.lower()])
        
        # 模型
        model_text = self.model_combo.currentText()
        model = model_text
        
        # 如果选择了带别名的模型显示项，则提取主模型名称
        if " (" in model_text and ")" in model_text:
            model = model_text.split(" (")[0]
        # 如果是别名，则转换为主名称
        elif model_text in self.model_map:
            model = self.model_map[model_text]
        
        # 添加模型扩展参数
        if self.gamma_checkbox.isChecked():
            model += f"+G{self.gamma_spinbox.value()}"
        
        if self.invar_checkbox.isChecked():
            model += "+I"
            
        if self.freerate_checkbox.isChecked():
            model += "+R"

        stfreq = self.state_freq_combo.currentText()
        model += {"Estimated": "", "Empirical (+F)": "+F", "ML-optimized (+FO)": "+FO", "Equal (+FQ)": "+FQ"}[stfreq]
        
        if model != 'auto':
            params.extend(["-m", model])
        
        params.extend(["-redo"])
        
        # 线程数
        threads = self.threads_spinbox.value()
        if threads > 1:
            params.extend(["-nt", str(threads)])
        
        return params
    
    def run_analysis(self):
        """运行分析"""
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "IQ-TREE executable file not found!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting ML distance calculation...", "info")
        
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
        
        # 在单独的线程中运行IQ-TREE
        self.analysis_thread = MLDistanceThread(
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
        
        # 保存输出文件
        self.current_output_files = output_files
        
        # 显示结果
        self.display_results(output_files)
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        # 添加控制台消息
        self.add_console_message(f"ML distance calculation completed successfully! Found {len(output_files)} result file(s)", "info")
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from == "YR_MPEA":
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        
        QMessageBox.information(self, "Completed", "ML distance calculation completed!")
    
    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 添加控制台消息
        self.add_console_message(f"ML distance calculation failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"ML distance calculation failed: {error_message}")
    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "ML distance calculation has been aborted.")
    
    def display_results(self, output_files):
        """显示结果为表格形式"""
        if not output_files:
            return
            
        try:
            with open(output_files[0], 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 解析距离矩阵
            lines = content.strip().split('\n')
            if not lines:
                return
            
            lines = lines[1:]
                
            # 获取序列名称（第一列）
            sequence_names = []
            for line in lines:
                parts = line.split()
                if len(parts) > 0:
                    sequence_names.append(parts[0])
            
            # 构建距离矩阵
            matrix_data = []
            for i, line in enumerate(lines):
                parts = line.split()
                row_data = []
                for j in range(1, len(parts)):
                    try:
                        row_data.append(float(parts[j]))
                    except ValueError:
                        row_data.append(0.0)
                matrix_data.append(row_data)
            
            # 设置表格大小
            num_sequences = len(sequence_names)
            self.distance_table.setRowCount(num_sequences)
            self.distance_table.setColumnCount(num_sequences)
            
            # 设置表头
            for i in range(num_sequences):
                self.distance_table.setHorizontalHeaderItem(i, QTableWidgetItem(str(i+1)))
                self.distance_table.setVerticalHeaderItem(i, QTableWidgetItem(sequence_names[i]))
            
            # 填充数据
            for i in range(num_sequences):
                for j in range(num_sequences):
                    if i == j:
                        self.distance_table.setItem(i, j, QTableWidgetItem("0.0000"))
                    else:
                        value = matrix_data[i][j]
                        item = QTableWidgetItem(f"{value:.4f}")
                        item.setTextAlignment(Qt.AlignRight)
                        self.distance_table.setItem(i, j, item)
                        
        except Exception as e:
            self.add_console_message(f"Error displaying results: {str(e)}", "error")
    
    def import_distance_matrices_to_platform(self, matrices):
        """将距离矩阵导入到平台的工作区"""
        # 发送信号将数据导入到平台
        self.export_distance_result_signal.emit({"type": "distance_matrix", "data": matrices})
    
    def import_to_platform(self):
        """将结果导入到当前平台"""
        if not hasattr(self, 'current_output_files') or not self.current_output_files:
            QMessageBox.warning(self, "Warning", "No ML distance results to import.")
            return
            
        try:
            # 读取ML距离矩阵文件内容
            distance_matrices = []
            for output_file in self.current_output_files:
                with open(output_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                    distance_matrices.append({
                        'filename': os.path.basename(output_file),
                        'content': content
                    })
            
            if not distance_matrices:
                QMessageBox.warning(self, "Warning", "No distance matrices found in results.")
                return
                
            # 发送信号将数据导入到平台
            self.import_distance_matrices_to_platform(distance_matrices)
            
            # 显示成功消息
            QMessageBox.information(self, "Success", f"Successfully imported {len(distance_matrices)} distance matrix(es) to the platform.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import distance matrices: {str(e)}")
    
    def import_alignment_to_platform(self, sequences):
        """将距离矩阵导入到平台的工作区"""
        # 发送信号将数据导入到平台
        self.import_alignment_signal.emit(sequences)
    
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
        elif isinstance(import_data, dict) and 'model' in import_data:
            # 处理模型导入
            imported_model = import_data['model']
            
            # 检查导入的模型是否为别名，并转换为主模型名
            model_entries = imported_model.split('+')
            if hasattr(self, 'model_map'):
                # 提取模型名（去除可能的参数）
                base_model = model_entries[0]
                if base_model in self.model_map:
                    # 替换为完整模型名
                    imported_model = imported_model.replace(base_model, self.model_map[base_model], 1)
            
            # 设置模型选择框的值
            self.model_combo.setCurrentText(imported_model)
        else:
            self.import_file = None
            self.imported_files = []


# 插件入口点
class MLDistancePluginEntry:
    """ML Distance插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return MLDistancePlugin(import_from=import_from, import_data=import_data)