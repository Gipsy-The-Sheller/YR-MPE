# uncorrected_distance_plugin.py
#
# Copyright (c) 2026 YRTools Team
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

import os
import sys
import tempfile
from typing import List, Optional, Set, Tuple

# 添加项目根目录到路径
plugin_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(plugin_dir))
sys.path.insert(0, project_root)

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QTextEdit,
                             QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView,
                             QDialog, QListWidget, QListWidgetItem, QFrame)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QTextCursor, QIcon
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..templates.base_plugin_ui import BasePlugin

# 相对导入距离计算工具
distance_utils_path = os.path.join(project_root, 'YR_MPE', 'plugins', 'components', 'methods')
sys.path.insert(0, distance_utils_path)
try:
    import distance_utils
    from distance_utils import (
        calculate_distance_matrix, 
        validate_sequences, 
        format_distance_matrix,
        get_substitution_filter,
        IUPAC_DNA_MAP,
        DNA_SUBSTITUTION_TYPES
    )
except ImportError:
    # 备用导入方式
    from .components.methods.distance_utils import (
        calculate_distance_matrix, 
        validate_sequences, 
        format_distance_matrix,
        get_substitution_filter,
        IUPAC_DNA_MAP,
        DNA_SUBSTITUTION_TYPES
    )


class AdvancedRulesDialog(QDialog):
    """高级替代规则配置对话框"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Advanced DNA Substitution Rules")
        self.setModal(True)
        self.resize(400, 300)
        self.selected_substitutions = set()
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        
        # 说明标签
        info_label = QLabel("Select DNA substitution types to include in distance calculation:")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # 替代类型列表
        self.sub_list = QListWidget()
        self.sub_list.setSelectionMode(QListWidget.MultiSelection)
        
        # 添加所有6种替代类型
        substitutions = [
            ("A ↔ G", "Transition (purine-purine)"),
            ("T ↔ C", "Transition (pyrimidine-pyrimidine)"), 
            ("A ↔ T", "Transversion"),
            ("A ↔ C", "Transversion"),
            ("G ↔ T", "Transversion"),
            ("G ↔ C", "Transversion")
        ]
        
        for i, (sub, desc) in enumerate(substitutions):
            item = QListWidgetItem(f"{sub} - {desc}")
            item.setData(Qt.UserRole, i)  # 存储索引
            self.sub_list.addItem(item)
        
        layout.addWidget(self.sub_list)
        
        # 按钮布局
        button_layout = QHBoxLayout()
        
        # 预设按钮
        all_button = QPushButton("Select All")
        all_button.clicked.connect(self.select_all)
        none_button = QPushButton("Select None") 
        none_button.clicked.connect(self.select_none)
        defaults_button = QPushButton("Defaults (All)")
        defaults_button.clicked.connect(self.select_defaults)
        
        button_layout.addWidget(all_button)
        button_layout.addWidget(none_button)
        button_layout.addWidget(defaults_button)
        button_layout.addStretch()
        
        # 确定/取消按钮
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        
        button_layout.addWidget(ok_button)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
        self.setLayout(layout)
        
        # 设置默认选择（全部选中）
        self.select_defaults()
    
    def select_all(self):
        """选择所有替代类型"""
        for i in range(self.sub_list.count()):
            self.sub_list.item(i).setSelected(True)
    
    def select_none(self):
        """取消选择所有替代类型"""
        for i in range(self.sub_list.count()):
            self.sub_list.item(i).setSelected(False)
    
    def select_defaults(self):
        """选择默认替代类型（全部）"""
        self.select_all()
    
    def get_selected_substitutions(self) -> Set[Tuple[str, str]]:
        """获取选中的替代类型"""
        selected = set()
        for i in range(self.sub_list.count()):
            item = self.sub_list.item(i)
            if item.isSelected():
                idx = item.data(Qt.UserRole)
                # 映射到实际的碱基对
                sub_map = [
                    ('A', 'G'), ('G', 'A'),  # A ↔ G
                    ('T', 'C'), ('C', 'T'),  # T ↔ C  
                    ('A', 'T'), ('T', 'A'),  # A ↔ T
                    ('A', 'C'), ('C', 'A'),  # A ↔ C
                    ('G', 'T'), ('T', 'G'),  # G ↔ T
                    ('G', 'C'), ('C', 'G')   # G ↔ C
                ]
                # 每个用户选择对应两个方向
                base_idx = idx * 2
                selected.add(sub_map[base_idx])
                selected.add(sub_map[base_idx + 1])
        
        return selected


class UncorrectedDistancePlugin(BasePlugin):
    """未校正遗传距离(p-distance)计算插件"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_distance_result_signal = pyqtSignal(dict)  # 导出距离矩阵信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化未校正距离插件"""
        super().__init__(import_from, import_data)
        
        # 初始化用于存储计算结果的属性
        self.calculated_distance_results = []
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
    
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "Uncorrected Distance (p-distance)"
        self.tool_name = "Built-in Python Implementation"
        self.citation = []  # 空的引用列表
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "NEXUS": ["nex", "nexus"]}
        self.output_types = {"Distance Matrix": ".dist"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    
    def config(self):
        """
        检查插件配置
        对于纯Python实现的插件，始终返回True
        """
        # 这是一个纯Python实现的插件，不需要外部工具依赖
        return True
    
    def show_config_guide(self):
        """
        显示配置指南
        对于纯Python插件，这个方法不应该被调用
        """
        pass
    
    def setup_ui(self):
        """设置用户界面"""
        # 创建标签页
        self.tab_widget = QTabWidget()
        self.layout().addWidget(self.tab_widget)
        
        # 输入/参数标签页
        self.input_tab = QWidget()
        self.setup_input_tab()
        self.tab_widget.addTab(self.input_tab, "Input & Parameters")
        
        # 输出标签页
        self.output_tab = QWidget()
        self.setup_output_tab()
        self.tab_widget.addTab(self.output_tab, "Output")
        
        # 控制面板
        self.setup_control_panel()
        
        # 设置默认值
        self.reset_to_defaults()
    
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
        
        # 参数组
        params_group = QGroupBox("Distance Calculation Parameters")
        params_form_layout = QFormLayout()
        params_group.setLayout(params_form_layout)
        layout.addWidget(params_group)
        
        # 序列类型
        self.seq_type_combo = QComboBox()
        self.seq_type_combo.addItems(["AUTO", "DNA", "Protein"])
        params_form_layout.addRow("Sequence Type:", self.seq_type_combo)
        
        # Gap处理方式
        self.gap_treatment_combo = QComboBox()
        self.gap_treatment_combo.addItems(["Pairwise Deletion", "Complete Deletion"])
        params_form_layout.addRow("Gap Treatment:", self.gap_treatment_combo)
        
        # DNA替代规则
        rule_layout = QHBoxLayout()
        self.dna_rule_combo = QComboBox()
        self.dna_rule_combo.addItems(["All", "Transition", "Transversion", "Advanced Rules"])
        self.advanced_rules_btn = QPushButton("Advanced...")
        self.advanced_rules_btn.clicked.connect(self.show_advanced_rules)
        rule_layout.addWidget(self.dna_rule_combo)
        rule_layout.addWidget(self.advanced_rules_btn)
        params_form_layout.addRow("DNA Substitution Rules:", rule_layout)
        
        layout.addStretch()
        
        if not hasattr(self, 'imported_files'):
            self.imported_files = []
        if not hasattr(self, 'file_tags'):
            self.file_tags = []
    
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
        
        # 添加到布局
        layout.addWidget(QLabel("Distance Matrix Preview:"))
        layout.addWidget(self.distance_table)
    
    def setup_control_panel(self):
        """设置控制面板"""
        # 先调用父类方法设置基本控件
        super().setup_control_panel()
        
        # 然后添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        self.control_layout.addWidget(self.import_to_platform_btn)
    
    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.fas *.fna *.fa *.fasta *.nex *.nexus);;All files (*)"
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select alignment files", "", file_filter)
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
    
    def add_file_tag(self, file_path):
        """添加文件标签"""
        if file_path not in self.imported_files:
            self.imported_files.append(file_path)
            
            # Create file tag widget (使用IQ-TREE plugin的灰底风格)
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
            
            # 文件名标签
            name_label = QLabel(os.path.basename(file_path))
            name_label.setStyleSheet("color: #495057;")
            tag_layout.addWidget(name_label)
            
            # 关闭按钮
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
            
            self.file_tags_layout.addWidget(tag_widget)
            self.file_tags.append((file_path, tag_widget))
            self.file_tags_container.setVisible(True)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def remove_file_tag(self, tag_widget, file_path):
        """移除文件标签"""
        # 从列表中移除
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        
        # 从UI中移除
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # 更新显示
        if not self.imported_files:
            self.file_tags_container.setVisible(False)
            self.file_path_edit.clear()
        else:
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def on_text_changed(self):
        """文本输入变化处理"""
        if self.sequence_text.toPlainText().strip():
            # 如果有文本输入，清空文件选择
            self.imported_files.clear()
            for tag in self.file_tags:
                self.file_tags_layout.removeWidget(tag)
                tag.deleteLater()
            self.file_tags.clear()
            self.file_tags_container.setVisible(False)
        elif not self.imported_files:
            # 如果既没有文本也没有文件，恢复初始状态
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
    
    def show_advanced_rules(self):
        """显示高级替代规则对话框"""
        dialog = AdvancedRulesDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            self.selected_substitutions = dialog.get_selected_substitutions()
            # 更新显示
            if self.dna_rule_combo.currentText() != "Advanced Rules":
                self.dna_rule_combo.setCurrentText("Advanced Rules")
    
    def reset_to_defaults(self):
        """重置参数为默认值"""
        self.seq_type_combo.setCurrentText("AUTO")
        self.gap_treatment_combo.setCurrentText("Pairwise Deletion")
        self.dna_rule_combo.setCurrentText("All")
        self.advanced_rules_btn.setEnabled(True)
    
    def prepare_input_files(self):
        """准备输入文件"""
        try:
            input_files = []
            if self.imported_files:
                for file_path in self.imported_files:
                    if os.path.exists(file_path):
                        input_files.append(file_path)
                return input_files
            elif self.sequence_text.toPlainText().strip():
                # 从文本创建临时文件
                temp_file = self.create_temp_file(suffix='.fas')
                with open(temp_file, 'w') as f:
                    f.write(self.sequence_text.toPlainText())
                self.temp_files.append(temp_file)
                return [temp_file]
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return None
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
    
    def get_parameters(self):
        """获取参数设置"""
        params = {}
        
        # 序列类型
        seq_type = self.seq_type_combo.currentText()
        params['seq_type'] = seq_type
        
        # Gap处理方式
        gap_treatment = self.gap_treatment_combo.currentText()
        params['gap_treatment'] = 'pairwise' if 'Pairwise' in gap_treatment else 'complete'
        
        # DNA替代规则
        rule_type = self.dna_rule_combo.currentText()
        if rule_type == "Advanced Rules" and hasattr(self, 'selected_substitutions'):
            params['substitution_filter'] = self.selected_substitutions
        else:
            params['substitution_filter'] = get_substitution_filter(rule_type.lower())
        
        return params
    
    def run_analysis(self):
        """运行距离分析"""
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip():
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting uncorrected distance calculation...", "info")
        
        # 准备输入文件
        input_files = self.prepare_input_files()
        if not input_files:
            return
        
        # 获取参数
        params = self.get_parameters()
        
        # 设置运行状态
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 未知进度
        self.import_to_platform_btn.setVisible(False)  # 分析期间隐藏导入按钮
        
        try:
            # 直接在内存中存储计算结果，而不是创建临时文件
            self.calculated_distance_results = []
            
            for i, input_file in enumerate(input_files):
                if not self.is_running:
                    break
                
                self.progress_bar.setFormat(f"Processing file {i+1}/{len(input_files)}...")
                self.add_console_message(f"Processing file {i+1}/{len(input_files)}: {os.path.basename(input_file)}", "info")
                
                # 读取序列
                sequences = list(SeqIO.parse(input_file, "fasta"))
                
                # 验证序列
                is_valid, msg = validate_sequences(sequences)
                if not is_valid:
                    self.add_console_message(f"Validation error: {msg}", "error")
                    continue
                
                # 计算距离矩阵
                distance_matrix = calculate_distance_matrix(
                    sequences, 
                    gap_treatment=params['gap_treatment'],
                    substitution_filter=params.get('substitution_filter')
                )
                
                # 直接保存计算结果到内存
                sequence_names = [seq.id for seq in sequences]
                formatted_matrix = format_distance_matrix(distance_matrix, sequence_names)
                
                result_data = {
                    'filename': f"{os.path.splitext(os.path.basename(input_file))[0]}_uncorrected_distance.dist",
                    'content': formatted_matrix,
                    'matrix': distance_matrix,
                    'sequence_names': sequence_names
                }
                
                self.calculated_distance_results.append(result_data)
                self.add_console_message(f"Calculated distance matrix for {len(sequences)} sequences", "info")
                
                # 显示结果预览
                self.display_results_preview(distance_matrix, sequence_names)
            
            # 分析完成，重置UI状态并显示导入按钮
            self.analysis_finished([], [])
            
        except Exception as e:
            self.analysis_error(str(e))
    
    def display_results_preview(self, distance_matrix, sequence_names):
        """显示结果预览"""
        try:
            n = len(sequence_names)
            self.distance_table.setRowCount(n)
            self.distance_table.setColumnCount(n)
            
            # 设置表头
            for i in range(n):
                self.distance_table.setHorizontalHeaderItem(i, QTableWidgetItem(str(i+1)))
                self.distance_table.setVerticalHeaderItem(i, QTableWidgetItem(sequence_names[i]))
            
            # 填充数据
            for i in range(n):
                for j in range(n):
                    if i == j:
                        self.distance_table.setItem(i, j, QTableWidgetItem("0.0000"))
                    else:
                        value = distance_matrix[i][j]
                        item = QTableWidgetItem(f"{value:.4f}")
                        item.setTextAlignment(Qt.AlignRight)
                        self.distance_table.setItem(i, j, item)
                        
        except Exception as e:
            self.add_console_message(f"Error displaying results: {str(e)}", "error")
    
    def analysis_finished(self, output_files, reports):
        """分析完成回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.progress_bar.setFormat("")
        
        # 显示完成消息
        QMessageBox.information(self, "Completed", "Distance calculation completed successfully!")
        
        # 分析完成后显示导入按钮（如果有计算结果）
        if hasattr(self, 'calculated_distance_results') and self.calculated_distance_results:
            self.import_to_platform_btn.setVisible(True)
            self.add_console_message("Analysis completed. Click 'Import to Current Platform' to export results.", "info")
    
    def analysis_error(self, error_msg):
        """分析错误回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.progress_bar.setFormat("")
        self.import_to_platform_btn.setVisible(False)  # 错误时隐藏导入按钮
        
        QMessageBox.critical(self, "Error", f"Analysis failed: {error_msg}")
    
    def prepare_platform_import(self):
        """准备将结果导入平台的数据"""
        
        if not hasattr(self, 'calculated_distance_results') or not self.calculated_distance_results:
            QMessageBox.warning(self, "Warning", "No distance results to import.")
            return False
        
        try:
            # 直接使用内存中的计算结果
            distance_matrices = []
            for result in self.calculated_distance_results:
                distance_matrices.append({
                    'filename': result['filename'],
                    'content': result['content']
                })
            
            if not distance_matrices:
                QMessageBox.warning(self, "Warning", "No distance matrices found in results.")
                return False
            
            # 准备要发送的数据
            export_data = {
                'type': 'distance_matrix',
                'data': distance_matrices,
                'source_plugin': 'Uncorrected Distance'
            }
            
            # 发送信号将数据导入到平台
            self.export_distance_result_signal.emit(export_data)
            
            # 显示成功消息
            QMessageBox.information(self, "Success", f"Successfully imported {len(distance_matrices)} distance matrix(es) to the platform.")
            return True
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import distance matrices: {str(e)}")
            return False

    def import_to_platform(self):
        """将结果导入到当前平台"""
        success = self.prepare_platform_import()
        if success:
            # 强制更新UI
            self.import_to_platform_btn.update()
            self.add_console_message("Results imported to platform successfully.", "info")
        return success
    
    def handle_import_data(self, import_data):
        """处理从YR-MPEA导入的数据"""
        if isinstance(import_data, dict) and 'type' in import_data and import_data['type'] == 'alignment':
            # 处理对齐序列数据
            sequences = import_data['data']
            if isinstance(sequences, list) and len(sequences) > 0:
                # 创建临时文件来存储导入的序列数据
                temp_file = self.create_temp_file(suffix='.fas')
                with open(temp_file, 'w') as f:
                    for seq in sequences:
                        f.write(f">{seq.id}\n{seq.seq}\n")
                self.temp_files.append(temp_file)
                self.import_file = temp_file
                self.imported_files = [temp_file]
                
                # 更新UI显示导入的文件
                if hasattr(self, 'file_path_edit') and self.file_path_edit:
                    self.file_path_edit.setText(temp_file)
                    
                # 添加控制台消息
                self.add_console_message(f"Imported alignment data with {len(sequences)} sequences from YR-MPEA", "info")
        elif isinstance(import_data, list):
            # 兼容旧版本的导入数据格式
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
                
            # 添加控制台消息
            self.add_console_message(f"Imported alignment data with {len(import_data)} sequences from YR-MPEA", "info")
        else:
            self.import_file = None
            self.imported_files = []
            self.add_console_message("Invalid import data format", "error")
            
        # 确保 calculated_distance_results 属性存在
        if not hasattr(self, 'calculated_distance_results'):
            self.calculated_distance_results = []
        
        # 根据导入来源设置导入按钮的可见性
        if self.import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"]:
            self.import_to_platform_btn.setVisible(False)
        else:
            self.import_to_platform_btn.setVisible(False)


# 插件入口点
class UncorrectedDistancePluginEntry:
    """未校正距离插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return UncorrectedDistancePlugin(import_from=import_from, import_data=import_data)