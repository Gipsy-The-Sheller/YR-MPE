# pdguide_ui.py
#
# Copyright (c) 2026 Zhi-Jie Xu
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

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, 
                            QLabel, QLineEdit, QPushButton, QComboBox, 
                            QTableWidget, QTableWidgetItem, QMessageBox,
                            QButtonGroup, QRadioButton, QDialog, QAbstractItemView,
                            QHeaderView, QGroupBox, QFileDialog, QSpinBox, QSizePolicy)
from PyQt5.QtCore import Qt
import os
import json
from typing import List, Dict
from .md_mrca import MDMRCA


class PdGuideUI(QWidget):
    """PD-Guide用户界面组件"""
    
    def __init__(self):
        super().__init__()
        self.calibration_points = []
        self.calibration_file = None
        self.all_param_widgets = {}  # 存储所有参数控件
        self.init_ui()
    
    def init_ui(self):
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # Dating methods
        dating_methods_layout = QHBoxLayout()
        
        self.methods_combo = QComboBox()
        self.methods_combo.addItems(['LSD2', 'PhyloBayes', 'MrBayes', 'MCMCTree'])
        self.methods_combo.currentTextChanged.connect(self.on_method_changed)

        dating_methods_label = QLabel('Dating methods:')
        dating_methods_label.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)

        self.methods_combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        dating_methods_layout.addWidget(dating_methods_label)
        dating_methods_layout.addWidget(self.methods_combo)

        main_layout.addLayout(dating_methods_layout)

        # 预创建所有参数控件
        self.create_all_parameter_widgets()
        
        # 初始显示LSD2参数
        self.show_parameters_for_method('LSD2')
    
    def create_all_parameter_widgets(self):
        """预创建所有参数控件"""
        main_layout = self.layout()
        
        # 输入文件管理
        input_files_group = QGroupBox("Input Files")
        input_files_layout = QVBoxLayout()
        
        # 树文件输入
        tree_file_layout = QHBoxLayout()
        tree_file_layout.addWidget(QLabel("Tree file (Newick):"))
        self.tree_file_path = QLineEdit()
        self.tree_file_path.setPlaceholderText("Select a tree file...")
        self.tree_file_path.setReadOnly(True)
        tree_file_layout.addWidget(self.tree_file_path)
        
        self.load_tree_button = QPushButton("Browse")
        self.load_tree_button.clicked.connect(self.load_tree_file)
        tree_file_layout.addWidget(self.load_tree_button)
        input_files_layout.addLayout(tree_file_layout)
        
        # 序列文件输入（可选，用于序列长度计算）
        seq_file_layout = QHBoxLayout()
        seq_file_layout.addWidget(QLabel("Sequence file (optional):"))
        self.seq_file_path = QLineEdit()
        self.seq_file_path.setPlaceholderText("Select a sequence file for sequence length...")
        self.seq_file_path.setReadOnly(True)
        seq_file_layout.addWidget(self.seq_file_path)
        
        self.load_seq_button = QPushButton("Browse")
        self.load_seq_button.clicked.connect(self.load_sequence_file)
        seq_file_layout.addWidget(self.load_seq_button)
        input_files_layout.addLayout(seq_file_layout)
        
        # 序列长度输入
        seq_length_layout = QHBoxLayout()
        seq_length_layout.addWidget(QLabel("Sequence length:"))
        
        self.seq_length_combo = QComboBox()
        self.seq_length_combo.addItems(['Disabled', 'Manual', 'From sequence file'])
        self.seq_length_combo.currentIndexChanged.connect(self.on_seq_length_mode_changed)
        seq_length_layout.addWidget(self.seq_length_combo)
        
        self.sequence_length_spin = QSpinBox()
        self.sequence_length_spin.setRange(0, 1000000)
        self.sequence_length_spin.setValue(0)
        self.sequence_length_spin.setEnabled(False)  # 初始禁用
        seq_length_layout.addWidget(self.sequence_length_spin)
        
        seq_length_layout.addStretch()
        input_files_layout.addLayout(seq_length_layout)
        
        input_files_group.setLayout(input_files_layout)
        main_layout.addWidget(input_files_group)
        
        # Infer tree with seqs 按钮
        self.infer_tree_button = QPushButton("Infer tree with seqs")
        self.infer_tree_button.clicked.connect(self.launch_iqtree_analysis)
        main_layout.addWidget(self.infer_tree_button)
        
        # Rooting options
        root_layout = QHBoxLayout()
        root_layout.addWidget(QLabel("Root the tree:"))
        
        self.choose_root_button = QPushButton("choose")
        root_layout.addWidget(self.choose_root_button)
        
        root_layout.addStretch()
        main_layout.addLayout(root_layout)
        
        # Rooting option radio buttons
        rooting_option_layout = QHBoxLayout()
        rooting_option_layout.addWidget(QLabel("Rooting option:"))
        
        self.rooting_button_group = QButtonGroup(self)
        self.use_chosen_root_radio = QRadioButton("use chosen root")
        self.estimate_root_radio = QRadioButton("estimate the root")
        self.estimate_root_radio.setChecked(True)  # 默认选中
        
        self.rooting_button_group.addButton(self.use_chosen_root_radio)
        self.rooting_button_group.addButton(self.estimate_root_radio)
        
        rooting_option_layout.addWidget(self.use_chosen_root_radio)
        rooting_option_layout.addWidget(self.estimate_root_radio)
        rooting_option_layout.addStretch()
        main_layout.addLayout(rooting_option_layout)
        
        # Calibration points
        calibration_layout = QHBoxLayout()
        calibration_layout.addWidget(QLabel("Calibration points:"))
        
        self.add_calibration_button = QPushButton("+")
        self.add_calibration_button.clicked.connect(self.add_calibration_point)
        calibration_layout.addWidget(self.add_calibration_button)
        calibration_layout.addStretch()
        main_layout.addLayout(calibration_layout)
        
        # Calibration table
        self.calibration_table = QTableWidget(0, 4)
        self.calibration_table.setHorizontalHeaderLabels(['Name', 'Taxa', 'Type', 'Options'])
        self.calibration_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.calibration_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        main_layout.addWidget(self.calibration_table)
        
        # Warning label for calibration types
        self.calibration_warning_label = QLabel(
            "* WARNING: Calibrations except for uniform / point / upper / lower "
            "will be automatically transferred to uniform(lower 95HPD, upper 95HPD)"
        )
        self.calibration_warning_label.setStyleSheet("color: red;")
        self.calibration_warning_label.setVisible(False)  # 初始隐藏
        self.calibration_warning_label.setWordWrap(True)
        main_layout.addWidget(self.calibration_warning_label)
        
        # 其他方法的参数组（简化示例）
        self.phylobayes_group = QGroupBox("PhyloBayes Parameters")
        phylobayes_layout = QVBoxLayout()
        phylobayes_layout.addWidget(QLabel("PhyloBayes parameters will be implemented here"))
        self.phylobayes_group.setLayout(phylobayes_layout)
        main_layout.addWidget(self.phylobayes_group)
        
        self.mrbayes_group = QGroupBox("MrBayes Parameters")
        mrbayes_layout = QVBoxLayout()
        mrbayes_layout.addWidget(QLabel("MrBayes parameters will be implemented here"))
        self.mrbayes_group.setLayout(mrbayes_layout)
        main_layout.addWidget(self.mrbayes_group)
        
        self.mcmctree_group = QGroupBox("MCMCTree Parameters")
        mcmctree_layout = QVBoxLayout()
        mcmctree_layout.addWidget(QLabel("MCMCTree parameters will be implemented here"))
        self.mcmctree_group.setLayout(mcmctree_layout)
        main_layout.addWidget(self.mcmctree_group)
        
        # 存储控件引用（LSD2内容现在直接在主布局中，不需要group）
        self.all_param_widgets = {
            'PhyloBayes': self.phylobayes_group,
            'MrBayes': self.mrbayes_group,
            'MCMCTree': self.mcmctree_group
        }
        
        # LSD2相关的控件（需要根据方法选择显示/隐藏）
        self.lsd2_widgets = [
            input_files_group,
            self.infer_tree_button,
            self.choose_root_button,
            self.use_chosen_root_radio,
            self.estimate_root_radio,
            self.add_calibration_button,
            self.calibration_table,
            self.calibration_warning_label
        ]
    
    def on_method_changed(self, method):
        """当方法选择改变时的回调"""
        self.show_parameters_for_method(method)
    
    def show_parameters_for_method(self, method):
        """显示指定方法的参数控件"""
        # 显示/隐藏LSD2相关控件
        is_lsd2 = (method == 'LSD2')
        for widget in self.lsd2_widgets:
            if widget:
                widget.setVisible(is_lsd2)
        
        # 显示/隐藏其他方法的参数组
        for method_name, widget in self.all_param_widgets.items():
            if method_name == method:
                widget.show()
            else:
                widget.hide()
    
    def setup_input_files(self):
        # 这个方法现在不需要了，因为我们在init_ui中已经处理了
        pass
    
    def setup_input_files_lsd2(self):
        # 这个方法现在不需要了
        pass
    
    def setup_input_files_phylobayes(self):
        # 这个方法现在不需要了
        pass
    
    def setup_input_files_mrbayes(self):
        # 这个方法现在不需要了
        pass
    
    def setup_input_files_mcmctree(self):
        # 这个方法现在不需要了
        pass
    
    def clear_layout(self):
        # 这个方法现在不需要了
        pass
    
    def launch_iqtree_analysis(self):
        """启动IQ-TREE分析"""
        try:
            from YR_MPE.plugins.iqtree_plugin import IQTreePluginEntry
            plugin_entry = IQTreePluginEntry()
            plugin_entry.run()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to launch IQ-TREE plugin: {str(e)}")
    
    def add_calibration_point(self):
        """添加校准点"""
        row = self.calibration_table.rowCount()
        self.calibration_table.insertRow(row)
        
        # Name列 - 可编辑
        name_item = QTableWidgetItem("")
        self.calibration_table.setItem(row, 0, name_item)
        
        # Taxa列 - 禁用编辑
        taxa_item = QTableWidgetItem("0 taxa")
        taxa_item.setFlags(taxa_item.flags() & ~Qt.ItemIsEditable)
        self.calibration_table.setItem(row, 1, taxa_item)
        
        # Type列 - 禁用编辑  
        type_item = QTableWidgetItem("")
        type_item.setFlags(type_item.flags() & ~Qt.ItemIsEditable)
        self.calibration_table.setItem(row, 2, type_item)
        
        # Options列 - 包含Edit和Discard按钮
        options_widget = QWidget()
        options_layout = QHBoxLayout(options_widget)
        options_layout.setContentsMargins(0, 0, 0, 0)
        
        edit_button = QPushButton("Edit")
        discard_button = QPushButton("Discard")
        
        options_layout.addWidget(edit_button)
        options_layout.addWidget(discard_button)
        
        self.calibration_table.setCellWidget(row, 3, options_widget)
        
        # 连接按钮信号 - 使用默认参数捕获当前行
        edit_button.clicked.connect(lambda checked, btn=edit_button: self.edit_calibration_point_from_button(btn))
        discard_button.clicked.connect(lambda checked, btn=discard_button: self.discard_calibration_point_from_button(btn))
        
        # 初始化校准数据存储
        if not hasattr(self, 'calibration_data'):
            self.calibration_data = {}
        self.calibration_data[row] = None  # 标记为未配置
    
    def edit_calibration_point_from_button(self, button):
        """通过Edit按钮编辑校准点"""
        # 找到按钮所在的Options widget
        options_widget = button.parent()
        if not options_widget:
            return
            
        # 找到该widget在表格中的行号
        row = self.calibration_table.indexAt(options_widget.pos()).row()
        if row < 0:
            return
            
        # 获取当前Newick树字符串
        tree_file = self.tree_file_path.text().strip()
        if not tree_file or not os.path.exists(tree_file):
            QMessageBox.warning(self, "Warning", "Please load a tree file first.")
            return
        
        try:
            with open(tree_file, 'r') as f:
                newick_text = f.read().strip()
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to read tree file: {str(e)}")
            return
        
        if not newick_text:
            QMessageBox.warning(self, "Warning", "The tree file is empty.")
            return
            
        # 创建并显示MD_MRCA对话框
        md_mrca_dialog = MDMRCA(parent=self)
        md_mrca_dialog.newick_string = newick_text
        md_mrca_dialog.annotated_newick_str = newick_text
        
        # 加载数据
        md_mrca_dialog.load_newick_data()
        
        # 以模态对话框方式显示
        if md_mrca_dialog.exec_() == QDialog.Accepted:
            # 用户点击了Apply按钮
            self.handle_md_mrca_result(row, md_mrca_dialog)
        
        # 对话框会自动销毁
        md_mrca_dialog.deleteLater()
        
    def handle_md_mrca_result(self, row, md_mrca_dialog):
        """处理MD_MRCA的结果"""
        if md_mrca_dialog:
            selected_taxa = [md_mrca_dialog.selected_taxa_list.item(i).text() 
                           for i in range(md_mrca_dialog.selected_taxa_list.count())]
            taxon_set_name = md_mrca_dialog.taxon_set_name.text().strip()
            
            if taxon_set_name and selected_taxa:
                # 获取校准类型
                cal_type = md_mrca_dialog.tmrca_type_combo.currentText()
                
                # 获取校准参数值
                cal_values = []
                try:
                    if cal_type == 'Point':
                        cal_values = [float(md_mrca_dialog.point_value.text())]
                    elif cal_type == 'Uniform':
                        cal_values = [
                            float(md_mrca_dialog.uniform_lower.text()),
                            float(md_mrca_dialog.uniform_upper.text())
                        ]
                    elif cal_type == 'Upper Boundary':
                        cal_values = [float(md_mrca_dialog.upper_boundary.text())]
                    elif cal_type == 'Lower Boundary':
                        cal_values = [float(md_mrca_dialog.lower_boundary.text())]
                    elif cal_type == 'Normal':
                        cal_values = [
                            float(md_mrca_dialog.normal_mean.text()),
                            float(md_mrca_dialog.normal_std.text())
                        ]
                    elif cal_type == 'Lognormal':
                        cal_values = [
                            float(md_mrca_dialog.lognormal_mean.text()),
                            float(md_mrca_dialog.lognormal_std.text())
                        ]
                except ValueError as e:
                    QMessageBox.warning(self, "Invalid Parameters", 
                                      f"Please check the calibration parameters: {str(e)}")
                    return
                
                # 映射到LSD2使用的类型名称
                type_mapping = {
                    'Point': 'fixed',
                    'Uniform': 'interval',
                    'Upper Boundary': 'upper',
                    'Lower Boundary': 'lower',
                    'Normal': 'uniform',  # LSD2不支持Normal，转换为uniform
                    'Lognormal': 'uniform'  # LSD2不支持Lognormal，转换为uniform
                }
                lsd2_type = type_mapping.get(cal_type, 'uniform')
                
                # 更新表格显示
                self.calibration_table.item(row, 0).setText(taxon_set_name)
                self.calibration_table.item(row, 1).setText(f"{len(selected_taxa)} taxa")
                self.calibration_table.item(row, 2).setText(cal_type)
                
                # 存储完整的校准点数据（使用行号作为键）
                if not hasattr(self, 'calibration_data'):
                    self.calibration_data = {}
                
                self.calibration_data[row] = {
                    'name': taxon_set_name,
                    'set': selected_taxa,
                    'type': lsd2_type,
                    'values': cal_values,
                    'display_type': cal_type  # 保存原始显示类型
                }
                
                # 触发校准类型检查
                self.check_calibration_types()
    
    def load_tree_file(self):
        """加载树文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Tree File",
            "",
            "Newick Files (*.nwk *.tree *.tre *.newick);;All Files (*.*)"
        )
        
        if file_path:
            # 验证文件格式
            if self.validate_tree_file(file_path):
                self.tree_file_path.setText(file_path)
                # 自动检测并设置序列长度为0（让LSD2自动处理）
                self.sequence_length_spin.setValue(0)
                # 提取OTU列表
                self.extract_otu_from_tree(file_path)
            else:
                QMessageBox.warning(self, "Invalid File", "The selected file is not a valid Newick tree file.")
    
    def load_sequence_file(self):
        """加载序列文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Sequence File",
            "",
            "FASTA Files (*.fas *.fna *.fa *.fasta);;All Files (*.*)"
        )
        
        if file_path:
            # 验证文件格式
            if self.validate_sequence_file(file_path):
                self.seq_file_path.setText(file_path)
                # 计算序列长度
                seq_length = self.calculate_sequence_length(file_path)
                if seq_length > 0:
                    self.sequence_length_spin.setValue(seq_length)
            else:
                QMessageBox.warning(self, "Invalid File", "The selected file is not a valid FASTA sequence file.")
    
    def validate_tree_file(self, file_path):
        """验证树文件格式"""
        try:
            with open(file_path, 'r') as f:
                content = f.read().strip()
            
            # 简单的Newick格式验证
            # 必须以分号结尾
            if not content.endswith(';'):
                return False
            
            # 检查括号匹配
            open_parens = content.count('(')
            close_parens = content.count(')')
            if open_parens != close_parens:
                return False
            
            # 检查是否有至少一个分类单元
            # 简单检查：应该有逗号分隔的分类单元
            if ',' not in content:
                return False
            
            return True
            
        except Exception as e:
            print(f"Tree file validation error: {e}")
            return False
    
    def validate_sequence_file(self, file_path):
        """验证序列文件格式"""
        try:
            from Bio import SeqIO
            # 尝试读取第一个序列
            record = next(SeqIO.parse(file_path, 'fasta'))
            if record:
                return True
            return False
        except Exception as e:
            print(f"Sequence file validation error: {e}")
            return False
    
    def calculate_sequence_length(self, file_path):
        """计算序列长度"""
        try:
            from Bio import SeqIO
            # 读取第一个序列的长度
            record = next(SeqIO.parse(file_path, 'fasta'))
            return len(record.seq)
        except Exception as e:
            print(f"Sequence length calculation error: {e}")
            return 0
    
    def extract_otu_from_tree(self, file_path):
        """从树文件中提取OTU列表"""
        try:
            with open(file_path, 'r') as f:
                newick_str = f.read().strip()
            
            # 简单的OTU提取：提取所有不在括号内的单词
            # 这里可以使用ETE3来更准确地提取
            try:
                from ete3 import Tree
                tree = Tree(newick_str)
                self.all_otu = [leaf.name for leaf in tree.get_leaves()]
                print(f"Extracted {len(self.all_otu)} OTUs from tree")
            except ImportError:
                # 如果ETE3不可用，使用简单的正则表达式
                import re
                # 移除分支长度
                clean_newick = re.sub(r':\d+\.?\d*', '', newick_str)
                # 提取分类单元名称
                self.all_otu = re.findall(r'([A-Za-z0-9_\-\.]+)', clean_newick)
                # 过滤掉空的和括号
                self.all_otu = [name for name in self.all_otu if name and name not in ['(', ')', ',']]
                print(f"Extracted {len(self.all_otu)} OTUs from tree (simple method)")
                
        except Exception as e:
            print(f"OTU extraction error: {e}")
            self.all_otu = []
    
    def discard_calibration_point_from_button(self, button):
        """通过Discard按钮删除校准点"""
        options_widget = button.parent()
        if not options_widget:
            return
            
        row = self.calibration_table.indexAt(options_widget.pos()).row()
        if row >= 0:
            self.calibration_table.removeRow(row)
            self.check_calibration_types()
    
    def check_calibration_types(self):
        """检查校准类型并显示警告"""
        complex_types = ['Normal', 'Lognormal']
        has_complex = False
        
        for row in range(self.calibration_table.rowCount()):
            cal_type = self.calibration_table.item(row, 2).text()
            if cal_type in complex_types:
                has_complex = True
                break
        
        self.calibration_warning_label.setVisible(has_complex)
    
    def get_calibration_data(self):
        """获取所有校准点数据（LSD2格式）"""
        calibration_data = []
        
        if not hasattr(self, 'calibration_data'):
            return calibration_data
        
        for row in range(self.calibration_table.rowCount()):
            # 只返回已配置的校准点
            if row in self.calibration_data and self.calibration_data[row] is not None:
                cal_data = self.calibration_data[row]
                calibration_data.append({
                    'name': cal_data['name'],
                    'set': cal_data['set'],
                    'type': cal_data['type'],
                    'values': cal_data['values']
                })
        
        return calibration_data
    
    def on_seq_length_mode_changed(self, index):
        """当序列长度模式改变时的回调"""
        mode = self.seq_length_combo.currentText()
        
        if mode == 'Disabled':
            self.sequence_length_spin.setEnabled(False)
            self.sequence_length_spin.setValue(0)
        elif mode == 'Manual':
            self.sequence_length_spin.setEnabled(True)
            if self.sequence_length_spin.value() == 0:
                self.sequence_length_spin.setValue(1000)  # 设置默认值
        elif mode == 'From sequence file':
            self.sequence_length_spin.setEnabled(False)
            # 从序列文件自动计算
            seq_file = self.seq_file_path.text().strip()
            if seq_file and os.path.exists(seq_file):
                seq_length = self.calculate_sequence_length(seq_file)
                if seq_length > 0:
                    self.sequence_length_spin.setValue(seq_length)
            else:
                self.sequence_length_spin.setValue(0)