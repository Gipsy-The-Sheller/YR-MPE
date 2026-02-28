# lsd2_ui.py
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
                            QHeaderView, QGroupBox, QFileDialog, QSpinBox, QSizePolicy, QCheckBox)
from PyQt5.QtCore import Qt
import os
from typing import List, Dict


class LSD2UI(QWidget):
    """LSD2用户界面组件"""
    
    def __init__(self):
        super().__init__()
        self.calibration_data = {}  # 校准点数据
        self.all_otu = []  # 从树中提取的所有OTU
        self.rooted_tree_file = None  # 置根后的树文件路径
        self.tip_calibrations = {}  # Tip标定数据
        self.init_ui()
    
    def init_ui(self):
        """初始化UI"""
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # Dating methods（虽然只有一个方法，但保留以备未来扩展）
        dating_methods_layout = QHBoxLayout()
        
        self.methods_combo = QComboBox()
        self.methods_combo.addItems(['LSD2'])
        self.methods_combo.setEnabled(False)  # 暂时禁用，因为只有一个方法

        dating_methods_label = QLabel('Dating methods:')
        dating_methods_label.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)

        self.methods_combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        dating_methods_layout.addWidget(dating_methods_label)
        dating_methods_layout.addWidget(self.methods_combo)

        main_layout.addLayout(dating_methods_layout)

        # 输入文件管理
        self.setup_input_files(main_layout)
        
        # 置根选项
        self.setup_rooting_options(main_layout)
        
        # Tip Dating 选项
        self.setup_tip_dating(main_layout)
        
        # 校准点管理
        self.setup_calibration_points(main_layout)
    
    def setup_input_files(self, parent_layout):
        """设置输入文件管理"""
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
        parent_layout.addWidget(input_files_group)
    
    def setup_rooting_options(self, parent_layout):
        """设置置根选项"""
        root_layout = QHBoxLayout()
        root_layout.addWidget(QLabel("Root the tree:"))
        
        self.choose_root_button = QPushButton("Choose root")
        self.choose_root_button.clicked.connect(self.on_choose_root_clicked)
        root_layout.addWidget(self.choose_root_button)
        
        # 显示置根状态的标签
        self.root_status_label = QLabel("(Unrooted)")
        self.root_status_label.setStyleSheet("color: gray;")
        root_layout.addWidget(self.root_status_label)
        
        root_layout.addStretch()
        parent_layout.addLayout(root_layout)
        
        # Rooting option radio buttons
        rooting_option_layout = QHBoxLayout()
        rooting_option_layout.addWidget(QLabel("Rooting option:"))
        
        self.rooting_button_group = QButtonGroup(self)
        self.use_chosen_root_radio = QRadioButton("Use chosen root")
        self.estimate_root_radio = QRadioButton("Estimate the root")
        self.estimate_root_radio.setChecked(True)  # 默认选中
        
        self.rooting_button_group.addButton(self.use_chosen_root_radio)
        self.rooting_button_group.addButton(self.estimate_root_radio)
        
        # Use chosen root 默认禁用，只有置根后才启用
        self.use_chosen_root_radio.setEnabled(False)
        
        rooting_option_layout.addWidget(self.use_chosen_root_radio)
        rooting_option_layout.addWidget(self.estimate_root_radio)
        rooting_option_layout.addStretch()
        parent_layout.addLayout(rooting_option_layout)
    
    def setup_tip_dating(self, parent_layout):
        """设置Tip Dating选项"""
        tip_dating_layout = QHBoxLayout()
        
        self.use_tip_dates_checkbox = QCheckBox("Use tip dates")
        self.use_tip_dates_checkbox.stateChanged.connect(self.on_tip_dates_toggled)
        tip_dating_layout.addWidget(self.use_tip_dates_checkbox)
        
        self.configure_tip_dates_button = QPushButton("Configure tip dates")
        self.configure_tip_dates_button.clicked.connect(self.on_configure_tip_dates)
        self.configure_tip_dates_button.setEnabled(False)  # 初始禁用
        tip_dating_layout.addWidget(self.configure_tip_dates_button)
        
        # 显示已配置tip数量的标签
        self.tip_count_label = QLabel("(0 tips configured)")
        self.tip_count_label.setStyleSheet("color: gray;")
        tip_dating_layout.addWidget(self.tip_count_label)
        
        tip_dating_layout.addStretch()
        parent_layout.addLayout(tip_dating_layout)
    
    def setup_calibration_points(self, parent_layout):
        """设置校准点管理"""
        calibration_layout = QHBoxLayout()
        calibration_layout.addWidget(QLabel("Calibration points:"))
        
        self.add_calibration_button = QPushButton("+")
        self.add_calibration_button.clicked.connect(self.add_calibration_point)
        calibration_layout.addWidget(self.add_calibration_button)
        calibration_layout.addStretch()
        parent_layout.addLayout(calibration_layout)
        
        # Calibration table
        self.calibration_table = QTableWidget(0, 4)
        self.calibration_table.setHorizontalHeaderLabels(['Name', 'Taxa', 'Type', 'Options'])
        self.calibration_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.calibration_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        parent_layout.addWidget(self.calibration_table)
        
        # Warning label for calibration types
        self.calibration_warning_label = QLabel(
            "* WARNING: Calibrations except for uniform / point / upper / lower "
            "will be automatically transferred to uniform(lower 95HPD, upper 95HPD)"
        )
        self.calibration_warning_label.setStyleSheet("color: red;")
        self.calibration_warning_label.setVisible(False)  # 初始隐藏
        self.calibration_warning_label.setWordWrap(True)
        parent_layout.addWidget(self.calibration_warning_label)
        
        # Confidence intervals settings
        ci_group = QGroupBox("Confidence Intervals")
        ci_layout = QVBoxLayout()
        
        # Calculate confidence intervals checkbox
        ci_checkbox_layout = QHBoxLayout()
        self.calculate_ci_checkbox = QCheckBox("Calculate confidence intervals")
        self.calculate_ci_checkbox.setChecked(False)
        ci_checkbox_layout.addWidget(self.calculate_ci_checkbox)
        ci_layout.addLayout(ci_checkbox_layout)
        
        # Lognormal distribution std
        lognormal_std_layout = QHBoxLayout()
        lognormal_std_layout.addWidget(QLabel("Lognormal distribution std (q):"))
        self.lognormal_std_edit = QLineEdit("0.2")
        self.lognormal_std_edit.setPlaceholderText("Default: 0.2")
        self.lognormal_std_edit.setEnabled(False)  # 初始禁用
        lognormal_std_layout.addWidget(self.lognormal_std_edit)
        ci_layout.addLayout(lognormal_std_layout)
        
        ci_group.setLayout(ci_layout)
        parent_layout.addWidget(ci_group)
        
        # 连接信号
        self.calculate_ci_checkbox.toggled.connect(self.on_ci_checkbox_toggled)
    
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
                # 提取OTU列表
                self.extract_otu_from_tree(file_path)
                # 重置置根状态
                self.reset_rooting_status()
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
            if not content.endswith(';'):
                return False
            
            # 检查括号匹配
            open_parens = content.count('(')
            close_parens = content.count(')')
            if open_parens != close_parens:
                return False
            
            # 检查是否有至少一个分类单元
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
            
            # 使用ETE3提取OTU
            try:
                from ete3 import Tree
                # 使用 format=1 来处理引号节点名称
                tree = Tree(newick_str, format=1)
                self.all_otu = [leaf.name for leaf in tree.get_leaves()]
                print(f"Extracted {len(self.all_otu)} OTUs from tree")
            except Exception as e:
                print(f"ETE3 parsing failed: {e}, falling back to simple method")
                # 如果ETE3解析失败，使用简单的正则表达式
                import re
                clean_newick = re.sub(r':\d+\.?\d*', '', newick_str)
                self.all_otu = re.findall(r'([A-Za-z0-9_\-\.]+)', clean_newick)
                self.all_otu = [name for name in self.all_otu if name and name not in ['(', ')', ',']]
                print(f"Extracted {len(self.all_otu)} OTUs from tree (simple method)")
                
        except Exception as e:
            print(f"OTU extraction error: {e}")
            self.all_otu = []
    
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
            seq_file = self.seq_file_path.text().strip()
            if seq_file and os.path.exists(seq_file):
                seq_length = self.calculate_sequence_length(seq_file)
                if seq_length > 0:
                    self.sequence_length_spin.setValue(seq_length)
            else:
                self.sequence_length_spin.setValue(0)
    
    def on_choose_root_clicked(self):
        """点击Choose root按钮的回调"""
        tree_file = self.tree_file_path.text().strip()
        if not tree_file or not os.path.exists(tree_file):
            QMessageBox.warning(self, "No Tree File", "Please load a tree file first.")
            return
        
        try:
            with open(tree_file, 'r') as f:
                newick_str = f.read().strip()
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to read tree file: {str(e)}")
            return
        
        # 导入并打开IcyTree Dialog（在主插件中实现）
        # 这里只是发出信号，由主插件处理
        if hasattr(self, 'open_icytree_dialog_callback'):
            self.open_icytree_dialog_callback(newick_str)
        else:
            QMessageBox.warning(self, "Not Implemented", "IcyTree Dialog not yet implemented.")
    
    def on_rooted_tree_loaded(self, rooted_tree_file):
        """当置根后的树加载完成时的回调"""
        self.rooted_tree_file = rooted_tree_file
        self.use_chosen_root_radio.setEnabled(True)
        self.use_chosen_root_radio.setChecked(True)
        self.root_status_label.setText("(Rooted)")
        self.root_status_label.setStyleSheet("color: green;")
    
    def reset_rooting_status(self):
        """重置置根状态"""
        self.rooted_tree_file = None
        self.use_chosen_root_radio.setEnabled(False)
        self.estimate_root_radio.setChecked(True)
        self.root_status_label.setText("(Unrooted)")
        self.root_status_label.setStyleSheet("color: gray;")
    
    def on_tip_dates_toggled(self, state):
        """当Use tip dates复选框状态改变时的回调"""
        is_checked = (state == Qt.Checked)
        self.configure_tip_dates_button.setEnabled(is_checked)
        
        if not is_checked:
            # 清空tip标定
            self.tip_calibrations.clear()
            self.update_tip_count_label()
    
    def on_configure_tip_dates(self):
        """点击Configure tip dates按钮的回调"""
        if not self.all_otu:
            QMessageBox.warning(self, "No OTUs", "Please load a tree file first.")
            return
        
        # 导入并打开Tip Dating Dialog（在主插件中实现）
        # 这里只是发出信号，由主插件处理
        if hasattr(self, 'open_tip_dating_dialog_callback'):
            self.open_tip_dating_dialog_callback(self.all_otu)
        else:
            QMessageBox.warning(self, "Not Implemented", "Tip Dating Dialog not yet implemented.")
    
    def on_tip_calibrations_updated(self, tip_calibrations):
        """当Tip标定更新时的回调"""
        self.tip_calibrations = tip_calibrations
        self.update_tip_count_label()
    
    def update_tip_count_label(self):
        """更新已配置tip数量的标签"""
        count = len([cal for cal in self.tip_calibrations.values() if cal is not None])
        if count > 0:
            self.tip_count_label.setText(f"({count} tips configured)")
            self.tip_count_label.setStyleSheet("color: green;")
        else:
            self.tip_count_label.setText("(0 tips configured)")
            self.tip_count_label.setStyleSheet("color: gray;")
    
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
        
        # 连接按钮信号
        edit_button.clicked.connect(lambda checked, btn=edit_button: self.edit_calibration_point_from_button(btn))
        discard_button.clicked.connect(lambda checked, btn=discard_button: self.discard_calibration_point_from_button(btn))
        
        # 初始化校准数据存储
        self.calibration_data[row] = None
    
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
            
        # 获取当前Newick树字符串（使用置根后的树，如果有的话）
        tree_file = self.rooted_tree_file if self.rooted_tree_file else self.tree_file_path.text().strip()
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
            
        # 创建并显示MD_MRCA对话框（由主插件处理）
        if hasattr(self, 'open_md_mrca_dialog_callback'):
            self.open_md_mrca_dialog_callback(row, newick_text)
        else:
            QMessageBox.warning(self, "Not Implemented", "MD-MRCA Dialog not yet implemented.")
    
    def on_md_mrca_result(self, row, mrca_result):
        """处理MD_MRCA的结果"""
        selected_taxa = mrca_result['selected_taxa']
        taxon_set_name = mrca_result['taxon_set_name']
        cal_type = mrca_result['cal_type']
        cal_values = mrca_result['cal_values']
        
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
        
        # 存储完整的校准点数据
        self.calibration_data[row] = {
            'name': taxon_set_name,
            'set': selected_taxa,
            'type': lsd2_type,
            'values': cal_values,
            'display_type': cal_type  # 保存原始显示类型
        }
        
        # 触发校准类型检查
        self.check_calibration_types()
    
    def discard_calibration_point_from_button(self, button):
        """通过Discard按钮删除校准点"""
        options_widget = button.parent()
        if not options_widget:
            return
            
        row = self.calibration_table.indexAt(options_widget.pos()).row()
        if row >= 0:
            self.calibration_table.removeRow(row)
            # 删除对应的校准数据
            if row in self.calibration_data:
                del self.calibration_data[row]
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
    
    def get_tip_calibrations(self):
        """获取Tip标定数据"""
        return self.tip_calibrations
    
    def get_rooted_tree_file(self):
        """获取置根后的树文件路径"""
        return self.rooted_tree_file
    
    def get_original_tree_file(self):
        """获取原始树文件路径"""
        return self.tree_file_path.text().strip()
    
    def get_sequence_length(self):
        """获取序列长度"""
        seq_length_mode = self.seq_length_combo.currentText()
        if seq_length_mode == 'Disabled':
            return 0
        else:
            return self.sequence_length_spin.value()
    
    def get_use_chosen_root(self):
        """是否使用用户选择的根"""
        return self.use_chosen_root_radio.isChecked() and self.rooted_tree_file is not None
    
    def on_ci_checkbox_toggled(self, checked):
        """当置信区间复选框状态改变时的回调"""
        self.lognormal_std_edit.setEnabled(checked)
    
    def get_calculate_ci(self):
        """是否计算置信区间"""
        return self.calculate_ci_checkbox.isChecked()
    
    def get_lognormal_std(self):
        """获取lognormal分布标准差"""
        try:
            return float(self.lognormal_std_edit.text())
        except ValueError:
            return 0.2  # 默认值