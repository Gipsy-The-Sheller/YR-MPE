# mrbayes_partition_ui.py
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

"""
MrBayes 分区模式UI组件
提供分区定义、模式选择和文件导入导出功能
"""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout,
                             QPushButton, QGroupBox, QLabel, QComboBox,
                             QTableWidget, QTableWidgetItem, QHeaderView,
                             QLineEdit, QFileDialog, QMessageBox, QCheckBox,
                             QSpinBox, QRadioButton, QButtonGroup, QFrame)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
from typing import List, Optional
import os

from .partition_mode import (
    PartitionMode, MrBayesPartitionDefinition, 
    MrBayesModelConverter, mrbayes_partition_to_dict, mrbayes_partition_from_dict
)


class MrBayesPartitionEditDialog(QDialog):
    """MrBayes分区编辑对话框"""
    
    def __init__(self, parent=None, partition: Optional[MrBayesPartitionDefinition] = None):
        super().__init__(parent)
        self.setWindowTitle("Edit Partition")
        self.setMinimumSize(400, 350)
        
        self.partition = partition or MrBayesPartitionDefinition(
            name="", range="", seq_type="DNA"
        )
        
        self.init_ui()
        self.load_data()
    
    def init_ui(self):
        """初始化UI"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        form_layout = QFormLayout()
        
        # 名称
        self.name_edit = QLineEdit()
        form_layout.addRow("Name:", self.name_edit)
        
        # 范围
        self.range_edit = QLineEdit()
        self.range_edit.setPlaceholderText("e.g., 1-675\\3 3784-4893\\3")
        form_layout.addRow("Range:", self.range_edit)
        
        # 序列类型
        self.seq_type_combo = QComboBox()
        self.seq_type_combo.addItems(["DNA", "Protein"])
        self.seq_type_combo.currentTextChanged.connect(self.on_seq_type_changed)
        form_layout.addRow("Seq Type:", self.seq_type_combo)
        
        # 分隔线
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        layout.addWidget(separator)
        
        # DNA模型参数（仅DNA显示）
        self.dna_group = QGroupBox("DNA Model Parameters")
        dna_layout = QFormLayout()
        self.dna_group.setLayout(dna_layout)
        
        self.nst_combo = QComboBox()
        self.nst_combo.addItem("JC69 (nst=1)", 1)
        self.nst_combo.addItem("HKY85 (nst=2)", 2)
        self.nst_combo.addItem("GTR (nst=6)", 6)
        self.nst_combo.setCurrentIndex(2)
        dna_layout.addRow("Substitution Model:", self.nst_combo)
        
        layout.addWidget(self.dna_group)
        
        # Protein模型参数（仅Protein显示）
        self.protein_group = QGroupBox("Protein Model Parameters")
        protein_layout = QFormLayout()
        self.protein_group.setLayout(protein_layout)
        
        self.aamodel_combo = QComboBox()
        self.aamodel_combo.addItems([
            "Blosum62", "Blosum", "Wag", "Lg", 
            "gtr", "jones", "mtrev", "Poisson", "mixed"
        ])
        protein_layout.addRow("Amino Acid Model:", self.aamodel_combo)
        
        layout.addWidget(self.protein_group)
        
        # 通用参数
        common_group = QGroupBox("Common Parameters")
        common_layout = QFormLayout()
        common_group.setLayout(common_layout)
        
        self.rates_combo = QComboBox()
        self.rates_combo.addItems([
            "Equal", "Gamma (+G)", "InvGamma (+G+I)", 
            "PropInv (+I)", "Lognormal", "Adgamma"
        ])
        common_layout.addRow("Rate Heterogenity:", self.rates_combo)
        
        self.ngammacat_spinbox = QSpinBox()
        self.ngammacat_spinbox.setRange(1, 100)
        self.ngammacat_spinbox.setValue(4)
        common_layout.addRow("Gamma Categories:", self.ngammacat_spinbox)
        
        layout.addWidget(common_group)
        
        # 按钮
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        button_layout.addWidget(ok_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
    
    def on_seq_type_changed(self, text: str):
        """序列类型改变时更新UI"""
        is_dna = text.upper() == "DNA"
        self.dna_group.setVisible(is_dna)
        self.protein_group.setVisible(not is_dna)
    
    def load_data(self):
        """加载分区数据"""
        self.name_edit.setText(self.partition.name)
        self.range_edit.setText(self.partition.range)
        self.seq_type_combo.setCurrentText(self.partition.seq_type)
        
        # DNA参数
        if self.partition.seq_type == "DNA":
            self.nst_combo.setCurrentIndex(
                self.nst_combo.findData(self.partition.nst) or 2
            )
        
        # Protein参数
        if self.partition.seq_type == "PROTEIN":
            self.aamodel_combo.setCurrentText(self.partition.aamodel)
        
        # 通用参数
        rates_map = {
            "equal": "Equal",
            "gamma": "Gamma (+G)",
            "invgamma": "InvGamma (+G+I)",
            "propinv": "PropInv (+I)",
            "lnorm": "Lognormal",
            "adgamma": "Adgamma"
        }
        self.rates_combo.setCurrentText(rates_map.get(self.partition.rates, "Gamma (+G)"))
        self.ngammacat_spinbox.setValue(self.partition.ngammacat)
    
    def get_partition(self) -> Optional[MrBayesPartitionDefinition]:
        """获取编辑后的分区定义"""
        name = self.name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Warning", "Partition name cannot be empty!")
            return None
        
        range_def = self.range_edit.text().strip()
        if not range_def:
            QMessageBox.warning(self, "Warning", "Partition range cannot be empty!")
            return None
        
        seq_type = self.seq_type_combo.currentText().upper()
        
        # DNA参数
        nst = 6
        if seq_type == "DNA":
            nst = self.nst_combo.currentData() or 6
        
        # Protein参数
        aamodel = "mixed"
        if seq_type == "PROTEIN":
            aamodel = self.aamodel_combo.currentText()
        
        # 通用参数
        rates_map = {
            "Equal": "equal",
            "Gamma (+G)": "gamma",
            "InvGamma (+G+I)": "invgamma",
            "PropInv (+I)": "propinv",
            "Lognormal": "lnorm",
            "Adgamma": "adgamma"
        }
        rates = rates_map.get(self.rates_combo.currentText(), "gamma")
        ngammacat = self.ngammacat_spinbox.value()
        
        partition = MrBayesPartitionDefinition(
            name=name,
            range=range_def,
            seq_type=seq_type,
            nst=nst,
            aamodel=aamodel,
            rates=rates,
            ngammacat=ngammacat
        )
        
        # 验证
        valid, error = partition.validate()
        if not valid:
            QMessageBox.warning(self, "Validation Error", error)
            return None
        
        return partition


class MrBayesPartitionDialog(QDialog):
    """MrBayes 分区模式对话框"""
    
    partitions_updated = pyqtSignal(list)
    
    def __init__(self, parent=None, partitions: Optional[List[MrBayesPartitionDefinition]] = None,
                 mode: PartitionMode = PartitionMode.EL):
        super().__init__(parent)
        self.setWindowTitle("MrBayes Partition Configuration")
        self.setMinimumSize(900, 700)
        
        self.partitions = partitions or []
        self.partition_mode = mode
        
        self.init_ui()
        self.update_table()
    
    def init_ui(self):
        """初始化UI"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # 分区模式选择组
        mode_group = QGroupBox("Partition Mode")
        mode_layout = QVBoxLayout()
        mode_group.setLayout(mode_layout)
        
        mode_desc = QLabel(
            "MrBayes支持两种分区模式:\n"
            "• Edge-linked: 共享拓扑和分支长度，每个分区有独立的模型参数\n"
            "• Edge-unlinked: 共享拓扑，但所有参数都独立"
        )
        mode_desc.setWordWrap(True)
        mode_desc.setStyleSheet("color: #666;")
        mode_layout.addWidget(mode_desc)
        
        mode_radio_layout = QHBoxLayout()
        self.el_radio = QRadioButton("Edge-linked")
        self.eul_radio = QRadioButton("Edge-unlinked")
        
        self.mode_group = QButtonGroup()
        self.mode_group.addButton(self.el_radio, 0)
        self.mode_group.addButton(self.eul_radio, 1)
        
        if self.partition_mode == PartitionMode.EL:
            self.el_radio.setChecked(True)
        else:
            self.eul_radio.setChecked(True)
        
        mode_radio_layout.addWidget(self.el_radio)
        mode_radio_layout.addWidget(self.eul_radio)
        mode_radio_layout.addStretch()
        mode_layout.addLayout(mode_radio_layout)
        
        layout.addWidget(mode_group)
        
        # 分区定义组
        partition_group = QGroupBox("Partition Definitions")
        partition_layout = QVBoxLayout()
        partition_group.setLayout(partition_layout)
        
        # 分区表格
        self.partition_table = QTableWidget()
        self.partition_table.setColumnCount(6)
        self.partition_table.setHorizontalHeaderLabels([
            "Name", "Range", "Seq Type", "Model", "Rates", "ngammacat"
        ])
        self.partition_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.partition_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.partition_table.setSelectionMode(QTableWidget.SingleSelection)
        partition_layout.addWidget(self.partition_table)
        
        # 表格按钮
        table_btn_layout = QHBoxLayout()
        self.add_partition_btn = QPushButton("+ Add")
        self.edit_partition_btn = QPushButton("Edit")
        self.remove_partition_btn = QPushButton("- Remove")
        self.import_partition_btn = QPushButton("Import File")
        self.export_partition_btn = QPushButton("Export File")
        self.clear_partitions_btn = QPushButton("Clear All")
        
        table_btn_layout.addWidget(self.add_partition_btn)
        table_btn_layout.addWidget(self.edit_partition_btn)
        table_btn_layout.addWidget(self.remove_partition_btn)
        table_btn_layout.addWidget(self.import_partition_btn)
        table_btn_layout.addWidget(self.export_partition_btn)
        table_btn_layout.addWidget(self.clear_partitions_btn)
        table_btn_layout.addStretch()
        
        partition_layout.addLayout(table_btn_layout)
        layout.addWidget(partition_group)
        
        # 预览组
        preview_group = QGroupBox("Preview MrBayes Data Block")
        preview_layout = QVBoxLayout()
        preview_group.setLayout(preview_layout)
        
        preview_btn_layout = QHBoxLayout()
        self.show_preview_btn = QPushButton("Show Preview")
        self.copy_preview_btn = QPushButton("Copy to Clipboard")
        preview_btn_layout.addWidget(self.show_preview_btn)
        preview_btn_layout.addWidget(self.copy_preview_btn)
        preview_btn_layout.addStretch()
        
        preview_layout.addLayout(preview_btn_layout)
        layout.addWidget(preview_group)
        
        # 控制按钮
        control_layout = QHBoxLayout()
        self.ok_btn = QPushButton("OK")
        self.ok_btn.setDefault(True)
        self.ok_btn.clicked.connect(self.accept)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        
        control_layout.addStretch()
        control_layout.addWidget(self.ok_btn)
        control_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(control_layout)
        
        # 连接信号
        self.add_partition_btn.clicked.connect(self.add_partition)
        self.edit_partition_btn.clicked.connect(self.edit_partition)
        self.remove_partition_btn.clicked.connect(self.remove_partition)
        self.import_partition_btn.clicked.connect(self.import_partition_file)
        self.export_partition_btn.clicked.connect(self.export_partition_file)
        self.clear_partitions_btn.clicked.connect(self.clear_partitions)
        self.show_preview_btn.clicked.connect(self.show_preview)
        self.copy_preview_btn.clicked.connect(self.copy_preview)
    
    def update_table(self):
        """更新分区表格"""
        self.partition_table.setRowCount(len(self.partitions))
        
        for row, partition in enumerate(self.partitions):
            self.partition_table.setItem(row, 0, QTableWidgetItem(partition.name))
            self.partition_table.setItem(row, 1, QTableWidgetItem(partition.range))
            self.partition_table.setItem(row, 2, QTableWidgetItem(partition.seq_type))
            self.partition_table.setItem(row, 3, QTableWidgetItem(partition.get_model_display()))
            self.partition_table.setItem(row, 4, QTableWidgetItem(partition.rates))
            self.partition_table.setItem(row, 5, QTableWidgetItem(str(partition.ngammacat)))
    
    def add_partition(self):
        """添加分区"""
        dialog = MrBayesPartitionEditDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            partition = dialog.get_partition()
            if partition:
                self.partitions.append(partition)
                self.update_table()
    
    def edit_partition(self):
        """编辑分区"""
        current_row = self.partition_table.currentRow()
        if current_row < 0:
            QMessageBox.warning(self, "Warning", "Please select a partition to edit!")
            return
        
        partition = self.partitions[current_row]
        dialog = MrBayesPartitionEditDialog(self, partition)
        if dialog.exec_() == QDialog.Accepted:
            new_partition = dialog.get_partition()
            if new_partition:
                self.partitions[current_row] = new_partition
                self.update_table()
    
    def remove_partition(self):
        """移除分区"""
        current_row = self.partition_table.currentRow()
        if current_row < 0:
            QMessageBox.warning(self, "Warning", "Please select a partition to remove!")
            return
        
        partition = self.partitions[current_row]
        reply = QMessageBox.question(
            self, "Confirm Removal",
            f"Are you sure you want to remove partition '{partition.name}'?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            del self.partitions[current_row]
            self.update_table()
    
    def clear_partitions(self):
        """清空所有分区"""
        if not self.partitions:
            return
        
        reply = QMessageBox.question(
            self, "Confirm Clear",
            f"Are you sure you want to clear all {len(self.partitions)} partitions?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            self.partitions.clear()
            self.update_table()
    
    def import_partition_file(self):
        """导入分区文件"""
        file_filter = "NEXUS files (*.nex *.nexus);;All files (*)"
        file_path, _ = QFileDialog.getOpenFileName(self, "Import Partition File", "", file_filter)
        
        if not file_path:
            return
        
        try:
            # 尝试从NEXUS文件解析分区
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            partitions, mode = MrBayesModelConverter.parse_mrbayes_block(content)
            
            if partitions:
                self.partitions = partitions
                if mode == PartitionMode.EUL:
                    self.eul_radio.setChecked(True)
                else:
                    self.el_radio.setChecked(True)
                self.update_table()
                QMessageBox.information(
                    self, "Import Successful",
                    f"Successfully imported {len(partitions)} partition(s)"
                )
            else:
                QMessageBox.warning(self, "Import Failed", "No partitions found in the file")
                
        except Exception as e:
            QMessageBox.critical(self, "Import Error", f"Failed to import partition file: {str(e)}")
    
    def export_partition_file(self):
        """导出分区文件"""
        if not self.partitions:
            QMessageBox.warning(self, "Warning", "No partitions to export!")
            return
        
        file_filter = "NEXUS files (*.nex *.nexus);;All files (*)"
        file_path, _ = QFileDialog.getSaveFileName(self, "Export Partition File", "", file_filter)
        
        if not file_path:
            return
        
        try:
            # 生成MrBayes分区命令
            mode = PartitionMode.EUL if self.eul_radio.isChecked() else PartitionMode.EL
            commands = MrBayesModelConverter.generate_partition_commands(self.partitions, mode)
            
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write("#nexus\n")
                f.write("begin mrbayes;\n")
                for cmd in commands:
                    f.write(f"    {cmd}\n")
                f.write("end;\n")
            
            QMessageBox.information(self, "Export Successful", f"Partition file saved to:\n{file_path}")
            
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to export partition file: {str(e)}")
    
    def show_preview(self):
        """显示MrBayes数据块预览"""
        if not self.partitions:
            QMessageBox.warning(self, "Warning", "No partitions to preview!")
            return
        
        preview = self.generate_mrbayes_block()
        
        dialog = QDialog(self)
        dialog.setWindowTitle("MrBayes Data Block Preview")
        dialog.setMinimumSize(700, 500)
        
        layout = QVBoxLayout()
        
        from PyQt5.QtWidgets import QTextEdit
        text_edit = QTextEdit()
        text_edit.setFont(QFont("Consolas", 10))
        text_edit.setStyleSheet("""
            QTextEdit {
                background-color: #272822;
                color: #f8f8f2;
                border: 1px solid #3c3c3c;
            }
        """)
        text_edit.setPlainText(preview)
        text_edit.setReadOnly(True)
        
        layout.addWidget(text_edit)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)
        
        dialog.exec_()
    
    def copy_preview(self):
        """复制预览到剪贴板"""
        if not self.partitions:
            QMessageBox.warning(self, "Warning", "No partitions to copy!")
            return
        
        preview = self.generate_mrbayes_block()
        
        from PyQt5.QtWidgets import QApplication
        clipboard = QApplication.clipboard()
        clipboard.setText(preview)
        
        QMessageBox.information(self, "Copied", "MrBayes data block copied to clipboard!")
    
    def generate_mrbayes_block(self) -> str:
        """生成MrBayes数据块"""
        mode = PartitionMode.EUL if self.eul_radio.isChecked() else PartitionMode.EL
        commands = MrBayesModelConverter.generate_partition_commands(self.partitions, mode)
        
        block = ["begin mrbayes;"]
        for cmd in commands:
            block.append(f"    {cmd}")
        block.append("end;")
        
        return "\n".join(block)
    
    def get_partitions(self) -> List[MrBayesPartitionDefinition]:
        """获取分区列表"""
        return self.partitions.copy()
    
    def get_mode(self) -> PartitionMode:
        """获取分区模式"""
        return PartitionMode.EUL if self.eul_radio.isChecked() else PartitionMode.EL
    
    def accept(self):
        """接受对话框"""
        if not self.partitions:
            reply = QMessageBox.question(
                self, "Confirm",
                "No partitions defined. Continue with single partition mode?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        super().accept()