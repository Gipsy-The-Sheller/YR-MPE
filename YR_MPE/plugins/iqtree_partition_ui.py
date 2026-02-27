# iqtree_partition_ui.py
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
IQ-TREE 分区模式UI组件
提供分区定义、文件导入导出和模式选择功能
"""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout,
                             QPushButton, QGroupBox, QLabel, QComboBox,
                             QTableWidget, QTableWidgetItem, QHeaderView,
                             QLineEdit, QFileDialog, QMessageBox, QCheckBox,
                             QSpinBox, QDoubleSpinBox)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
from typing import List, Optional
import os

from .partition_mode import PartitionDefinition, PartitionMode


class IQTreePartitionDialog(QDialog):
    """IQ-TREE 分区模式对话框"""
    
    # 定义信号
    partitions_selected = pyqtSignal(list)  # 分区定义已选择
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("IQ-TREE Partition Configuration")
        self.setMinimumSize(900, 700)
        
        # 分区定义列表
        self.partitions: List[PartitionDefinition] = []
        
        # 分区模式
        self.partition_mode = PartitionMode.EL
        
        # RCluster 设置
        self.rcluster_enabled = False
        self.rcluster_percent = 10
        
        self.init_ui()
    
    def init_ui(self):
        """初始化UI"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # 分区模式选择组
        mode_group = QGroupBox("Partition Mode")
        mode_layout = QFormLayout()
        mode_group.setLayout(mode_layout)
        
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Edge-linked Proportional (-p)",
            "Edge-linked Equal (-q)",
            "Edge-unlinked (-Q)",
            "Separate Tree (-S)"
        ])
        self.mode_combo.setCurrentIndex(0)  # 默认选择-p
        mode_layout.addRow("Mode:", self.mode_combo)
        
        # Bootstrap 参数
        bootstrap_layout = QHBoxLayout()
        self.enable_bootstrap = QCheckBox("Enable Bootstrap")
        self.enable_bootstrap.setChecked(True)
        self.enable_bootstrap.stateChanged.connect(self.on_bootstrap_toggled)
        bootstrap_layout.addWidget(self.enable_bootstrap)
        
        self.bootstrap_replicates = QSpinBox()
        self.bootstrap_replicates.setRange(100, 100000)
        self.bootstrap_replicates.setValue(1000)
        self.bootstrap_replicates.setEnabled(True)
        bootstrap_layout.addWidget(QLabel("Replicates:"))
        bootstrap_layout.addWidget(self.bootstrap_replicates)
        bootstrap_layout.addStretch()
        mode_layout.addRow("Bootstrap:", bootstrap_layout)
        
        # 抽样类型选择
        sampling_layout = QHBoxLayout()
        self.sampling_type_combo = QComboBox()
        self.sampling_type_combo.addItems([
            "SITE (Site-wise resampling)",  # 默认
            "GENE (Gene-wise resampling)",
            "GENESITE (Gene then Site resampling)"
        ])
        self.sampling_type_combo.setCurrentIndex(0)  # 默认SITE
        self.sampling_type_combo.setEnabled(True)
        bootstrap_layout.addWidget(QLabel("Sampling:"))
        bootstrap_layout.addWidget(self.sampling_type_combo)
        mode_layout.addRow("", sampling_layout)
        
        layout.addWidget(mode_group)
        
        # 分区定义组
        partition_group = QGroupBox("Partition Definitions")
        partition_layout = QVBoxLayout()
        partition_group.setLayout(partition_layout)
        
        # 分区表格
        self.partition_table = QTableWidget()
        self.partition_table.setColumnCount(5)
        self.partition_table.setHorizontalHeaderLabels([
            "Name", "File", "Seq Type", "Range", "Model"
        ])
        self.partition_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        partition_layout.addWidget(self.partition_table)
        
        # 表格按钮
        table_btn_layout = QHBoxLayout()
        self.add_partition_btn = QPushButton("+ Add Partition")
        self.remove_partition_btn = QPushButton("- Remove")
        self.import_partition_btn = QPushButton("Import File")
        self.export_partition_btn = QPushButton("Export File")
        self.clear_partitions_btn = QPushButton("Clear All")
        
        table_btn_layout.addWidget(self.add_partition_btn)
        table_btn_layout.addWidget(self.remove_partition_btn)
        table_btn_layout.addWidget(self.import_partition_btn)
        table_btn_layout.addWidget(self.export_partition_btn)
        table_btn_layout.addWidget(self.clear_partitions_btn)
        table_btn_layout.addStretch()
        
        partition_layout.addLayout(table_btn_layout)
        layout.addWidget(partition_group)
        
        # 高级选项组
        advanced_group = QGroupBox("Advanced Options")
        advanced_layout = QFormLayout()
        advanced_group.setLayout(advanced_layout)
        
        # RCluster 选项
        rcluster_layout = QHBoxLayout()
        self.use_rcluster = QCheckBox("Use RCluster for accelerated partition merging")
        self.use_rcluster.setChecked(False)
        self.rcluster_percent_spin = QSpinBox()
        self.rcluster_percent_spin.setRange(1, 100)
        self.rcluster_percent_spin.setValue(10)
        self.rcluster_percent_spin.setSuffix("%")
        self.rcluster_percent_spin.setEnabled(False)
        
        self.use_rcluster.stateChanged.connect(
            lambda state: self.rcluster_percent_spin.setEnabled(state == Qt.Checked)
        )
        
        rcluster_layout.addWidget(self.use_rcluster)
        rcluster_layout.addWidget(QLabel("Check"))
        rcluster_layout.addWidget(self.rcluster_percent_spin)
        rcluster_layout.addWidget(QLabel("of partition combinations"))
        rcluster_layout.addStretch()
        advanced_layout.addRow("RCluster:", rcluster_layout)
        
        # 输出前缀
        self.output_prefix_edit = QLineEdit()
        self.output_prefix_edit.setPlaceholderText("auto")
        advanced_layout.addRow("Output Prefix:", self.output_prefix_edit)
        
        layout.addWidget(advanced_group)
        
        # 控制按钮
        control_layout = QHBoxLayout()
        self.ok_btn = QPushButton("OK")
        self.ok_btn.setDefault(True)
        self.ok_btn.clicked.connect(self.accept_partitions)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        
        control_layout.addStretch()
        control_layout.addWidget(self.ok_btn)
        control_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(control_layout)
        
        # 连接信号
        self.add_partition_btn.clicked.connect(self.add_partition)
        self.remove_partition_btn.clicked.connect(self.remove_partition)
        self.import_partition_btn.clicked.connect(self.import_partition_file)
        self.export_partition_btn.clicked.connect(self.export_partition_file)
        self.clear_partitions_btn.clicked.connect(self.clear_partitions)
        
        # 初始化表格
        self.update_table()
    
    def set_partitions(self, partitions: List[PartitionDefinition]):
        """设置分区定义"""
        self.partitions = partitions.copy()
        self.update_table()
    
    def set_partition_mode(self, mode: PartitionMode):
        """设置分区模式"""
        self.partition_mode = mode
        
        # 更新下拉框
        mode_map = {
            PartitionMode.EL: 0,
            PartitionMode.EUL: 2,
            PartitionMode.TUL: 3
        }
        index = mode_map.get(mode, 0)
        self.mode_combo.setCurrentIndex(index)
    
    def get_partitions(self) -> List[PartitionDefinition]:
        """获取分区定义"""
        return self.partitions.copy()
    
    def get_partition_mode(self) -> PartitionMode:
        """获取分区模式"""
        mode_index = self.mode_combo.currentIndex()
        mode_map = {
            0: PartitionMode.EL,  # -p
            1: PartitionMode.EL,  # -q
            2: PartitionMode.EUL, # -Q
            3: PartitionMode.TUL  # -S
        }
        return mode_map.get(mode_index, PartitionMode.EL)
    
    def get_rcluster_enabled(self) -> bool:
        """获取是否启用RCluster"""
        return self.use_rcluster.isChecked()
    
    def get_rcluster_percent(self) -> Optional[int]:
        """获取RCluster百分比"""
        if self.use_rcluster.isChecked():
            return self.rcluster_percent_spin.value()
        return None
    
    def update_table(self):
        """更新分区表格"""
        self.partition_table.setRowCount(len(self.partitions))
        
        for row, partition in enumerate(self.partitions):
            # 分区名称
            name_item = QTableWidgetItem(partition.name)
            self.partition_table.setItem(row, 0, name_item)
            
            # 文件路径
            file_item = QTableWidgetItem(os.path.basename(partition.file_path) if partition.file_path else "")
            self.partition_table.setItem(row, 1, file_item)
            
            # 序列类型
            type_item = QTableWidgetItem(partition.seq_type)
            self.partition_table.setItem(row, 2, type_item)
            
            # 模型范围
            range_item = QTableWidgetItem(partition.model_range)
            self.partition_table.setItem(row, 3, range_item)
            
            # 选定模型
            model_item = QTableWidgetItem(partition.selected_model or "Not selected")
            self.partition_table.setItem(row, 4, model_item)
    
    def add_partition(self):
        """添加新分区"""
        # 简单实现：添加一个默认分区
        # 实际应用中应该弹出一个对话框让用户输入详细信息
        partition = PartitionDefinition(
            name=f"Partition{len(self.partitions) + 1}",
            file_path="",
            seq_type="DNA",
            model_range=f"{len(self.partitions) * 1000 + 1}-{(len(self.partitions) + 1) * 1000}",
            selected_model=None
        )
        self.partitions.append(partition)
        self.update_table()
    
    def remove_partition(self):
        """删除选中的分区"""
        current_row = self.partition_table.currentRow()
        if current_row >= 0 and current_row < len(self.partitions):
            del self.partitions[current_row]
            self.update_table()
    
    def clear_partitions(self):
        """清空所有分区"""
        self.partitions.clear()
        self.update_table()
    
    def import_partition_file(self):
        """导入分区文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Import Partition File",
            "",
            "Partition files (*.nex *.nexus *.partition);;All files (*)"
        )
        
        if not file_path:
            return
        
        try:
            from .partition_mode import PartitionFileIO
            
            # 根据文件扩展名选择解析器
            file_ext = os.path.splitext(file_path)[1].lower()
            
            if file_ext in ['.nex', '.nexus']:
                partitions = PartitionFileIO.parse_nexus_partition(file_path)
            else:
                partitions = PartitionFileIO.parse_raxml_partition(file_path)
            
            if partitions:
                self.partitions = partitions
                self.update_table()
                QMessageBox.information(self, "Success", f"Imported {len(partitions)} partitions from {os.path.basename(file_path)}")
            else:
                QMessageBox.warning(self, "Warning", "No partitions found in the file.")
        
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import partition file: {str(e)}")
    
    def export_partition_file(self):
        """导出分区文件"""
        if not self.partitions:
            QMessageBox.warning(self, "Warning", "No partitions to export.")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Partition File",
            "",
            "NEXUS format (*.nex);;RAxML format (*.partition)"
        )
        
        if not file_path:
            return
        
        try:
            from .partition_mode import PartitionFileIO
            
            # 根据文件扩展名选择导出格式
            file_ext = os.path.splitext(file_path)[1].lower()
            
            if file_ext in ['.nex', '.nexus']:
                PartitionFileIO.export_nexus_partition(self.partitions, file_path)
            else:
                PartitionFileIO.export_raxml_partition(self.partitions, file_path)
            
            QMessageBox.information(self, "Success", f"Exported {len(self.partitions)} partitions to {os.path.basename(file_path)}")
        
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export partition file: {str(e)}")
    
    def get_bootstrap_enabled(self) -> bool:
        """获取是否启用Bootstrap"""
        return self.enable_bootstrap.isChecked()
    
    def get_bootstrap_replicates(self) -> int:
        """获取Bootstrap重复次数"""
        return self.bootstrap_replicates.value()
    
    def get_sampling_type(self) -> str:
        """获取抽样类型"""
        return self.sampling_type_combo.currentText().split()[0]  # 返回GENE/SITE/GENESITE
    
    def on_bootstrap_toggled(self, state):
        """Bootstrap复选框切换"""
        enabled = (state == Qt.Checked)
        self.bootstrap_replicates.setEnabled(enabled)
        self.sampling_type_combo.setEnabled(enabled)
    
    def accept_partitions(self):
        """接受分区配置"""
        # 验证分区定义
        if not self.partitions:
            reply = QMessageBox.question(
                self,
                "Confirm",
                "No partitions defined. Continue with partition mode disabled?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        # 发送信号
        self.partitions_selected.emit(self.partitions)
        
        # 接受对话框
        self.accept()