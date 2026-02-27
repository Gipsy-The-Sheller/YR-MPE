# model_finder_partition_ui.py
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
ModelFinder分区模式UI组件
提供分区定义管理、模式选择和文件导入导出功能
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout,
    QPushButton, QLabel, QComboBox, QTableWidget, QTableWidgetItem,
    QHeaderView, QMessageBox, QFileDialog, QGroupBox, QSpinBox,
    QCheckBox, QLineEdit, QAbstractItemView, QMenu, QAction
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor
import os
from typing import List, Optional

from .partition_mode import (
    PartitionDefinition, PartitionMode, PartitionModeConverter,
    PartitionFileIO, partition_to_dict, create_partition_from_dict
)


class PartitionTableModel:
    """分区表格模型"""
    
    def __init__(self):
        self.partitions: List[PartitionDefinition] = []
        self.headers = ["Name", "File", "Model Range", "Selected Model"]
        self.validation_status = {}  # partition_name -> (is_valid, error_message)
    
    def add_partition(self, partition: PartitionDefinition) -> bool:
        """添加分区"""
        # 验证分区
        valid, error = partition.validate()
        self.validation_status[partition.name] = (valid, error)
        
        if not valid:
            return False
        
        self.partitions.append(partition)
        return True
    
    def remove_partition(self, index: int):
        """移除分区"""
        if 0 <= index < len(self.partitions):
            partition = self.partitions.pop(index)
            if partition.name in self.validation_status:
                del self.validation_status[partition.name]
    
    def update_partition(self, index: int, partition: PartitionDefinition) -> bool:
        """更新分区"""
        if 0 <= index < len(self.partitions):
            # 验证分区
            valid, error = partition.validate()
            self.validation_status[partition.name] = (valid, error)
            
            if not valid:
                return False
            
            self.partitions[index] = partition
            return True
        return False
    
    def get_partition(self, index: int) -> Optional[PartitionDefinition]:
        """获取分区"""
        if 0 <= index < len(self.partitions):
            return self.partitions[index]
        return None
    
    def clear(self):
        """清空所有分区"""
        self.partitions.clear()
        self.validation_status.clear()
    
    def row_count(self) -> int:
        """返回行数"""
        return len(self.partitions)
    
    def column_count(self) -> int:
        """返回列数"""
        return len(self.headers)


class ModelFinderPartitionDialog(QDialog):
    """ModelFinder分区模式对话框"""
    
    # 定义信号
    partitions_selected = pyqtSignal(list)  # 分区定义已选择信号
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("ModelFinder Partition Mode")
        self.setMinimumSize(900, 700)
        
        # 初始化变量
        self.table_model = PartitionTableModel()
        self.current_partition_mode = PartitionMode.EL  # 默认EL模式
        
        # 初始化UI
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
            "Topo-unlinked (-S)"
        ])
        self.mode_combo.setCurrentIndex(0)  # Default EL mode
        self.mode_combo.currentIndexChanged.connect(self.on_mode_changed)
        mode_layout.addRow("Mode:", self.mode_combo)
        
        # 分区速率选项
        self.rcluster_checkbox = QCheckBox("Use hierarchical clustering (--rcluster)")
        self.rcluster_checkbox.setChecked(False)
        self.rcluster_percent_spinbox = QSpinBox()
        self.rcluster_percent_spinbox.setRange(1, 100)
        self.rcluster_percent_spinbox.setValue(10)
        self.rcluster_percent_spinbox.setSuffix("%")
        self.rcluster_percent_spinbox.setEnabled(False)
        self.rcluster_checkbox.stateChanged.connect(
            lambda state: self.rcluster_percent_spinbox.setEnabled(state == Qt.Checked)
        )
        
        rcluster_layout = QHBoxLayout()
        rcluster_layout.addWidget(self.rcluster_checkbox)
        rcluster_layout.addWidget(self.rcluster_percent_spinbox)
        mode_layout.addRow("", rcluster_layout)
        
        layout.addWidget(mode_group)
        
        # 分区定义组
        partition_group = QGroupBox("Partition Definitions")
        partition_layout = QVBoxLayout()
        partition_group.setLayout(partition_layout)
        
        # 分区表格
        self.partition_table = QTableWidget()
        self.partition_table.setColumnCount(4)
        self.partition_table.setHorizontalHeaderLabels(self.table_model.headers)
        self.partition_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.partition_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.partition_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.partition_table.setAlternatingRowColors(True)
        partition_layout.addWidget(self.partition_table)
        
        # 表格按钮
        btn_layout = QHBoxLayout()
        
        self.add_partition_btn = QPushButton("+ Add")
        self.add_partition_btn.clicked.connect(self.add_partition)
        
        self.edit_partition_btn = QPushButton("Edit")
        self.edit_partition_btn.clicked.connect(self.edit_partition)
        self.edit_partition_btn.setEnabled(False)
        
        self.remove_partition_btn = QPushButton("- Remove")
        self.remove_partition_btn.clicked.connect(self.remove_partition)
        self.remove_partition_btn.setEnabled(False)
        
        self.import_partition_btn = QPushButton("Import")
        self.import_partition_btn.clicked.connect(self.import_partition_file)
        
        self.export_partition_btn = QPushButton("Export")
        self.export_partition_btn.clicked.connect(self.export_partition_file)
        self.export_partition_btn.setEnabled(False)
        
        self.clear_partitions_btn = QPushButton("Clear")
        self.clear_partitions_btn.clicked.connect(self.clear_partitions)
        self.clear_partitions_btn.setEnabled(False)
        
        btn_layout.addWidget(self.add_partition_btn)
        btn_layout.addWidget(self.edit_partition_btn)
        btn_layout.addWidget(self.remove_partition_btn)
        btn_layout.addWidget(self.import_partition_btn)
        btn_layout.addWidget(self.export_partition_btn)
        btn_layout.addWidget(self.clear_partitions_btn)
        btn_layout.addStretch()
        
        partition_layout.addLayout(btn_layout)
        layout.addWidget(partition_group)
        
        # 连接表格选择信号
        self.partition_table.itemSelectionChanged.connect(self.on_selection_changed)
        
        # 底部按钮
        control_layout = QHBoxLayout()
        control_layout.addStretch()
        
        self.ok_btn = QPushButton("OK")
        self.ok_btn.clicked.connect(self.accept_partitions)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        
        control_layout.addWidget(self.ok_btn)
        control_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(control_layout)
        
        # 刷新表格
        self.refresh_table()
    
    def set_partitions(self, partitions: List[PartitionDefinition]):
        """设置分区定义"""
        self.table_model.clear()
        for partition in partitions:
            self.table_model.add_partition(partition)
        self.refresh_table()
    
    def set_partition_mode(self, mode: PartitionMode):
        """设置分区模式"""
        self.current_partition_mode = mode
        
        # Update dropdown
        mode_map = {
            PartitionMode.EL: 0,
            PartitionMode.TL: 1,
            PartitionMode.EUL: 2,
            PartitionMode.TUL: 3
        }
        index = mode_map.get(mode, 0)
        self.mode_combo.setCurrentIndex(index)
    
    def get_partitions(self) -> List[PartitionDefinition]:
        """获取分区定义"""
        return self.table_model.partitions.copy()
    
    def get_partition_mode(self) -> PartitionMode:
        """获取分区模式"""
        return self.current_partition_mode
    
    def get_rcluster_enabled(self) -> bool:
        """获取是否启用RCluster"""
        return self.rcluster_checkbox.isChecked()
    
    def get_rcluster_percent(self) -> Optional[int]:
        """获取RCluster百分比"""
        if self.rcluster_checkbox.isChecked():
            return self.rcluster_percent_spinbox.value()
        return None
    
    def accept_partitions(self):
        """接受分区配置"""
        # 验证分区定义
        if self.table_model.row_count() == 0:
            reply = QMessageBox.question(
                self,
                "Confirm",
                "No partitions defined. Continue with partition mode disabled?",
                QMessageBox.Yes | QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        # 发送信号
        self.partitions_selected.emit(self.table_model.partitions.copy())
        
        # 接受对话框
        self.accept()
    
    def on_mode_changed(self, index: int):
        """Mode change handler"""
        mode_names = [PartitionMode.EL, PartitionMode.TL, PartitionMode.EUL, PartitionMode.TUL]
        self.current_partition_mode = mode_names[index]
    
    def on_selection_changed(self):
        """选择改变处理"""
        has_selection = len(self.partition_table.selectedItems()) > 0
        self.edit_partition_btn.setEnabled(has_selection)
        self.remove_partition_btn.setEnabled(has_selection)
    
    def refresh_table(self):
        """刷新表格显示"""
        self.partition_table.setRowCount(self.table_model.row_count())
        
        for row in range(self.table_model.row_count()):
            partition = self.table_model.get_partition(row)
            
            # Name
            name_item = QTableWidgetItem(partition.name)
            self.partition_table.setItem(row, 0, name_item)
            
            # File
            file_item = QTableWidgetItem(partition.file_path if partition.file_path else "Same file")
            self.partition_table.setItem(row, 1, file_item)
            
            # Model Range
            range_item = QTableWidgetItem(partition.get_display_range())
            self.partition_table.setItem(row, 2, range_item)
            
            # Selected Model
            model_item = QTableWidgetItem(partition.selected_model if partition.selected_model else "Auto")
            self.partition_table.setItem(row, 3, model_item)
        
        # 更新按钮状态
        has_partitions = self.table_model.row_count() > 0
        self.export_partition_btn.setEnabled(has_partitions)
        self.clear_partitions_btn.setEnabled(has_partitions)
    
    def add_partition(self):
        """添加分区"""
        dialog = PartitionEditDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            partition = dialog.get_partition()
            if self.table_model.add_partition(partition):
                self.refresh_table()
            else:
                valid, error = self.table_model.validation_status.get(partition.name, (False, "Unknown error"))
                QMessageBox.warning(self, "Validation Error", f"Invalid partition: {error}")
    
    def edit_partition(self):
        """编辑分区"""
        selected_rows = self.partition_table.selectionModel().selectedRows()
        if not selected_rows:
            return
        
        row = selected_rows[0].row()
        partition = self.table_model.get_partition(row)
        
        if partition:
            dialog = PartitionEditDialog(self, partition)
            if dialog.exec_() == QDialog.Accepted:
                updated_partition = dialog.get_partition()
                if self.table_model.update_partition(row, updated_partition):
                    self.refresh_table()
                else:
                    valid, error = self.table_model.validation_status.get(updated_partition.name, (False, "Unknown error"))
                    QMessageBox.warning(self, "Validation Error", f"Invalid partition: {error}")
    
    def remove_partition(self):
        """移除分区"""
        selected_rows = self.partition_table.selectionModel().selectedRows()
        if not selected_rows:
            return
        
        row = selected_rows[0].row()
        reply = QMessageBox.question(
            self, "Confirm Removal",
            f"Are you sure you want to remove partition '{self.table_model.get_partition(row).name}'?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            self.table_model.remove_partition(row)
            self.refresh_table()
    
    def import_partition_file(self):
        """导入分区文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Import Partition File", "",
            "Partition files (*.nex *.nexus *.partitions);;All files (*)"
        )
        
        if file_path:
            try:
                valid, error, partitions = PartitionFileIO.validate_partition_file(file_path)
                
                if not valid:
                    QMessageBox.critical(self, "Import Error", error)
                    return
                
                # 询问是否清空现有分区
                if self.table_model.row_count() > 0:
                    reply = QMessageBox.question(
                        self, "Clear Existing Partitions",
                        "Do you want to clear existing partitions before importing?",
                        QMessageBox.Yes | QMessageBox.No
                    )
                    if reply == QMessageBox.Yes:
                        self.table_model.clear()
                
                # 添加导入的分区
                failed_partitions = []
                for partition in partitions:
                    if not self.table_model.add_partition(partition):
                        valid, error = self.table_model.validation_status.get(partition.name, (False, "Unknown error"))
                        failed_partitions.append(f"{partition.name}: {error}")
                
                self.refresh_table()
                
                if failed_partitions:
                    warning_msg = "Some partitions failed to import:\n" + "\n".join(failed_partitions)
                    QMessageBox.warning(self, "Import Warning", warning_msg)
                else:
                    QMessageBox.information(self, "Import Success", f"Successfully imported {len(partitions)} partitions")
                    
            except Exception as e:
                QMessageBox.critical(self, "Import Error", f"Failed to import partition file: {str(e)}")
    
    def export_partition_file(self):
        """导出分区文件"""
        if self.table_model.row_count() == 0:
            QMessageBox.warning(self, "Export Error", "No partitions to export")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Partition File", "",
            "NEXUS format (*.nex *.nexus);;RAxML format (*.partitions);;All files (*)"
        )
        
        if file_path:
            try:
                # 根据文件扩展名选择格式
                file_ext = os.path.splitext(file_path)[1].lower()
                
                if file_ext in ['.nex', '.nexus']:
                    PartitionFileIO.export_nexus_partition(self.table_model.partitions, file_path)
                else:
                    PartitionFileIO.export_raxml_partition(self.table_model.partitions, file_path)
                
                QMessageBox.information(self, "Export Success", f"Successfully exported {self.table_model.row_count()} partitions")
                
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to export partition file: {str(e)}")
    
    def clear_partitions(self):
        """清空所有分区"""
        reply = QMessageBox.question(
            self, "Confirm Clear",
            "Are you sure you want to clear all partitions?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            self.table_model.clear()
            self.refresh_table()


class PartitionEditDialog(QDialog):
    """分区编辑对话框"""
    
    def __init__(self, parent=None, partition: Optional[PartitionDefinition] = None):
        super().__init__(parent)
        self.setWindowTitle("Edit Partition" if partition else "Add Partition")
        self.setMinimumSize(500, 300)
        
        self.partition = partition
        self.init_ui()
    
    def init_ui(self):
        """初始化UI"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        form_layout = QFormLayout()
        
        # 分区名称
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("Enter partition name...")
        form_layout.addRow("Name:", self.name_edit)
        
        # 文件路径
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Leave empty for same file...")
        self.file_path_edit.setReadOnly(True)
        self.browse_file_btn = QPushButton("Browse...")
        self.browse_file_btn.clicked.connect(self.browse_file)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.browse_file_btn)
        form_layout.addRow("File:", file_layout)
        
        # 模型范围
        self.model_range_edit = QLineEdit()
        self.model_range_edit.setPlaceholderText("e.g., 1-1000 or gene1.fas:1-1000")
        form_layout.addRow("Model Range:", self.model_range_edit)
        
        # 选定模型
        self.selected_model_edit = QLineEdit()
        self.selected_model_edit.setPlaceholderText("Leave empty for auto model selection...")
        form_layout.addRow("Selected Model:", self.selected_model_edit)
        
        layout.addLayout(form_layout)
        
        # 如果是编辑模式，填充现有数据
        if self.partition:
            self.name_edit.setText(self.partition.name)
            self.file_path_edit.setText(self.partition.file_path)
            self.model_range_edit.setText(self.partition.model_range)
            if self.partition.selected_model:
                self.selected_model_edit.setText(self.partition.selected_model)
        
        # 底部按钮
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        
        self.ok_btn = QPushButton("OK")
        self.ok_btn.clicked.connect(self.validate_and_accept)
        
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        
        btn_layout.addWidget(self.ok_btn)
        btn_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(btn_layout)
    
    def browse_file(self):
        """浏览文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Sequence File", "",
            "Sequence files (*.fas *.fasta *.fa *.phy);;All files (*)"
        )
        
        if file_path:
            self.file_path_edit.setText(file_path)
    
    def validate_and_accept(self):
        """验证并接受"""
        # 创建分区对象
        partition = PartitionDefinition(
            name=self.name_edit.text().strip(),
            file_path=self.file_path_edit.text().strip(),
            seq_type="DNA",  # 默认使用DNA类型，seq_type由全局参数控制
            model_range=self.model_range_edit.text().strip(),
            selected_model=self.selected_model_edit.text().strip() or None
        )
        
        # 验证
        valid, error = partition.validate()
        
        if not valid:
            QMessageBox.warning(self, "Validation Error", error)
            return
        
        self.partition = partition
        self.accept()
    
    def get_partition(self) -> PartitionDefinition:
        """获取分区"""
        return self.partition
