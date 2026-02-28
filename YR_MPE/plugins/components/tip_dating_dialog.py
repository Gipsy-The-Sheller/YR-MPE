# tip_dating_dialog.py
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

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QTableWidget, 
                            QTableWidgetItem, QComboBox, QLineEdit, QPushButton,
                            QMessageBox, QHeaderView, QLabel, QWidget, QAbstractItemView)
from PyQt5.QtCore import Qt
from typing import Dict, List, Tuple, Optional


class TipDatingDialog(QDialog):
    """Tip Dating配置对话框"""
    
    def __init__(self, otu_list: List[str], tip_calibrations: Dict = None, parent=None):
        super().__init__(parent)
        self.otu_list = otu_list
        # 如果传入了已有的tip标定数据，则使用它；否则初始化为空
        self.tip_calibrations = tip_calibrations.copy() if tip_calibrations else {}
        self.init_ui()
        self.setup_table()
    
    def init_ui(self):
        """初始化UI"""
        self.setWindowTitle("Configure Tip Dates")
        self.setMinimumSize(800, 600)
        
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # 说明标签
        instruction_label = QLabel(
            "Configure tip dates for each OTU. By default, all tips are uncalibrated. "
            "Select a calibration type and set parameters for each tip you want to calibrate."
        )
        instruction_label.setWordWrap(True)
        instruction_label.setStyleSheet("color: #666; font-style: italic;")
        main_layout.addWidget(instruction_label)
        
        # Tip表格
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(['OTU Name', 'Type', 'Value', 'Modify'])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)  # 禁用直接编辑
        main_layout.addWidget(self.table)
        
        # 底部信息标签
        info_layout = QHBoxLayout()
        info_layout.addStretch()
        
        self.info_label = QLabel(f"Tip Count: {len(self.otu_list)}")
        info_layout.addWidget(self.info_label)
        
        main_layout.addLayout(info_layout)
        
        # 按钮布局
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        self.ok_button = QPushButton("OK")
        self.ok_button.setMinimumWidth(100)
        self.ok_button.clicked.connect(self.accept)
        button_layout.addWidget(self.ok_button)
        
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setMinimumWidth(100)
        self.cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(self.cancel_button)
        
        main_layout.addLayout(button_layout)
    
    def setup_table(self):
        """设置表格"""
        self.table.setRowCount(len(self.otu_list))
        
        for row, otu in enumerate(self.otu_list):
            # OTU Name列（只读）
            name_item = QTableWidgetItem(otu)
            self.table.setItem(row, 0, name_item)
            
            # Type列（ComboBox）
            type_combo = QComboBox()
            type_combo.addItems(['Uncalibrated', 'Point', 'Uniform', 'Upper Boundary', 'Lower Boundary'])
            type_combo.currentTextChanged.connect(lambda text, r=row: self.on_type_changed(r, text))
            self.table.setCellWidget(row, 1, type_combo)
            
            # Value列（LineEdit或只读）
            value_edit = QLineEdit()
            value_edit.setPlaceholderText("Enter value...")
            value_edit.editingFinished.connect(lambda r=row: self.on_value_changed(r))
            self.table.setCellWidget(row, 2, value_edit)
            
            # Modify列（按钮）
            modify_button = QPushButton("Modify...")
            modify_button.clicked.connect(lambda checked, r=row: self.on_modify_clicked(r))
            modify_button.setEnabled(False)  # 初始禁用
            self.table.setCellWidget(row, 3, modify_button)
            
            # 初始化状态
            if otu in self.tip_calibrations and self.tip_calibrations[otu] is not None:
                # 已有标定
                cal_data = self.tip_calibrations[otu]
                display_type = self._lsd2_to_display_type(cal_data['type'])
                type_combo.setCurrentText(display_type)
                
                if display_type == 'Point':
                    value_edit.setText(str(cal_data['values'][0]))
                else:
                    value_edit.setText(self._format_lsd2_value(cal_data['type'], cal_data['values']))
                    value_edit.setEnabled(False)
                    modify_button.setEnabled(True)
            else:
                # 无标定
                type_combo.setCurrentText('Uncalibrated')
                value_edit.setEnabled(False)
                modify_button.setEnabled(False)
    
    def _lsd2_to_display_type(self, lsd2_type: str) -> str:
        """将LSD2类型转换为显示类型"""
        type_mapping = {
            'fixed': 'Point',
            'interval': 'Uniform',
            'upper': 'Upper Boundary',
            'lower': 'Lower Boundary'
        }
        return type_mapping.get(lsd2_type, 'Uncalibrated')
    
    def _display_to_lsd2_type(self, display_type: str) -> str:
        """将显示类型转换为LSD2类型"""
        type_mapping = {
            'Point': 'fixed',
            'Uniform': 'interval',
            'Upper Boundary': 'upper',
            'Lower Boundary': 'lower'
        }
        return type_mapping.get(display_type, None)
    
    def _format_lsd2_value(self, lsd2_type: str, values: List[float]) -> str:
        """格式化LSD2值用于显示"""
        if lsd2_type == 'fixed':
            return str(values[0])
        elif lsd2_type == 'interval':
            return f'b({values[0]},{values[1]})'
        elif lsd2_type == 'upper':
            return f'u({values[0]})'
        elif lsd2_type == 'lower':
            return f'l({values[0]})'
        else:
            return str(values)
    
    def on_type_changed(self, row: int, display_type: str):
        """当类型改变时的回调"""
        value_edit = self.table.cellWidget(row, 2)
        modify_button = self.table.cellWidget(row, 3)
        otu = self.table.item(row, 0).text()
        
        if display_type == 'Uncalibrated':
            # 清除标定
            value_edit.setText("")
            value_edit.setEnabled(False)
            modify_button.setEnabled(False)
            if otu in self.tip_calibrations:
                del self.tip_calibrations[otu]
        elif display_type == 'Point':
            # Point类型，启用Value编辑框
            value_edit.setEnabled(True)
            value_edit.setPlaceholderText("Enter age value...")
            modify_button.setEnabled(False)
            
            # 如果之前有值，保留它
            if otu in self.tip_calibrations and self.tip_calibrations[otu]:
                value_edit.setText(str(self.tip_calibrations[otu]['values'][0]))
        else:
            # 其他类型，禁用Value编辑框，启用Modify按钮
            value_edit.setText("")
            value_edit.setEnabled(False)
            value_edit.setPlaceholderText("Click Modify to configure...")
            modify_button.setEnabled(True)
            
            # 如果之前有值，保留它
            if otu in self.tip_calibrations and self.tip_calibrations[otu]:
                lsd2_type = self.tip_calibrations[otu]['type']
                values = self.tip_calibrations[otu]['values']
                value_edit.setText(self._format_lsd2_value(lsd2_type, values))
    
    def on_value_changed(self, row: int):
        """当值改变时的回调"""
        otu = self.table.item(row, 0).text()
        type_combo = self.table.cellWidget(row, 1)
        value_edit = self.table.cellWidget(row, 2)
        
        display_type = type_combo.currentText()
        
        if display_type == 'Point' and value_edit.text().strip():
            try:
                value = float(value_edit.text())
                self.tip_calibrations[otu] = {
                    'type': 'fixed',
                    'values': [value]
                }
            except ValueError:
                QMessageBox.warning(self, "Invalid Value", "Please enter a valid numeric value.")
                value_edit.setText("")
                if otu in self.tip_calibrations:
                    del self.tip_calibrations[otu]
        elif display_type == 'Point':
            # 清空值
            if otu in self.tip_calibrations:
                del self.tip_calibrations[otu]
    
    def on_modify_clicked(self, row: int):
        """点击Modify按钮时的回调"""
        otu = self.table.item(row, 0).text()
        type_combo = self.table.cellWidget(row, 1)
        
        display_type = type_combo.currentText()
        
        # 打开参数编辑对话框
        dialog = CalibrationParamDialog(display_type, self)
        if dialog.exec_() == QDialog.Accepted:
            values = dialog.get_values()
            
            # 转换为LSD2格式
            lsd2_type = self._display_to_lsd2_type(display_type)
            
            # 更新表格
            value_edit = self.table.cellWidget(row, 2)
            value_edit.setText(self._format_lsd2_value(lsd2_type, values))
            
            # 更新标定数据
            self.tip_calibrations[otu] = {
                'type': lsd2_type,
                'values': values
            }
    
    def get_tip_calibrations(self) -> Dict:
        """获取Tip标定数据"""
        return self.tip_calibrations


class CalibrationParamDialog(QDialog):
    """校准参数编辑对话框"""
    
    def __init__(self, cal_type: str, parent=None):
        super().__init__(parent)
        self.cal_type = cal_type
        self.values = []
        self.init_ui()
    
    def init_ui(self):
        """初始化UI"""
        self.setWindowTitle(f"Configure {self.cal_type} Calibration")
        
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        from PyQt5.QtWidgets import QFormLayout
        
        form_layout = QFormLayout()
        
        if self.cal_type == 'Point':
            self.value_edit = QLineEdit()
            self.value_edit.setPlaceholderText("e.g., 1.5")
            form_layout.addRow("Age (time units):", self.value_edit)
            
        elif self.cal_type == 'Uniform':
            self.lower_edit = QLineEdit()
            self.lower_edit.setPlaceholderText("e.g., 0.5")
            form_layout.addRow("Minimum age:", self.lower_edit)
            
            self.upper_edit = QLineEdit()
            self.upper_edit.setPlaceholderText("e.g., 2.0")
            form_layout.addRow("Maximum age:", self.upper_edit)
            
        elif self.cal_type == 'Upper Boundary':
            self.upper_edit = QLineEdit()
            self.upper_edit.setPlaceholderText("e.g., 3.0")
            form_layout.addRow("Maximum age:", self.upper_edit)
            
        elif self.cal_type == 'Lower Boundary':
            self.lower_edit = QLineEdit()
            self.lower_edit.setPlaceholderText("e.g., 0.2")
            form_layout.addRow("Minimum age:", self.lower_edit)
        
        main_layout.addLayout(form_layout)
        
        # 按钮
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        self.ok_button = QPushButton("OK")
        self.ok_button.clicked.connect(self.validate_and_accept)
        button_layout.addWidget(self.ok_button)
        
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(self.cancel_button)
        
        main_layout.addLayout(button_layout)
    
    def validate_and_accept(self):
        """验证输入并接受"""
        try:
            if self.cal_type == 'Point':
                value = float(self.value_edit.text())
                self.values = [value]
                
            elif self.cal_type == 'Uniform':
                lower = float(self.lower_edit.text())
                upper = float(self.upper_edit.text())
                if lower >= upper:
                    QMessageBox.warning(self, "Invalid Range", "Minimum age must be less than maximum age.")
                    return
                self.values = [lower, upper]
                
            elif self.cal_type == 'Upper Boundary':
                upper = float(self.upper_edit.text())
                self.values = [upper]
                
            elif self.cal_type == 'Lower Boundary':
                lower = float(self.lower_edit.text())
                self.values = [lower]
            
            self.accept()
            
        except ValueError as e:
            QMessageBox.warning(self, "Invalid Input", f"Please enter valid numeric values: {str(e)}")
    
    def get_values(self) -> List[float]:
        """获取参数值"""
        return self.values


def main():
    """测试函数"""
    from PyQt5.QtWidgets import QApplication
    import sys
    
    app = QApplication(sys.argv)
    
    # 测试数据
    otu_list = [f"Taxon_{i}" for i in range(1, 11)]
    
    dialog = TipDatingDialog(otu_list)
    if dialog.exec_() == QDialog.Accepted:
        print("Tip calibrations:")
        for otu, cal in dialog.get_tip_calibrations().items():
            print(f"  {otu}: {cal}")
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()