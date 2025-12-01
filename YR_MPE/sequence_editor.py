# sequence_editor.py
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

import os
import re
import json
from typing import Dict, List, Optional, Tuple, Union
from collections import OrderedDict
from dataclasses import dataclass
from enum import Enum

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

try:
    from .sequence_models import SequenceModel, SequenceItem
    from .sequence_commands import *
    from .sequence_highlighter import SequenceHighlighter
    from .sequence_utils import SequenceUtils, SequenceValidator
except ImportError:
    # Fallback for direct execution
    from sequence_models import SequenceModel, SequenceItem
    from sequence_commands import *
    from sequence_highlighter import SequenceHighlighter
    from sequence_utils import SequenceUtils, SequenceValidator


class ToolProgressDialog(QDialog):
    """工具运行进度对话框"""
    def __init__(self, tool_name, parent=None):
        super().__init__(parent)
        self.tool_name = tool_name
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle(f'{self.tool_name} Progress')
        self.setModal(True)
        self.resize(600, 400)
        
        layout = QVBoxLayout()
        
        # 工具名称标签
        self.tool_label = QLabel(f'Running {self.tool_name}...')
        self.tool_label.setAlignment(Qt.AlignCenter)
        font = self.tool_label.font()
        font.setPointSize(12)
        font.setBold(True)
        self.tool_label.setFont(font)
        layout.addWidget(self.tool_label)
        
        # 进度条
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)  # Indeterminate progress
        layout.addWidget(self.progress_bar)
        
        # 控制台输出文本框
        self.console_output = QTextEdit()
        self.console_output.setReadOnly(True)
        self.console_output.setFont(QFont("Consolas", 9))
        layout.addWidget(self.console_output)
        
        # 按钮布局
        button_layout = QHBoxLayout()
        self.close_button = QPushButton('Close')
        self.close_button.clicked.connect(self.accept)
        self.close_button.setEnabled(False)  # 运行完成前禁用
        button_layout.addStretch()
        button_layout.addWidget(self.close_button)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
        
    def append_output(self, text):
        """添加输出文本"""
        self.console_output.append(text)
        self.console_output.verticalScrollBar().setValue(
            self.console_output.verticalScrollBar().maximum()
        )
        
    def finish_run(self, success=True):
        """完成运行"""
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)
        if success:
            self.tool_label.setText(f'{self.tool_name} completed successfully!')
        else:
            self.tool_label.setText(f'{self.tool_name} failed!')
        self.close_button.setEnabled(True)


class SequenceAlignmentViewer(QMainWindow):
    # 定义信号，用于将保存的序列数据发送回主平台
    sequences_saved = pyqtSignal(list, str)  # sequences, file_path
    # 定义信号，用于将比对结果导入到主平台
    import_alignment_signal = pyqtSignal(list)
    
    def __init__(self, sequences = []):
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        super().__init__()
        self.initUI()
        self.sequences = sequences  # 存储序列数据
        if sequences:
            self.display_alignment()
        # 初始化时更新序列列表
        self.update_sequences_list()
        
    def initUI(self):
        # 设置窗口属性
        self.setWindowTitle('Sequence Alignment Viewer')
        self.setGeometry(100, 100, 1000, 700)
        
        # 创建中央部件和布局
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)
        
        # 创建菜单栏
        self.create_menu_bar()
        
        # 创建工具栏
        self.create_toolbar()
        
        # 创建显示区域 - 使用分割器包含两个表格
        self.splitter = QSplitter(Qt.Horizontal)
        
        # 创建序列名称表格（固定在左侧）
        self.font = QFont("Courier New", 10)

        self.headers_table = QTableWidget(0, 1)
        self.headers_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.headers_table.horizontalHeader().setVisible(False)
        self.headers_table.verticalHeader().setVisible(False)
        self.headers_table.setFont(self.font)
        self.headers_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        
        # 创建序列内容表格
        self.alignment_table = QTableWidget(0, 0)
        self.alignment_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.alignment_table.horizontalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.alignment_table.verticalHeader().setVisible(False)
        self.alignment_table.setFont(self.font)
        self.alignment_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        
        # 添加表格到分割器
        self.splitter.addWidget(self.headers_table)
        self.splitter.addWidget(self.alignment_table)
        
        # 设置分割器初始大小
        self.splitter.setSizes([200, 800])
        
        # 添加到布局
        self.layout.addWidget(self.splitter)
        
        # 连接两个表格的垂直滚动条
        self.headers_table.verticalScrollBar().valueChanged.connect(
            self.alignment_table.verticalScrollBar().setValue)
        self.alignment_table.verticalScrollBar().valueChanged.connect(
            self.headers_table.verticalScrollBar().setValue)
        
        # 创建状态栏
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_label = QLabel("Ready")
        self.status_bar.addPermanentWidget(self.status_label)
        
        # 添加序列选择和位置选择控件
        self.create_selection_controls()
        
        # 初始化变量
        self.current_file = None
        
    def create_selection_controls(self):
        """创建序列和位置选择控件"""
        # 序列选择组合框
        self.sequence_combo = QComboBox()
        self.sequence_combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.sequence_combo.setMinimumWidth(100)
        self.sequence_combo.setMaximumWidth(200)
        
        # 列索引输入框
        self.column_input = QLineEdit()
        self.column_input.setPlaceholderText("Column")
        self.column_input.setMaximumWidth(80)
        
        # 行索引输入框
        self.row_input = QLineEdit()
        self.row_input.setPlaceholderText("Row")
        self.row_input.setMaximumWidth(80)
        
        # 连接信号
        self.sequence_combo.currentIndexChanged.connect(self.on_sequence_changed)
        self.column_input.editingFinished.connect(self.on_position_changed)
        self.row_input.editingFinished.connect(self.on_position_changed)
        
        # 添加到状态栏
        self.status_bar.addPermanentWidget(QLabel("Sequence:"))
        self.status_bar.addPermanentWidget(self.sequence_combo)
        self.status_bar.addPermanentWidget(QLabel("Column:"))
        self.status_bar.addPermanentWidget(self.column_input)
        self.status_bar.addPermanentWidget(QLabel("Row:"))
        self.status_bar.addPermanentWidget(self.row_input)
        
    def on_sequence_changed(self, index):
        """当序列选择改变时"""
        if index >= 0:
            # 获取当前选中的行列
            current_row = self.alignment_table.currentRow()
            current_col = self.alignment_table.currentColumn()
            
            # 如果当前没有选中单元格，则默认选中第一列
            if current_row < 0:
                current_row = 0
            if current_col < 0:
                current_col = 0
                
            # 更新表格选择
            self.alignment_table.setCurrentCell(index, current_col)
            
    def on_position_changed(self):
        """当位置输入改变时"""
        try:
            # 获取输入的行列值
            row_text = self.row_input.text()
            col_text = self.column_input.text()
            
            if row_text and col_text:
                # 转换为整数并调整为0基索引
                row = int(row_text) - 1
                col = int(col_text) - 1
                
                # 检查范围有效性
                if 0 <= row < self.alignment_table.rowCount() and 0 <= col < self.alignment_table.columnCount():
                    # 选择对应的单元格
                    self.alignment_table.setCurrentCell(row, col)
                    # 确保单元格可见
                    self.alignment_table.scrollToItem(self.alignment_table.item(row, col), QAbstractItemView.PositionAtCenter)
        except ValueError:
            # 输入不是有效数字，忽略
            pass
        
    def update_sequences_list(self):
        """更新序列列表"""
        self.sequence_combo.clear()
        for i, seq_data in enumerate(self.sequences):
            self.sequence_combo.addItem(f"{i+1}. {seq_data['header']}", i)
        
    def create_menu_bar(self):
        # 文件菜单
        file_menu = self.menuBar().addMenu('&File')
        
        open_action = QAction('&Open', self)
        open_action.setShortcut('Ctrl+O')
        open_action.setStatusTip('Open sequence alignment file')
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)
        
        # 添加 Insert sequences from file 操作
        self.insert_action = QAction('&Insert Sequences from File', self)
        self.insert_action.setStatusTip('Insert sequences from another alignment file')
        self.insert_action.triggered.connect(self.insert_sequences)
        self.insert_action.setEnabled(False)  # 初始时禁用，直到加载文件
        file_menu.addAction(self.insert_action)
        
        # 添加 Save changes 操作
        self.save_action = QAction('&Save Changes', self)
        self.save_action.setShortcut('Ctrl+S')
        self.save_action.setStatusTip('Save changes to the current file')
        self.save_action.triggered.connect(self.save_changes)
        self.save_action.setEnabled(False)  # 初始时禁用，直到加载文件
        file_menu.addAction(self.save_action)
        
        # 添加 Save as 操作
        self.save_as_action = QAction('Save &As...', self)
        self.save_as_action.setShortcut('Ctrl+Shift+S')
        self.save_as_action.setStatusTip('Save changes to a new file')
        self.save_as_action.triggered.connect(self.save_as)
        self.save_as_action.setEnabled(False)  # 初始时禁用，直到加载文件
        file_menu.addAction(self.save_as_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction('&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # 对齐菜单
        alignment_menu = self.menuBar().addMenu('&Alignment')
        
        mafft_action = QAction('MAFFT', self)
        mafft_action.setStatusTip('Run MAFFT alignment')
        mafft_action.triggered.connect(lambda: self.run_alignment_tool('MAFFT'))
        alignment_menu.addAction(mafft_action)
        
        clustal_omega_action = QAction('Clustal Omega', self)
        clustal_omega_action.setStatusTip('Run Clustal Omega alignment')
        clustal_omega_action.triggered.connect(lambda: self.run_alignment_tool('Clustal Omega'))
        alignment_menu.addAction(clustal_omega_action)
        
        muscle5_action = QAction('Muscle5', self)
        muscle5_action.setStatusTip('Run Muscle5 alignment')
        muscle5_action.triggered.connect(lambda: self.run_alignment_tool('Muscle5'))
        alignment_menu.addAction(muscle5_action)
        
        # 修剪菜单
        trim_menu = self.menuBar().addMenu('&Trim')
        
        trimal_action = QAction('TrimAl', self)
        trimal_action.setStatusTip('Run TrimAl trimming')
        trimal_action.triggered.connect(lambda: self.run_trim_tool('TrimAl'))
        trim_menu.addAction(trimal_action)
        
        gblocks_action = QAction('GBlocks', self)
        gblocks_action.setStatusTip('Run GBlocks trimming')
        gblocks_action.triggered.connect(lambda: self.run_trim_tool('GBlocks'))
        trim_menu.addAction(gblocks_action)
        
        # 视图菜单
        view_menu = self.menuBar().addMenu('&View')
        
        font_action = QAction('&Font', self)
        font_action.setStatusTip('Change font')
        font_action.triggered.connect(self.change_font)
        view_menu.addAction(font_action)
        
    def create_toolbar(self):
        toolbar = QToolBar('Main Toolbar')
        self.addToolBar(toolbar)

        # get default action height and apply to icon        
        open_action = QAction('Open', self)
        font = open_action.font()
        height = QFontMetrics(font).height()
        toolbar.setIconSize(QSize(height, height))


        open_action.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/open.svg')))
        open_action.setStatusTip('Open sequence alignment file')
        open_action.triggered.connect(self.open_file)

        toolbar.addAction(open_action)
        
        # 添加插入序列按钮
        self.insert_action_toolbar = QAction('Insert', self)
        self.insert_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/insert.svg')))
        self.insert_action_toolbar.setStatusTip('Insert sequences from another alignment file')
        self.insert_action_toolbar.triggered.connect(self.insert_sequences)
        self.insert_action_toolbar.setEnabled(False)
        toolbar.addAction(self.insert_action_toolbar)
        
        # 添加保存按钮
        self.save_action_toolbar = QAction('Save', self)
        self.save_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/save.svg')))
        self.save_action_toolbar.setStatusTip('Save changes to the current file')
        self.save_action_toolbar.triggered.connect(self.save_changes)
        self.save_action_toolbar.setEnabled(False)
        toolbar.addAction(self.save_action_toolbar)
        
        self.save_as_action_toolbar = QAction('Save As', self)
        self.save_as_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/export.svg')))
        self.save_as_action_toolbar.setStatusTip('Save changes to a new file')
        self.save_as_action_toolbar.triggered.connect(self.save_as)
        self.save_as_action_toolbar.setEnabled(False)
        toolbar.addAction(self.save_as_action_toolbar)
        
        toolbar.addSeparator()
        
        # 添加对齐工具按钮组
        alignment_toolbar = QToolBar('Alignment Tools')
        self.addToolBar(Qt.TopToolBarArea, alignment_toolbar)
        alignment_toolbar.setWindowTitle('Alignment')
        
        mafft_action_toolbar = QAction('MAFFT', self)
        mafft_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/software/mafft.svg')))
        mafft_action_toolbar.setStatusTip('Run MAFFT alignment')
        mafft_action_toolbar.triggered.connect(lambda: self.run_alignment_tool('MAFFT'))
        alignment_toolbar.addAction(mafft_action_toolbar)
        
        clustal_omega_action_toolbar = QAction('Clustal Omega', self)
        clustal_omega_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/software/clustalo.svg')))
        clustal_omega_action_toolbar.setStatusTip('Run Clustal Omega alignment')
        clustal_omega_action_toolbar.triggered.connect(lambda: self.run_alignment_tool('Clustal Omega'))
        alignment_toolbar.addAction(clustal_omega_action_toolbar)
        
        muscle5_action_toolbar = QAction('Muscle5', self)
        muscle5_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/software/muscle.svg')))
        muscle5_action_toolbar.setStatusTip('Run Muscle5 alignment')
        muscle5_action_toolbar.triggered.connect(lambda: self.run_alignment_tool('Muscle5'))
        alignment_toolbar.addAction(muscle5_action_toolbar)
        
        # 添加修剪工具按钮组
        trim_toolbar = QToolBar('Trim Tools')
        self.addToolBar(Qt.TopToolBarArea, trim_toolbar)
        trim_toolbar.setWindowTitle('Trim')
        
        trimal_action_toolbar = QAction('TrimAl', self)
        trimal_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/software/trimal.svg')))
        trimal_action_toolbar.setStatusTip('Run TrimAl trimming')
        trimal_action_toolbar.triggered.connect(lambda: self.run_trim_tool('TrimAl'))
        trim_toolbar.addAction(trimal_action_toolbar)
        
        gblocks_action_toolbar = QAction('GBlocks', self)
        gblocks_action_toolbar.setIcon(QIcon(os.path.join(self.plugin_path, 'icons/software/gblocks.svg')))
        gblocks_action_toolbar.setStatusTip('Run GBlocks trimming')
        gblocks_action_toolbar.triggered.connect(lambda: self.run_trim_tool('GBlocks'))
        trim_toolbar.addAction(gblocks_action_toolbar)
        
    def open_file(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(
            self, 'Open Sequence Alignment File', '', 
            'FASTA Files (*.fasta *.fa *.aln *.fas *.fna);;All Files (*)', 
            options=options)
        
        if file_name:
            self.load_alignment(file_name)
            
    def load_alignment(self, file_path):
        """加载序列比对文件"""
        try:
            with open(file_path, 'r') as file:
                content = file.read()
                
            # 解析FASTA格式文件
            self.parse_fasta(content)
            self.display_alignment()
            self.current_file = file_path
            self.status_label.setText(f"Loaded: {file_path[0:20]+'...' if len(file_path) > 20 else file_path}")
            
            # 启用保存和插入操作
            self.save_action.setEnabled(True)
            self.save_as_action.setEnabled(True)
            self.save_action_toolbar.setEnabled(True)
            self.save_as_action_toolbar.setEnabled(True)
            self.insert_action.setEnabled(True)
            self.insert_action_toolbar.setEnabled(True)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file:\n{str(e)}")
            
    def parse_fasta(self, content):
        """解析FASTA格式的内容"""
        self.sequences = []
        lines = content.strip().split('\n')
        
        current_header = ""
        current_sequence = ""
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # 如果已有上一个序列，则保存它
                if current_header and current_sequence:
                    self.sequences.append({
                        'header': current_header,
                        'sequence': current_sequence.upper()
                    })
                # 开始新的序列
                current_header = line[1:]  # 去掉 '>' 符号
                current_sequence = ""
            else:
                # 继续累积序列数据
                current_sequence += line.replace(' ', '').replace('\t', '')
                
        # 不要忘记最后一个序列
        if current_header and current_sequence:
            self.sequences.append({
                'header': current_header,
                'sequence': current_sequence.upper()
            })
            
    def display_alignment(self):
        """在表格中显示序列比对"""
        if not self.sequences:
            return
            
        # 清空表格
        self.headers_table.clear()
        self.alignment_table.clear()
        
        # 获取最长序列的长度作为列数
        max_length = max(len(seq['sequence']) for seq in self.sequences) if self.sequences else 0
        
        # 设置行列数 (左侧表格增加一行用于标题)
        self.headers_table.setRowCount(len(self.sequences) + 1)  # +1 用于标题行
        self.headers_table.setColumnCount(1)
        
        self.alignment_table.setRowCount(len(self.sequences))
        self.alignment_table.setColumnCount(max_length)
        
        # 检查每个位点的保守性
        conservation_marks = self.check_conservation()
        
        # 设置序列内容表格的表头 (包含保守性标记)
        headers = []
        for i in range(max_length):
            if conservation_marks[i] == 2:
                mark = "*"
            elif conservation_marks[i] == 1:
                mark = "."
            else:
                mark = " "
            headers.append(mark)
        self.alignment_table.setHorizontalHeaderLabels(headers)
        
        # 添加左上角的标题单元格
        index_header_item = QTableWidgetItem("Index")
        index_header_item.setFlags(index_header_item.flags() ^ Qt.ItemIsEditable)
        index_header_item.setTextAlignment(Qt.AlignCenter)
        index_header_item.setBackground(QColor(220, 220, 220))  # 灰色背景突出显示
        self.headers_table.setItem(0, 0, index_header_item)
        
        # 填充数据
        for row, seq_data in enumerate(self.sequences):
            # 在左侧表格显示序列名称 (row+1 因为第0行是标题)
            header_item = QTableWidgetItem(f"{row+1}. {seq_data['header']}")
            header_item.setFlags(header_item.flags() ^ Qt.ItemIsEditable)
            header_item.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
            self.headers_table.setItem(row + 1, 0, header_item)
            
            # 在右侧表格显示序列内容
            sequence = seq_data['sequence']
            for col, nucleotide in enumerate(sequence):
                item = QTableWidgetItem(nucleotide)
                item.setTextAlignment(Qt.AlignCenter)
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                
                # 根据核苷酸类型着色
                if nucleotide.upper() == 'A':
                    item.setBackground(QColor(255, 150, 150))  # 浅红色
                elif nucleotide.upper() == 'T' or nucleotide.upper() == 'U':
                    item.setBackground(QColor(150, 150, 255))  # 浅蓝色
                elif nucleotide.upper() == 'G':
                    item.setBackground(QColor(150, 255, 150))  # 浅绿色
                elif nucleotide.upper() == 'C':
                    item.setBackground(QColor(255, 255, 150))  # 浅黄色
                elif nucleotide.upper() == '-':
                    item.setBackground(QColor(200, 200, 200))  # 浅灰色，表示缺口
                    
                self.alignment_table.setItem(row, col, item)
                
            # 对于较短的序列，在末尾填充空白
            for col in range(len(sequence), max_length):
                item = QTableWidgetItem("")
                item.setFlags(item.flags() ^ Qt.ItemIsEditable)
                self.alignment_table.setItem(row, col, item)
                
        # 调整序列名称列宽度
        self.headers_table.resizeColumnToContents(0)
        self.headers_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        
        # 设置默认列宽，避免一次性设置过多列宽导致卡顿
        default_column_width = 0
        self.alignment_table.horizontalHeader().setDefaultSectionSize(default_column_width)
        
        # 设置默认行高，让显示更紧凑
        default_row_height = 0
        self.headers_table.verticalHeader().setDefaultSectionSize(default_row_height)
        self.alignment_table.verticalHeader().setDefaultSectionSize(default_row_height)
            
        # 确保标题行高度与内容行一致
        self.headers_table.resizeRowToContents(0)
        
        # 更新序列列表
        self.update_sequences_list()
        
        # 连接表格点击事件
        self.alignment_table.cellClicked.connect(self.on_alignment_cell_clicked)

        # 设置插入 保存 另存为按钮可用
        self.insert_action.setEnabled(True)
        self.save_action.setEnabled(True)
        self.save_as_action.setEnabled(True)
        self.insert_action_toolbar.setEnabled(True)
        self.save_as_action_toolbar.setEnabled(True)
        self.save_action_toolbar.setEnabled(True)
        
    def on_alignment_cell_clicked(self, row, column):
        """处理比对表格单元格点击事件"""
        # 更新序列组合框，阻止触发on_sequence_changed事件
        self.sequence_combo.blockSignals(True)
        self.sequence_combo.setCurrentIndex(row)
        self.sequence_combo.blockSignals(False)
        
        # 更新行列输入框
        self.row_input.setText(str(row + 1))  # 行索引从1开始显示
        self.column_input.setText(str(column + 1))  # 列索引从1开始显示
        
    def on_sequence_changed(self, index):
        """当序列选择改变时"""
        if index >= 0:
            # 获取当前选中的行列
            current_row = self.alignment_table.currentRow()
            current_col = self.alignment_table.currentColumn()
            
            # 如果当前没有选中单元格，则默认选中第一列
            if current_row < 0:
                current_row = 0
            if current_col < 0:
                current_col = 0
                
            # 更新表格选择
            self.alignment_table.setCurrentCell(index, current_col)
            
    def check_conservation(self):
        """检查每个位点的保守性，如果所有序列在该位点碱基一致则返回True"""
        if not self.sequences:
            return []
            
        max_length = max(len(seq['sequence']) for seq in self.sequences)
        conservation_marks = [True] * max_length
        
        # 对每个位点检查保守性
        for col in range(max_length):
            bases = []
            gap = [False for _ in range(max_length)]
            for seq in self.sequences:
                if col < len(seq['sequence']):
                    base = seq['sequence'][col]
                    # 忽略空白
                    if base.strip() and base != '-':
                        bases.append(base)
                    if base == '-':
                        gap[col] = True
            
            # 如果该列所有碱基相同且非空，则标记为保守
            if bases and all(base == bases[0] for base in bases):
                if not gap[col]:
                    conservation_marks[col] = 2
                else:
                    conservation_marks[col] = 1
            else:
                conservation_marks[col] = 0
                
        return conservation_marks
            
    def change_font(self):
        font, ok = QFontDialog.getFont(self.alignment_table.font(), self)
        if ok:
            self.font = font
            self.headers_table.setFont(font)
            self.alignment_table.setFont(font)
            self.status_label.setText("Font changed")

    def save_changes(self):
        """保存更改到当前文件或通过信号发送回平台"""
        if self.current_file:
            # 单独打开时保存到源文件
            self.save_to_file(self.current_file)
        else:
            # 通过YR-MPEA打开时，发送信号将更改输送回YR-MPEA
            if self.sequences:
                # 发送序列数据和一个标识符回主平台
                self.sequences_saved.emit(self.sequences, "modified_sequences")
                self.status_label.setText("Changes sent back to YR-MPEA")
                QMessageBox.information(self, "Success", "Changes have been sent back to YR-MPEA.")
            else:
                QMessageBox.warning(self, "Warning", "No sequences to save.")

    def run_alignment_tool(self, tool_name):
        """运行指定的比对工具"""
        if not self.sequences:
            QMessageBox.warning(self, "Warning", "No sequences loaded. Please load sequences first.")
            return
            
        # 将当前序列转换为FASTA格式
        fasta_content = ""
        for seq_data in self.sequences:
            fasta_content += f">{seq_data['header']}\n{seq_data['sequence']}\n"
            
        # 根据工具名称导入相应的插件并显示对话框
        self.show_alignment_plugin_dialog(tool_name, fasta_content)
        
    def run_trim_tool(self, tool_name):
        """运行指定的修剪工具"""
        if not self.sequences:
            QMessageBox.warning(self, "Warning", "No sequences loaded. Please load sequences first.")
            return
            
        # 将当前序列转换为FASTA格式
        fasta_content = ""
        for seq_data in self.sequences:
            fasta_content += f">{seq_data['header']}\n{seq_data['sequence']}\n"
            
        # 根据工具名称导入相应的插件并显示对话框
        self.show_trim_plugin_dialog(tool_name, fasta_content)
        
    def show_alignment_plugin_dialog(self, tool_name, fasta_content):
        """显示比对插件对话框"""
        # 创建对话框
        dialog = QDialog(self)
        dialog.setMinimumSize(800, 600)
        
        # 设置对话框标题和图标
        if tool_name == 'MAFFT':
            dialog.setWindowTitle("MAFFT - Sequence Viewer")
            from .plugins.mafft_plugin import MAFFTPluginEntry
            plugin_entry = MAFFTPluginEntry()
            dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/mafft.svg")))
        elif tool_name == 'Clustal Omega':
            dialog.setWindowTitle("Clustal Omega - Sequence Viewer")
            from .plugins.clustal_omega_plugin import ClustalOmegaPluginEntry
            plugin_entry = ClustalOmegaPluginEntry()
            dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/clustalo.svg")))
        elif tool_name == 'Muscle5':
            dialog.setWindowTitle("Muscle5 - Sequence Viewer")
            from .plugins.muscle5_plugin import Muscle5PluginEntry
            plugin_entry = Muscle5PluginEntry()
            dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/muscle.svg")))
        else:
            return
            
        # 设置布局
        dialog.setLayout(QVBoxLayout())
        
        # 将字典格式的序列数据转换为BioPython SeqRecord对象列表
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        records = []
        for seq_data in self.sequences:
            record = SeqRecord(
                Seq(seq_data['sequence']),
                id=seq_data['header'],
                name=seq_data['header'],
                description=""
            )
            records.append(record)
        
        # 创建插件实例
        plugin_instance = plugin_entry.run(import_from="seq_viewer", import_data=records)
        
        # 连接信号，以便接收比对结果
        plugin_instance.import_alignment_signal.connect(self.receive_alignment_result)
        
        # 添加插件到对话框
        dialog.layout().addWidget(plugin_instance)
        
        # 显示对话框
        dialog.exec_()
        
    def show_trim_plugin_dialog(self, tool_name, fasta_content):
        """显示修剪插件对话框"""
        # 创建对话框
        dialog = QDialog(self)
        dialog.setMinimumSize(800, 600)
        
        # 设置对话框标题和图标
        if tool_name == 'TrimAl':
            dialog.setWindowTitle("TrimAl - Sequence Viewer")
            from .plugins.trimal_plugin import TrimAlPluginEntry
            plugin_entry = TrimAlPluginEntry()
            dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/trimal.svg")))
        elif tool_name == 'GBlocks':
            dialog.setWindowTitle("GBlocks - Sequence Viewer")
            from .plugins.gblocks_plugin import GBlocksPluginEntry
            plugin_entry = GBlocksPluginEntry()
            dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/gblocks.svg")))
        else:
            return
            
        # 设置布局
        dialog.setLayout(QVBoxLayout())
        
        # 将字典格式的序列数据转换为BioPython SeqRecord对象列表
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        records = []
        for seq_data in self.sequences:
            record = SeqRecord(
                Seq(seq_data['sequence']),
                id=seq_data['header'],
                name=seq_data['header'],
                description=""
            )
            records.append(record)
        
        # 创建插件实例
        plugin_instance = plugin_entry.run(import_from="seq_viewer", import_data=records)
        
        # 连接信号，以便接收修剪结果
        plugin_instance.import_alignment_signal.connect(self.receive_alignment_result)
        
        # 添加插件到对话框
        dialog.layout().addWidget(plugin_instance)
        
        # 显示对话框
        dialog.exec_()
        
    def receive_alignment_result(self, alignment_data):
        """接收比对结果"""
        self.status_label.setText("Alignment result received")
        # alignment_data是从插件返回的序列列表
        if isinstance(alignment_data, list) and len(alignment_data) > 0:
            # 更新序列数据
            self.sequences = []
            for seq_record in alignment_data:
                self.sequences.append({
                    'header': seq_record.id,
                    'sequence': str(seq_record.seq)
                })
            
            # 刷新显示
            self.display_alignment()
            
            # 更新状态
            self.status_label.setText(f"Loaded {len(self.sequences)} sequences from alignment result")
            
            # 启用保存和插入操作
            self.save_action.setEnabled(True)
            self.save_as_action.setEnabled(True)
            self.save_action_toolbar.setEnabled(True)
            self.save_as_action_toolbar.setEnabled(True)
            self.insert_action.setEnabled(True)
            self.insert_action_toolbar.setEnabled(True)
            
    def save_as(self):
        """另存为新文件"""
        if not self.sequences:
            QMessageBox.warning(self, "Warning", "No sequences to save.")
            return
            
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getSaveFileName(
            self, 'Save Sequence Alignment File', '', 
            'FASTA Files (*.fasta *.fa *.aln *.fas *.fna);;All Files (*)', 
            options=options)
        
        if file_name:
            self.save_to_file(file_name)
            self.current_file = file_name
            self.status_label.setText(f"Saved: {file_name[0:20]+'...' if len(file_name) > 20 else file_name}")
            
            # 启用保存操作
            self.save_action.setEnabled(True)
            self.save_as_action.setEnabled(True)
            self.save_action_toolbar.setEnabled(True)
            self.save_as_action_toolbar.setEnabled(True)
            
    def save_to_file(self, file_path):
        """将序列保存到文件"""
        try:
            with open(file_path, 'w') as file:
                for seq_data in self.sequences:
                    file.write(f">{seq_data['header']}\n")
                    # 按每行80个字符格式化序列
                    sequence = seq_data['sequence']
                    for i in range(0, len(sequence), 80):
                        file.write(sequence[i:i+80] + '\n')
            QMessageBox.information(self, "Success", f"File saved successfully to:\n{file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save file:\n{str(e)}")

    def insert_sequences(self):
        """从文件插入序列到当前比对中"""
        if not self.sequences:
            QMessageBox.warning(self, "Warning", "Please load a sequence alignment file first.")
            return
            
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(
            self, 'Insert Sequences from File', '', 
            'FASTA Files (*.fasta *.fa *.aln *.fas *.fna);;All Files (*)', 
            options=options)
        
        if file_name:
            try:
                with open(file_name, 'r') as file:
                    content = file.read()
                    
                # 解析FASTA格式文件
                new_sequences = []
                lines = content.strip().split('\n')
                
                current_header = ""
                current_sequence = ""
                
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'):
                        # 如果已有上一个序列，则保存它
                        if current_header and current_sequence:
                            new_sequences.append({
                                'header': current_header,
                                'sequence': current_sequence.upper()
                            })
                        # 开始新的序列
                        current_header = line[1:]  # 去掉 '>' 符号
                        current_sequence = ""
                    else:
                        # 继续累积序列数据
                        current_sequence += line.replace(' ', '').replace('\t', '')
                        
                # 不要忘记最后一个序列
                if current_header and current_sequence:
                    new_sequences.append({
                        'header': current_header,
                        'sequence': current_sequence.upper()
                    })
                    
                if new_sequences:
                    # 检查序列长度是否一致
                    lengths = [len(seq['sequence']) for seq in self.sequences]
                    lengths.extend([len(seq['sequence']) for seq in new_sequences])
                    max_len = max(lengths)
                    min_len = min(lengths)
                    
                    # 如果长度不一致，提示用户
                    if max_len != min_len:
                        reply = QMessageBox.question(
                            self, 
                            "Sequence Length Mismatch", 
                            f"Sequence lengths are inconsistent.\n"
                            f"Current alignment: {min_len}-{max_len} characters\n"
                            f"Do you want to pad shorter sequences with gaps '-' to match the longest?",
                            QMessageBox.Yes | QMessageBox.No, 
                            QMessageBox.No
                        )
                        
                        if reply == QMessageBox.Yes:
                            # 用'-'填充所有序列到最大长度
                            for seq in self.sequences:
                                seq['sequence'] = seq['sequence'].ljust(max_len, '-')
                                
                            for seq in new_sequences:
                                seq['sequence'] = seq['sequence'].ljust(max_len, '-')
                    
                    # 将新序列添加到现有序列列表中
                    self.sequences.extend(new_sequences)
                    
                    # 重新显示比对
                    self.display_alignment()
                    
                    self.status_label.setText(f"Inserted {len(new_sequences)} sequences from file.")
                else:
                    QMessageBox.warning(self, "Warning", "No sequences found in the selected file.")
                    
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to insert sequences from file:\n{str(e)}")

class SequenceEditorEntry:
    """序列编辑器入口类"""
    def run(self):
        return SequenceAlignmentViewer()
