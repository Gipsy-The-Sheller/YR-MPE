#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Advanced Sequence Editor for YR-MPE
A modern, feature-rich sequence editor for molecular phylogenetics analysis
"""

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


class SequenceAlignmentViewer(QMainWindow):
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
        
        exit_action = QAction('&Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.setStatusTip('Exit application')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
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


class SequenceEditorEntry:
    """序列编辑器入口类"""
    def run(self):
        return SequenceAlignmentViewer()


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    editor = SequenceEditor()
    editor.show()
    sys.exit(app.exec_())
