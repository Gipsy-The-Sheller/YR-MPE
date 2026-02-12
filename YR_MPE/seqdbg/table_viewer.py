#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SeqDBG表格查看器

提供基于表格的基因邻接关系可视化，支持：
- 基因统计表格（gene_name, left_nodes, right_nodes, taxon）
- 颜色着色（left/right_nodes越少越红，taxon越少越蓝）
- 点击绘制桑基图（左邻接 - taxon - 右邻接）
- Taxon管理组件
"""

import logging
import copy
import os
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Tuple

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem,
    QHeaderView, QSplitter, QLabel, QPushButton, QGroupBox, QMessageBox,
    QListWidget, QFileDialog, QTabWidget, QComboBox, QLineEdit, QStyledItemDelegate,
    QRadioButton, QCheckBox, QDialog, QDialogButtonBox, QAbstractItemView, QListWidgetItem
)
from PyQt5.QtCore import Qt, pyqtSignal, QSize
from PyQt5.QtGui import QColor, QFont, QBrush

from .models import GeneStats, GeneInfo
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

def hex_to_rgb(hex_color: str, alpha: float = 1.0) -> Tuple[int, int, int, float]:
    """
    将十六进制颜色转换为RGBA
    
    Args:
        hex_color: 十六进制颜色字符串，如 "#4e699a"
        alpha: 透明度 (0.0-1.0)
        
    Returns:
        RGBA元组 (r, g, b, a)
    """
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return (r, g, b, alpha)


def interpolate_color(
    value: float,
    min_value: float,
    max_value: float
) -> QColor:
    """
    根据值在颜色渐变中插值
    
    颜色渐变：[#ba3e45, #ffffff, #4e699a]
    最小值 -> 红色 (#ba3e45)
    中间值 -> 白色 (#ffffff)
    最大值 -> 蓝色 (#4e699a)
    
    Args:
        value: 当前值
        min_value: 最小值
        max_value: 最大值
        
    Returns:
        QColor对象
    """
    # 如果min和max相等，返回中间色（白色）
    if max_value == min_value:
        return QColor(255, 255, 255, 255)
    
    # 确保min_value <= max_value，如果不满足则交换
    actual_min = min(min_value, max_value)
    actual_max = max(min_value, max_value)
    is_reversed = min_value > max_value
    
    # 归一化到0-1范围
    ratio = (value - actual_min) / (actual_max - actual_min)
    ratio = max(0.0, min(1.0, ratio))
    
    # 如果是反转的，反转ratio
    if is_reversed:
        ratio = 1.0 - ratio
    
    # 定义颜色节点
    red_color = hex_to_rgb("#ba3e45", 1.0)[:3]    # 最小值颜色
    white_color = (255, 255, 255)                  # 中间值颜色
    blue_color = hex_to_rgb("#4e699a", 1.0)[:3]   # 最大值颜色
    
    # 根据比例插值
    if ratio < 0.5:
        # 红色到白色
        local_ratio = ratio * 2  # 0-1
        r = int(red_color[0] + (white_color[0] - red_color[0]) * local_ratio)
        g = int(red_color[1] + (white_color[1] - red_color[1]) * local_ratio)
        b = int(red_color[2] + (white_color[2] - red_color[2]) * local_ratio)
    else:
        # 白色到蓝色
        local_ratio = (ratio - 0.5) * 2  # 0-1
        r = int(white_color[0] + (blue_color[0] - white_color[0]) * local_ratio)
        g = int(white_color[1] + (blue_color[1] - white_color[1]) * local_ratio)
        b = int(white_color[2] + (blue_color[2] - white_color[2]) * local_ratio)
    
    return QColor(r, g, b, 255)


class GeneStatsTable(QTableWidget):
    """
    基因统计表格
    
    显示基因的邻接关系和taxon统计信息，支持颜色着色
    紧凑式设计，参考MiniTracer的风格
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.stats_data: List[GeneStats] = []
        self.max_ln = 0
        self.max_rn = 0
        self.max_lt = 0
        self.max_rt = 0
        self.max_taxon = 0
        
        self.init_ui()
    
    def init_ui(self):
        """初始化表格（紧凑式设计）"""
        # 应用紧凑式样式
        self.setStyleSheet("QTableWidgetItem { padding: 0px; } QTableWidget { gridline-color: lightgray; }")
        
        # 设置列（缩写列名）
        self.setColumnCount(6)
        self.setHorizontalHeaderLabels(["Gene", "LN", "RN", "LT", "RT", "Taxon"])
        
        # 设置列宽
        header = self.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeToContents)
        
        # 设置选择模式
        self.setSelectionBehavior(QTableWidget.SelectRows)
        self.setSelectionMode(QTableWidget.SingleSelection)
        self.setAlternatingRowColors(True)
        self.setSortingEnabled(False)
        self.setEditTriggers(QTableWidget.NoEditTriggers)
        
        # 应用紧凑式字体
        font = QFont()
        font.setPointSize(9)
        self.setFont(font)
        self.horizontalHeader().setFont(font)
        self.verticalHeader().setFont(font)
        self.verticalHeader().setVisible(False)  # 隐藏垂直表头
    
    def set_data(self, stats: List[GeneStats]):
        """
        设置数据
        
        Args:
            stats: GeneStats列表
        """
        self.stats_data = stats
        
        # 计算最大值
        if stats:
            self.max_ln = max(s.left_nodes for s in stats)
            self.max_rn = max(s.right_nodes for s in stats)
            self.max_lt = max(s.left for s in stats)
            self.max_rt = max(s.right for s in stats)
            self.max_taxon = max(s.taxon for s in stats)
        else:
            self.max_ln = 0
            self.max_rn = 0
            self.max_lt = 0
            self.max_rt = 0
            self.max_taxon = 0
        
        # 清空表格
        self.setRowCount(len(stats))
        
        # 填充数据
        for row, stat in enumerate(stats):
            # Gene Name - 根据基因类型着色
            name_item = QTableWidgetItem(stat.gene_name)
            name_item.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
            name_item.setFont(QFont("Arial", 9, QFont.Bold))
            
            # 根据基因类型应用不同的背景颜色
            gene_type = getattr(stat, 'gene_type', 'gene')
            if gene_type == "CDS":
                # CDS: 浅绿色
                name_item.setBackground(QBrush(QColor(180, 220, 180, 100)))
            elif gene_type == "rRNA":
                # rRNA: 浅橙色
                name_item.setBackground(QBrush(QColor(255, 200, 150, 100)))
            elif gene_type == "tRNA":
                # tRNA: 浅紫色
                name_item.setBackground(QBrush(QColor(200, 180, 255, 100)))
            else:
                # gene/other: 浅灰色
                name_item.setBackground(QBrush(QColor(240, 240, 240, 80)))
            
            self.setItem(row, 0, name_item)
            
            # LN (Left Nodes) - 越大越红（反转插值参数）
            ln_item = QTableWidgetItem(str(stat.left_nodes))
            ln_item.setTextAlignment(Qt.AlignCenter)
            if self.max_ln > 0:
                # 反转：最大值传给min_value，0传给max_value，使值越大越红
                bg_color = interpolate_color(stat.left_nodes, self.max_ln, 0)
                ln_item.setBackground(QBrush(bg_color))
            self.setItem(row, 1, ln_item)
            
            # RN (Right Nodes) - 越大越红（反转插值参数）
            rn_item = QTableWidgetItem(str(stat.right_nodes))
            rn_item.setTextAlignment(Qt.AlignCenter)
            if self.max_rn > 0:
                # 反转：最大值传给min_value，0传给max_value，使值越大越红
                bg_color = interpolate_color(stat.right_nodes, self.max_rn, 0)
                rn_item.setBackground(QBrush(bg_color))
            self.setItem(row, 2, rn_item)
            
            # LT (Left Taxon) - 越大越蓝
            lt_item = QTableWidgetItem(str(stat.left))
            lt_item.setTextAlignment(Qt.AlignCenter)
            if self.max_lt > 0:
                bg_color = interpolate_color(stat.left, 0, self.max_lt)
                lt_item.setBackground(QBrush(bg_color))
            self.setItem(row, 3, lt_item)
            
            # RT (Right Taxon) - 越大越蓝
            rt_item = QTableWidgetItem(str(stat.right))
            rt_item.setTextAlignment(Qt.AlignCenter)
            if self.max_rt > 0:
                bg_color = interpolate_color(stat.right, 0, self.max_rt)
                rt_item.setBackground(QBrush(bg_color))
            self.setItem(row, 4, rt_item)
            
            # Taxon - 越大越蓝
            taxon_item = QTableWidgetItem(str(stat.taxon))
            taxon_item.setTextAlignment(Qt.AlignCenter)
            if self.max_taxon > 0:
                bg_color = interpolate_color(stat.taxon, 0, self.max_taxon)
                taxon_item.setBackground(QBrush(bg_color))
            self.setItem(row, 5, taxon_item)
            
            # 设置固定行高（紧凑式设计）
            self.setRowHeight(row, 25)
        
        logger.info(f"表格更新完成: {len(stats)} 行")


class ComboboxEditDelegate(QStyledItemDelegate):
    """
    ComboboxEdit 委托
    
    允许单元格既可以输入文本，也可以从下拉列表中选择
    """
    
    def __init__(self, items: List[str] = None, parent=None):
        super().__init__(parent)
        self.items = items or []
    
    def createEditor(self, parent, option, index):
        """创建编辑器"""
        editor = QLineEdit(parent)
        return editor
    
    def setEditorData(self, editor: QLineEdit, index):
        """设置编辑器数据"""
        text = index.model().data(index, Qt.EditRole)
        if text:
            editor.setText(text)
    
    def setModelData(self, editor: QLineEdit, model, index):
        """从编辑器设置模型数据"""
        text = editor.text()
        model.setData(index, text, Qt.EditRole)
    
    def updateEditorGeometry(self, editor, option, index):
        """更新编辑器几何形状"""
        editor.setGeometry(option.rect)


class StandardizationTable(QTableWidget):
    """
    标准化规则表格
    
    紧凑布局，支持增删行，每个单元格都是 ComboboxEdit
    """
    
    # 信号定义
    rulesChanged = pyqtSignal()  # 规则变化时发出
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.available_genes: List[str] = []  # 可用的基因名称列表
        self.init_ui()
    
    def init_ui(self):
        """初始化用户界面"""
        # 应用紧凑式样式
        self.setStyleSheet("QTableWidgetItem { padding: 0px; } QTableWidget { gridline-color: lightgray; }")
        
        # 设置列
        self.setColumnCount(2)
        self.setHorizontalHeaderLabels(["Old Name", "New Name"])
        
        # 设置列宽
        header = self.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        
        # 设置选择模式
        self.setSelectionBehavior(QTableWidget.SelectRows)
        self.setSelectionMode(QTableWidget.SingleSelection)
        self.setAlternatingRowColors(True)
        self.setEditTriggers(QTableWidget.DoubleClicked | QTableWidget.EditKeyPressed)
        
        # 应用紧凑式字体
        font = QFont()
        font.setPointSize(9)
        self.setFont(font)
        self.horizontalHeader().setFont(font)
        self.verticalHeader().setVisible(False)
        
        # 设置行高
        self.verticalHeader().setDefaultSectionSize(25)
        
        # 设置委托
        self.setItemDelegate(ComboboxEditDelegate(parent=self))
    
    def set_available_genes(self, genes: List[str]):
        """
        设置可用的基因名称列表
        
        Args:
            genes: 基因名称列表
        """
        self.available_genes = sorted(set(genes))
        # 更新委托的可选项
        for row in range(self.rowCount()):
            self.update_cell_options(row, 0)
            self.update_cell_options(row, 1)
    
    def update_cell_options(self, row: int, column: int):
        """
        更新单元格的可选项
        
        Args:
            row: 行索引
            column: 列索引
        """
        item = self.item(row, column)
        if item:
            # 这里可以添加下拉功能，目前简化为文本编辑
            pass
    
    def add_rule(self, old_name: str = "", new_name: str = ""):
        """
        添加一条规则
        
        Args:
            old_name: 旧名称
            new_name: 新名称
        """
        row = self.rowCount()
        self.insertRow(row)
        
        old_item = QTableWidgetItem(old_name)
        old_item.setFlags(Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        self.setItem(row, 0, old_item)
        
        new_item = QTableWidgetItem(new_name)
        new_item.setFlags(Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
        self.setItem(row, 1, new_item)
    
    def add_rules(self, rules: List[Tuple[str, str]]):
        """
        添加多条规则
        
        Args:
            rules: [(old_name, new_name), ...]
        """
        for old_name, new_name in rules:
            self.add_rule(old_name, new_name)
    
    def get_rules(self) -> Dict[str, str]:
        """
        获取所有规则
        
        Returns:
            {old_name: new_name} 字典
        """
        rules = {}
        for row in range(self.rowCount()):
            old_item = self.item(row, 0)
            new_item = self.item(row, 1)
            if old_item and new_item:
                old_name = old_item.text().strip()
                new_name = new_item.text().strip()
                if old_name and new_name and old_name != new_name:
                    rules[old_name] = new_name
        return rules
    
    def remove_selected_rule(self):
        """移除选中的规则"""
        current_row = self.currentRow()
        if current_row >= 0:
            self.removeRow(current_row)
            self.rulesChanged.emit()
    
    def clear_rules(self):
        """清空所有规则"""
        self.setRowCount(0)
        self.rulesChanged.emit()
    
    def on_cell_changed(self, item: QTableWidgetItem):
        """
        单元格内容变化事件
        
        Args:
            item: 变化的 QTableWidgetItem 对象
        """
        if item:
            self.rulesChanged.emit()


class TaxonManager(QWidget):
    """
    Taxon管理组件
    
    支持多选的QListWidget和'+' '-'按钮，用于导入GenBank文件
    """
    
    # 信号定义
    filesAdded = pyqtSignal(list)  # [file_paths]
    entriesRemoved = pyqtSignal(list)  # [entry_ids]
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.all_entries: List[str] = []  # 所有可用的Entry
        
        self.init_ui()
    
    def init_ui(self):
        """初始化用户界面"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # 标题和按钮
        header_layout = QHBoxLayout()
        
        title_label = QLabel("Taxon Management")
        title_label.setFont(QFont("Arial", 10, QFont.Bold))
        header_layout.addWidget(title_label)
        
        header_layout.addStretch()
        
        # 添加按钮
        self.add_btn = QPushButton("+")
        self.add_btn.setMaximumWidth(30)
        self.add_btn.setToolTip("Add GenBank Files")
        self.add_btn.clicked.connect(self.add_files)
        header_layout.addWidget(self.add_btn)
        
        # 移除按钮
        self.remove_btn = QPushButton("-")
        self.remove_btn.setMaximumWidth(30)
        self.remove_btn.setToolTip("Remove Selected Entries")
        self.remove_btn.clicked.connect(self.remove_selected_entries)
        header_layout.addWidget(self.remove_btn)
        
        layout.addLayout(header_layout)
        
        # Entry列表
        self.entry_list = QListWidget()
        self.entry_list.setSelectionMode(QListWidget.MultiSelection)
        layout.addWidget(self.entry_list)
    
    def add_files(self):
        """添加GenBank文件"""
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select GenBank Files",
            "",
            "GenBank Files (*.gb *.gbk);;All Files (*)"
        )
        
        if files:
            # 发送信号，让父窗口处理文件加载
            self.filesAdded.emit(files)
            logger.info(f"选择文件: {len(files)} 个")
    
    def add_entry(self, entry_id: str):
        """
        添加单个Entry
        
        Args:
            entry_id: Entry ID
        """
        self.entry_list.addItem(entry_id)
        logger.info(f"添加Entry: {entry_id}")
    
    def add_entries(self, entry_ids: List[str]):
        """
        添加多个Entry
        
        Args:
            entry_ids: Entry ID列表
        """
        for entry_id in entry_ids:
            self.entry_list.addItem(entry_id)
        logger.info(f"添加 {len(entry_ids)} 个Entry")
    
    def remove_selected_entries(self):
        """移除选中的Entry"""
        selected_items = self.entry_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "Warning", "Please select entries to remove")
            return
        
        # 获取要移除的Entry ID
        entry_ids = [item.text() for item in selected_items]
        
        # 移除选中的项
        for item in selected_items:
            row = self.entry_list.row(item)
            self.entry_list.takeItem(row)
            logger.info(f"移除Entry: {item.text()}")
        
        # 发送信号
        self.entriesRemoved.emit(entry_ids)
    
    def get_selected_entries(self) -> List[str]:
        """
        获取当前选中的Entry
        
        Returns:
            Entry ID列表
        """
        selected_items = self.entry_list.selectedItems()
        return [item.text() for item in selected_items]
    
    def get_all_entries(self) -> List[str]:
        """
        获取所有当前显示的Entry
        
        Returns:
            Entry ID列表
        """
        return [self.entry_list.item(i).text() for i in range(self.entry_list.count())]
    
    def clear(self):
        """清空Entry列表"""
        self.entry_list.clear()


class SankeyPlotCanvas(FigureCanvas):
    """
    桑基图画布
    
    绘制三层桑基图：左邻接节点 - taxon - 右邻接节点
    """
    
    def __init__(self, parent=None):
        self.fig, self.ax = plt.subplots(figsize=(14, 8), facecolor='white')
        super().__init__(self.fig)
        self.setParent(parent)
        
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.set_xlim(-0.5, 2.5)
        self.ax.set_ylim(0, 1)
        
        self.fig.tight_layout()
    
    def plot_sankey(
        self,
        left_data: List[Tuple[str, int]],    # [(node_name, count), ...]
        taxon_data: List[Tuple[str, int]],   # [(taxon_name, count), ...] 已废弃
        right_data: List[Tuple[str, int]],   # [(node_name, count), ...]
        left_matrix: List[List[int]],         # [left][taxon] 已废弃
        right_matrix: List[List[int]],        # [taxon][right] 已废弃
        gene_name: str = "",
        gene_type_map: Optional[Dict[str, str]] = None  # {node_name: gene_type}
    ):
        """
        绘制桑基图
        
        新结构：左邻接节点 - 当前基因 - 右邻接节点
        
        Args:
            left_data: 左邻接节点数据
            taxon_data: 已废弃，保留参数以兼容旧调用
            right_data: 右邻接节点数据
            left_matrix: 已废弃，保留参数以兼容旧调用
            right_matrix: 已废弃，保留参数以兼容旧调用
            gene_name: 当前基因名称（用于标题和中间节点）
            gene_type_map: 基因类型映射字典 {node_name: gene_type}
        """
        self.ax.clear()
        
        # 计算总高度
        left_values = [v for _, v in left_data]
        right_values = [v for _, v in right_data]
        
        # 中间节点：当前基因，高度等于所有连接的总数
        current_gene_value = sum(left_values) if left_values else sum(right_values)
        
        total_height = max(
            sum(left_values) if left_values else 1,
            current_gene_value if current_gene_value else 1,
            sum(right_values) if right_values else 1
        )
        
        # 归一化高度
        left_heights = np.array(left_values) / total_height
        center_height = current_gene_value / total_height
        right_heights = np.array(right_values) / total_height
        
        # 计算起始位置
        left_starts = np.concatenate(([0], np.cumsum(left_heights)[:-1]))
        center_start = 0.5 - center_height / 2  # 垂直居中
        right_starts = np.concatenate(([0], np.cumsum(right_heights)[:-1]))
        
        # 颜色映射
        num_left = len(left_data)
        num_right = len(right_data)
        
        # 基因类型颜色映射（与表格保持一致）
        def get_gene_type_color(gene_type: str) -> Tuple[float, float, float]:
            """根据基因类型返回RGB颜色"""
            if gene_type == "CDS":
                return (180/255, 220/255, 180/255)  # 浅绿色
            elif gene_type == "rRNA":
                return (255/255, 200/255, 150/255)  # 浅橙色
            elif gene_type == "tRNA":
                return (200/255, 180/255, 255/255)  # 浅紫色
            else:
                return (240/255, 240/255, 240/255)  # 浅灰色
        
        # 左节点颜色（根据基因类型）
        left_colors = []
        for name, _ in left_data:
            gene_type = gene_type_map.get(name, "gene") if gene_type_map else "gene"
            left_colors.append(get_gene_type_color(gene_type))
        
        # 中间节点颜色（根据当前基因类型）
        center_gene_type = gene_type_map.get(gene_name, "gene") if gene_type_map else "gene"
        center_color = get_gene_type_color(center_gene_type)
        
        # 右节点颜色（根据基因类型）
        right_colors = []
        for name, _ in right_data:
            gene_type = gene_type_map.get(name, "gene") if gene_type_map else "gene"
            right_colors.append(get_gene_type_color(gene_type))
        
        # 绘制连接（流线）
        # 左 -> 中间（当前基因）
        if left_data:
            for i, (name, count) in enumerate(left_data):
                y_start = left_starts[i]
                h_start = left_heights[i]
                
                # 左节点直接连接到中间节点
                flow = count
                
                # 流线宽度
                line_width = max(0.5, 3.0 * flow / total_height * 100)
                
                # 贝塞尔曲线
                x_start = 0
                x_end = 1
                cp_x = 0.5
                path_data = [
                    (Path.MOVETO, (x_start, y_start + h_start / 2)),
                    (Path.CURVE4, (cp_x, y_start + h_start / 2)),
                    (Path.CURVE4, (cp_x, center_start + center_height / 2)),
                    (Path.CURVE4, (x_end, center_start + center_height / 2)),
                ]
                codes, verts = zip(*path_data)
                path = Path(verts, codes)
                patch = PathPatch(
                    path,
                    facecolor='none',
                    edgecolor=left_colors[i],
                    alpha=0.6,
                    linewidth=line_width,
                    zorder=1
                )
                self.ax.add_patch(patch)
        
        # 中间（当前基因）-> 右
        if right_data:
            for j, (name, count) in enumerate(right_data):
                y_end = right_starts[j]
                h_end = right_heights[j]
                
                # 中间节点直接连接到右节点
                flow = count
                
                # 流线宽度
                line_width = max(0.5, 3.0 * flow / total_height * 100)
                
                # 贝塞尔曲线
                x_start = 1
                x_end = 2
                cp_x = 1.5
                path_data = [
                    (Path.MOVETO, (x_start, center_start + center_height / 2)),
                    (Path.CURVE4, (cp_x, center_start + center_height / 2)),
                    (Path.CURVE4, (cp_x, y_end + h_end / 2)),
                    (Path.CURVE4, (x_end, y_end + h_end / 2)),
                ]
                codes, verts = zip(*path_data)
                path = Path(verts, codes)
                patch = PathPatch(
                    path,
                    facecolor='none',
                    edgecolor=right_colors[j],
                    alpha=0.6,
                    linewidth=line_width,
                    zorder=1
                )
                self.ax.add_patch(patch)
        
        # 绘制节点（柱子）
        bar_width = 0.3
        
        # 左节点
        for i, (name, count) in enumerate(left_data):
            y_start = left_starts[i]
            h = left_heights[i]
            rect = FancyBboxPatch(
                (0 - bar_width / 2, y_start),
                bar_width,
                h,
                boxstyle="round,pad=0,rounding_size=0.01",
                facecolor=left_colors[i],
                edgecolor='black',
                linewidth=1.5,
                zorder=2
            )
            self.ax.add_patch(rect)
            self.ax.text(
                0 - bar_width / 2 - 0.02,
                y_start + h / 2,
                name,
                ha='right',
                va='center',
                fontsize=9,
                rotation=0
            )
            # 根据亮度判断文字颜色（使用亮度公式：0.299*R + 0.587*G + 0.114*B）
            brightness = 0.299 * left_colors[i][0] + 0.587 * left_colors[i][1] + 0.114 * left_colors[i][2]
            text_color = 'white' if brightness < 0.6 else 'black'
            self.ax.text(
                0,
                y_start + h / 2,
                str(count),
                ha='center',
                va='center',
                fontsize=8,
                fontweight='bold',
                color=text_color,
                zorder=3
            )
        
        # 中间节点（当前基因）
        if gene_name:
            rect = FancyBboxPatch(
                (1 - bar_width / 2, center_start),
                bar_width,
                center_height,
                boxstyle="round,pad=0,rounding_size=0.01",
                facecolor=center_color,
                edgecolor='black',
                linewidth=2.0,  # 加粗边框突出显示
                zorder=2
            )
            self.ax.add_patch(rect)
            self.ax.text(
                1 - bar_width / 2 - 0.02,
                center_start + center_height / 2,
                gene_name,
                ha='right',
                va='center',
                fontsize=10,
                fontweight='bold',
                rotation=0
            )
            # 根据亮度判断文字颜色
            brightness = 0.299 * center_color[0] + 0.587 * center_color[1] + 0.114 * center_color[2]
            text_color = 'white' if brightness < 0.6 else 'black'
            self.ax.text(
                1,
                center_start + center_height / 2,
                str(int(current_gene_value)),
                ha='center',
                va='center',
                fontsize=9,
                fontweight='bold',
                color=text_color,
                zorder=3
            )
        
        # 右节点
        for i, (name, count) in enumerate(right_data):
            y_start = right_starts[i]
            h = right_heights[i]
            rect = FancyBboxPatch(
                (2 - bar_width / 2, y_start),
                bar_width,
                h,
                boxstyle="round,pad=0,rounding_size=0.01",
                facecolor=right_colors[i],
                edgecolor='black',
                linewidth=1.5,
                zorder=2
            )
            self.ax.add_patch(rect)
            self.ax.text(
                2 + bar_width / 2 + 0.02,
                y_start + h / 2,
                name,
                ha='left',
                va='center',
                fontsize=9,
                rotation=0
            )
            # 根据亮度判断文字颜色（使用亮度公式：0.299*R + 0.587*G + 0.114*B）
            brightness = 0.299 * right_colors[i][0] + 0.587 * right_colors[i][1] + 0.114 * right_colors[i][2]
            text_color = 'white' if brightness < 0.6 else 'black'
            self.ax.text(
                2,
                y_start + h / 2,
                str(count),
                ha='center',
                va='center',
                fontsize=8,
                fontweight='bold',
                color=text_color,
                zorder=3
            )
        
        # 设置标题
        title = f"Gene Adjacency Sankey: {gene_name}" if gene_name else "Gene Adjacency Sankey"
        self.ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # 添加层标签
        self.ax.text(0, 1.02, "Left Adjacent", ha='center', fontsize=11, fontweight='bold', transform=self.ax.transAxes)
        self.ax.text(0.5, 1.02, "Current Gene", ha='center', fontsize=11, fontweight='bold', transform=self.ax.transAxes)
        self.ax.text(1, 1.02, "Right Adjacent", ha='center', fontsize=11, fontweight='bold', transform=self.ax.transAxes)
        
        # 调整视图范围
        self.ax.set_xlim(-0.6, 2.6)
        self.ax.set_ylim(0, 1)
        
        self.draw()
    
    def clear_plot(self):
        """清空图表"""
        self.ax.clear()
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.set_xlim(-0.5, 2.5)
        self.ax.set_ylim(0, 1)
        self.ax.set_title("Select a gene to view adjacency relationships", fontsize=14)
        self.draw()


class SeqDBGTableViewer(QWidget):
    """
    SeqDBG表格查看器
    
    提供基因统计表格和桑基图可视化
    新布局：
    - 主布局：QSplitter，竖直方向
    - 左侧：QSplitter，水平方向
      - 左上：taxon管理
      - 左下：主体表格
    - 右侧：多标签页
      - Plot 标签页：桑基图
      - Standardization 标签页：名称标准化规则表格
    """
    
    # 信号定义
    geneSelected = pyqtSignal(str)  # gene_name
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.stats_data: List[GeneStats] = []
        self.adjacency_data: Dict[str, Dict] = {}
        
        # 标准化相关
        self.original_stats_data: List[GeneStats] = []  # 保留原始数据副本
        self.original_adjacency_data: Dict[str, Dict] = {}  # 保留原始邻接数据副本
        
        self.init_ui()
    
    def init_ui(self):
        """初始化用户界面"""
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        main_layout.setContentsMargins(0, 0, 0, 0)
        
        # 创建主分割器（水平，分为左右）
        main_splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(main_splitter)
        
        # 左侧：分割器（竖直，分为左上和左下）
        left_splitter = QSplitter(Qt.Vertical)
        main_splitter.addWidget(left_splitter)
        
        # 左上：taxon管理（无GroupBox）
        taxon_container = QWidget()
        taxon_layout = QVBoxLayout()
        taxon_layout.setContentsMargins(0, 0, 0, 0)
        taxon_container.setLayout(taxon_layout)
        
        self.taxon_manager = TaxonManager()
        taxon_layout.addWidget(self.taxon_manager)
        
        left_splitter.addWidget(taxon_container)
        
        # 左下：主体表格（无GroupBox）
        table_container = QWidget()
        table_layout = QVBoxLayout()
        table_layout.setContentsMargins(0, 0, 0, 0)
        table_container.setLayout(table_layout)
        
        # # 说明标签
        # info_label = QLabel(
        #     "<b>Color Legend:</b> "
        #     "<span style='background-color:#f0b8b8; padding:2px 8px;'>Red</span> "
        #     "LN/RN (larger = redder) | "
        #     "<span style='background-color:#b8d4f3; padding:2px 8px;'>Blue</span> "
        #     "LT/RT/Taxon (larger = bluer)"
        # )
        # info_label.setWordWrap(True)
        # info_label.setFont(QFont("Arial", 8))
        # table_layout.addWidget(info_label)
        
        # 表格
        self.table = GeneStatsTable()
        self.table.cellClicked.connect(self.on_table_cell_clicked)
        table_layout.addWidget(self.table)
        
        left_splitter.addWidget(table_container)
        
        # 右侧：多标签页
        self.tab_widget = QTabWidget()
        
        # Plot 标签页：桑基图
        plot_container = QWidget()
        plot_layout = QVBoxLayout()
        plot_layout.setContentsMargins(0, 0, 0, 0)
        plot_container.setLayout(plot_layout)
        
        self.sankey_canvas = SankeyPlotCanvas()
        plot_layout.addWidget(self.sankey_canvas)
        
        self.tab_widget.addTab(plot_container, "Plot")
        
        # Standardization 标签页：名称标准化
        std_container = QWidget()
        std_layout = QVBoxLayout()
        std_layout.setContentsMargins(5, 5, 5, 5)
        std_container.setLayout(std_layout)
        
        # 标题
        std_title = QLabel("Name Standardization Rules")
        std_title.setFont(QFont("Arial", 10, QFont.Bold))
        std_layout.addWidget(std_title)
        
        # 标准化规则表格
        self.std_table = StandardizationTable()
        self.std_table.itemChanged.connect(self.std_table.on_cell_changed)
        std_layout.addWidget(self.std_table)
        
        # 按钮区域
        std_button_layout = QHBoxLayout()
        
        # 添加规则按钮
        add_rule_btn = QPushButton("+ Add Rule")
        add_rule_btn.setMaximumWidth(100)
        add_rule_btn.clicked.connect(self.on_add_standardization_rule)
        std_button_layout.addWidget(add_rule_btn)
        
        # 删除规则按钮
        remove_rule_btn = QPushButton("- Remove Rule")
        remove_rule_btn.setMaximumWidth(120)
        remove_rule_btn.clicked.connect(self.on_remove_standardization_rule)
        std_button_layout.addWidget(remove_rule_btn)
        
        # Import rules 按钮
        import_btn = QPushButton("Import rules")
        import_btn.setMaximumWidth(110)
        import_btn.clicked.connect(self.on_import_rules)
        std_button_layout.addWidget(import_btn)
        
        # Export rules 按钮
        export_btn = QPushButton("Export rules")
        export_btn.setMaximumWidth(110)
        export_btn.clicked.connect(self.on_export_rules)
        std_button_layout.addWidget(export_btn)
        
        std_button_layout.addStretch()
        
        # Reset 按钮
        reset_btn = QPushButton("Reset")
        reset_btn.setMaximumWidth(80)
        reset_btn.clicked.connect(self.reset_to_original)
        std_button_layout.addWidget(reset_btn)
        
        # Refresh 按钮
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.setMaximumWidth(100)
        self.refresh_btn.clicked.connect(self.on_refresh_standardization)
        std_button_layout.addWidget(self.refresh_btn)
        
        std_layout.addLayout(std_button_layout)
        
        self.tab_widget.addTab(std_container, "Standardization")
        
        # Dataset 标签页：Dataset 导出
        self.dataset_widget = DatasetExportWidget()
        self.tab_widget.addTab(self.dataset_widget, "Dataset")
        
        main_splitter.addWidget(self.tab_widget)
        
        # 设置分割器比例
        # # 主分割器：左侧占40%，右侧占60% (4:6)
        main_splitter.setStretchFactor(0, 6)
        main_splitter.setStretchFactor(1, 6)
        
        # 左侧分割器：taxon管理占40%，表格占60% (4:6)
        left_splitter.setStretchFactor(0, 4)
        left_splitter.setStretchFactor(1, 6)
    
    def set_stats_data(self, stats: List[GeneStats], loaded_files: Optional[List[str]] = None,
                        gene_info_data: Optional[Dict[str, Dict[str, GeneInfo]]] = None,
                        sequence_data: Optional[Dict[str, Dict[str, SeqRecord]]] = None,
                        gene_type_map: Optional[Dict[str, str]] = None):
        """
        设置统计数据
        
        Args:
            stats: GeneStats列表
            loaded_files: 加载的文件列表（可选）
            gene_info_data: 基因信息数据 {gene_name: {genome_id: GeneInfo}}
            sequence_data: 序列数据 {gene_name: {genome_id: SeqRecord}}
            gene_type_map: 基因类型映射 {gene_name: gene_type}
        """
        # 保存原始数据副本（如果还没有）
        if not self.original_stats_data:
            self.original_stats_data = copy.deepcopy(stats)
            logger.info(f"保存原始数据副本: {len(self.original_stats_data)} 个基因")
        
        self.stats_data = stats
        self.table.set_data(stats)
        
        # 更新标准化表格的可选项
        all_genes = [stat.gene_name for stat in stats]
        self.std_table.set_available_genes(all_genes)
        
        # 更新 dataset widget 的数据
        self.dataset_widget.set_data(stats, loaded_files or [], gene_info_data, sequence_data, gene_type_map)
    
    def set_adjacency_data(self, adjacency_data: Dict[str, Dict]):
        """
        设置邻接数据
        
        Args:
            adjacency_data: {gene_name: adjacency_info} 字典
        """
        # 保存原始数据副本（如果还没有）
        if not self.original_adjacency_data:
            self.original_adjacency_data = copy.deepcopy(adjacency_data)
            logger.info(f"保存原始邻接数据副本: {len(self.original_adjacency_data)} 个基因")
        
        self.adjacency_data = adjacency_data
    
    def set_all_entries(self, entries: List[str]):
        """
        设置所有可用的Entry（已废弃，不再需要）
        
        Args:
            entries: Entry ID列表
        """
        # 不再需要维护all_entries列表
        pass
    
    def on_table_cell_clicked(self, row: int, column: int):
        """
        表格单元格点击事件
        
        Args:
            row: 行索引
            column: 列索引
        """
        if row >= len(self.stats_data):
            return
        
        stat = self.stats_data[row]
        gene_name = stat.gene_name
        
        # 更新桑基图
        self.update_sankey_plot(gene_name)
        
        # 发送信号
        self.geneSelected.emit(gene_name)
    
    def update_sankey_plot(self, gene_name: str):
        """
        更新桑基图
        
        Args:
            gene_name: 基因名称
        """
        if gene_name not in self.adjacency_data:
            logger.warning(f"未找到基因 {gene_name} 的邻接数据")
            return
        
        adj_info = self.adjacency_data[gene_name]
        
        # 准备数据
        left_data = adj_info.get("left_adjacent", [])
        right_data = adj_info.get("right_adjacent", [])
        
        # 收集基因类型信息
        gene_type_map = {}
        for stat in self.stats_data:
            gene_type_map[stat.gene_name] = getattr(stat, 'gene_type', 'gene')
        
        # 绘制桑基图（不再需要 taxon_data 和连接矩阵）
        self.sankey_canvas.plot_sankey(
            left_data=left_data,
            taxon_data=[],  # 不再使用
            right_data=right_data,
            left_matrix=[],  # 不再使用
            right_matrix=[],  # 不再使用
            gene_name=gene_name,
            gene_type_map=gene_type_map
        )
    
    def clear_all(self):
        """清空所有数据"""
        self.stats_data.clear()
        self.adjacency_data.clear()
        self.original_stats_data.clear()
        self.original_adjacency_data.clear()
        self.table.setRowCount(0)
        self.taxon_manager.clear()
        self.sankey_canvas.clear_plot()
        self.std_table.clear_rules()
    
    def on_add_standardization_rule(self):
        """添加标准化规则"""
        self.std_table.add_rule()
    
    def on_import_rules(self):
        """从 tab 文件导入规则"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Import Standardization Rules",
            "",
            "Tab-separated files (*.tab *.txt);;All Files (*)"
        )
        
        if not file_path:
            return
        
        try:
            import_count = 0
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # 解析 tab 分隔的行
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        old_name = parts[0].strip()
                        new_name = parts[1].strip()
                        if old_name and new_name and old_name != new_name:
                            self.std_table.add_rule(old_name, new_name)
                            import_count += 1
            
            logger.info(f"从 {file_path} 导入了 {import_count} 条规则")
            QMessageBox.information(
                self,
                "Import Complete",
                f"Successfully imported {import_count} rules from:\n{file_path}"
            )
            
        except Exception as e:
            logger.error(f"导入规则失败: {e}")
            QMessageBox.critical(
                self,
                "Import Error",
                f"Failed to import rules:\n{e}"
            )
    
    def on_export_rules(self):
        """导出规则到 tab 文件"""
        default_filename = "standardization.tab"
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Standardization Rules",
            default_filename,
            "Tab-separated files (*.tab *.txt);;All Files (*)"
        )
        
        if not file_path:
            return
        
        try:
            rules = self.std_table.get_rules()
            
            with open(file_path, 'w', encoding='utf-8') as f:
                # 写入头部注释
                f.write("# SeqDBG Standardization Rules\n")
                f.write("# Format: Old Name<TAB>New Name\n")
                f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("# Number of rules: " + str(len(rules)) + "\n\n")
                
                # 写入规则
                for old_name, new_name in rules.items():
                    f.write(f"{old_name}\t{new_name}\n")
            
            logger.info(f"导出 {len(rules)} 条规则到 {file_path}")
            QMessageBox.information(
                self,
                "Export Complete",
                f"Successfully exported {len(rules)} rules to:\n{file_path}"
            )
            
        except Exception as e:
            logger.error(f"导出规则失败: {e}")
            QMessageBox.critical(
                self,
                "Export Error",
                f"Failed to export rules:\n{e}"
            )
    
    def on_remove_standardization_rule(self):
        """移除标准化规则"""
        self.std_table.remove_selected_rule()
    
    def on_refresh_standardization(self):
        """
        应用标准化规则并刷新
        
        使用标准化规则将原始数据中的旧名称替换为新名称，
        然后更新表格和图表
        
        同时会自动将预定义的规则（来自 normalizer）添加到表格中
        """
        # 如果没有原始数据，先保存
        if not self.original_stats_data and self.stats_data:
            self.original_stats_data = copy.deepcopy(self.stats_data)
            self.original_adjacency_data = copy.deepcopy(self.adjacency_data)
            logger.info(f"保存原始数据副本: {len(self.original_stats_data)} 个基因")
        
        # 如果没有原始数据，无法继续
        if not self.original_stats_data:
            QMessageBox.warning(self, "Warning", "No original data to refresh from")
            return
        
        # 获取标准化规则（表格中的规则）
        table_rules = self.std_table.get_rules()
        
        # 获取所有原始基因名称
        all_original_names = set(stat.gene_name for stat in self.original_stats_data)
        
        # 应用表格中的规则
        logger.info(f"应用表格规则 {len(table_rules)} 条: {table_rules}")
        
        # 应用规则：替换原始数据中的基因名称
        new_stats_data = []
        name_mapping = {}  # {original_name: new_name}
        
        for stat in self.original_stats_data:
            gene_name = stat.gene_name
            new_name = table_rules.get(gene_name, gene_name)
            
            if gene_name != new_name:
                name_mapping[gene_name] = new_name
                logger.debug(f"标准化: {gene_name} -> {new_name}")
            
            # 创建新的统计对象
            new_stat = GeneStats(
                gene_name=new_name,
                left=stat.left,
                right=stat.right,
                left_nodes=stat.left_nodes,
                right_nodes=stat.right_nodes,
                taxon=stat.taxon,
                gene_type=stat.gene_type
            )
            new_stats_data.append(new_stat)
        
        # 更新邻接数据中的基因名称
        new_adjacency_data = {}
        for gene_name, adj_info in self.original_adjacency_data.items():
            new_gene_name = table_rules.get(gene_name, gene_name)
            
            # 转换左邻接数据
            new_left_adjacent = [(table_rules.get(name, name), count) for name, count in adj_info.get("left_adjacent", [])]
            
            # 转换右邻接数据
            new_right_adjacent = [(table_rules.get(name, name), count) for name, count in adj_info.get("right_adjacent", [])]
            
            new_adjacency_data[new_gene_name] = {
                "gene_name": new_gene_name,
                "left_adjacent": new_left_adjacent,
                "right_adjacent": new_right_adjacent,
                "taxon_info": adj_info.get("taxon_info", [])
            }
        
        # 合并相同基因的统计信息（因为多个旧名称可能映射到同一个新名称）
        merged_stats = {}
        for stat in new_stats_data:
            if stat.gene_name not in merged_stats:
                merged_stats[stat.gene_name] = stat
            else:
                # 合并统计信息
                existing = merged_stats[stat.gene_name]
                existing.left += stat.left
                existing.right += stat.right
                existing.left_nodes = max(existing.left_nodes, stat.left_nodes)
                existing.right_nodes = max(existing.right_nodes, stat.right_nodes)
                existing.taxon += stat.taxon
        
        self.stats_data = list(merged_stats.values())
        
        # 合并邻接数据
        for gene_name, adj_info in new_adjacency_data.items():
            if gene_name not in self.adjacency_data:
                self.adjacency_data[gene_name] = adj_info
            else:
                # 合并邻接信息
                existing = self.adjacency_data[gene_name]
                # 合并左邻接
                existing_left = dict(existing.get("left_adjacent", []))
                for name, count in adj_info.get("left_adjacent", []):
                    existing_left[name] = existing_left.get(name, 0) + count
                existing["left_adjacent"] = sorted(existing_left.items(), key=lambda x: -x[1])
                
                # 合并右邻接
                existing_right = dict(existing.get("right_adjacent", []))
                for name, count in adj_info.get("right_adjacent", []):
                    existing_right[name] = existing_right.get(name, 0) + count
                existing["right_adjacent"] = sorted(existing_right.items(), key=lambda x: -x[1])
        
        # 更新表格和图表
        self.table.set_data(self.stats_data)
        
        # 更新标准化表格的可选项
        all_genes = [stat.gene_name for stat in self.stats_data]
        self.std_table.set_available_genes(all_genes)
        
        # 清空桑基图
        self.sankey_canvas.clear_plot()
        
        logger.info(f"标准化完成: {len(name_mapping)} 个基因名称被修改")
        
        QMessageBox.information(
            self,
            "Standardization Complete",
            f"Applied {len(table_rules)} rules.\n{len(name_mapping)} gene names changed."
        )
    
    def load_preset_rules(self, rules: Dict[str, str]):
        """
        加载预定义的标准化规则到表格中
        
        Args:
            rules: {old_name: new_name} 字典
        """
        for old_name, new_name in rules.items():
            if old_name != new_name:
                self.std_table.add_rule(old_name, new_name)
        logger.info(f"加载 {len(rules)} 条预定义规则")
    
    def reset_to_original(self):
        """
        重置到原始数据
        
        清除所有标准化规则，恢复到原始数据
        """
        if not self.original_stats_data:
            QMessageBox.warning(self, "Warning", "No original data to restore")
            return
        
        import copy
        
        # 恢复原始数据
        self.stats_data = copy.deepcopy(self.original_stats_data)
        self.adjacency_data = copy.deepcopy(self.original_adjacency_data)
        
        # 清空标准化规则
        self.std_table.clear_rules()
        
        # 更新表格
        self.table.set_data(self.stats_data)
        
        # 更新标准化表格的可选项
        all_genes = [stat.gene_name for stat in self.stats_data]
        self.std_table.set_available_genes(all_genes)
        
        # 清空桑基图
        self.sankey_canvas.clear_plot()
        
        logger.info("重置到原始数据")
        QMessageBox.information(self, "Reset Complete", "Reset to original data.")


class CustomizeLociDialog(QDialog):
    """
    自定义 Loci 选择对话框
    
    参考 MD-MRCA 的 taxon set 设计
    """
    
    def __init__(self, available_loci: List[str], selected_loci: Optional[List[str]] = None, parent=None):
        """
        初始化对话框
        
        Args:
            available_loci: 所有可用的 loci 列表
            selected_loci: 已选择的 loci 列表（可选）
            parent: 父窗口
        """
        super().__init__(parent)
        self.setWindowTitle("Customize Loci Selection")
        self.setMinimumSize(600, 400)
        
        self.available_loci = available_loci
        self.selected_loci = selected_loci or []
        
        self.init_ui()
        self.load_data()
    
    def init_ui(self):
        """初始化用户界面"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # 标题
        title_label = QLabel("Select loci to export:")
        title_label.setFont(QFont("Arial", 10, QFont.Bold))
        layout.addWidget(title_label)
        
        # 主布局：两个列表 + 按钮
        main_layout = QHBoxLayout()
        layout.addLayout(main_layout)
        
        # 左侧：可用 loci
        left_layout = QVBoxLayout()
        main_layout.addLayout(left_layout)
        
        left_label = QLabel("Available Loci:")
        left_layout.addWidget(left_label)
        
        self.available_list = QListWidget()
        self.available_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        left_layout.addWidget(self.available_list)
        
        # 中间：按钮
        button_layout = QVBoxLayout()
        main_layout.addLayout(button_layout)
        
        add_btn = QPushButton("Add >")
        add_btn.clicked.connect(self.add_loci)
        button_layout.addWidget(add_btn)
        
        remove_btn = QPushButton("< Remove")
        remove_btn.clicked.connect(self.remove_loci)
        button_layout.addWidget(remove_btn)
        
        button_layout.addStretch()
        
        select_all_btn = QPushButton("Select All")
        select_all_btn.clicked.connect(self.select_all)
        button_layout.addWidget(select_all_btn)
        
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_selected)
        button_layout.addWidget(clear_btn)
        
        # 右侧：已选择 loci
        right_layout = QVBoxLayout()
        main_layout.addLayout(right_layout)
        
        right_label = QLabel("Selected Loci:")
        right_layout.addWidget(right_label)
        
        self.selected_list = QListWidget()
        self.selected_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        right_layout.addWidget(self.selected_list)
        
        # 底部：按钮
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
    
    def load_data(self):
        """加载数据到列表"""
        # 加载可用 loci
        self.available_list.clear()
        for locus in self.available_loci:
            item = QListWidgetItem(locus)
            self.available_list.addItem(item)
        
        # 加载已选择 loci
        self.selected_list.clear()
        for locus in self.selected_loci:
            item = QListWidgetItem(locus)
            self.selected_list.addItem(item)
    
    def add_loci(self):
        """添加选中的 loci"""
        selected_items = self.available_list.selectedItems()
        for item in selected_items:
            locus = item.text()
            # 检查是否已存在
            if not self.selected_list.findItems(locus, Qt.MatchExactly):
                new_item = QListWidgetItem(locus)
                self.selected_list.addItem(new_item)
        self.available_list.clearSelection()
    
    def remove_loci(self):
        """移除选中的 loci"""
        selected_items = self.selected_list.selectedItems()
        for item in selected_items:
            row = self.selected_list.row(item)
            self.selected_list.takeItem(row)
    
    def select_all(self):
        """选择所有可用 loci"""
        self.available_list.selectAll()
    
    def clear_selected(self):
        """清空已选择 loci"""
        self.selected_list.clear()
    
    def get_selected_loci(self) -> List[str]:
        """
        获取已选择的 loci 列表
        
        Returns:
            loci 名称列表
        """
        return [self.selected_list.item(i).text() for i in range(self.selected_list.count())]


class DatasetExportWidget(QWidget):
    """
    Dataset 导出组件
    
    支持将基因组导出为 FASTA 格式的 dataset
    """
    
    # 信号定义
    exportToDataset = pyqtSignal(list)  # 发送 dataset 到 YR-MPEA
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.stats_data: List[GeneStats] = []
        self.available_loci: List[str] = []
        self.selected_loci: List[str] = []
        self.loaded_files: List[str] = []
        
        # 基因信息和序列数据 {gene_name: {genome_id: GeneInfo}}
        self.gene_info_data: Dict[str, Dict[str, GeneInfo]] = {}
        # 序列数据 {gene_name: {genome_id: SeqRecord}}
        self.sequence_data: Dict[str, Dict[str, SeqRecord]] = {}
        
        # 基因类型映射 {gene_name: gene_type}
        self.gene_type_map: Dict[str, str] = {}
        
        # 基因类型颜色映射（与表格保持一致）
        self.gene_type_colors = {
            "CDS": QColor(180, 220, 180),    # 浅绿色
            "rRNA": QColor(255, 200, 150),   # 浅橙色
            "tRNA": QColor(200, 180, 255),   # 浅紫色
            "gene": QColor(240, 240, 240),   # 浅灰色
            "other": QColor(240, 240, 240)   # 浅灰色
        }
        
        self.init_ui()
    
    def init_ui(self):
        """初始化用户界面"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # 标题
        title_label = QLabel("Dataset Export")
        title_label.setFont(QFont("Arial", 12, QFont.Bold))
        layout.addWidget(title_label)
        
        # Dataset type
        type_group = QGroupBox("Dataset Type")
        type_layout = QVBoxLayout()
        type_group.setLayout(type_layout)
        
        self.nucleotide_radio = QRadioButton("Nucleotide (DNA / RNA)")
        self.nucleotide_radio.setChecked(True)
        self.nucleotide_radio.toggled.connect(self.on_dataset_type_changed)
        type_layout.addWidget(self.nucleotide_radio)
        
        # Nucleotide 选项
        nucleotide_options_layout = QHBoxLayout()
        nucleotide_options_layout.setContentsMargins(20, 5, 5, 5)
        type_layout.addLayout(nucleotide_options_layout)
        
        self.cds_checkbox = QCheckBox("CDS")
        self.cds_checkbox.setChecked(True)
        self.cds_checkbox.setStyleSheet(f"background-color: {self.gene_type_colors['CDS'].name()}; padding: 2px 5px;")
        nucleotide_options_layout.addWidget(self.cds_checkbox)
        
        self.trna_checkbox = QCheckBox("tRNA")
        self.trna_checkbox.setChecked(True)
        self.trna_checkbox.setStyleSheet(f"background-color: {self.gene_type_colors['tRNA'].name()}; padding: 2px 5px;")
        nucleotide_options_layout.addWidget(self.trna_checkbox)
        
        self.rrna_checkbox = QCheckBox("rRNA")
        self.rrna_checkbox.setChecked(True)
        self.rrna_checkbox.setStyleSheet(f"background-color: {self.gene_type_colors['rRNA'].name()}; padding: 2px 5px;")
        nucleotide_options_layout.addWidget(self.rrna_checkbox)
        
        self.other_nucleotide_checkbox = QCheckBox("Others")
        nucleotide_options_layout.addWidget(self.other_nucleotide_checkbox)
        
        self.amino_acid_radio = QRadioButton("Amino Acid (Proteins)")
        self.amino_acid_radio.toggled.connect(self.on_dataset_type_changed)
        type_layout.addWidget(self.amino_acid_radio)
        
        # Amino Acid 选项
        amino_options_layout = QHBoxLayout()
        amino_options_layout.setContentsMargins(20, 5, 5, 5)
        type_layout.addLayout(amino_options_layout)
        
        self.amino_cds_checkbox = QCheckBox("CDS (Protein sequences)")
        self.amino_cds_checkbox.setChecked(True)
        amino_options_layout.addWidget(self.amino_cds_checkbox)
        
        layout.addWidget(type_group)
        
        # Loci 选项
        loci_group = QGroupBox("Loci Selection")
        loci_layout = QVBoxLayout()
        loci_group.setLayout(loci_layout)
        
        loci_radio_layout = QHBoxLayout()
        loci_layout.addLayout(loci_radio_layout)
        
        self.all_loci_radio = QRadioButton("All")
        self.all_loci_radio.setChecked(True)
        self.all_loci_radio.toggled.connect(self.on_loci_selection_changed)
        loci_radio_layout.addWidget(self.all_loci_radio)
        
        self.customize_loci_radio = QRadioButton("Customize")
        self.customize_loci_radio.toggled.connect(self.on_loci_selection_changed)
        loci_radio_layout.addWidget(self.customize_loci_radio)
        
        self.customize_loci_btn = QPushButton("Customize loci")
        self.customize_loci_btn.clicked.connect(self.on_customize_loci)
        self.customize_loci_btn.setEnabled(False)
        loci_layout.addWidget(self.customize_loci_btn)
        
        # 显示选中的 loci 数量
        self.loci_count_label = QLabel("Selected: All available loci")
        loci_layout.addWidget(self.loci_count_label)
        
        layout.addWidget(loci_group)
        
        layout.addStretch()
        
        # Export 按钮
        button_layout = QHBoxLayout()
        layout.addLayout(button_layout)
        
        self.export_button = QPushButton("Export (to multiple FASTA files)")
        self.export_button.setMinimumHeight(40)
        self.export_button.clicked.connect(self.on_export)
        button_layout.addWidget(self.export_button)
        
        self.export_to_dataset_button = QPushButton("Export to YR-MPEA")
        self.export_to_dataset_button.setMinimumHeight(40)
        self.export_to_dataset_button.clicked.connect(self.on_export_to_dataset)
        button_layout.addWidget(self.export_to_dataset_button)
    
    def set_data(self, stats: List[GeneStats], loaded_files: List[str],
                gene_info_data: Optional[Dict[str, Dict[str, GeneInfo]]] = None,
                sequence_data: Optional[Dict[str, Dict[str, SeqRecord]]] = None,
                gene_type_map: Optional[Dict[str, str]] = None):
        """
        设置数据
        
        Args:
            stats: GeneStats 列表
            loaded_files: 加载的文件列表
            gene_info_data: 基因信息数据 {gene_name: {genome_id: GeneInfo}}
            sequence_data: 序列数据 {gene_name: {genome_id: SeqRecord}}
            gene_type_map: 基因类型映射 {gene_name: gene_type}
        """
        self.stats_data = stats
        self.loaded_files = loaded_files
        self.gene_info_data = gene_info_data or {}
        self.sequence_data = sequence_data or {}
        self.gene_type_map = gene_type_map or {}
        
        # 更新可用的 loci
        self.available_loci = [stat.gene_name for stat in stats]
        self.update_loci_count_label()
    
    def on_dataset_type_changed(self):
        """Dataset type 改变时处理"""
        is_nucleotide = self.nucleotide_radio.isChecked()
        
        # 启用/禁用对应的 checkbox
        self.cds_checkbox.setEnabled(is_nucleotide)
        self.trna_checkbox.setEnabled(is_nucleotide)
        self.rrna_checkbox.setEnabled(is_nucleotide)
        self.other_nucleotide_checkbox.setEnabled(is_nucleotide)
        
        self.amino_cds_checkbox.setEnabled(not is_nucleotide)
        
        # TODO: 需要触发序列数据重新获取（需要从 parent 获取 builder）
        # 暂时只是更新可用 loci
        self.update_loci_count_label()
    
    def on_loci_selection_changed(self):
        """Loci 选择方式改变时处理"""
        self.customize_loci_btn.setEnabled(self.customize_loci_radio.isChecked())
        self.update_loci_count_label()
    
    def on_customize_loci(self):
        """打开自定义 loci 对话框"""
        dialog = CustomizeLociDialog(self.available_loci, self.selected_loci, self)
        if dialog.exec_() == QDialog.Accepted:
            self.selected_loci = dialog.get_selected_loci()
            self.update_loci_count_label()
            logger.info(f"选择了 {len(self.selected_loci)} 个 loci")
    
    def update_loci_count_label(self):
        """更新 loci 数量标签"""
        if self.all_loci_radio.isChecked():
            self.loci_count_label.setText(f"Selected: All ({len(self.available_loci)}) available loci")
        else:
            self.loci_count_label.setText(f"Selected: {len(self.selected_loci)} / {len(self.available_loci)} loci")
    
    def get_selected_gene_types(self) -> List[str]:
        """
        获取选中的基因类型
        
        Returns:
            基因类型列表
        """
        if self.nucleotide_radio.isChecked():
            types = []
            if self.cds_checkbox.isChecked():
                types.append("CDS")
            if self.trna_checkbox.isChecked():
                types.append("tRNA")
            if self.rrna_checkbox.isChecked():
                types.append("rRNA")
            if self.other_nucleotide_checkbox.isChecked():
                types.append("gene")
                types.append("other")
            return types
        else:
            if self.amino_cds_checkbox.isChecked():
                return ["CDS"]
            return []
    
    def get_selected_loci_for_export(self) -> List[str]:
        """
        获取用于导出的 loci 列表（根据基因类型过滤）
        
        Returns:
            loci 名称列表
        """
        if self.all_loci_radio.isChecked():
            selected = self.available_loci.copy()
        else:
            selected = self.selected_loci.copy()
        
        # 根据选中的基因类型过滤
        gene_types = self.get_selected_gene_types()
        
        if not gene_types:
            return selected
        
        filtered_loci = []
        for locus in selected:
            locus_type = self.gene_type_map.get(locus, "gene")
            if locus_type in gene_types:
                filtered_loci.append(locus)
        
        return filtered_loci
    
    def on_export(self):
        """导出到 FASTA 文件"""
        if not self.stats_data:
            QMessageBox.warning(self, "Warning", "No data to export")
            return
        
        # 选择导出目录
        output_dir = QFileDialog.getExistingDirectory(
            self,
            "Select Output Directory"
        )
        
        if not output_dir:
            return
        
        try:
            # 获取选中的基因类型和 loci
            gene_types = self.get_selected_gene_types()
            loci = self.get_selected_loci_for_export()
            
            if not gene_types:
                QMessageBox.warning(self, "Warning", "Please select at least one gene type")
                return
            
            if not loci:
                QMessageBox.warning(self, "Warning", f"No loci match the selected gene types: {', '.join(gene_types)}")
                return
            
            is_amino_acid = self.amino_acid_radio.isChecked()
            exported_count = 0
            skipped_count = 0
            
            for locus in loci:
                sequences = []
                
                # 收集该 locus 的所有序列
                if locus in self.sequence_data:
                    for genome_id, seq_records in self.sequence_data[locus].items():
                        # seq_records 现在是一个列表，遍历所有序列
                        for seq_record in seq_records:
                            # 使用 seq_record.id 作为序列名（包含 entry_id）
                            seq_name = seq_record.id
                            seq_str = str(seq_record.seq).upper()
                            
                            # 如果是 Amino Acid 模式，使用蛋白质序列
                            if is_amino_acid:
                                # 检查是否是蛋白质序列（带有 _protein 后缀）
                                if "_protein" in seq_name:
                                    sequences.append(f">{seq_name}\n{seq_str}\n")
                            else:
                                # Nucleotide 模式，使用 DNA 序列（跳过蛋白质序列）
                                if "_protein" not in seq_name:
                                    sequences.append(f">{seq_name}\n{seq_str}\n")
                
                # 如果没有找到序列，跳过
                if not sequences:
                    logger.warning(f"未找到 locus {locus} 的序列（模式: {'Amino Acid' if is_amino_acid else 'Nucleotide'}）")
                    skipped_count += 1
                    continue
                
                # 写入文件
                file_path = os.path.join(output_dir, f"{locus}.fasta")
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.writelines(sequences)
                
                exported_count += 1
                logger.info(f"导出 {locus}.fasta: {len(sequences)} 条序列 ({'protein' if is_amino_acid else 'DNA'})")
            
            if exported_count == 0:
                QMessageBox.warning(
                    self,
                    "Export Warning",
                    "No sequences found for the selected loci.\n\n"
                    "Possible reasons:\n"
                    "1. Make sure 'extract_sequences' is enabled when loading GenBank files\n"
                    "2. For Amino Acid mode, ensure CDS features have /translation fields"
                )
            else:
                message = f"Successfully exported {exported_count} FASTA files to:\n{output_dir}\n\n"
                message += f"Mode: {'Amino Acid (Protein sequences)' if is_amino_acid else 'Nucleotide (DNA/RNA sequences)'}"
                if skipped_count > 0:
                    message += f"\nSkipped: {skipped_count} loci (no sequences available)"
                
                QMessageBox.information(
                    self,
                    "Export Complete",
                    message
                )
            
        except Exception as e:
            logger.error(f"导出失败: {e}")
            QMessageBox.critical(
                self,
                "Export Error",
                f"Failed to export:\n{e}"
            )
    
    def on_export_to_dataset(self):
        """导出到 YR-MPEA dataset"""
        if not self.stats_data:
            QMessageBox.warning(self, "Warning", "No data to export")
            return
        
        try:
            # TODO: 实现实际的导出到 YR-MPEA 的逻辑
            # 这里发送信号到 YR-MPEA
            dataset_items = []  # TODO: 构建 DatasetItem 列表
            self.exportToDataset.emit(dataset_items)
            
            QMessageBox.information(
                self,
                "Export Complete",
                "Dataset exported to YR-MPEA"
            )
            logger.info("导出 dataset 到 YR-MPEA")
            
        except Exception as e:
            logger.error(f"导出到 YR-MPEA 失败: {e}")
            QMessageBox.critical(
                self,
                "Export Error",
                f"Failed to export to YR-MPEA:\n{e}"
            )