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
from typing import List, Dict, Optional, Tuple

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem,
    QHeaderView, QSplitter, QLabel, QPushButton, QGroupBox, QMessageBox,
    QListWidget, QFileDialog
)
from PyQt5.QtCore import Qt, pyqtSignal, QSize
from PyQt5.QtGui import QColor, QFont, QBrush

from .models import GeneStats
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path

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
            # Gene Name
            name_item = QTableWidgetItem(stat.gene_name)
            name_item.setTextAlignment(Qt.AlignLeft | Qt.AlignVCenter)
            name_item.setFont(QFont("Arial", 9, QFont.Bold))
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
        taxon_data: List[Tuple[str, int]],   # [(taxon_name, count), ...]
        right_data: List[Tuple[str, int]],   # [(node_name, count), ...]
        left_matrix: List[List[int]],         # [left][taxon]
        right_matrix: List[List[int]],        # [taxon][right]
        gene_name: str = ""
    ):
        """
        绘制桑基图
        
        Args:
            left_data: 左邻接节点数据
            taxon_data: Taxon数据
            right_data: 右邻接节点数据
            left_matrix: 左到taxon的连接矩阵
            right_matrix: Taxon到右的连接矩阵
            gene_name: 当前基因名称（用于标题）
        """
        self.ax.clear()
        
        # 计算总高度
        left_values = [v for _, v in left_data]
        taxon_values = [v for _, v in taxon_data]
        right_values = [v for _, v in right_data]
        
        total_height = max(
            sum(left_values) if left_values else 1,
            sum(taxon_values) if taxon_values else 1,
            sum(right_values) if right_values else 1
        )
        
        # 归一化高度
        left_heights = np.array(left_values) / total_height
        taxon_heights = np.array(taxon_values) / total_height
        right_heights = np.array(right_values) / total_height
        
        # 计算起始位置
        left_starts = np.concatenate(([0], np.cumsum(left_heights)[:-1]))
        taxon_starts = np.concatenate(([0], np.cumsum(taxon_heights)[:-1]))
        right_starts = np.concatenate(([0], np.cumsum(right_heights)[:-1]))
        
        # 颜色映射
        num_left = len(left_data)
        num_taxon = len(taxon_data)
        num_right = len(right_data)
        
        left_colors = [plt.cm.Blues(0.3 + 0.5 * i / max(1, num_left - 1)) for i in range(num_left)]
        taxon_colors = [plt.cm.Reds(0.3 + 0.5 * i / max(1, num_taxon - 1)) for i in range(num_taxon)]
        right_colors = [plt.cm.Greens(0.3 + 0.5 * i / max(1, num_right - 1)) for i in range(num_right)]
        
        # 绘制连接（流线）
        # 左 -> Taxon
        if left_matrix and taxon_data:
            left_matrix = np.array(left_matrix, dtype=float)
            for i in range(len(left_data)):
                y_start = left_starts[i]
                h_start = left_heights[i]
                for j in range(len(taxon_data)):
                    flow = left_matrix[i, j]
                    if flow <= 0:
                        continue
                    
                    y_end = taxon_starts[j]
                    h_end = taxon_heights[j]
                    
                    # 流线宽度
                    line_width = max(0.5, 3.0 * flow / total_height * 100)
                    
                    # 贝塞尔曲线
                    x_start = 0
                    x_end = 1
                    cp_x = 0.5
                    path_data = [
                        (Path.MOVETO, (x_start, y_start + h_start / 2)),
                        (Path.CURVE4, (cp_x, y_start + h_start / 2)),
                        (Path.CURVE4, (cp_x, y_end + h_end / 2)),
                        (Path.CURVE4, (x_end, y_end + h_end / 2)),
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
        
        # Taxon -> 右
        if right_matrix and right_data:
            right_matrix = np.array(right_matrix, dtype=float)
            for i in range(len(taxon_data)):
                y_start = taxon_starts[i]
                h_start = taxon_heights[i]
                for j in range(len(right_data)):
                    flow = right_matrix[i, j]
                    if flow <= 0:
                        continue
                    
                    y_end = right_starts[j]
                    h_end = right_heights[j]
                    
                    # 流线宽度
                    line_width = max(0.5, 3.0 * flow / total_height * 100)
                    
                    # 贝塞尔曲线
                    x_start = 1
                    x_end = 2
                    cp_x = 1.5
                    path_data = [
                        (Path.MOVETO, (x_start, y_start + h_start / 2)),
                        (Path.CURVE4, (cp_x, y_start + h_start / 2)),
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
            self.ax.text(
                0,
                y_start + h / 2,
                str(count),
                ha='center',
                va='center',
                fontsize=8,
                fontweight='bold',
                color='white' if left_colors[i][0] < 0.5 else 'black',
                zorder=3
            )
        
        # Taxon节点
        for i, (name, count) in enumerate(taxon_data):
            y_start = taxon_starts[i]
            h = taxon_heights[i]
            rect = FancyBboxPatch(
                (1 - bar_width / 2, y_start),
                bar_width,
                h,
                boxstyle="round,pad=0,rounding_size=0.01",
                facecolor=taxon_colors[i],
                edgecolor='black',
                linewidth=1.5,
                zorder=2
            )
            self.ax.add_patch(rect)
            self.ax.text(
                1 - bar_width / 2 - 0.02,
                y_start + h / 2,
                name,
                ha='right',
                va='center',
                fontsize=9,
                rotation=0
            )
            self.ax.text(
                1,
                y_start + h / 2,
                str(count),
                ha='center',
                va='center',
                fontsize=8,
                fontweight='bold',
                color='white' if taxon_colors[i][0] < 0.5 else 'black',
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
            self.ax.text(
                2,
                y_start + h / 2,
                str(count),
                ha='center',
                va='center',
                fontsize=8,
                fontweight='bold',
                color='white' if right_colors[i][0] < 0.5 else 'black',
                zorder=3
            )
        
        # 设置标题
        title = f"Gene Adjacency Sankey: {gene_name}" if gene_name else "Gene Adjacency Sankey"
        self.ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # 添加层标签
        self.ax.text(0, 1.02, "Left Adjacent", ha='center', fontsize=11, fontweight='bold', transform=self.ax.transAxes)
        self.ax.text(0.5, 1.02, "Taxon (Entry)", ha='center', fontsize=11, fontweight='bold', transform=self.ax.transAxes)
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
    - 右侧：Sankey图
    """
    
    # 信号定义
    geneSelected = pyqtSignal(str)  # gene_name
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.stats_data: List[GeneStats] = []
        self.adjacency_data: Dict[str, Dict] = {}
        
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
        
        # 右侧：Sankey图（无GroupBox）
        sankey_container = QWidget()
        sankey_layout = QVBoxLayout()
        sankey_layout.setContentsMargins(0, 0, 0, 0)
        sankey_container.setLayout(sankey_layout)
        
        self.sankey_canvas = SankeyPlotCanvas()
        sankey_layout.addWidget(self.sankey_canvas)
        
        main_splitter.addWidget(sankey_container)
        
        # 设置分割器比例
        # 主分割器：左侧占40%，右侧占60%
        main_splitter.setStretchFactor(0, 40)
        main_splitter.setStretchFactor(1, 60)
        
        # 左侧分割器：taxon管理占30%，表格占70%
        left_splitter.setStretchFactor(0, 30)
        left_splitter.setStretchFactor(1, 70)
    
    def set_stats_data(self, stats: List[GeneStats]):
        """
        设置统计数据
        
        Args:
            stats: GeneStats列表
        """
        self.stats_data = stats
        self.table.set_data(stats)
    
    def set_adjacency_data(self, adjacency_data: Dict[str, Dict]):
        """
        设置邻接数据
        
        Args:
            adjacency_data: {gene_name: adjacency_info} 字典
        """
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
        taxon_data = adj_info.get("taxon_info", [])
        right_data = adj_info.get("right_adjacent", [])
        
        # 构建连接矩阵
        # 左 -> Taxon 矩阵
        left_matrix = []
        if left_data and taxon_data:
            # 简化处理：假设每个左节点平均分配到所有taxon
            num_left = len(left_data)
            num_taxon = len(taxon_data)
            left_matrix = [[0] * num_taxon for _ in range(num_left)]
            for i, (name, count) in enumerate(left_data):
                for j in range(num_taxon):
                    left_matrix[i][j] = count // num_taxon if num_taxon > 0 else 0
        
        # Taxon -> 右 矩阵
        right_matrix = []
        if taxon_data and right_data:
            # 简化处理：假设每个taxon平均分配到所有右节点
            num_taxon = len(taxon_data)
            num_right = len(right_data)
            right_matrix = [[0] * num_right for _ in range(num_taxon)]
            for i, (name, count) in enumerate(taxon_data):
                for j in range(num_right):
                    right_matrix[i][j] = count // num_right if num_right > 0 else 0
        
        # 绘制桑基图
        self.sankey_canvas.plot_sankey(
            left_data=left_data,
            taxon_data=taxon_data,
            right_data=right_data,
            left_matrix=left_matrix,
            right_matrix=right_matrix,
            gene_name=gene_name
        )
    
    def clear_all(self):
        """清空所有数据"""
        self.stats_data.clear()
        self.adjacency_data.clear()
        self.table.setRowCount(0)
        self.taxon_manager.clear()
        self.sankey_canvas.clear_plot()