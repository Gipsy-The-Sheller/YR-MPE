#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SeqDBG可视化查看器

提供交互式的注释邻接图可视化界面。
"""

import logging
from typing import Dict, List, Optional, Set, Tuple
from math import sqrt

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGraphicsView, QGraphicsScene,
    QGraphicsItem, QGraphicsEllipseItem, QGraphicsLineItem, QToolButton,
    QLabel, QPushButton, QGroupBox, QTextEdit, QFileDialog, QSplitter,
    QCheckBox, QSlider, QListWidget, QTabWidget, QComboBox, QInputDialog,
    QMessageBox
)
from PyQt5.QtCore import Qt, QPointF, QRectF, QTimer, pyqtSignal, QObject
from PyQt5.QtGui import (
    QPainter, QPen, QBrush, QColor, QFont, QPainterPath,
    QCursor, QWheelEvent, QMouseEvent
)

from .models import AnnotationGraph, AnnotationNode
from .layout import ForceDirectedLayout, NodePosition

logger = logging.getLogger(__name__)


class ViewerSignals(QObject):
    """查看器信号类"""
    nodeClicked = pyqtSignal(str, bool)  # (node_name, selected)
    nodeDoubleClicked = pyqtSignal(str)  # node_name
    selectionChanged = pyqtSignal(list)  # [node_names]
    exportRequested = pyqtSignal(list, str)  # (selected_nodes, export_type)


class GraphNodeItem(QGraphicsEllipseItem):
    """
    图形节点项
    
    表示图中的一个注释节点
    """
    
    def __init__(
        self, 
        name: str, 
        node: AnnotationNode, 
        radius: float = 20.0,
        parent=None
    ):
        super().__init__(-radius, -radius, radius * 2, radius * 2, parent)
        
        self.node_name = name
        self.node = node
        self.radius = radius
        self.base_radius = radius
        
        # 状态
        self.selected = False
        self.hovered = False
        
        # 回调函数
        self.on_clicked = None  # (node_name, event) -> None
        self.on_double_clicked = None  # (node_name) -> None
        
        # 设置属性
        self.setAcceptHoverEvents(True)
        self.setFlag(QGraphicsItem.ItemIsSelectable)
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges)
        
        # 更新样式
        self.update_appearance()
    
    def update_appearance(self):
        """更新节点外观"""
        # 计算颜色（基于覆盖率和选择状态）
        coverage = self.node.count
        total_genomes = len(self.node.sequences) if hasattr(self.node, 'sequences') else 1
        
        if self.selected:
            # 选中状态：金黄色
            color = QColor(255, 215, 0)
            border_color = QColor(255, 140, 0)
            border_width = 3
        elif self.hovered:
            # 悬停状态：浅蓝色
            color = QColor(173, 216, 230)
            border_color = QColor(70, 130, 180)
            border_width = 2
        else:
            # 正常状态：根据覆盖率变化
            intensity = min(255, 100 + coverage * 30)
            color = QColor(100, 150, intensity)
            border_color = QColor(50, 50, 50)
            border_width = 1
        
        # 设置画刷和画笔
        self.setBrush(QBrush(color))
        self.setPen(QPen(border_color, border_width))
        
        # 更新大小
        scale_factor = 1.5 if self.selected or self.hovered else 1.0
        new_radius = self.base_radius * scale_factor
        
        if abs(self.radius - new_radius) > 0.1:
            self.radius = new_radius
            self.setRect(-self.radius, -self.radius, self.radius * 2, self.radius * 2)
    
    def hoverEnterEvent(self, event):
        """Mouse hover enter"""
        self.hovered = True
        self.update_appearance()
        self.setToolTip(f"{self.node_name}\nGenomes: {self.node.count}\nOccurrences: {self.node.occurrences}")
        super().hoverEnterEvent(event)
    
    def hoverLeaveEvent(self, event):
        """鼠标悬停离开"""
        self.hovered = False
        self.update_appearance()
        super().hoverLeaveEvent(event)
    
    def mousePressEvent(self, event):
        """鼠标点击事件"""
        modifiers = event.modifiers()
        
        if modifiers & Qt.ControlModifier:
            # Ctrl+点击：切换选择
            self.selected = not self.selected
        elif modifiers & Qt.ShiftModifier:
            # Shift+点击：添加选择
            self.selected = True
        else:
            # 普通点击：单选
            self.selected = True
        
        self.update_appearance()
        
        # 调用回调函数
        if self.on_clicked:
            self.on_clicked(self.node_name, event)
        
        super().mousePressEvent(event)
    
    def mouseDoubleClickEvent(self, event):
        """鼠标双击事件"""
        if self.on_double_clicked:
            self.on_double_clicked(self.node_name)
        super().mouseDoubleClickEvent(event)


class GraphEdgeItem(QGraphicsLineItem):
    """
    图形边项
    
    表示图中的连接关系
    """
    
    def __init__(
        self, 
        source_item: GraphNodeItem, 
        target_item: GraphNodeItem,
        count: int = 1,
        parent=None
    ):
        super().__init__(parent)
        
        self.source_item = source_item
        self.target_item = target_item
        self.count = count
        
        # 设置属性
        self.setZValue(-1)  # 在节点下方
        self.update_appearance()
    
    def update_appearance(self):
        """更新边外观"""
        # 根据计数设置粗细
        width = max(0.5, min(4, 0.5 + self.count * 0.3))
        
        # 颜色基于计数
        intensity = min(255, 150 + self.count * 20)
        color = QColor(intensity, intensity, intensity)
        
        self.setPen(QPen(color, width))
    
    def update_position(self):
        """更新边的位置"""
        line = self.line()
        line.setP1(self.source_item.pos())
        line.setP2(self.target_item.pos())
        self.setLine(line)


class SeqDBGGraphicsView(QGraphicsView):
    """
    图形视图
    
    支持缩放、拖拽等交互
    """
    
    def __init__(self, scene: QGraphicsScene, parent=None):
        super().__init__(scene, parent)
        
        self.setRenderHint(QPainter.Antialiasing)
        self.setDragMode(QGraphicsView.RubberBandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorUnderMouse)
        
        # 背景色
        self.setBackgroundBrush(QBrush(QColor(245, 245, 245)))
    
    def wheelEvent(self, event: QWheelEvent):
        """鼠标滚轮缩放"""
        zoom_in_factor = 1.15
        zoom_out_factor = 1 / zoom_in_factor
        
        if event.angleDelta().y() > 0:
            zoom_factor = zoom_in_factor
        else:
            zoom_factor = zoom_out_factor
        
        self.scale(zoom_factor, zoom_factor)


class SeqDBGViewer(QWidget):
    """
    SeqDBG查看器主窗口
    
    提供完整的可视化界面，包括：
    - 图形显示区域
    - 节点信息面板
    - 控制面板
    - 导出功能
    """
    
    def __init__(self, graph: AnnotationGraph, parent=None):
        super().__init__(parent)
        
        self.graph = graph
        self.layout_engine = ForceDirectedLayout()
        self.signals = ViewerSignals()
        
        # 可视化项
        self.node_items: Dict[str, GraphNodeItem] = {}
        self.edge_items: List[GraphEdgeItem] = []
        self.selected_nodes: Set[str] = set()
        
        # 初始化UI
        self.init_ui()
        
        # 构建图形
        self.build_graph()
    
    def init_ui(self):
        """初始化用户界面"""
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # 创建分割器
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # 左侧：图形视图
        self.scene = QGraphicsScene()
        self.view = SeqDBGGraphicsView(self.scene)
        splitter.addWidget(self.view)
        
        # 右侧：信息面板
        info_panel = self.create_info_panel()
        splitter.addWidget(info_panel)
        
        # 设置分割器比例
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)
        
        # 底部：控制面板
        control_panel = self.create_control_panel()
        main_layout.addWidget(control_panel)
    
    def create_info_panel(self) -> QWidget:
        """Create information panel"""
        panel = QGroupBox("Node Information")
        layout = QVBoxLayout()
        panel.setLayout(layout)
        
        # Basic information
        self.node_info_label = QLabel("Select a node to view details")
        self.node_info_label.setWordWrap(True)
        layout.addWidget(self.node_info_label)
        
        # Genome information
        layout.addWidget(QLabel("Genomes containing this node:"))
        self.genomes_list = QListWidget()
        self.genomes_list.setMaximumHeight(150)
        layout.addWidget(self.genomes_list)
        
        # Connection information
        layout.addWidget(QLabel("Adjacent nodes:"))
        self.connections_list = QListWidget()
        self.connections_list.setMaximumHeight(150)
        layout.addWidget(self.connections_list)
        
        # Sequence information
        self.has_sequences_label = QLabel("")
        layout.addWidget(self.has_sequences_label)
        
        layout.addStretch()
        return panel
    
    def create_control_panel(self) -> QWidget:
        """Create control panel"""
        panel = QWidget()
        layout = QHBoxLayout()
        panel.setLayout(layout)
        
        # Selection controls
        select_group = QGroupBox("Selection")
        select_layout = QHBoxLayout()
        select_group.setLayout(select_layout)
        
        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.clicked.connect(self.select_all_nodes)
        select_layout.addWidget(self.select_all_btn)
        
        self.clear_selection_btn = QPushButton("Clear Selection")
        self.clear_selection_btn.clicked.connect(self.clear_selection)
        select_layout.addWidget(self.clear_selection_btn)
        
        layout.addWidget(select_group)
        
        # Export controls
        export_group = QGroupBox("Export")
        export_layout = QHBoxLayout()
        export_group.setLayout(export_layout)
        
        self.export_dataset_btn = QPushButton("Export to Dataset")
        self.export_dataset_btn.clicked.connect(lambda: self.export_nodes("dataset"))
        export_layout.addWidget(self.export_dataset_btn)
        
        self.export_fasta_btn = QPushButton("Export to FASTA")
        self.export_fasta_btn.clicked.connect(lambda: self.export_nodes("fasta"))
        export_layout.addWidget(self.export_fasta_btn)
        
        layout.addWidget(export_group)
        
        # Layout controls
        layout_group = QGroupBox("Layout")
        layout_layout = QHBoxLayout()
        layout_group.setLayout(layout_layout)
        
        self.relayout_btn = QPushButton("Relayout")
        self.relayout_btn.clicked.connect(self.rebuild_graph)
        layout_layout.addWidget(self.relayout_btn)
        
        layout.addWidget(layout_group)
        
        layout.addStretch()
        return panel
    
    def build_graph(self):
        """构建图形"""
        # 清空场景
        self.scene.clear()
        self.node_items.clear()
        self.edge_items.clear()
        self.selected_nodes.clear()
        
        if not self.graph.nodes:
            logger.warning("图中没有节点")
            return
        
        # 获取节点和边
        node_names = list(self.graph.nodes.keys())
        edges = self.graph.get_edges()
        
        # 执行布局
        positions = self.layout_engine.layout(node_names, edges)
        
        # 居中布局
        self.layout_engine.center_layout()
        
        # 创建节点项
        for name, pos in positions.items():
            node = self.graph.nodes[name]
            
            # 计算节点大小（基于覆盖率）
            base_radius = max(15, min(40, 20 + node.count * 2))
            
            node_item = GraphNodeItem(name, node, base_radius)
            node_item.setPos(pos.x, pos.y)
            
            # 设置回调函数
            node_item.on_clicked = self.on_node_clicked
            node_item.on_double_clicked = self.on_node_double_clicked
            
            self.scene.addItem(node_item)
            self.node_items[name] = node_item
        
        # 创建边项
        for source, target, count in edges:
            if source in self.node_items and target in self.node_items:
                edge_item = GraphEdgeItem(
                    self.node_items[source],
                    self.node_items[target],
                    count
                )
                edge_item.update_position()
                self.scene.addItem(edge_item)
                self.edge_items.append(edge_item)
        
        # 调整场景大小
        self.adjust_scene_size()
        
        logger.info(f"图形构建完成: {len(self.node_items)} 个节点, {len(self.edge_items)} 条边")
    
    def adjust_scene_size(self):
        """调整场景大小"""
        if not self.node_items:
            return
        
        # 计算边界
        min_x = min(item.pos().x() - item.radius for item in self.node_items.values())
        max_x = max(item.pos().x() + item.radius for item in self.node_items.values())
        min_y = min(item.pos().y() - item.radius for item in self.node_items.values())
        max_y = max(item.pos().y() + item.radius for item in self.node_items.values())
        
        # 添加边距
        margin = 100
        self.scene.setSceneRect(
            min_x - margin,
            min_y - margin,
            max_x - min_x + 2 * margin,
            max_y - min_y + 2 * margin
        )
    
    def on_node_clicked(self, node_name: str, event):
        """节点点击事件"""
        # 清除其他节点的选择（如果不是多选模式）
        modifiers = event.modifiers()
        if not (modifiers & Qt.ControlModifier) and not (modifiers & Qt.ShiftModifier):
            for name, item in self.node_items.items():
                if name != node_name:
                    item.selected = False
                    item.update_appearance()
        
        # 更新选择状态
        self.selected_nodes.clear()
        for name, item in self.node_items.items():
            if item.selected:
                self.selected_nodes.add(name)
        
        # 更新信息面板
        if self.selected_nodes:
            self.update_info_panel(node_name)
        else:
            self.clear_info_panel()
        
        # 发送信号
        is_selected = self.node_items[node_name].selected
        self.signals.nodeClicked.emit(node_name, is_selected)
        self.signals.selectionChanged.emit(list(self.selected_nodes))
    
    def on_node_double_clicked(self, node_name: str):
        """节点双击事件"""
        self.signals.nodeDoubleClicked.emit(node_name)
    
    def update_info_panel(self, node_name: str):
        """更新信息面板"""
        if node_name not in self.graph.nodes:
            return
        
        node = self.graph.nodes[node_name]
        
        # 基本信息
        self.node_info_label.setText(
            f"<b>注释名称:</b> {node_name}<br>"
            f"<b>基因组数:</b> {node.count}<br>"
            f"<b>总出现次数:</b> {node.occurrences}"
        )
        
        # 基因组列表
        self.genomes_list.clear()
        for genome_id in node.get_unique_genomes():
            occurrences = len(node.genomes.get(genome_id, []))
            self.genomes_list.addItem(f"{genome_id} ({occurrences} 次)")
        
        # 连接列表
        self.connections_list.clear()
        for target_name, count in sorted(node.connections.items(), key=lambda x: -x[1]):
            self.connections_list.addItem(f"{target_name}: {count} 次")
        
        # Sequence information
        if node.sequences:
            total_seqs = sum(len(seqs) for seqs in node.sequences.values())
            self.has_sequences_label.setText(f"✓ Contains {total_seqs} sequences")
            self.has_sequences_label.setStyleSheet("color: green; font-weight: bold;")
        else:
            self.has_sequences_label.setText("✗ No sequence data")
            self.has_sequences_label.setStyleSheet("color: red;")
    
    def clear_info_panel(self):
        """Clear information panel"""
        self.node_info_label.setText("Select a node to view details")
        self.genomes_list.clear()
        self.connections_list.clear()
        self.has_sequences_label.setText("")
    
    def select_all_nodes(self):
        """选择所有节点"""
        for node_item in self.node_items.values():
            node_item.selected = True
            node_item.update_appearance()
        
        self.selected_nodes = set(self.node_items.keys())
        self.signals.selectionChanged.emit(list(self.selected_nodes))
    
    def clear_selection(self):
        """清除选择"""
        for node_item in self.node_items.values():
            node_item.selected = False
            node_item.update_appearance()
        
        self.selected_nodes.clear()
        self.clear_info_panel()
        self.signals.selectionChanged.emit([])
    
    def rebuild_graph(self):
        """重新构建图形"""
        self.build_graph()
    
    def export_nodes(self, export_type: str):
        """Export selected nodes"""
        if not self.selected_nodes:
            QMessageBox.warning(self, "Warning", "No nodes selected")
            return
        
        selected_list = list(self.selected_nodes)
        self.signals.exportRequested.emit(selected_list, export_type)
