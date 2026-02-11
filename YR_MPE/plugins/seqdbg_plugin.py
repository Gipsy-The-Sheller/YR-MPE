#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SeqDBG插件 - 序列数据库可视化工具

用于分析多基因组的注释邻接关系，发现结构变异和注释不统一。
"""

import logging
from typing import Optional, List, Dict
from pathlib import Path

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QMessageBox
)
from PyQt5.QtCore import pyqtSignal

from ..seqdbg import SeqDBGraphBuilder, AnnotationNormalizer
from ..seqdbg.table_viewer import SeqDBGTableViewer

logger = logging.getLogger(__name__)


class SeqDBGWindow(QMainWindow):
    """
    SeqDBG主窗口
    
    采用SeqMatrix的设计模式，独立运行或从YR-MPEA调用
    """
    
    # 信号定义
    import_dataset_signal = pyqtSignal(list)  # 发送DatasetItem列表
    
    def __init__(self, import_from=None):
        """初始化SeqDBG窗口"""
        super().__init__()
        
        # Track if opened from YR-MPEA
        self.opened_from_yrmpea = import_from == "YR_MPEA"
        
        # 初始化核心组件
        self.builder = SeqDBGraphBuilder()
        self.table_viewer: Optional[SeqDBGTableViewer] = None
        
        # 文件和entry追踪
        self.loaded_files: List[str] = []  # 已加载的文件列表
        self.entry_to_file_map: Dict[str, str] = {}  # {entry_id: file_path}
        
        # 初始化UI
        self.init_ui()
        
        # 添加一些标准化规则
        self.add_default_rules()
    
    def get_window_title(self) -> str:
        """返回窗口标题"""
        if self.opened_from_yrmpea:
            return "SeqDBG - Annotation Graph [YR-MPEA]"
        return "SeqDBG - Annotation Graph"
    
    def init_ui(self):
        """初始化用户界面"""
        self.setWindowTitle(self.get_window_title())
        self.setGeometry(100, 100, 1400, 900)
        
        main_widget = QWidget()
        main_layout = QVBoxLayout()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)
        
        # 创建表格查看器
        self.table_viewer = SeqDBGTableViewer()
        self.viewer_layout = QVBoxLayout()
        self.viewer_layout.setContentsMargins(0, 0, 0, 0)
        self.viewer_layout.addWidget(self.table_viewer)
        main_layout.addLayout(self.viewer_layout)
        
        # 连接信号
        self.table_viewer.taxon_manager.filesAdded.connect(self.on_files_added)
        self.table_viewer.taxon_manager.entriesRemoved.connect(self.on_entries_removed)
    
    
    
    def add_default_rules(self):
        """添加默认的常见规则"""
        common_rules = {
            "cox1": "COX1",
            "cox2": "COX2",
            "cox3": "COX3",
            "cob": "COB",
            "nad1": "NAD1",
            "nad2": "NAD2",
            "nad3": "NAD3",
            "nad4": "NAD4",
            "nad4L": "NAD4L",
            "nad5": "NAD5",
            "nad6": "NAD6",
            "atp6": "ATP6",
            "atp8": "ATP8",
        }
        
        for alias, canonical in common_rules.items():
            self.builder.normalizer.add_synonym_rule(alias, canonical)
    
    def on_files_added(self, files: List[str]):
        """
        处理添加的GenBank文件
        
        Args:
            files: 文件路径列表
        """
        try:
            # 存储文件
            self.loaded_files.extend(files)
            
            # 构建图
            self.builder.reset()
            graph = self.builder.build_multi_genome_graph(
                self.loaded_files,
                extract_sequences=False  # 不提取序列以加快速度
            )
            
            # 统计信息
            stats = self.builder.get_stats()
            logger.info(f"Graph built: {stats}")
            
            # 更新查看器
            self.update_viewer(graph)
            
            # 添加entries到列表和映射
            for genome_id, entry_ids in graph.genome_entries.items():
                file_path = next((f for f in files if Path(f).stem == genome_id), files[0])
                for entry_id in entry_ids:
                    self.entry_to_file_map[entry_id] = file_path
                    self.table_viewer.taxon_manager.add_entry(entry_id)
            
            logger.info(f"成功加载 {len(files)} 个文件")
            
        except Exception as e:
            logger.error(f"文件加载失败: {e}")
            QMessageBox.critical(self, "Error", f"Failed to load files:\n{e}")
    
    def on_entries_removed(self, entry_ids: List[str]):
        """
        处理移除的Entries
        
        Args:
            entry_ids: Entry ID列表
        """
        logger.info(f"移除 {len(entry_ids)} 个entries")
        
        # 从映射中移除
        for entry_id in entry_ids:
            if entry_id in self.entry_to_file_map:
                del self.entry_to_file_map[entry_id]
        
        # 获取当前选中的entries
        current_entries = self.table_viewer.taxon_manager.get_all_entries()
        
        if not current_entries:
            # 没有entries了，清空表格
            self.table_viewer.table.setRowCount(0)
            self.table_viewer.sankey_canvas.clear_plot()
            return
        
        # 重新构建图（只包含选中的entries）
        try:
            self.builder.reset()
            graph = self.builder.build_multi_genome_graph(
                self.loaded_files,
                extract_sequences=False
            )
            
            # 更新查看器，只显示选中的entries
            self.update_viewer(graph, filter_entries=current_entries)
            
            logger.info(f"表格已更新，当前 {len(current_entries)} 个entries")
            
        except Exception as e:
            logger.error(f"更新表格失败: {e}")
            QMessageBox.critical(self, "Error", f"Failed to update table:\n{e}")
    
    def update_viewer(self, graph, filter_entries: Optional[List[str]] = None):
        """
        更新查看器
        
        Args:
            graph: 注释图
            filter_entries: 可选，只显示这些entries中的基因
        """
        # 获取统计数据
        stats = self.builder.get_gene_stats(filter_entries=filter_entries)
        self.table_viewer.set_stats_data(stats)
        
        # 获取邻接数据
        adjacency_data = {}
        for node_name in graph.nodes.keys():
            adjacency_data[node_name] = self.builder.get_adjacency_info(node_name)
        self.table_viewer.set_adjacency_data(adjacency_data)
        
        # 设置entries（不再需要维护all_entries列表）
        self.table_viewer.set_all_entries([])
        
        logger.info(f"查看器更新完成: {len(stats)} 个基因")
    
    
    
    
    
    
    
    def closeEvent(self, event):
        """处理关闭事件"""
        # 对于新的表格查看器，暂时不处理选择导出
        # 正常关闭
        super().closeEvent(event)
    
    def cleanup(self):
        """清理资源"""
        if self.builder:
            self.builder.clear_cache()
        if self.table_viewer:
            self.table_viewer.clear_all()
        logger.info("SeqDBG cleanup completed")


class SeqDBGPluginEntry:
    """
    SeqDBG插件入口
    
    采用SeqMatrix的设计模式
    """
    
    def __init__(self):
        self.window = None
        
    def run(self, import_from=None):
        """运行SeqDBG插件"""
        return SeqDBGWindow(import_from=import_from)


if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    
    app = QApplication(sys.argv)
    window = SeqDBGWindow()
    window.show()
    sys.exit(app.exec_())
