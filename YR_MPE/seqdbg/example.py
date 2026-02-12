#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SeqDBG使用示例

演示如何使用SeqDBG的核心功能。
"""

import sys
from pathlib import Path

# 添加项目路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QPushButton
from PyQt5.QtCore import Qt

from seqdbg import SeqDBGraphBuilder, AnnotationNormalizer
from seqdbg.viewer import SeqDBGViewer


class SeqDBGExampleWindow(QMainWindow):
    """SeqDBG示例窗口"""
    
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("SeqDBG 使用示例")
        self.setGeometry(100, 100, 1200, 800)
        
        # 创建构建器
        self.builder = SeqDBGraphBuilder()
        
        # 添加一些标准化规则
        self.builder.normalizer.add_synonym_rule("cox1", "COX1")
        self.builder.normalizer.add_synonym_rule("cob", "COB")
        self.builder.normalizer.add_synonym_rule("cox2", "COX2")
        
        # 创建示例数据
        self.create_example_data()
        
        # 创建UI
        self.init_ui()
    
    def create_example_data(self):
            """Create example data (simulating GeneBank parsing results)"""
            from seqdbg.models import GeneInfo
    
            # Simulate genome 1
            genes1 = [
                GeneInfo(name="cox1", start=100, end=200, strand=1),
                GeneInfo(name="cob", start=300, end=400, strand=1),
                GeneInfo(name="cox2", start=500, end=600, strand=-1),
            ]
    
            # Simulate genome 2 (with different annotations)
            genes2 = [
                GeneInfo(name="COX1", start=100, end=200, strand=1),
                GeneInfo(name="COB", start=300, end=400, strand=1),
                GeneInfo(name="nad5", start=500, end=600, strand=1),
            ]
    
            # Simulate genome 3 (sharing some annotations)
            genes3 = [
                GeneInfo(name="cox1", start=100, end=200, strand=1),
                GeneInfo(name="nad5", start=300, end=400, strand=1),
            ]
    
            # Add to graph
            for gene in genes1:
                gene.normalized_name = self.builder.normalizer.normalize(gene.name)
                self.builder.graph.add_node_occurrence("genome1", gene)
    
            for i, gene in enumerate(genes1[:-1]):
                next_gene = genes1[i + 1]
                self.builder.graph.add_edge(
                    gene.normalized_name,
                    self.builder.normalizer.normalize(next_gene.name)
                )
    
            for gene in genes2:
                gene.normalized_name = self.builder.normalizer.normalize(gene.name)
                self.builder.graph.add_node_occurrence("genome2", gene)
    
            for i, gene in enumerate(genes2[:-1]):
                next_gene = genes2[i + 1]
                self.builder.graph.add_edge(
                    gene.normalized_name,
                    self.builder.normalizer.normalize(next_gene.name)
                )
    
            for gene in genes3:
                gene.normalized_name = self.builder.normalizer.normalize(gene.name)
                self.builder.graph.add_node_occurrence("genome3", gene)
    
            for i, gene in enumerate(genes3[:-1]):
                next_gene = genes3[i + 1]
                self.builder.graph.add_edge(
                    gene.normalized_name,
                    self.builder.normalizer.normalize(next_gene.name)
                )
    
            # Apply normalization mapping
            self.builder._apply_normalization_mapping()
    
            print(f"Example data created:")
            print(f"  Nodes: {len(self.builder.graph.nodes)}")
            print(f"  Genomes: {len(self.builder.graph.genomes)}")
            print(f"  Edges: {len(self.builder.graph.get_edges())}")    
    def init_ui(self):
        """初始化UI"""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        layout = QVBoxLayout()
        central_widget.setLayout(layout)
        
        # 创建查看器（已经内置了完整的控制面板）
        self.viewer = SeqDBGViewer(self.builder.graph)
        layout.addWidget(self.viewer)
        
        # 添加一个简单的调试按钮
        button_layout = QHBoxLayout()
        self.print_info_btn = QPushButton("打印图信息到控制台")
        self.print_info_btn.clicked.connect(self.print_graph_info)
        button_layout.addWidget(self.print_info_btn)
        button_layout.addStretch()
        layout.addLayout(button_layout)
    
    def print_graph_info(self):
            """Print graph information"""
            print("\n" + "=" * 50)
            print("Graph Information:")
            print(f"  Nodes: {len(self.builder.graph.nodes)}")
            print(f"  Genomes: {len(self.builder.graph.genomes)}")
    
            print("\nNode List (sorted by genome count):")
            for node in self.builder.graph.get_nodes_sorted_by_count():
                print(f"  {node.name}: {node.count} genomes, {node.occurrences} occurrences")
    
            print("\nEdge List (sorted by count):")
            for source, target, count in self.builder.graph.get_edges_sorted_by_count():
                print(f"  {source} -> {target}: {count} times")
    
            print("=" * 50 + "\n")

def main():
    """主函数"""
    app = QApplication(sys.argv)
    
    window = SeqDBGExampleWindow()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()