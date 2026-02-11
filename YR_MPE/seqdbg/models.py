#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SeqDBG核心数据模型

定义注释节点、邻接图和基因信息等核心数据结构。
所有模型都是纯Python类，不依赖Qt，便于测试和复用。
"""

from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from Bio.SeqRecord import SeqRecord


@dataclass
class GeneInfo:
    """
    基因信息数据类
    
    存储从GeneBank文件中提取的基因信息
    """
    name: str                          # 基因名称（原始）
    start: int                         # 起始位置（0-based）
    end: int                           # 结束位置（0-based, exclusive）
    strand: int                        # 链方向: 1 (正链) 或 -1 (负链)
    entry_id: str = ""                 # 所属的Entry ID
    normalized_name: str = ""          # 标准化后的名称
    sequence: Optional[SeqRecord] = None  # 基因序列（可选，延迟加载）
    
    def __post_init__(self):
        """初始化后处理"""
        if not self.normalized_name:
            self.normalized_name = self.name
    
    @property
    def length(self) -> int:
        """基因长度"""
        return self.end - self.start
    
    def __str__(self) -> str:
        strand_str = "+" if self.strand == 1 else "-"
        return f"{self.name}[{self.start}-{self.end}:{strand_str}]"
    
    def __repr__(self) -> str:
        return (f"GeneInfo(name='{self.name}', start={self.start}, "
                f"end={self.end}, strand={self.strand}, entry_id='{self.entry_id}')")


class AnnotationNode:
    """
    注释节点
    
    代表一个标准化的基因注释，可能出现在多个基因组中
    """
    
    def __init__(self, name: str):
        """
        初始化注释节点
        
        Args:
            name: 标准化的注释名称
        """
        self.name = name
        self.count = 0                          # 在多少个基因组中出现过
        self.occurrences = 0                    # 总出现次数（一个基因组可能多次）
        self.genomes: Dict[str, List[GeneInfo]] = {}  # {genome_id: [GeneInfo]}
        self.connections: Dict[str, int] = {}   # {target_name: edge_count}
        self.sequences: Dict[str, List[SeqRecord]] = {}  # {genome_id: [SeqRecord]}
    
    def add_occurrence(self, genome_id: str, gene_info: GeneInfo):
        """
        添加一次出现
        
        Args:
            genome_id: 基因组ID
            gene_info: 基因信息
        """
        if genome_id not in self.genomes:
            self.genomes[genome_id] = []
            self.count += 1
        self.genomes[genome_id].append(gene_info)
        self.occurrences += 1
    
    def add_sequence(self, genome_id: str, sequence: SeqRecord):
        """
        添加序列
        
        Args:
            genome_id: 基因组ID
            sequence: 序列记录
        """
        if genome_id not in self.sequences:
            self.sequences[genome_id] = []
        self.sequences[genome_id].append(sequence)
    
    def add_connection(self, target_name: str):
        """
        添加连接（边）
        
        Args:
            target_name: 目标节点名称
        """
        if target_name not in self.connections:
            self.connections[target_name] = 0
        self.connections[target_name] += 1
    
    def merge_with(self, other: 'AnnotationNode'):
        """
        合并另一个节点
        
        Args:
            other: 要合并的节点
        """
        self.count += other.count
        self.occurrences += other.occurrences
        
        # 合并基因组信息
        for genome_id, gene_infos in other.genomes.items():
            if genome_id not in self.genomes:
                self.genomes[genome_id] = []
            self.genomes[genome_id].extend(gene_infos)
        
        # 合并序列信息
        for genome_id, sequences in other.sequences.items():
            if genome_id not in self.sequences:
                self.sequences[genome_id] = []
            self.sequences[genome_id].extend(sequences)
        
        # 合并连接
        for target_name, count in other.connections.items():
            if target_name not in self.connections:
                self.connections[target_name] = 0
            self.connections[target_name] += count
    
    def get_all_sequences(self) -> List[SeqRecord]:
        """
        获取所有序列
        
        Returns:
            所有序列记录的列表
        """
        all_seqs = []
        for sequences in self.sequences.values():
            all_seqs.extend(sequences)
        return all_seqs
    
    def get_unique_genomes(self) -> Set[str]:
        """
        获取包含此节点的唯一基因组集合
        
        Returns:
            基因组ID集合
        """
        return set(self.genomes.keys())
    
    def __str__(self) -> str:
        return f"AnnotationNode('{self.name}', count={self.count}, occurrences={self.occurrences})"
    
    def __repr__(self) -> str:
        return f"AnnotationNode(name='{self.name}', count={self.count})"
    
    def __hash__(self) -> int:
        return hash(self.name)
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, AnnotationNode):
            return False
        return self.name == other.name


class AnnotationGraph:
    """
    注释邻接图
    
    管理所有注释节点和它们之间的邻接关系
    """
    
    def __init__(self):
        """初始化图"""
        self.nodes: Dict[str, AnnotationNode] = {}  # {name: AnnotationNode}
        self.genomes: List[str] = []                # 基因组ID列表
        self.genome_entries: Dict[str, List[str]] = {}  # {genome_id: [entry_ids]}
    
    def add_genome(self, genome_id: str):
        """
        添加基因组
        
        Args:
            genome_id: 基因组ID
        """
        if genome_id not in self.genomes:
            self.genomes.append(genome_id)
            self.genome_entries[genome_id] = []
    
    def add_entry(self, genome_id: str, entry_id: str):
        """
        添加条目到基因组
        
        Args:
            genome_id: 基因组ID
            entry_id: 条目ID
        """
        if genome_id not in self.genome_entries:
            self.genome_entries[genome_id] = []
        if entry_id not in self.genome_entries[genome_id]:
            self.genome_entries[genome_id].append(entry_id)
    
    def get_or_create_node(self, name: str) -> AnnotationNode:
        """
        获取或创建节点
        
        Args:
            name: 节点名称
            
        Returns:
            AnnotationNode实例
        """
        if name not in self.nodes:
            self.nodes[name] = AnnotationNode(name)
        return self.nodes[name]
    
    def add_node_occurrence(self, genome_id: str, gene_info: GeneInfo):
        """
        添加节点出现
        
        Args:
            genome_id: 基因组ID
            gene_info: 基因信息
        """
        node = self.get_or_create_node(gene_info.normalized_name)
        node.add_occurrence(genome_id, gene_info)
    
    def add_edge(self, source_name: str, target_name: str):
        """
        添加边（邻接关系）
        
        Args:
            source_name: 源节点名称
            target_name: 目标节点名称
        """
        source_node = self.get_or_create_node(source_name)
        target_node = self.get_or_create_node(target_name)
        source_node.add_connection(target_name)
    
    def merge_nodes(self, from_name: str, to_name: str):
        """
        合并两个节点（用于注释标准化）
        
        Args:
            from_name: 源节点名称
            to_name: 目标节点名称
        """
        if from_name not in self.nodes or to_name not in self.nodes:
            return
        
        from_node = self.nodes[from_name]
        to_node = self.nodes[to_name]
        
        # 合并节点数据
        to_node.merge_with(from_node)
        
        # 更新所有指向from_node的连接
        for node in self.nodes.values():
            if from_name in node.connections:
                count = node.connections[from_name]
                del node.connections[from_name]
                if to_name not in node.connections:
                    node.connections[to_name] = 0
                node.connections[to_name] += count
        
        # 删除旧节点
        del self.nodes[from_name]
    
    def merge_graph(self, other: 'AnnotationGraph'):
        """
        合并另一个图（用于多基因组叠加）
        
        Args:
            other: 要合并的图
        """
        # 合并基因组信息
        for genome_id in other.genomes:
            if genome_id not in self.genomes:
                self.genomes.append(genome_id)
        
        for genome_id, entries in other.genome_entries.items():
            if genome_id not in self.genome_entries:
                self.genome_entries[genome_id] = []
            self.genome_entries[genome_id].extend(entries)
        
        # 合并节点
        for name, other_node in other.nodes.items():
            if name in self.nodes:
                self.nodes[name].merge_with(other_node)
            else:
                self.nodes[name] = other_node
    
    def get_nodes_sorted_by_count(self) -> List[AnnotationNode]:
        """
        按基因组数量排序返回节点
        
        Returns:
            排序后的节点列表
        """
        return sorted(self.nodes.values(), key=lambda x: x.count, reverse=True)
    
    def get_edges(self) -> List[Tuple[str, str, int]]:
        """
        获取所有边
        
        Returns:
            [(source_name, target_name, count), ...]
        """
        edges = []
        for source_name, node in self.nodes.items():
            for target_name, count in node.connections.items():
                edges.append((source_name, target_name, count))
        return edges
    
    def get_edges_sorted_by_count(self) -> List[Tuple[str, str, int]]:
        """
        按计数排序返回边
        
        Returns:
            排序后的边列表
        """
        return sorted(self.get_edges(), key=lambda x: x[2], reverse=True)
    
    def clear(self):
        """清空图"""
        self.nodes.clear()
        self.genomes.clear()
        self.genome_entries.clear()
    
    def __len__(self) -> int:
        """返回节点数量"""
        return len(self.nodes)
    
    def __str__(self) -> str:
        return f"AnnotationGraph(nodes={len(self.nodes)}, genomes={len(self.genomes)})"
    
    def __repr__(self) -> str:
        return f"AnnotationGraph({len(self.nodes)} nodes, {len(self.genomes)} genomes)"


@dataclass
class GeneStats:
    """
    基因统计信息数据类
    
    用于表格显示，存储基因的邻接关系和taxon数量统计
    """
    gene_name: str          # 基因名称
    left: int               # 左邻接边数量（指向该节点的总边数）
    right: int              # 右邻接边数量（该节点指向其他节点的总边数）
    left_nodes: int         # 左邻接节点数量（唯一的不同节点名）
    right_nodes: int        # 右邻接节点数量（唯一的不同节点名）
    taxon: int              # 包含此节点的Entry数量（taxa数量）
    
    def __str__(self) -> str:
        return f"GeneStats({self.gene_name}, L:{self.left}, R:{self.right}, LNodes:{self.left_nodes}, RNodes:{self.right_nodes}, T:{self.taxon})"
    
    def __repr__(self) -> str:
        return f"GeneStats('{self.gene_name}', left={self.left}, right={self.right}, left_nodes={self.left_nodes}, right_nodes={self.right_nodes}, taxon={self.taxon})"