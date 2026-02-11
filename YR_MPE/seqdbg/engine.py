#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SeqDBG图构建引擎

负责从解析的GeneBank数据构建注释邻接图，
支持多基因组叠加和注释标准化。
"""

import logging
from typing import List, Dict, Optional, Tuple
from pathlib import Path

from .models import AnnotationGraph, GeneInfo, AnnotationNode, GeneStats
from .parser import GeneBankParser
from .normalizer import AnnotationNormalizer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class SeqDBGraphBuilder:
    """
    SeqDBG图构建器
    
    主要功能：
    1. 解析GeneBank文件
    2. 应用注释标准化规则
    3. 构建注释邻接图
    4. 支持多基因组叠加
    5. 提取基因序列
    """
    
    def __init__(
        self, 
        parser: Optional[GeneBankParser] = None,
        normalizer: Optional[AnnotationNormalizer] = None
    ):
        """
        初始化图构建器
        
        Args:
            parser: GeneBank解析器（可选，默认创建新实例）
            normalizer: 注释标准化器（可选，默认创建新实例）
        """
        self.parser = parser or GeneBankParser()
        self.normalizer = normalizer or AnnotationNormalizer()
        self.graph = AnnotationGraph()
        self.file_cache: Dict[str, List[Tuple[str, List[GeneInfo]]]] = {}
        
        # 统计信息
        self.stats = {
            "total_files": 0,
            "total_entries": 0,
            "total_genes": 0,
            "normalized_genes": 0,
            "unique_annotations": 0
        }
    
    def load_file(
        self, 
        file_path: str,
        extract_sequences: bool = True
    ) -> bool:
        """
        加载单个GeneBank文件并构建图
        
        Args:
            file_path: GeneBank文件路径
            extract_sequences: 是否提取序列
            
        Returns:
            是否成功
        """
        try:
            # 解析文件
            entries = self.parser.parse_file(file_path)
            self.file_cache[file_path] = entries
            
            # 添加基因组
            genome_id = Path(file_path).stem
            self.graph.add_genome(genome_id)
            
            # 处理每个条目
            for entry_id, genes in entries:
                self.graph.add_entry(genome_id, entry_id)
                self._process_entry(genome_id, entry_id, genes, file_path, extract_sequences)
            
            self.stats["total_files"] += 1
            self.stats["total_entries"] += len(entries)
            self.stats["total_genes"] += sum(len(genes) for _, genes in entries)
            
            logger.info(f"成功加载 {file_path}: {len(entries)} 个条目, "
                       f"{sum(len(genes) for _, genes in entries)} 个基因")
            
            return True
            
        except Exception as e:
            logger.error(f"加载文件失败 {file_path}: {e}")
            return False
    
    def load_files(
        self, 
        file_paths: List[str],
        extract_sequences: bool = True,
        progress_callback=None
    ) -> Dict[str, bool]:
        """
        批量加载GeneBank文件
        
        Args:
            file_paths: GeneBank文件路径列表
            extract_sequences: 是否提取序列
            progress_callback: 可选，进度回调函数 (message, current, total)
            
        Returns:
            {file_path: success} 字典
        """
        results = {}
        total_files = len(file_paths)
        
        for i, file_path in enumerate(file_paths):
            # 调用进度回调
            if progress_callback:
                progress_callback(f"正在解析: {Path(file_path).name}", i + 1, total_files)
            
            results[file_path] = self.load_file(file_path, extract_sequences)
        
        return results
    
    def _process_entry(
        self,
        genome_id: str,
        entry_id: str,
        genes: List[GeneInfo],
        file_path: str,
        extract_sequences: bool
    ):
        """
        处理单个条目
        
        Args:
            genome_id: 基因组ID
            entry_id: 条目ID
            genes: 基因信息列表
            file_path: 文件路径
            extract_sequences: 是否提取序列
        """
        # 标准化基因名称
        normalized_genes = []
        for gene in genes:
            original_name = gene.name
            canonical_name = self.normalizer.normalize(gene.name)
            
            # 设置entry_id
            gene.entry_id = entry_id
            
            if original_name != canonical_name:
                gene.normalized_name = canonical_name
                self.stats["normalized_genes"] += 1
                logger.debug(f"标准化: {original_name} -> {canonical_name}")
            
            normalized_genes.append(gene)
        
        # 应用注释映射（合并节点）
        self._apply_normalization_mapping()
        
        # 构建图的节点和边
        for i, gene in enumerate(normalized_genes):
            # 添加节点
            self.graph.add_node_occurrence(genome_id, gene)
            
            # 提取序列
            if extract_sequences:
                # 提取 DNA 序列
                seq_record = self.parser.extract_sequence(file_path, entry_id, gene)
                if seq_record:
                    node = self.graph.get_or_create_node(gene.normalized_name)
                    node.add_sequence(genome_id, seq_record)
                
                # 如果是 CDS 类型，也提取蛋白质序列
                if gene.gene_type == "CDS":
                    protein_record = self.parser.extract_protein_sequence(file_path, entry_id, gene)
                    if protein_record:
                        node = self.graph.get_or_create_node(gene.normalized_name)
                        # 蛋白质序列ID使用 entry_id_protein 后缀，可通过后缀区分
                        node.add_sequence(genome_id, protein_record)
        
        # 构建邻接关系（基于基因组坐标顺序）
        # 按基因组坐标排序后的相邻基因即为邻接关系
        # 邻接关系是有向的：前面基因 -> 后面基因
        for i in range(len(normalized_genes) - 1):
            gene1 = normalized_genes[i]
            gene2 = normalized_genes[i + 1]
            
            # 添加有向边：gene1 -> gene2
            self.graph.add_edge(
                gene1.normalized_name,
                gene2.normalized_name
            )
            logger.debug(f"邻接: {gene1.normalized_name}({gene1.start}-{gene1.end}) -> {gene2.normalized_name}({gene2.start}-{gene2.end})")
    
    def _apply_normalization_mapping(self):
        """
        应用标准化映射，合并节点
        """
        mapping = self.normalizer.get_mapping()
        
        for original, canonical in mapping.items():
            if original != canonical:
                # 检查是否已经合并过
                if original in self.graph.nodes:
                    self.graph.merge_nodes(original, canonical)
                    logger.debug(f"合并节点: {original} -> {canonical}")
    
    def build_multi_genome_graph(
        self, 
        file_paths: List[str],
        extract_sequences: bool = True,
        progress_callback=None
    ) -> AnnotationGraph:
        """
        构建多基因组邻接图
        
        Args:
            file_paths: GeneBank文件路径列表
            extract_sequences: 是否提取序列
            progress_callback: 可选，进度回调函数 (message, current, total)
            
        Returns:
            构建好的图
        """
        # 清空现有图
        self.graph = AnnotationGraph()
        self.stats = {
            "total_files": 0,
            "total_entries": 0,
            "total_genes": 0,
            "normalized_genes": 0,
            "unique_annotations": 0
        }
        
        # 加载所有文件
        self.load_files(file_paths, extract_sequences, progress_callback)
        
        # 统计唯一注释
        self.stats["unique_annotations"] = len(self.graph.nodes)
        
        logger.info(f"多基因组图构建完成: {len(self.graph.nodes)} 个唯一注释")
        
        return self.graph
    
    def export_to_dataset_items(
        self, 
        selected_nodes: List[str],
        dataset_id: str = ""
    ) -> List:
        """
        将选中的节点导出为DatasetItem
        
        Args:
            selected_nodes: 选中的节点名称列表
            dataset_id: 数据集ID
            
        Returns:
            DatasetItem列表
        """
        from ..platforms.methods.dataset_models import DatasetItem, ITEM_TYPE_SEQUENCE
        
        items = []
        
        for node_name in selected_nodes:
            if node_name not in self.graph.nodes:
                logger.warning(f"节点不存在: {node_name}")
                continue
            
            node = self.graph.nodes[node_name]
            all_sequences = node.get_all_sequences()
            
            if not all_sequences:
                logger.warning(f"节点 {node_name} 没有序列数据")
                continue
            
            # 创建DatasetItem
            item = DatasetItem(item_type=ITEM_TYPE_SEQUENCE)
            item.dataset_id = dataset_id
            item.loci_name = node_name
            item.sequences = all_sequences
            item.sequence_count = len(all_sequences)
            
            if all_sequences:
                item.length = len(all_sequences[0].seq)
            
            item.data = {
                "annotation": node_name,
                "genome_count": node.count,
                "total_occurrences": node.occurrences,
                "genomes": list(node.get_unique_genomes())
            }
            
            items.append(item)
            logger.info(f"导出节点: {node_name}, {len(all_sequences)} 条序列")
        
        return items
    
    def export_to_fasta(
        self, 
        selected_nodes: List[str],
        output_dir: str,
        combine: bool = False
    ) -> Dict[str, str]:
        """
        将选中的节点导出为FASTA文件
        
        Args:
            selected_nodes: 选中的节点名称列表
            output_dir: 输出目录
            combine: 是否合并为一个文件
            
        Returns:
            {node_name: file_path} 字典
        """
        from Bio import SeqIO
        from pathlib import Path
        
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        results = {}
        
        if combine:
            # 合并为一个文件
            output_path = Path(output_dir) / "selected_nodes.fasta"
            all_sequences = []
            
            for node_name in selected_nodes:
                if node_name not in self.graph.nodes:
                    continue
                
                node = self.graph.nodes[node_name]
                sequences = node.get_all_sequences()
                
                for seq in sequences:
                    # 修改序列ID以包含节点名称
                    seq.id = f"{node_name}_{seq.id}"
                    seq.description = f"{node_name} - {seq.description}"
                    all_sequences.append(seq)
            
            if all_sequences:
                SeqIO.write(all_sequences, output_path, "fasta")
                results["combined"] = str(output_path)
                logger.info(f"导出合并FASTA: {output_path}")
        
        else:
            # 每个节点一个文件
            for node_name in selected_nodes:
                if node_name not in self.graph.nodes:
                    continue
                
                node = self.graph.nodes[node_name]
                sequences = node.get_all_sequences()
                
                if sequences:
                    output_path = Path(output_dir) / f"{node_name}.fasta"
                    SeqIO.write(sequences, output_path, "fasta")
                    results[node_name] = str(output_path)
                    logger.info(f"导出节点 {node_name}: {output_path}")
        
        return results
    
    def get_stats(self) -> Dict:
        """
        获取统计信息
        
        Returns:
            统计信息字典
        """
        return self.stats.copy()
    
    def clear_cache(self):
        """清空缓存"""
        self.file_cache.clear()
        self.parser.clear_cache()
        self.normalizer.mapping_history.clear()
        logger.info("缓存已清空")
    
    def reset(self):
        """重置图构建器"""
        self.graph = AnnotationGraph()
        self.file_cache.clear()
        self.stats = {
            "total_files": 0,
            "total_entries": 0,
            "total_genes": 0,
            "normalized_genes": 0,
            "unique_annotations": 0
        }
        logger.info("图构建器已重置")
    
    def get_gene_stats(self, filter_entries: Optional[List[str]] = None) -> List[GeneStats]:
        """
        获取基因统计信息
        
        计算每个基因的左邻接边数量、右邻接边数量、左邻接节点数量、右邻接节点数量和taxon数量
        
        Args:
            filter_entries: 可选，只统计这些entries中的基因。如果为None，统计所有entries。
        
        Returns:
            GeneStats列表
        """
        stats_list = []
        
        if not filter_entries:
            # 没有过滤，使用原始图的统计数据
            left_counts: Dict[str, int] = {}      # {node_name: 左侧边数}
            right_counts: Dict[str, int] = {}     # {node_name: 右侧边数}
            left_node_counts: Dict[str, int] = {}  # {node_name: 左侧唯一节点数}
            right_node_counts: Dict[str, int] = {} # {node_name: 右侧唯一节点数}
            
            for source_name, node in self.graph.nodes.items():
                right_counts[source_name] = sum(node.connections.values())
                right_node_counts[source_name] = len(node.connections)
                
            for source_name, node in self.graph.nodes.items():
                for target_name, count in node.connections.items():
                    if target_name not in left_counts:
                        left_counts[target_name] = 0
                    left_counts[target_name] += count
                    
                    if target_name not in left_node_counts:
                        left_node_counts[target_name] = 0
                    left_node_counts[target_name] += 1
            
            taxon_counts: Dict[str, int] = {}
            gene_type_map: Dict[str, Dict[str, int]] = {}  # {node_name: {gene_type: count}}
            for node_name, node in self.graph.nodes.items():
                taxon_counts[node_name] = node.occurrences
                
                # 收集基因类型信息
                gene_type_map[node_name] = {}
                for genome_id, gene_infos in node.genomes.items():
                    for gene_info in gene_infos:
                        gene_type = gene_info.gene_type or "gene"
                        if gene_type not in gene_type_map[node_name]:
                            gene_type_map[node_name][gene_type] = 0
                        gene_type_map[node_name][gene_type] += 1
            
            for node_name in self.graph.nodes.keys():
                left = left_counts.get(node_name, 0)
                right = right_counts.get(node_name, 0)
                left_nodes = left_node_counts.get(node_name, 0)
                right_nodes = right_node_counts.get(node_name, 0)
                taxon = taxon_counts.get(node_name, 0)
                
                # 确定主要基因类型（出现次数最多）
                gene_type = "gene"
                if node_name in gene_type_map and gene_type_map[node_name]:
                    gene_type = max(gene_type_map[node_name].items(), key=lambda x: x[1])[0]
                
                stats = GeneStats(
                    gene_name=node_name,
                    left=left,
                    right=right,
                    left_nodes=left_nodes,
                    right_nodes=right_nodes,
                    taxon=taxon,
                    gene_type=gene_type
                )
                stats_list.append(stats)
        else:
            # 基于entries过滤，重新计算所有统计数据
            # 收集过滤后的所有entry中的基因信息
            filter_entries_set = set(filter_entries)
            gene_in_entries: Dict[str, List[Tuple[str, int]]] = {}  # {gene_name: [(entry_id, index), ...]}
            entry_genes: Dict[str, List[Tuple[str, int]]] = {}  # {entry_id: [(gene_name, index), ...]}
            
            # 从所有节点中收集在过滤entries中的基因
            gene_type_map: Dict[str, Dict[str, int]] = {}  # {node_name: {gene_type: count}}
            entry_genes_full: Dict[str, List[GeneInfo]] = {}  # {entry_id: [GeneInfo, ...]} 存储完整基因信息
            for node_name, node in self.graph.nodes.items():
                for genome_id, gene_infos in node.genomes.items():
                    for gene_info in gene_infos:
                        if gene_info.entry_id in filter_entries_set:
                            if node_name not in gene_in_entries:
                                gene_in_entries[node_name] = []
                            gene_in_entries[node_name].append((gene_info.entry_id, gene_info.start))
                            
                            if gene_info.entry_id not in entry_genes_full:
                                entry_genes_full[gene_info.entry_id] = []
                            entry_genes_full[gene_info.entry_id].append(gene_info)
                            
                            # 收集基因类型信息
                            gene_type = gene_info.gene_type or "gene"
                            if node_name not in gene_type_map:
                                gene_type_map[node_name] = {}
                            if gene_type not in gene_type_map[node_name]:
                                gene_type_map[node_name][gene_type] = 0
                            gene_type_map[node_name][gene_type] += 1
            
            # 重新计算统计（基于基因组坐标顺序）
            taxon_counts: Dict[str, int] = {}
            left_counts: Dict[str, int] = {}
            right_counts: Dict[str, int] = {}
            left_node_counts: Dict[str, int] = {}
            right_node_counts: Dict[str, int] = {}
            
            for gene_name, occurrences in gene_in_entries.items():
                taxon_counts[gene_name] = len(occurrences)
            
            # 对于每个entry中的基因，按坐标顺序构建邻接关系
            for entry_id, genes in entry_genes_full.items():
                # 按起始位置排序
                genes.sort(key=lambda x: x.start)
                
                # 构建邻接关系（基于基因组坐标顺序）
                for i in range(len(genes) - 1):
                    gene1 = genes[i]
                    gene2 = genes[i + 1]
                    
                    # 右邻接：gene1 -> gene2
                    if gene1.normalized_name not in right_counts:
                        right_counts[gene1.normalized_name] = 0
                    right_counts[gene1.normalized_name] += 1
                    
                    if gene1.normalized_name not in right_node_counts:
                        right_node_counts[gene1.normalized_name] = set()
                    right_node_counts[gene1.normalized_name].add(gene2.normalized_name)
                    
                    # 左邻接：gene2 -> gene1
                    if gene2.normalized_name not in left_counts:
                        left_counts[gene2.normalized_name] = 0
                    left_counts[gene2.normalized_name] += 1
                    
                    if gene2.normalized_name not in left_node_counts:
                        left_node_counts[gene2.normalized_name] = set()
                    left_node_counts[gene2.normalized_name].add(gene1.normalized_name)
            
            # 转换集合为计数
            left_node_counts_final = {k: len(v) if isinstance(v, set) else v for k, v in left_node_counts.items()}
            right_node_counts_final = {k: len(v) if isinstance(v, set) else v for k, v in right_node_counts.items()}
            
            # 生成统计列表
            for gene_name in gene_in_entries.keys():
                left = left_counts.get(gene_name, 0)
                right = right_counts.get(gene_name, 0)
                left_nodes = left_node_counts_final.get(gene_name, 0)
                right_nodes = right_node_counts_final.get(gene_name, 0)
                taxon = taxon_counts.get(gene_name, 0)
                
                # 确定主要基因类型（出现次数最多）
                gene_type = "gene"
                if gene_name in gene_type_map and gene_type_map[gene_name]:
                    gene_type = max(gene_type_map[gene_name].items(), key=lambda x: x[1])[0]
                
                stats = GeneStats(
                    gene_name=gene_name,
                    left=left,
                    right=right,
                    left_nodes=left_nodes,
                    right_nodes=right_nodes,
                    taxon=taxon,
                    gene_type=gene_type
                )
                stats_list.append(stats)
        
        # 按taxon数量排序（降序）
        stats_list.sort(key=lambda x: x.taxon, reverse=True)
        
        logger.info(f"生成统计信息: {len(stats_list)} 个基因")
        return stats_list
    
    def get_adjacency_info(self, gene_name: str) -> Dict:
        """
        获取指定基因的邻接信息
        
        Args:
            gene_name: 基因名称
            
        Returns:
            包含左邻接节点、右邻接节点和taxon信息的字典
        """
        if gene_name not in self.graph.nodes:
            return {
                "gene_name": gene_name,
                "left_adjacent": [],
                "right_adjacent": [],
                "taxon_info": {},
                "error": "Gene not found"
            }
        
        node = self.graph.nodes[gene_name]
        
        # 左邻接节点（指向该节点的节点）
        left_adjacent: Dict[str, int] = {}
        for source_name, source_node in self.graph.nodes.items():
            if gene_name in source_node.connections:
                left_adjacent[source_name] = source_node.connections[gene_name]
        
        # 右邻接节点（该节点指向的节点）
        right_adjacent: Dict[str, int] = node.connections.copy()
        
        # Taxon信息
        taxon_info: Dict[str, int] = {}
        for genome_id, gene_infos in node.genomes.items():
            taxon_info[genome_id] = len(gene_infos)
        
        return {
            "gene_name": gene_name,
            "left_adjacent": sorted(left_adjacent.items(), key=lambda x: -x[1]),
            "right_adjacent": sorted(right_adjacent.items(), key=lambda x: -x[1]),
            "taxon_info": sorted(taxon_info.items(), key=lambda x: -x[1])
        }
    
    def get_gene_info_and_sequences(self, use_protein: bool = False) -> Tuple[Dict[str, Dict[str, GeneInfo]], Dict[str, Dict[str, List[SeqRecord]]], Dict[str, str]]:
        """
        获取所有基因的信息和序列数据
        
        Args:
            use_protein: 已废弃，保留参数以兼容旧调用。现在总是返回所有序列。
        
        Returns:
            (gene_info_data, sequence_data, gene_type_map)
            - gene_info_data: {gene_name: {genome_id: GeneInfo}}
            - sequence_data: {gene_name: {genome_id: [SeqRecord]}}  # 返回所有序列列表（DNA和蛋白质）
            - gene_type_map: {gene_name: gene_type}
        """
        gene_info_data = {}
        sequence_data = {}
        gene_type_map = {}
        
        for node_name, node in self.graph.nodes.items():
            # 收集基因信息
            gene_info_data[node_name] = {}
            for genome_id, gene_infos in node.genomes.items():
                for gene_info in gene_infos:
                    # 只保存第一个基因信息（避免重复）
                    if genome_id not in gene_info_data[node_name]:
                        gene_info_data[node_name][genome_id] = gene_info
            
            # 收集所有序列数据（包括DNA和蛋白质）
            # 序列ID中的 _protein 后缀可用于区分序列类型
            sequence_data[node_name] = {}
            for genome_id, seq_records in node.sequences.items():
                # 直接返回所有序列，让调用者根据需要过滤
                sequence_data[node_name][genome_id] = seq_records.copy()
            
            # 收集基因类型
            if node_name not in gene_type_map:
                gene_type_map[node_name] = "gene"
            
            for genome_id, gene_infos in node.genomes.items():
                for gene_info in gene_infos:
                    if gene_info.gene_type and gene_info.gene_type != "gene":
                        gene_type_map[node_name] = gene_info.gene_type
                        break
        
        logger.info(f"收集了 {len(gene_info_data)} 个基因的信息和所有序列")
        return gene_info_data, sequence_data, gene_type_map