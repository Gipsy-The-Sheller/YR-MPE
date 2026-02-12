#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GeneBank文件解析器

解析GeneBank格式的基因组文件，提取基因信息和序列。
支持多条目文件，并提供完整的错误处理。
"""

import re
import os
import logging
from typing import Dict, List, Optional, Tuple
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from .models import GeneInfo

# 配置日志
logger = logging.getLogger(__name__)


class GeneBankParseError(Exception):
    """GeneBank解析错误"""
    pass


class GeneBankParser:
    """
    GeneBank文件解析器
    
    支持Biopython和正则表达式两种解析方式
    """
    
    def __init__(self, use_biopython: bool = True):
        """
        初始化解析器
        
        Args:
            use_biopython: 优先使用Biopython解析，否则使用正则表达式
        """
        self.use_biopython = use_biopython and self._check_biopython()
        self.records_cache: Dict[str, Dict] = {}  # 缓存解析结果
        self.seq_record_cache: Dict[Tuple[str, str], SeqRecord] = {}  # 缓存SeqRecord对象: {(file_path, entry_id): SeqRecord}
        
    @staticmethod
    def _check_biopython() -> bool:
        """检查Biopython是否可用"""
        try:
            import Bio
            return True
        except ImportError:
            return False
    
    def parse_file(self, file_path: str) -> List[Tuple[str, List[GeneInfo]]]:
        """
        解析GeneBank文件
        
        Args:
            file_path: GeneBank文件路径
            
        Returns:
            [(entry_id, [GeneInfo, ...]), ...] 条目列表
            
        Raises:
            GeneBankParseError: 解析失败时抛出
        """
        if not os.path.exists(file_path):
            raise GeneBankParseError(f"文件不存在: {file_path}")
        
        # 检查缓存
        if file_path in self.records_cache:
            logger.info(f"从缓存加载: {file_path}")
            return self.records_cache[file_path]
        
        entries = []
        
        try:
            if self.use_biopython:
                entries = self._parse_with_biopython(file_path)
            else:
                entries = self._parse_with_regex(file_path)
            
            # 缓存结果
            self.records_cache[file_path] = entries
            logger.info(f"成功解析 {file_path}: {len(entries)} 个条目")
            
        except Exception as e:
            logger.error(f"解析失败 {file_path}: {e}")
            raise GeneBankParseError(f"解析文件失败: {e}") from e
        
        return entries
    
    def _parse_with_biopython(self, file_path: str) -> List[Tuple[str, List[GeneInfo]]]:
        """
        使用Biopython解析GeneBank文件
        
        Args:
            file_path: 文件路径
            
        Returns:
            条目列表
        """
        entries = []
        
        try:
            for record_index, record in enumerate(SeqIO.parse(file_path, "genbank")):
                # 生成条目ID - 直接使用LOCUS字段(record.id)
                entry_id = record.id if record.id else f"entry_{record_index}"
                
                # 首先遍历所有特征，建立位置到类型的映射
                # GeneBank中CDS/rRNA/tRNA特征通常紧跟着gene特征，但基因名称可能不一致
                # 使用位置作为主键，避免基因名称不一致导致的匹配失败
                gene_feature_map: Dict[str, List[Tuple[str, int, int]]] = {}  # {gene_name: [(feature_type, start, end), ...]}
                feature_by_position: List[Tuple[int, int, str, str]] = []  # [(start, end, feature_type, gene_name)]
                
                for feature in record.features:
                    if feature.type not in ("CDS", "rRNA", "tRNA", "gene"):
                        continue
                    
                    # 获取基因名称
                    gene_name = None
                    if "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                    elif "locus_tag" in feature.qualifiers:
                        gene_name = feature.qualifiers["locus_tag"][0]
                    else:
                        gene_name = f"gene_{feature.location.start}_{feature.location.end}"
                    
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    
                    # 添加到位置列表
                    feature_by_position.append((start, end, feature.type, gene_name))
                    
                    # 同时维护按基因名称的映射（用于调试）
                    if gene_name not in gene_feature_map:
                        gene_feature_map[gene_name] = []
                    gene_feature_map[gene_name].append((feature.type, start, end))
                
                # 按起始位置排序，方便后续匹配
                feature_by_position.sort(key=lambda x: x[0])
                
                genes = []
                
                for feature in record.features:
                    if feature.type == "gene":
                        gene_info = self._extract_gene_info_biopython(
                            feature, 
                            record,
                            gene_feature_map=gene_feature_map,
                            feature_by_position=feature_by_position
                        )
                        if gene_info:
                            genes.append(gene_info)
                
                # 按起始位置排序
                genes.sort(key=lambda x: x.start)
                entries.append((entry_id, genes))
                
        except Exception as e:
            logger.error(f"Biopython解析失败: {e}")
            raise
        
        return entries
    
    def _extract_gene_info_biopython(
        self, 
        feature, 
        record: SeqRecord,
        gene_feature_map: Optional[Dict] = None,
        feature_by_position: Optional[List] = None
    ) -> Optional[GeneInfo]:
        """
        从Biopython特征中提取基因信息
        
        Args:
            feature: Biopython特征对象
            record: SeqRecord对象
            gene_feature_map: 基因名称到特征类型的映射（可选）
            feature_by_position: 按位置排序的特征列表（可选）
            
        Returns:
            GeneInfo对象或None
        """
        try:
            # 获取基因名称
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
            elif "locus_tag" in feature.qualifiers:
                gene_name = feature.qualifiers["locus_tag"][0]
            else:
                gene_name = f"gene_{feature.location.start}_{feature.location.end}"
            
            # 获取位置
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand if feature.location.strand else 1
            
            # 获取基因类型
            gene_type = "gene"  # 默认类型
            if gene_feature_map and feature_by_position:
                # 直接根据位置从 feature_by_position 中查找匹配的特征类型
                # 优先精确坐标匹配（距离=0），因为在线粒体基因组中gene和CDS坐标通常完全一致
                exact_match_type = None
                exact_match_name = None
                best_match_type = None
                min_distance = float('inf')
                matched_gene_name = None
                
                for f_start, f_end, f_type, f_name in feature_by_position:
                    # 跳过当前的 gene feature 自身
                    if f_type == "gene" and f_start == start and f_end == end:
                        continue
                    
                    # 计算位置差异（考虑起点和终点的差异）
                    distance = abs(f_start - start) + abs(f_end - end)
                    
                    # 优先检查精确匹配（距离=0）
                    if distance == 0 and f_type in ("CDS", "rRNA", "tRNA"):
                        exact_match_type = f_type
                        exact_match_name = f_name
                        break  # 找到精确匹配，立即停止
                    
                    # 如果没有精确匹配，记录最佳匹配（20bp容差）
                    if distance <= 20 and distance < min_distance and f_type in ("CDS", "rRNA", "tRNA"):
                        min_distance = distance
                        best_match_type = f_type
                        matched_gene_name = f_name
                
                # 优先使用精确匹配
                if exact_match_type:
                    gene_type = exact_match_type
                    if exact_match_name and exact_match_name != gene_name:
                        logger.debug(f"基因 {gene_name} ({start}-{end}) 精确匹配到 {exact_match_name} 的类型: {gene_type}")
                    else:
                        logger.debug(f"基因 {gene_name} ({start}-{end}) 精确匹配到类型: {gene_type}")
                # 如果没有精确匹配，使用最佳近似匹配
                elif best_match_type:
                    gene_type = best_match_type
                    if matched_gene_name and matched_gene_name != gene_name:
                        logger.debug(f"基因 {gene_name} ({start}-{end}) 近似匹配到 {matched_gene_name} 的类型: {gene_type} (距离: {min_distance})")
                    else:
                        logger.debug(f"基因 {gene_name} ({start}-{end}) 近似匹配到类型: {gene_type} (距离: {min_distance})")
                else:
                    logger.debug(f"基因 {gene_name} ({start}-{end}) 使用默认类型: {gene_type}")
            
            return GeneInfo(
                name=gene_name,
                start=start,
                end=end,
                strand=strand,
                gene_type=gene_type
            )
            
        except Exception as e:
            logger.warning(f"提取基因信息失败: {e}")
            return None
    
    def _parse_with_regex(self, file_path: str) -> List[Tuple[str, List[GeneInfo]]]:
        """
        使用正则表达式解析GeneBank文件（回退方案）
        
        Args:
            file_path: 文件路径
            
        Returns:
            条目列表
        """
        entries = []
        
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # 分割条目
            entry_pattern = r'(LOCUS\s+.+?)(?=LOCUS\s+|$)'
            entry_matches = re.finditer(entry_pattern, content, re.DOTALL)
            
            for entry_index, entry_match in enumerate(entry_matches):
                entry_content = entry_match.group(1)
                
                # 提取条目ID
                locus_match = re.search(r'LOCUS\s+(\S+)', entry_content)
                entry_id = locus_match.group(1) if locus_match else f"entry_{entry_index}"
                entry_id = f"{Path(file_path).stem}_{entry_id}"
                
                # 提取基因
                genes = self._extract_genes_regex(entry_content)
                
                # 按起始位置排序
                genes.sort(key=lambda x: x.start)
                entries.append((entry_id, genes))
                
        except Exception as e:
            logger.error(f"正则表达式解析失败: {e}")
            raise
        
        return entries
    
    def _extract_genes_regex(self, entry_content: str) -> List[GeneInfo]:
        """
        从条目内容中提取基因（正则表达式方法）
        
        Args:
            entry_content: 条目内容
            
        Returns:
            基因信息列表
        """
        genes = []
        
        # 匹配基因特征
        gene_pattern = r'gene\s+(complement\()?(\d+)\.\.(\d+)\)?\s*.*?/gene="([^"]*)"'
        matches = re.finditer(gene_pattern, entry_content, re.DOTALL)
        
        for match in matches:
            try:
                complement = match.group(1) is not None
                start = int(match.group(2))
                end = int(match.group(3))
                gene_name = match.group(4)
                strand = -1 if complement else 1
                
                genes.append(GeneInfo(
                    name=gene_name,
                    start=start,
                    end=end,
                    strand=strand
                ))
            except (ValueError, IndexError) as e:
                logger.warning(f"解析基因特征失败: {e}")
                continue
        
        return genes
    
    def get_full_record(self, file_path: str, entry_id: str) -> Optional[SeqRecord]:
        """
        获取完整的SeqRecord对象（包含序列）
        
        使用缓存避免重复读取文件，提升性能。
        
        Args:
            file_path: 文件路径
            entry_id: 条目ID
            
        Returns:
            SeqRecord对象或None
        """
        # 检查缓存
        cache_key = (file_path, entry_id)
        if cache_key in self.seq_record_cache:
            logger.debug(f"从缓存获取记录: {entry_id}")
            return self.seq_record_cache[cache_key]
        
        try:
            if not self.use_biopython:
                logger.warning("Biopython不可用，无法提取序列")
                return None
            
            for record in SeqIO.parse(file_path, "genbank"):
                # entry_id是LOCUS字段，对应record.name
                # record.id是DEFINITION字段（通常是accession）
                # 所以使用record.name来匹配entry_id
                if record.name == entry_id or record.id == entry_id:
                    # 缓存记录以避免重复读取
                    self.seq_record_cache[cache_key] = record
                    return record
                    
        except Exception as e:
            logger.error(f"获取记录失败: {e}")
        
        return None
    
    def extract_sequence(
        self, 
        file_path: str, 
        entry_id: str, 
        gene_info: GeneInfo
    ) -> Optional[SeqRecord]:
        """
        提取基因序列
        
        Args:
            file_path: 文件路径
            entry_id: 条目ID
            gene_info: 基因信息
            
        Returns:
            SeqRecord对象或None
        """
        record = self.get_full_record(file_path, entry_id)
        if not record:
            return None
        
        try:
            # 提取序列片段
            seq_str = str(record.seq[gene_info.start:gene_info.end])
            
            # 如果是负链，进行反向互补
            if gene_info.strand == -1:
                seq_obj = Seq(seq_str)
                seq_str = str(seq_obj.reverse_complement())
            
            # 创建SeqRecord - 使用entry_id作为序列ID，方便按entry合并
            seq_record = SeqRecord(
                Seq(seq_str),
                id=entry_id,
                description=f"{gene_info.name} from {entry_id}"
            )
            
            return seq_record
            
        except Exception as e:
            logger.error(f"提取序列失败: {e}")
            return None
    
    def extract_protein_sequence(
        self, 
        file_path: str, 
        entry_id: str, 
        gene_info: GeneInfo
    ) -> Optional[SeqRecord]:
        """
        提取蛋白质序列（从 CDS 特征的 /translation 字段）
        
        Args:
            file_path: 文件路径
            entry_id: 条目ID
            gene_info: 基因信息
            
        Returns:
            SeqRecord对象或None
        """
        record = self.get_full_record(file_path, entry_id)
        if not record:
            return None
        
        try:
            # 查找对应的 CDS 特征
            cds_feature = None
            
            # 首先尝试精确位置匹配
            for feature in record.features:
                if feature.type == "CDS":
                    f_start = int(feature.location.start)
                    f_end = int(feature.location.end)
                    
                    # 检查位置是否精确匹配（考虑可能的1bp差异，因为有些特征使用<或>符号）
                    if abs(f_start - gene_info.start) <= 1 and abs(f_end - gene_info.end) <= 1:
                        cds_feature = feature
                        break
            
            # 如果没有找到精确位置匹配，尝试基因名称匹配（不区分大小写）
            if not cds_feature:
                for feature in record.features:
                    if feature.type == "CDS":
                        # 获取CDS的基因名称
                        cds_gene_name = None
                        if "gene" in feature.qualifiers:
                            cds_gene_name = feature.qualifiers["gene"][0]
                        
                        # 匹配基因名称（不区分大小写）
                        if cds_gene_name and cds_gene_name.lower() == gene_info.name.lower():
                            cds_feature = feature
                            break
            
            # 如果仍然没有找到，尝试更宽松的位置匹配（20bp容差）
            if not cds_feature:
                for feature in record.features:
                    if feature.type == "CDS":
                        f_start = int(feature.location.start)
                        f_end = int(feature.location.end)
                        
                        if abs(f_start - gene_info.start) <= 20 and abs(f_end - gene_info.end) <= 20:
                            cds_feature = feature
                            break
            
            if not cds_feature:
                logger.warning(f"未找到基因 {gene_info.name} 的 CDS 特征")
                return None
            
            # 从 /translation 字段获取蛋白质序列
            if "translation" in cds_feature.qualifiers:
                protein_seq = cds_feature.qualifiers["translation"][0]
            else:
                logger.warning(f"CDS 特征没有 /translation 字段: {gene_info.name}")
                return None
            
            # 创建SeqRecord - 添加 _protein 后缀以区分蛋白质序列
            seq_record = SeqRecord(
                Seq(protein_seq),
                id=f"{entry_id}_protein",
                description=f"Protein sequence of {gene_info.name} from {entry_id}"
            )
            
            return seq_record
            
        except Exception as e:
            logger.error(f"提取蛋白质序列失败: {e}")
            return None
    
    def clear_cache(self):
        """清空缓存"""
        self.records_cache.clear()
        self.seq_record_cache.clear()
        logger.info("缓存已清空")
    
    def get_cached_files(self) -> List[str]:
        """
        获取缓存的文件列表
        
        Returns:
            文件路径列表
        """
        # 从 seq_record_cache 中提取唯一的文件路径
        files = set(file_path for file_path, _ in self.seq_record_cache.keys())
        # 合并 records_cache 中的文件路径
        files.update(self.records_cache.keys())
        return list(files)