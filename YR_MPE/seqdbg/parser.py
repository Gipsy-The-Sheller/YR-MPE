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
                # 生成条目ID
                entry_id = record.id if record.id else f"entry_{record_index}"
                if not entry_id.startswith(f"{Path(file_path).stem}_"):
                    entry_id = f"{Path(file_path).stem}_{entry_id}"
                
                genes = []
                
                for feature in record.features:
                    if feature.type == "gene":
                        gene_info = self._extract_gene_info_biopython(feature, record)
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
        record: SeqRecord
    ) -> Optional[GeneInfo]:
        """
        从Biopython特征中提取基因信息
        
        Args:
            feature: Biopython特征对象
            record: SeqRecord对象
            
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
            
            return GeneInfo(
                name=gene_name,
                start=start,
                end=end,
                strand=strand
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
        
        Args:
            file_path: 文件路径
            entry_id: 条目ID
            
        Returns:
            SeqRecord对象或None
        """
        try:
            if not self.use_biopython:
                logger.warning("Biopython不可用，无法提取序列")
                return None
            
            for record in SeqIO.parse(file_path, "genbank"):
                record_id = f"{Path(file_path).stem}_{record.id}"
                if record_id == entry_id:
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
            
            # 创建SeqRecord
            seq_record = SeqRecord(
                Seq(seq_str),
                id=f"{gene_info.name}_{entry_id}",
                description=f"{gene_info.name} from {entry_id}"
            )
            
            return seq_record
            
        except Exception as e:
            logger.error(f"提取序列失败: {e}")
            return None
    
    def clear_cache(self):
        """清空缓存"""
        self.records_cache.clear()
        logger.info("缓存已清空")
    
    def get_cached_files(self) -> List[str]:
        """
        获取缓存的文件列表
        
        Returns:
            文件路径列表
        """
        return list(self.records_cache.keys())