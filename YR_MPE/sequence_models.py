#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Sequence Models for YR-MPE Sequence Editor
Data models and structures for sequence management
"""

import re
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum
from collections import OrderedDict

from PyQt5.QtCore import *
from PyQt5.QtGui import *


class SequenceType(Enum):
    """序列类型"""
    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "PROTEIN"
    UNKNOWN = "UNKNOWN"


class SequenceFormat(Enum):
    """序列格式"""
    FASTA = "fasta"
    PHYLIP = "phylip"
    NEXUS = "nexus"
    GENBANK = "genbank"


@dataclass
class SequenceItem:
    """序列项数据类"""
    name: str
    sequence: str
    description: str = ""
    sequence_type: SequenceType = SequenceType.UNKNOWN
    quality_scores: Optional[List[float]] = None
    features: Dict[str, str] = field(default_factory=dict)
    
    def __post_init__(self):
        """初始化后处理"""
        if not self.sequence_type or self.sequence_type == SequenceType.UNKNOWN:
            self.sequence_type = self.detect_sequence_type()
            
    def detect_sequence_type(self) -> SequenceType:
        """检测序列类型"""
        seq = self.sequence.upper().replace('-', '').replace('N', '')
        if not seq:
            return SequenceType.UNKNOWN
            
        # 检查DNA特征
        dna_chars = set('ATCG')
        if all(c in dna_chars for c in seq):
            return SequenceType.DNA
            
        # 检查RNA特征
        rna_chars = set('AUCG')
        if all(c in rna_chars for c in seq):
            return SequenceType.RNA
            
        # 检查蛋白质特征
        protein_chars = set('ACDEFGHIKLMNPQRSTVWY')
        if all(c in protein_chars for c in seq):
            return SequenceType.PROTEIN
            
        return SequenceType.UNKNOWN
        
    def get_length(self) -> int:
        """获取序列长度"""
        return len(self.sequence)
        
    def get_gc_content(self) -> float:
        """获取GC含量"""
        if self.sequence_type not in [SequenceType.DNA, SequenceType.RNA]:
            return 0.0
            
        seq = self.sequence.upper().replace('-', '')
        if not seq:
            return 0.0
            
        gc_count = seq.count('G') + seq.count('C')
        return (gc_count / len(seq)) * 100
        
    def get_ambiguous_bases(self) -> int:
        """获取模糊碱基数量"""
        if self.sequence_type not in [SequenceType.DNA, SequenceType.RNA]:
            return 0
            
        ambiguous_chars = 'NMRWSYKVHDBN'
        return sum(1 for c in self.sequence.upper() if c in ambiguous_chars)
        
    def reverse(self) -> 'SequenceItem':
        """反向序列"""
        return SequenceItem(
            name=self.name,
            sequence=self.sequence[::-1],
            description=self.description,
            sequence_type=self.sequence_type,
            quality_scores=self.quality_scores[::-1] if self.quality_scores else None,
            features=self.features.copy()
        )
        
    def complement(self) -> 'SequenceItem':
        """互补序列"""
        if self.sequence_type == SequenceType.DNA:
            complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        elif self.sequence_type == SequenceType.RNA:
            complement_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        else:
            return self
            
        complemented = ''.join(complement_map.get(c.upper(), c) for c in self.sequence)
        return SequenceItem(
            name=self.name,
            sequence=complemented,
            description=self.description,
            sequence_type=self.sequence_type,
            quality_scores=self.quality_scores.copy() if self.quality_scores else None,
            features=self.features.copy()
        )
        
    def reverse_complement(self) -> 'SequenceItem':
        """反向互补序列"""
        return self.reverse().complement()
        
    def translate(self, genetic_code: int = 1) -> 'SequenceItem':
        """翻译为蛋白质序列"""
        if self.sequence_type != SequenceType.DNA:
            return self
            
        # 简化的遗传密码表
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        # 确保序列长度是3的倍数
        seq = self.sequence.upper().replace('-', '')
        if len(seq) % 3 != 0:
            seq = seq[:-(len(seq) % 3)]
            
        protein = ''
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            protein += codon_table.get(codon, 'X')
            
        return SequenceItem(
            name=f"{self.name}_protein",
            sequence=protein,
            description=f"Translated from {self.name}",
            sequence_type=SequenceType.PROTEIN,
            features=self.features.copy()
        )


class SequenceModel(QAbstractTableModel):
    """序列数据模型"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequences: List[SequenceItem] = []
        self.headers = ["Name", "Length", "Type", "GC%", "Description"]
        
    def rowCount(self, parent=QModelIndex()):
        """返回行数"""
        return len(self.sequences)
        
    def columnCount(self, parent=QModelIndex()):
        """返回列数"""
        return len(self.headers)
        
    def data(self, index, role=Qt.DisplayRole):
        """返回数据"""
        if not index.isValid() or index.row() >= len(self.sequences):
            return QVariant()
            
        sequence = self.sequences[index.row()]
        
        if role == Qt.DisplayRole:
            if index.column() == 0:  # Name
                return sequence.name
            elif index.column() == 1:  # Length
                return str(sequence.get_length())
            elif index.column() == 2:  # Type
                return sequence.sequence_type.value
            elif index.column() == 3:  # GC%
                return f"{sequence.get_gc_content():.1f}%"
            elif index.column() == 4:  # Description
                return sequence.description
                
        elif role == Qt.ToolTipRole:
            return f"Sequence: {sequence.name}\nLength: {sequence.get_length()}\nType: {sequence.sequence_type.value}"
            
        return QVariant()
        
    def headerData(self, section, orientation, role=Qt.DisplayRole):
        """返回表头数据"""
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.headers[section]
        return QVariant()
        
    def load_sequences(self, sequences: List[SequenceItem]):
        """加载序列"""
        self.beginResetModel()
        self.sequences = sequences
        self.endResetModel()
        
    def add_sequence(self, sequence: SequenceItem):
        """添加序列"""
        self.beginInsertRows(QModelIndex(), len(self.sequences), len(self.sequences))
        self.sequences.append(sequence)
        self.endInsertRows()
        
    def remove_sequence(self, row: int):
        """删除序列"""
        if 0 <= row < len(self.sequences):
            self.beginRemoveRows(QModelIndex(), row, row)
            del self.sequences[row]
            self.endRemoveRows()
            
    def get_sequence(self, row: int) -> Optional[SequenceItem]:
        """获取序列"""
        if 0 <= row < len(self.sequences):
            return self.sequences[row]
        return None
        
    def get_sequences(self) -> List[SequenceItem]:
        """获取所有序列"""
        return self.sequences.copy()
        
    def clear(self):
        """清空序列"""
        self.beginResetModel()
        self.sequences.clear()
        self.endResetModel()


class SequenceTableModel(QAbstractTableModel):
    """序列表格模型 - 用于显示序列矩阵"""
    
    def __init__(self, sequences: List[SequenceItem] = None, parent=None):
        super().__init__(parent)
        self.sequences = sequences or []
        self.max_length = self._calculate_max_length()
        
    def _calculate_max_length(self) -> int:
        """计算最大序列长度"""
        if not self.sequences:
            return 0
        return max(seq.get_length() for seq in self.sequences)
        
    def rowCount(self, parent=QModelIndex()):
        """返回行数（序列数）"""
        return len(self.sequences)
        
    def columnCount(self, parent=QModelIndex()):
        """返回列数（最大序列长度 + 1 for name）"""
        return self.max_length + 1
        
    def data(self, index, role=Qt.DisplayRole):
        """返回数据"""
        if not index.isValid():
            return QVariant()
            
        row = index.row()
        col = index.column()
        
        if row >= len(self.sequences):
            return QVariant()
            
        sequence = self.sequences[row]
        
        if role == Qt.DisplayRole:
            if col == 0:  # 序列名称
                return sequence.name
            else:  # 序列字符
                pos = col - 1
                if pos < len(sequence.sequence):
                    return sequence.sequence[pos]
                else:
                    return '-'
                    
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter
            
        return QVariant()
        
    def headerData(self, section, orientation, role=Qt.DisplayRole):
        """返回表头数据"""
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                if section == 0:
                    return "Name"
                else:
                    return str(section)
            else:
                return str(section + 1)
        return QVariant()
        
    def setData(self, index, value, role=Qt.EditRole):
        """设置数据"""
        if not index.isValid() or role != Qt.EditRole:
            return False
            
        row = index.row()
        col = index.column()
        
        if row >= len(self.sequences) or col == 0:
            return False
            
        sequence = self.sequences[row]
        pos = col - 1
        
        # 确保序列长度足够
        while len(sequence.sequence) <= pos:
            sequence.sequence += '-'
            
        # 更新字符
        sequence.sequence = sequence.sequence[:pos] + str(value) + sequence.sequence[pos+1:]
        
        self.dataChanged.emit(index, index)
        return True
        
    def flags(self, index):
        """返回项目标志"""
        if not index.isValid():
            return Qt.NoItemFlags
            
        if index.column() == 0:  # 名称列不可编辑
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
            
    def load_sequences(self, sequences: List[SequenceItem]):
        """加载序列"""
        self.beginResetModel()
        self.sequences = sequences
        self.max_length = self._calculate_max_length()
        self.endResetModel()
        
    def add_sequence(self, sequence: SequenceItem):
        """添加序列"""
        self.beginInsertRows(QModelIndex(), len(self.sequences), len(self.sequences))
        self.sequences.append(sequence)
        self.max_length = max(self.max_length, sequence.get_length())
        self.endInsertRows()
        
    def remove_sequence(self, row: int):
        """删除序列"""
        if 0 <= row < len(self.sequences):
            self.beginRemoveRows(QModelIndex(), row, row)
            del self.sequences[row]
            self.max_length = self._calculate_max_length()
            self.endRemoveRows()
            
    def get_sequences(self) -> List[SequenceItem]:
        """获取所有序列"""
        return self.sequences.copy()




















