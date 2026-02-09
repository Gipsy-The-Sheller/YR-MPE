"""
Distance calculation utilities for uncorrected genetic distance (p-distance)
Supports DNA and protein sequences with IUPAC degenerate base handling
"""

import numpy as np
from typing import List, Tuple, Dict, Set, Optional
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# IUPAC简并碱基映射
IUPAC_DNA_MAP = {
    'R': ['A', 'G'],      # 嘌呤
    'Y': ['C', 'T'],      # 嘧啶
    'M': ['A', 'C'],      # 氨基
    'K': ['G', 'T'],      # 酮基
    'S': ['C', 'G'],      # 强氢键
    'W': ['A', 'T'],      # 弱氢键
    'H': ['A', 'C', 'T'], # 非G
    'B': ['C', 'G', 'T'], # 非A
    'V': ['A', 'C', 'G'], # 非T
    'D': ['A', 'G', 'T'], # 非C
    'N': ['A', 'C', 'G', 'T']  # 任意
}

# DNA替代类型定义
DNA_SUBSTITUTION_TYPES = {
    'transition': [('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')],
    'transversion': [('A', 'T'), ('A', 'C'), ('G', 'T'), ('G', 'C'), 
                     ('T', 'A'), ('T', 'G'), ('C', 'A'), ('C', 'G')]
}

def expand_degenerate_base(base: str) -> List[str]:
    """
    展开简并碱基为可能的基础碱基
    
    Args:
        base: 单个碱基字符（可能是简并碱基）
        
    Returns:
        基础碱基列表
    """
    base = base.upper()
    if base in IUPAC_DNA_MAP:
        return IUPAC_DNA_MAP[base]
    elif base in ['A', 'T', 'C', 'G']:
        return [base]
    else:
        # 对于非标准字符，返回空列表
        return []

def calculate_base_difference(base1: str, base2: str, substitution_filter: Optional[Set[Tuple[str, str]]] = None) -> float:
    """
    计算两个碱基之间的差异度
    
    Args:
        base1: 第一个碱基
        base2: 第二个碱基
        substitution_filter: 可选的替代类型过滤器
        
    Returns:
        差异度 (0.0-1.0)，0表示相同，1表示完全不同
    """
    # 展开简并碱基
    bases1 = expand_degenerate_base(base1)
    bases2 = expand_degenerate_base(base2)
    
    if not bases1 or not bases2:
        # 无法解析的字符视为完全差异
        return 1.0
    
    # 计算所有可能组合的最小差异
    min_diff = 1.0
    for b1 in bases1:
        for b2 in bases2:
            # 检查替代类型过滤
            if substitution_filter and (b1, b2) not in substitution_filter:
                continue
                
            diff = 0.0 if b1 == b2 else 1.0
            min_diff = min(min_diff, diff)
    
    return min_diff

def pairwise_deletion_compare(seq1: str, seq2: str,
                            substitution_filter: Optional[Set[Tuple[str, str]]] = None) -> Tuple[int, int]:
    """
    使用成对删除法比较两个序列
    
    Args:
        seq1: 第一个序列
        seq2: 第二个序列
        substitution_filter: 可选的DNA替代类型过滤器
        
    Returns:
        (差异位点数, 总比较位点数)
    """
    differences = 0.0
    total_sites = 0
    
    # 确保序列等长
    min_len = min(len(seq1), len(seq2))
    
    for i in range(min_len):
        base1 = seq1[i].upper()
        base2 = seq2[i].upper()
        
        # 跳过gap字符
        if base1 == '-' or base2 == '-':
            continue
            
        total_sites += 1
        diff = calculate_base_difference(base1, base2, substitution_filter)
        differences += diff
    
    return differences, total_sites

def complete_deletion_compare(seq1: str, seq2: str,
                            substitution_filter: Optional[Set[Tuple[str, str]]] = None) -> Tuple[int, int]:
    """
    使用完全删除法比较两个序列（Complete Deletion）
    排除所有含有gap的位点
    
    Args:
        seq1: 第一个序列
        seq2: 第二个序列
        substitution_filter: 可选的DNA替代类型过滤器
        
    Returns:
        (差异位点数, 总比较位点数)
    """
    differences = 0.0
    total_sites = 0
    
    # 确保序列等长
    min_len = min(len(seq1), len(seq2))
    
    # 首先找出所有不含gap的位点
    valid_positions = []
    for i in range(min_len):
        base1 = seq1[i].upper()
        base2 = seq2[i].upper()
        
        # 只有当两个序列在该位置都不是gap时才保留
        if base1 != '-' and base2 != '-':
            valid_positions.append(i)
    
    # 在有效位点上进行比较
    for i in valid_positions:
        base1 = seq1[i].upper()
        base2 = seq2[i].upper()
        
        diff = calculate_base_difference(base1, base2, substitution_filter)
        differences += diff
        total_sites += 1
    
    return differences, total_sites

def calculate_p_distance(seq1: str, seq2: str, gap_treatment: str = 'pairwise',
                        substitution_filter: Optional[Set[Tuple[str, str]]] = None) -> float:
    """
    计算两个序列间的未校正遗传距离(p-distance)
    
    Args:
        seq1: 第一个序列
        seq2: 第二个序列
        gap_treatment: gap处理方式 ('pairwise' 或 'complete')
        substitution_filter: 可选的DNA替代类型过滤器
        
    Returns:
        p-distance值 (0.0-1.0)
    """
    if gap_treatment.lower() == 'pairwise':
        differences, total_sites = pairwise_deletion_compare(seq1, seq2, substitution_filter)
    else:  # complete deletion (原来是partial)
        differences, total_sites = complete_deletion_compare(seq1, seq2, substitution_filter)
    
    if total_sites == 0:
        return 0.0  # 无有效比较位点时返回0
    
    return differences / total_sites

def calculate_distance_matrix(sequences: List[SeqRecord], gap_treatment: str = 'pairwise',
                            substitution_filter: Optional[Set[Tuple[str, str]]] = None) -> np.ndarray:
    """
    计算序列集的距离矩阵
    
    Args:
        sequences: 序列记录列表
        gap_treatment: gap处理方式
        substitution_filter: 可选的DNA替代类型过滤器
        
    Returns:
        距离矩阵 (numpy数组)
    """
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    
    # 计算上三角矩阵（对称矩阵）
    for i in range(n):
        for j in range(i + 1, n):
            dist = calculate_p_distance(
                str(sequences[i].seq), 
                str(sequences[j].seq),
                gap_treatment,
                substitution_filter
            )
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist  # 对称性
    
    return distance_matrix

def get_substitution_filter(rule_type: str) -> Optional[Set[Tuple[str, str]]]:
    """
    根据规则类型获取替代过滤器
    
    Args:
        rule_type: 规则类型 ('all', 'transition', 'transversion', 'custom')
        
    Returns:
        替代类型集合，None表示不过滤
    """
    if rule_type == 'all':
        return None  # 不过滤，包含所有替代
    elif rule_type == 'transition':
        return set(DNA_SUBSTITUTION_TYPES['transition'])
    elif rule_type == 'transversion':
        return set(DNA_SUBSTITUTION_TYPES['transversion'])
    else:
        return None  # custom需要额外参数

def validate_sequences(sequences: List[SeqRecord]) -> Tuple[bool, str]:
    """
    验证序列集的有效性
    
    Args:
        sequences: 序列记录列表
        
    Returns:
        (是否有效, 错误信息)
    """
    if not sequences:
        return False, "No sequences provided"
    
    if len(sequences) < 2:
        return False, "At least 2 sequences required"
    
    # 检查序列长度一致性
    first_length = len(sequences[0].seq)
    for seq_record in sequences[1:]:
        if len(seq_record.seq) != first_length:
            return False, f"Inconsistent sequence lengths: expected {first_length}, got {len(seq_record.seq)}"
    
    return True, ""

def format_distance_matrix(distance_matrix: np.ndarray, sequence_names: List[str]) -> str:
    """
    格式化距离矩阵为文本输出 Format: phylip matrix
    
    Args:
        distance_matrix: 距离矩阵
        sequence_names: 序列名称列表
        
    Returns:
        格式化的距离矩阵字符串
    """
    n = len(sequence_names)
    lines = []
    
    # # 表头
    # header = " " * 15 + "".join(f"{name:>12}" for name in sequence_names)
    # lines.append(header)
    # lines.append("-" * len(header))

    lines.append(f"{n}")

    name_len_max = max(len(name) for name in sequence_names)
    
    # 矩阵数据
    for i in range(n):
        row_data = f"{sequence_names[i]}{' ' * (name_len_max - len(sequence_names[i]))}"
        for j in range(n):
            row_data = f"{row_data} {distance_matrix[i, j]:>12.6f}"
        lines.append(row_data)
    
    return "\n".join(lines)