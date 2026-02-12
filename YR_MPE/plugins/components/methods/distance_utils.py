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
    
    # 检查是否为核苷酸序列（替代饱和度分析仅支持核苷酸序列）
    
    # 定义核苷酸有效字符
    valid_nucleotides = set('ACGTURYKMSWBDHVNacgturykmswbdhvn-?!')  # 标准IUPAC核苷酸代码和gap
    
    for seq_record in sequences:
        seq_str = str(seq_record.seq)
        for char in seq_str:
            char_upper = char.upper()
            # 检查是否为有效核苷酸字符
            if char_upper not in valid_nucleotides:
                return False, f"Invalid character '{char}' found in sequence '{seq_record.id}'. Substitution saturation analysis only supports nucleotide sequences (A, C, G, T and IUPAC degenerate bases)."
    
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

def calculate_transition_transversion_counts(seq1: str, seq2: str, 
                                           gap_treatment: str = 'pairwise') -> Tuple[int, int, int]:
    """
    计算两个序列间的转换(Transition, Ti)、颠换(Transversion, Tv)和总比较位点数
    
    Args:
        seq1: 第一个序列字符串
        seq2: 第二个序列字符串
        gap_treatment: gap处理方式 ('pairwise' 或 'complete')
        
    Returns:
        (Ti数量, Tv数量, 总比较位点数)
    """
    ti_count = 0
    tv_count = 0
    total_sites = 0
    
    # 转换定义：A↔G, T↔C (嘌呤↔嘌呤, 嘧啶↔嘧啶)
    transitions = {('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')}
    
    # 颠换定义：嘌呤↔嘧啶
    transversions = {('A', 'T'), ('A', 'C'), ('G', 'T'), ('G', 'C'),
                     ('T', 'A'), ('T', 'G'), ('C', 'A'), ('C', 'G')}
    
    # 确保序列等长
    min_len = min(len(seq1), len(seq2))
    
    if gap_treatment.lower() == 'pairwise':
        # 成对删除：只跳过gap位点
        for i in range(min_len):
            base1 = seq1[i].upper()
            base2 = seq2[i].upper()
            
            # 跳过gap字符
            if base1 == '-' or base2 == '-':
                continue
            
            # 只处理标准碱基
            if base1 not in ['A', 'T', 'C', 'G'] or base2 not in ['A', 'T', 'C', 'G']:
                continue
            
            total_sites += 1
            
            if base1 == base2:
                # 相同位点，不计数
                continue
            
            # 检查是否为转换
            if (base1, base2) in transitions:
                ti_count += 1
            # 检查是否为颠换
            elif (base1, base2) in transversions:
                tv_count += 1
    
    else:  # complete deletion
        # 完全删除：排除所有含有gap的位点
        valid_positions = []
        for i in range(min_len):
            base1 = seq1[i].upper()
            base2 = seq2[i].upper()
            
            # 只有两个序列在该位置都不是gap时才保留
            if base1 != '-' and base2 != '-':
                if base1 in ['A', 'T', 'C', 'G'] and base2 in ['A', 'T', 'C', 'G']:
                    valid_positions.append(i)
        
        # 在有效位点上进行比较
        for i in valid_positions:
            base1 = seq1[i].upper()
            base2 = seq2[i].upper()
            
            total_sites += 1
            
            if base1 == base2:
                # 相同位点，不计数
                continue
            
            # 检查是否为转换
            if (base1, base2) in transitions:
                ti_count += 1
            # 检查是否为颠换
            elif (base1, base2) in transversions:
                tv_count += 1
    
    return ti_count, tv_count, total_sites

def calculate_saturation_metrics(sequences: List[SeqRecord], 
                                gap_treatment: str = 'pairwise') -> Dict:
    """
    计算替代饱和度指标
    
    Args:
        sequences: 序列记录列表
        gap_treatment: gap处理方式 ('pairwise' 或 'complete')
        
    Returns:
        包含以下键的字典：
        - p_values: 所有序列对的p值列表
        - ti_values: 所有序列对的Ti值列表
        - tv_values: 所有序列对的Tv值列表
        - ti_tv_ratios: 所有序列对的Ti/Tv比率列表
        - c_value: C值 (std(Ti/Tv) / std(p))
        - std_ti_tv: Ti/Tv比率的标准差
        - std_p: p值的标准差
        - stats: 其他统计信息
    """
    n = len(sequences)
    if n < 2:
        return {
            'p_values': [],
            'ti_values': [],
            'tv_values': [],
            'ti_tv_ratios': [],
            'c_value': 0.0,
            'std_ti_tv': 0.0,
            'std_p': 0.0,
            'n_pairs': 0
        }
    
    # 如果是complete deletion，需要先在整个数据集级别删除含有gap的位点
    if gap_treatment.lower() == 'complete':
        sequences = _apply_complete_deletion(sequences)
        if len(sequences) < 2:
            return {
                'p_values': [],
                'ti_values': [],
                'tv_values': [],
                'ti_tv_ratios': [],
                'c_value': 0.0,
                'std_ti_tv': 0.0,
                'std_p': 0.0,
                'n_pairs': 0
            }
    
    p_values = []
    ti_values = []
    tv_values = []
    ti_tv_ratios = []
    
    for i in range(n):
        for j in range(i+1, n):
            ti, tv, total = calculate_transition_transversion_counts(
                str(sequences[i].seq), 
                str(sequences[j].seq),
                'pairwise'  # 完成complete deletion后，使用pairwise方式计算
            )
            
            if total > 0:
                p = (ti + tv) / total
                
                # 计算Ti/Tv比率，避免除零
                if tv > 0:
                    ti_tv_ratio = ti / tv
                else:
                    # 当Tv=0时，设置为一个较大的值表示饱和
                    ti_tv_ratio = 10.0  # 或使用float('inf')
                
                p_values.append(p)
                ti_values.append(ti)
                tv_values.append(tv)
                ti_tv_ratios.append(ti_tv_ratio)
    
    # 计算标准差
    if len(p_values) >= 2:
        p_array = np.array(p_values)
        ti_tv_array = np.array(ti_tv_ratios)
        
        std_p = np.std(p_array, ddof=1)  # 样本标准差（ddof=1）
        std_ti_tv = np.std(ti_tv_array, ddof=1)
        
        # 计算C值：C = std(Ti/Tv) / std(p)
        c_value = std_ti_tv / std_p if std_p > 0 else 0.0
    else:
        std_p = 0.0
        std_ti_tv = 0.0
        c_value = 0.0
    
    return {
        'p_values': p_values,
        'ti_values': ti_values,
        'tv_values': tv_values,
        'ti_tv_ratios': ti_tv_ratios,
        'c_value': c_value,
        'std_ti_tv': std_ti_tv,
        'std_p': std_p,
        'n_pairs': len(p_values)
    }

def _apply_complete_deletion(sequences: List[SeqRecord]) -> List[SeqRecord]:
    """
    应用complete deletion：排除所有含有gap的位点
    
    Args:
        sequences: 序列记录列表
        
    Returns:
        处理后的序列列表
    """
    if not sequences:
        return []
    
    # 检查序列长度，找到最大长度
    seq_lengths = [seq.seq.__len__() for seq in sequences]
    max_length = max(seq_lengths)
    min_length = min(seq_lengths)
    
    # 如果序列长度差异太大，使用最小长度
    seq_length = min_length
    
    # 找出所有不含gap的列索引
    keep_columns = []
    for col_idx in range(seq_length):
        column_has_gap = False
        for seq in sequences:
            char = str(seq.seq[col_idx]).upper()
            if char in ['-', '.']:  # 检查gap字符
                column_has_gap = True
                break
        if not column_has_gap:
            keep_columns.append(col_idx)
    
    # 如果没有保留的列，返回空列表
    if not keep_columns:
        return []
    
    # 构建新序列
    cleaned_sequences = []
    for seq in sequences:
        # 只保留无gap的列
        new_seq = ''.join(str(seq.seq[i]) for i in keep_columns if i < seq.seq.__len__())
        cleaned_sequences.append(SeqRecord(
            Seq(new_seq),
            id=seq.id,
            description=seq.description
        ))
    
    return cleaned_sequences