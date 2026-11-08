"""
分区工具模块
提供超级矩阵合并和分区定义创建的通用函数
"""

from typing import List, Dict, Tuple, Optional
from Bio.SeqRecord import SeqRecord


def create_supermatrix_from_alignments(
    alignment_items: List,
    output_file: Optional[str] = None
) -> Tuple[str, List[Dict]]:
    """
    从多个比对文件创建超级矩阵
    
    Args:
        alignment_items: DatasetItem 列表，每个代表一个分区比对
        output_file: 输出文件路径（可选），如果不提供则返回文件内容
    
    Returns:
        tuple: (supermatrix_file_path: str,
                partition_definitions: List[Dict])
    
    Raises:
        ValueError: 如果没有提供任何比对项目
    """
    if not alignment_items:
        raise ValueError("No alignment items provided for supermatrix creation")
    
    # 合并所有选中的 partition 为 supermatrix
    supermatrix_sequences = {}
    partition_definitions = []
    current_pos = 1
    
    # 获取所有序列名称（假设所有 partition 有相同的序列名称集合）
    if alignment_items and alignment_items[0].sequences:
        all_seq_names = [seq.id for seq in alignment_items[0].sequences]
    else:
        all_seq_names = []
    
    # 合并序列
    for seq_name in all_seq_names:
        supermatrix_seq = ""
        for item in alignment_items:
            # 找到对应的序列
            seq = next((s for s in item.sequences if s.id == seq_name), None)
            if seq:
                supermatrix_seq += str(seq.seq)
            else:
                # 关键：如果某个 partition 缺少该序列，用 ? 填充
                supermatrix_seq += "?" * item.length
        supermatrix_sequences[seq_name] = supermatrix_seq
    
    # 创建超级矩阵文件
    if output_file is None:
        import tempfile
        output_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False).name
    
    with open(output_file, 'w') as f:
        for seq_name, seq_content in supermatrix_sequences.items():
            f.write(f">{seq_name}\n{seq_content}\n")
    
    # 计算 partition 坐标并创建 partition definitions
    for item in alignment_items:
        end_pos = current_pos + item.length - 1
        partition_def = {
            'name': item.loci_name,
            'range': f"{current_pos}-{end_pos}",
            'seq_type': detect_sequence_type(item.sequences[0]) if item.sequences else "DNA",
            'model': item.data.get('model', None) if isinstance(item.data, dict) else None,
            'length': item.length
        }
        partition_definitions.append(partition_def)
        current_pos = end_pos + 1
    
    return output_file, partition_definitions


def detect_sequence_type(seq_record: SeqRecord) -> str:
    """
    检测序列类型（DNA 或 Protein）
    
    Args:
        seq_record: BioPython SeqRecord 对象
    
    Returns:
        str: "DNA" 或 "PROTEIN"
    """
    if not seq_record:
        return "DNA"
    
    seq_str = str(seq_record.seq).upper()
    
    # 统计 DNA 字符和 Protein 字符
    dna_chars = set('ACGTURYMKSWHBVDN')
    protein_chars = set('ACDEFGHIKLMNPQRSTVWY*')
    
    seq_chars = set(seq_str)
    
    # 计算 DNA 字符和 Protein 字符的比例
    dna_count = sum(1 for c in seq_str if c in dna_chars)
    protein_count = sum(1 for c in seq_str if c in protein_chars)
    
    total_chars = len(seq_str.replace('-', '').replace('?', ''))
    
    if total_chars == 0:
        return "DNA"
    
    dna_ratio = dna_count / total_chars
    protein_ratio = protein_count / total_chars
    
    # 如果 DNA 字符比例超过 0.9，认为是 DNA
    if dna_ratio > 0.9:
        return "DNA"
    else:
        return "PROTEIN"


def detect_sequence_type_conflicts(
    alignment_items: List
) -> Tuple[bool, str]:
    """
    检测多个比对之间的序列类型冲突
    
    Args:
        alignment_items: DatasetItem 列表
    
    Returns:
        tuple: (has_conflict: bool, message: str)
    """
    if not alignment_items:
        return False, ""
    
    seq_types = []
    for item in alignment_items:
        if item.sequences:
            seq_type = detect_sequence_type(item.sequences[0])
            seq_types.append((item.loci_name, seq_type))
    
    unique_types = set(st[1] for st in seq_types)
    
    if len(unique_types) > 1:
        # 有冲突
        type_details = ", ".join([f"{name} ({typ})" for name, typ in seq_types])
        message = f"Sequence type conflict detected: {type_details}"
        return True, message
    
    return False, ""


def create_partition_definitions_from_items(
    alignment_items: List,
    supermatrix_start_pos: int = 1
) -> List[Dict]:
    """
    从比对项目创建分区定义
    
    Args:
        alignment_items: DatasetItem 列表
        supermatrix_start_pos: 超级矩阵的起始位置（1-based）
    
    Returns:
        分区定义列表
    """
    partition_definitions = []
    current_pos = supermatrix_start_pos
    
    for item in alignment_items:
        end_pos = current_pos + item.length - 1
        partition_def = {
            'name': item.loci_name,
            'range': f"{current_pos}-{end_pos}",
            'seq_type': detect_sequence_type(item.sequences[0]) if item.sequences else "DNA",
            'model': item.data.get('model', None) if isinstance(item.data, dict) else None,
            'length': item.length
        }
        partition_definitions.append(partition_def)
        current_pos = end_pos + 1
    
    return partition_definitions


if __name__ == "__main__":
    # 测试代码
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from ..plugins.templates.base_plugin_ui import BasePlugin
    
    # 创建测试数据
    class TestItem:
        def __init__(self, name, length, sequences):
            self.loci_name = name
            self.length = length
            self.sequences = sequences
            self.data = {}
    
    # 分区1: 4个序列，每个10 bp
    item1 = TestItem(
        "partition1",
        10,
        [
            SeqRecord(Seq("ATCGATCGAT"), id="seq1"),
            SeqRecord(Seq("ATCGATCGAT"), id="seq2"),
            SeqRecord(Seq("ATCGATCGAT"), id="seq3"),
        ]
    )
    
    # 分区2: 3个序列，每个8 bp（seq3缺失）
    item2 = TestItem(
        "partition2",
        8,
        [
            SeqRecord(Seq("TGCATGCA"), id="seq1"),
            SeqRecord(Seq("TGCATGCA"), id="seq2"),
        ]
    )
    
    # 创建超级矩阵
    import tempfile
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False).name
    supermatrix_file, partitions = create_supermatrix_from_alignments([item1, item2], temp_file)
    
    print("Supermatrix file:", supermatrix_file)
    print("\nPartitions:")
    for part in partitions:
        print(f"  {part['name']}: {part['range']} ({part['length']} bp, {part['seq_type']})")
    
    # 检测序列类型冲突
    has_conflict, message = detect_sequence_type_conflicts([item1, item2])
    if has_conflict:
        print(f"\nConflict detected: {message}")
    else:
        print("\nNo sequence type conflicts detected")
    
    # 打印超级矩阵内容
    print("\nSupermatrix content:")
    with open(supermatrix_file, 'r') as f:
        print(f.read())