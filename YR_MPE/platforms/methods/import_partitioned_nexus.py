"""
导入分区 NEXUS 文件
支持从分区 NEXUS 文件中导入序列数据和分区设置
"""
import re
import os
from typing import List, Dict, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# 创建临时 DatasetItem 类以避免循环导入
class TempDatasetItem:
    """临时数据项类 - 用于导入分区 NEXUS"""
    def __init__(self):
        self.selected = False
        self.loci_name = ""
        self.length = 0
        self.sequence_count = 0
        self.is_aligned = False
        self.file_path = ""
        self.sequences = []
        self.name = ""


def _parse_nexus_with_regex(content: str) -> List[SeqRecord]:
    """
    使用正则表达式解析 NEXUS 文件（BioPython 失败时的回退方案）
    
    Args:
        content: NEXUS 文件内容
    
    Returns:
        SeqRecord 对象列表
    """
    records = []
    
    # 查找 matrix 块
    matrix_pattern = r'matrix\s*([\s\S]*?)(?=\s*end\s*;)'
    matrix_match = re.search(matrix_pattern, content)
    
    if not matrix_match:
        return records
    
    matrix_content = matrix_match.group(1)
    
    # 解析序列 - 支持多行序列（interleave 格式）
    sequences = {}
    current_name = None
    current_seq = []
    
    for line in matrix_content.split('\n'):
        line = line.strip()
        
        # 跳过空行和注释
        if not line or line.startswith(';'):
            continue
        
        # 检查是否是新的序列名
        # 如果行以字母开头且后续有空格，可能是序列名
        if re.match(r'^[A-Za-z_]', line) and re.search(r'\s', line):
            parts = line.split(maxsplit=1)
            if len(parts) >= 2:
                # 保存前一个序列
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                
                # 开始新序列
                current_name = parts[0]
                current_seq = [parts[1]]
                continue
        
        # 否则是序列的延续
        if current_name:
            # 移除空格和非序列字符
            seq_part = line.replace(' ', '').replace('\t', '')
            current_seq.append(seq_part)
    
    # 保存最后一个序列
    if current_name:
        sequences[current_name] = ''.join(current_seq)
    
    # 转换为 SeqRecord 对象
    for seq_name, seq_content in sequences.items():
        if seq_content:  # 只添加非空序列
            record = SeqRecord(
                Seq(seq_content),
                id=seq_name,
                description=""
            )
            records.append(record)
    
    return records


def parse_charset_range(range_str: str, total_length: int) -> List[int]:
    """
    解析 charset 范围定义，返回位点索引列表
    
    支持的格式:
    - 1-1000
    - 1-1000\3 (位点 1, 4, 7, ...)
    - 1-1000, 2000-3000 (多个范围)
    - 1, 5, 10 (单个位点)
    
    Args:
        range_str: 范围定义字符串
        total_length: 序列总长度
    
    Returns:
        位点索引列表（0-based）
    """
    positions = []
    
    # 移除空格
    range_str = range_str.strip()
    
    # 处理多个范围（用逗号分隔）
    parts = range_str.split(',')
    
    for part in parts:
        part = part.strip()
        
        # 检查是否有步长 (格式: start-end\step)
        step = 1
        if '\\' in part:
            range_part, step_part = part.split('\\', 1)
            step = int(step_part.strip())
            part = range_part.strip()
        
        # 处理单个位点
        if '-' not in part:
            pos = int(part) - 1  # 转换为 0-based
            if 0 <= pos < total_length:
                positions.append(pos)
            continue
        
        # 处理范围
        start, end = part.split('-', 1)
        start = int(start.strip()) - 1  # 转换为 0-based
        end = int(end.strip()) - 1      # 转换为 0-based
        
        # 修正范围
        if start < 0:
            start = 0
        if end >= total_length:
            end = total_length - 1
        
        # 生成位置列表
        for pos in range(start, end + 1, step):
            if 0 <= pos < total_length:
                positions.append(pos)
    
    return sorted(set(positions))  # 去重并排序


def import_partitioned_nexus(file_path: str) -> Tuple[List[TempDatasetItem], str, Dict]:
    """
    导入分区 NEXUS 文件

    Args:
        file_path: NEXUS 文件路径

    Returns:
        tuple: (dataset_items: List[TempDatasetItem],
                partition_scheme: str,
                summary: Dict)

    Raises:
        FileNotFoundError: 文件不存在
        ValueError: 文件格式错误或解析失败
    """
    # 1. 检查文件存在性
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件不存在: {file_path}")

    # 2. 读取文件内容（用于回退方案）
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 3. 使用 BioPython 读取序列（支持 interleave 格式）
    # 如果 BioPython 失败（如跨行 format 定义），回退到正则表达式方法
    try:
        records = list(SeqIO.parse(file_path, "nexus"))
    except Exception as e:
        print(f"Warning: BioPython 解析失败，使用正则表达式方法: {e}")
        # 回退：使用正则表达式读取序列
        records = _parse_nexus_with_regex(content)
    
    if not records:
        raise ValueError("NEXUS 文件中未找到序列")

    # 4. 提取序列信息
    ntax = len(records)
    nchar = len(records[0].seq) if records else 0

    # 验证所有序列长度一致
    for i, record in enumerate(records):
        if len(record.seq) != nchar:
            import warnings
            warnings.warn(f"序列 '{record.id}' 长度不一致: {len(record.seq)} (期望 {nchar})")

    # 4. 解析 partition 信息（手动读取 sets 块）
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    sets_pattern = r'begin\s+sets\s*;([\s\S]*?)end\s*;'
    sets_match = re.search(sets_pattern, content, re.IGNORECASE)

    partitions = []
    partition_names = []

    if sets_match:
        sets_content = sets_match.group(1)

        # 解析 charset 定义
        charset_pattern = r'charset\s+(\S+)\s*=\s*([^;]+);'
        charset_matches = re.findall(charset_pattern, sets_content, re.IGNORECASE)

        for name, range_def in charset_matches:
            name = name.strip()
            range_def = range_def.strip()

            # 跳过外部文件引用（暂时不支持）
            if ':' in range_def and not re.match(r'^\d', range_def):
                continue

            # 解析位点范围
            positions = parse_charset_range(range_def, nchar)

            if positions:
                partitions.append({
                    'name': name,
                    'positions': positions,
                    'range_def': range_def
                })
                partition_names.append(name)

    # 如果没有找到分区，创建单个分区
    if not partitions:
        partitions.append({
            'name': 'full',
            'positions': list(range(nchar)),
            'range_def': f'1-{nchar}'
        })
        partition_names.append('full')
    
    # 6. 拆分 BioPython 读取的序列为多个 DatasetItem
        dataset_items = []
    
        for partition in partitions:
            partition_name = partition['name']
            positions = partition['positions']
    
            if not positions:
                continue
    
            # 创建 DatasetItem
            dataset_item = TempDatasetItem()
            dataset_item.loci_name = partition_name
            dataset_item.file_path = file_path
            dataset_item.format = 'nexus'
    
            # 提取该分区的序列（从 BioPython 的 SeqRecord 对象）
            item_sequences = []
            for record in records:
                # 提取分区的序列片段
                full_seq = str(record.seq)
                partition_seq = ''.join([full_seq[pos] for pos in positions])
    
                # 创建新的 SeqRecord
                seq_record = SeqRecord(
                    Seq(partition_seq),
                    id=record.id,
                    description=f"{partition_name} partition"
                )
                item_sequences.append(seq_record)
    
            dataset_item.sequences = item_sequences
            dataset_item.length = len(positions)
            dataset_item.sequence_count = len(item_sequences)
    
            # 检查是否已比对（简单判断）
            all_chars = ''.join(str(record.seq) for record in item_sequences)
            has_missing = '?' in all_chars or '-' in all_chars
            dataset_item.is_aligned = not has_missing or (has_missing and len(set(all_chars) - set('?-\n')) > 1)
    
            dataset_items.append(dataset_item)
    
        # 7. 生成分区方案字符串
        partition_scheme = "begin sets;\n"
        current_pos = 1
    
        for i, dataset_item in enumerate(dataset_items):
            end_pos = current_pos + dataset_item.length - 1
            partition_scheme += f"    charset {dataset_item.loci_name} = {current_pos}-{end_pos};\n"
            current_pos = end_pos + 1
    
        partition_scheme += "end;"
    
        # 8. 生成摘要信息
        summary = {
            'source_file': file_path,
            'ntax': ntax,
            'nchar': nchar,
            'datatype': records[0].annotations.get('molecule_type', 'DNA') if records else 'DNA',
            'partition_count': len(dataset_items),
            'partition_names': partition_names,
            'total_taxa': dataset_items[0].sequence_count if dataset_items else 0
        }
    
        return dataset_items, partition_scheme, summary

if __name__ == "__main__":
    # 测试代码
    import sys
    
    if len(sys.argv) < 2:
        print("用法: python import_partitioned_nexus.py <nexus_file>")
        sys.exit(1)
    
    nexus_file = sys.argv[1]
    
    try:
        dataset_items, partition_scheme, summary = import_partitioned_nexus(nexus_file)
        
        print("=" * 60)
        print("导入成功！")
        print("=" * 60)
        print(f"\n摘要信息:")
        for key, value in summary.items():
            print(f"  {key}: {value}")
        
        print(f"\n分区方案:")
        print(partition_scheme)
        
        print(f"\n数据项:")
        for i, item in enumerate(dataset_items, 1):
            print(f"  {i}. {item.loci_name}")
            print(f"     长度: {item.length}")
            print(f"     序列数: {item.sequence_count}")
            print(f"     已比对: {item.is_aligned}")
        
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)