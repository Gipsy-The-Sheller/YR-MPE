#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Sequence Utilities for YR-MPE Sequence Editor
Utility functions for sequence processing and file I/O
"""

import os
import re
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path

try:
    from .sequence_models import SequenceItem, SequenceType, SequenceFormat
except ImportError:
    from sequence_models import SequenceItem, SequenceType, SequenceFormat


class SequenceValidator:
    """序列验证器"""
    
    @staticmethod
    def validate_dna_sequence(sequence: str) -> Tuple[bool, str]:
        """验证DNA序列"""
        sequence = sequence.upper().replace('-', '').replace('N', '')
        valid_chars = set('ATCG')
        
        if not sequence:
            return False, "Empty sequence"
            
        invalid_chars = set(sequence) - valid_chars
        if invalid_chars:
            return False, f"Invalid characters: {', '.join(invalid_chars)}"
            
        return True, "Valid DNA sequence"
        
    @staticmethod
    def validate_rna_sequence(sequence: str) -> Tuple[bool, str]:
        """验证RNA序列"""
        sequence = sequence.upper().replace('-', '').replace('N', '')
        valid_chars = set('AUCG')
        
        if not sequence:
            return False, "Empty sequence"
            
        invalid_chars = set(sequence) - valid_chars
        if invalid_chars:
            return False, f"Invalid characters: {', '.join(invalid_chars)}"
            
        return True, "Valid RNA sequence"
        
    @staticmethod
    def validate_protein_sequence(sequence: str) -> Tuple[bool, str]:
        """验证蛋白质序列"""
        sequence = sequence.upper().replace('-', '')
        valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
        
        if not sequence:
            return False, "Empty sequence"
            
        invalid_chars = set(sequence) - valid_chars
        if invalid_chars:
            return False, f"Invalid characters: {', '.join(invalid_chars)}"
            
        return True, "Valid protein sequence"
        
    @staticmethod
    def validate_sequence(sequence: str, sequence_type: SequenceType) -> Tuple[bool, str]:
        """验证序列"""
        if sequence_type == SequenceType.DNA:
            return SequenceValidator.validate_dna_sequence(sequence)
        elif sequence_type == SequenceType.RNA:
            return SequenceValidator.validate_rna_sequence(sequence)
        elif sequence_type == SequenceType.PROTEIN:
            return SequenceValidator.validate_protein_sequence(sequence)
        else:
            return True, "Unknown sequence type"


class SequenceUtils:
    """序列工具类"""
    
    @staticmethod
    def load_sequences(file_path: str) -> List[SequenceItem]:
        """加载序列文件"""
        file_path = Path(file_path)
        file_format = SequenceUtils._detect_format(file_path)
        
        if file_format == SequenceFormat.FASTA:
            return SequenceUtils._load_fasta(file_path)
        elif file_format == SequenceFormat.PHYLIP:
            return SequenceUtils._load_phylip(file_path)
        elif file_format == SequenceFormat.NEXUS:
            return SequenceUtils._load_nexus(file_path)
        elif file_format == SequenceFormat.GENBANK:
            return SequenceUtils._load_genbank(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")
            
    @staticmethod
    def save_sequences(sequences: List[SequenceItem], file_path: str, 
                      file_format: Optional[SequenceFormat] = None):
        """保存序列文件"""
        file_path = Path(file_path)
        
        if file_format is None:
            file_format = SequenceUtils._detect_format(file_path)
            
        if file_format == SequenceFormat.FASTA:
            SequenceUtils._save_fasta(sequences, file_path)
        elif file_format == SequenceFormat.PHYLIP:
            SequenceUtils._save_phylip(sequences, file_path)
        elif file_format == SequenceFormat.NEXUS:
            SequenceUtils._save_nexus(sequences, file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_format}")
            
    @staticmethod
    def _detect_format(file_path: Path) -> SequenceFormat:
        """检测文件格式"""
        suffix = file_path.suffix.lower()
        
        if suffix in ['.fasta', '.fas', '.fa']:
            return SequenceFormat.FASTA
        elif suffix in ['.phy', '.phylip']:
            return SequenceFormat.PHYLIP
        elif suffix in ['.nex', '.nexus']:
            return SequenceFormat.NEXUS
        elif suffix in ['.gb', '.genbank']:
            return SequenceFormat.GENBANK
        else:
            # 尝试从文件内容检测
            return SequenceUtils._detect_format_from_content(file_path)
            
    @staticmethod
    def _detect_format_from_content(file_path: Path) -> SequenceFormat:
        """从文件内容检测格式"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                first_line = f.readline().strip()
                
            if first_line.startswith('>'):
                return SequenceFormat.FASTA
            elif first_line.startswith('#NEXUS') or first_line.startswith('#nexus'):
                return SequenceFormat.NEXUS
            elif first_line.startswith('LOCUS'):
                return SequenceFormat.GENBANK
            else:
                # 尝试解析为Phylip格式
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        lines = f.readlines()
                    if len(lines) >= 2:
                        first_line_parts = lines[0].strip().split()
                        if len(first_line_parts) == 2 and first_line_parts[0].isdigit():
                            return SequenceFormat.PHYLIP
                except:
                    pass
                    
            return SequenceFormat.FASTA  # 默认格式
        except:
            return SequenceFormat.FASTA
            
    @staticmethod
    def _load_fasta(file_path: Path) -> List[SequenceItem]:
        """加载FASTA格式文件"""
        sequences = []
        current_name = ""
        current_sequence = ""
        current_description = ""
        
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # 保存前一个序列
                    if current_name:
                        sequences.append(SequenceItem(
                            name=current_name,
                            sequence=current_sequence,
                            description=current_description
                        ))
                    
                    # 开始新序列
                    header = line[1:].strip()
                    if ' ' in header:
                        current_name, current_description = header.split(' ', 1)
                    else:
                        current_name = header
                        current_description = ""
                    current_sequence = ""
                else:
                    current_sequence += line
                    
        # 保存最后一个序列
        if current_name:
            sequences.append(SequenceItem(
                name=current_name,
                sequence=current_sequence,
                description=current_description
            ))
            
        return sequences
        
    @staticmethod
    def _load_phylip(file_path: Path) -> List[SequenceItem]:
        """加载Phylip格式文件"""
        sequences = []
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        if len(lines) < 2:
            return sequences
            
        # 解析第一行
        first_line = lines[0].strip().split()
        num_sequences = int(first_line[0])
        sequence_length = int(first_line[1])
        
        # 解析序列
        current_sequence = ""
        current_name = ""
        
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
                
            # 检查是否是新的序列开始
            if len(current_sequence) == 0:
                # 新序列开始
                parts = line.split(None, 1)
                current_name = parts[0]
                current_sequence = parts[1] if len(parts) > 1 else ""
            else:
                # 继续当前序列
                current_sequence += line
                
            # 如果序列长度达到预期，保存序列
            if len(current_sequence.replace('-', '')) >= sequence_length:
                sequences.append(SequenceItem(
                    name=current_name,
                    sequence=current_sequence[:sequence_length],
                    description=""
                ))
                current_sequence = ""
                current_name = ""
                
        return sequences
        
    @staticmethod
    def _load_nexus(file_path: Path) -> List[SequenceItem]:
        """加载Nexus格式文件"""
        sequences = []
        
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
        # 查找DATA块
        data_match = re.search(r'BEGIN\s+DATA;.*?END;', content, re.DOTALL | re.IGNORECASE)
        if not data_match:
            return sequences
            
        data_block = data_match.group(0)
        
        # 解析序列
        sequence_pattern = r'(\w+)\s+([ATCGN-]+)'
        matches = re.findall(sequence_pattern, data_block, re.IGNORECASE)
        
        for name, sequence in matches:
            sequences.append(SequenceItem(
                name=name,
                sequence=sequence.upper(),
                description=""
            ))
            
        return sequences
        
    @staticmethod
    def _load_genbank(file_path: Path) -> List[SequenceItem]:
        """加载GenBank格式文件"""
        sequences = []
        
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
        # 分割记录
        records = content.split('//\n')
        
        for record in records:
            if not record.strip():
                continue
                
            # 提取序列名称
            locus_match = re.search(r'LOCUS\s+(\w+)', record)
            if not locus_match:
                continue
                
            name = locus_match.group(1)
            
            # 提取序列
            origin_match = re.search(r'ORIGIN\s+(.*?)(?=\n//|\Z)', record, re.DOTALL)
            if not origin_match:
                continue
                
            sequence_text = origin_match.group(1)
            sequence = re.sub(r'[^ATCGN]', '', sequence_text.upper())
            
            sequences.append(SequenceItem(
                name=name,
                sequence=sequence,
                description=""
            ))
            
        return sequences
        
    @staticmethod
    def _save_fasta(sequences: List[SequenceItem], file_path: Path):
        """保存为FASTA格式"""
        with open(file_path, 'w', encoding='utf-8') as f:
            for seq in sequences:
                f.write(f">{seq.name}")
                if seq.description:
                    f.write(f" {seq.description}")
                f.write("\n")
                
                # 每行80个字符
                sequence = seq.sequence
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")
                    
    @staticmethod
    def _save_phylip(sequences: List[SequenceItem], file_path: Path):
        """保存为Phylip格式"""
        if not sequences:
            return
            
        # 计算最大序列长度
        max_length = max(seq.get_length() for seq in sequences)
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(f"{len(sequences)} {max_length}\n")
            
            for seq in sequences:
                # 截断或填充序列名称到10个字符
                name = seq.name[:10].ljust(10)
                sequence = seq.sequence.ljust(max_length, '-')
                f.write(f"{name}{sequence}\n")
                
    @staticmethod
    def _save_nexus(sequences: List[SequenceItem], file_path: Path):
        """保存为Nexus格式"""
        if not sequences:
            return
            
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write("#NEXUS\n\n")
            f.write("BEGIN DATA;\n")
            f.write(f"DIMENSIONS NTAX={len(sequences)} NCHAR={max(seq.get_length() for seq in sequences)};\n")
            f.write("FORMAT DATATYPE=DNA GAP=- MISSING=?;\n")
            f.write("MATRIX\n")
            
            for seq in sequences:
                f.write(f"{seq.name} {seq.sequence}\n")
                
            f.write(";\n")
            f.write("END;\n")
            
    @staticmethod
    def find_motifs(sequence: str, pattern: str) -> List[Tuple[int, int]]:
        """查找序列中的模式"""
        matches = []
        for match in re.finditer(pattern, sequence, re.IGNORECASE):
            matches.append((match.start(), match.end()))
        return matches
        
    @staticmethod
    def calculate_consensus(sequences: List[SequenceItem], threshold: float = 0.5) -> SequenceItem:
        """计算一致性序列"""
        if not sequences:
            return SequenceItem("consensus", "", "Consensus sequence")
            
        # 找到最大长度
        max_length = max(seq.get_length() for seq in sequences)
        
        consensus = ""
        for pos in range(max_length):
            # 统计每个位置的碱基
            base_counts = {}
            total_count = 0
            
            for seq in sequences:
                if pos < len(seq.sequence):
                    base = seq.sequence[pos].upper()
                    if base != '-':
                        base_counts[base] = base_counts.get(base, 0) + 1
                        total_count += 1
                        
            if total_count == 0:
                consensus += '-'
            else:
                # 找到最频繁的碱基
                most_common_base = max(base_counts.items(), key=lambda x: x[1])
                if most_common_base[1] / total_count >= threshold:
                    consensus += most_common_base[0]
                else:
                    consensus += 'N'  # 模糊碱基
                    
        return SequenceItem(
            name="consensus",
            sequence=consensus,
            description=f"Consensus sequence (threshold: {threshold})"
        )
        
    @staticmethod
    def calculate_sequence_identity(seq1: str, seq2: str) -> float:
        """计算序列一致性"""
        if len(seq1) != len(seq2):
            return 0.0
            
        matches = sum(1 for a, b in zip(seq1, seq2) if a.upper() == b.upper())
        return matches / len(seq1) * 100
        
    @staticmethod
    def find_orfs(sequence: str, min_length: int = 100) -> List[Tuple[int, int, str]]:
        """查找开放阅读框"""
        if len(sequence) < min_length:
            return []
            
        orfs = []
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        # 三个阅读框
        for frame in range(3):
            for i in range(frame, len(sequence) - 2, 3):
                codon = sequence[i:i+3].upper()
                if codon in start_codons:
                    # 寻找终止密码子
                    for j in range(i + 3, len(sequence) - 2, 3):
                        stop_codon = sequence[j:j+3].upper()
                        if stop_codon in stop_codons:
                            orf_length = j - i + 3
                            if orf_length >= min_length:
                                orf_sequence = sequence[i:j+3]
                                orfs.append((i, j + 3, orf_sequence))
                            break
                            
        return orfs
