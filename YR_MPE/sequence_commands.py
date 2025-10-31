#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Sequence Commands for YR-MPE Sequence Editor
Undo/Redo commands for sequence operations
"""

from typing import List, Optional
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

try:
    from .sequence_models import SequenceItem, SequenceType
except ImportError:
    from sequence_models import SequenceItem, SequenceType


class SequenceCommand(QUndoCommand):
    """序列命令基类"""
    
    def __init__(self, description: str = ""):
        super().__init__(description)
        self.timestamp = QDateTime.currentDateTime()
        
    def get_timestamp(self) -> QDateTime:
        """获取时间戳"""
        return self.timestamp


class InsertSequenceCommand(SequenceCommand):
    """插入序列命令"""
    
    def __init__(self, sequences: List[SequenceItem], position: int, 
                 model, description: str = "Insert sequence"):
        super().__init__(description)
        self.sequences = sequences
        self.position = position
        self.model = model
        
    def redo(self):
        """重做操作"""
        for i, seq in enumerate(self.sequences):
            self.model.add_sequence_at_position(seq, self.position + i)
            
    def undo(self):
        """撤销操作"""
        for i in range(len(self.sequences)):
            self.model.remove_sequence_at_position(self.position)


class DeleteSequenceCommand(SequenceCommand):
    """删除序列命令"""
    
    def __init__(self, position: int, model, description: str = "Delete sequence"):
        super().__init__(description)
        self.position = position
        self.model = model
        self.deleted_sequence = None
        
    def redo(self):
        """重做操作"""
        self.deleted_sequence = self.model.get_sequence(self.position)
        self.model.remove_sequence_at_position(self.position)
        
    def undo(self):
        """撤销操作"""
        if self.deleted_sequence:
            self.model.insert_sequence_at_position(self.deleted_sequence, self.position)


class EditSequenceCommand(SequenceCommand):
    """编辑序列命令"""
    
    def __init__(self, position: int, old_sequence: SequenceItem, 
                 new_sequence: SequenceItem, model, description: str = "Edit sequence"):
        super().__init__(description)
        self.position = position
        self.old_sequence = old_sequence
        self.new_sequence = new_sequence
        self.model = model
        
    def redo(self):
        """重做操作"""
        self.model.set_sequence_at_position(self.new_sequence, self.position)
        
    def undo(self):
        """撤销操作"""
        self.model.set_sequence_at_position(self.old_sequence, self.position)


class EditSequenceTextCommand(SequenceCommand):
    """编辑序列文本命令"""
    
    def __init__(self, sequence_index: int, position: int, old_text: str, 
                 new_text: str, model, description: str = "Edit sequence text"):
        super().__init__(description)
        self.sequence_index = sequence_index
        self.position = position
        self.old_text = old_text
        self.new_text = new_text
        self.model = model
        
    def redo(self):
        """重做操作"""
        sequence = self.model.get_sequence(self.sequence_index)
        if sequence:
            old_sequence = sequence.sequence
            new_sequence = (old_sequence[:self.position] + 
                          self.new_text + 
                          old_sequence[self.position + len(self.old_text):])
            sequence.sequence = new_sequence
            self.model.sequence_changed.emit(self.sequence_index)
            
    def undo(self):
        """撤销操作"""
        sequence = self.model.get_sequence(self.sequence_index)
        if sequence:
            old_sequence = sequence.sequence
            new_sequence = (old_sequence[:self.position] + 
                          self.old_text + 
                          old_sequence[self.position + len(self.new_text):])
            sequence.sequence = new_sequence
            self.model.sequence_changed.emit(self.sequence_index)


class ReverseSequenceCommand(SequenceCommand):
    """反向序列命令"""
    
    def __init__(self, sequence_index: int, model, description: str = "Reverse sequence"):
        super().__init__(description)
        self.sequence_index = sequence_index
        self.model = model
        self.original_sequence = None
        
    def redo(self):
        """重做操作"""
        sequence = self.model.get_sequence(self.sequence_index)
        if sequence:
            self.original_sequence = sequence.sequence
            sequence.sequence = sequence.sequence[::-1]
            self.model.sequence_changed.emit(self.sequence_index)
            
    def undo(self):
        """撤销操作"""
        if self.original_sequence:
            sequence = self.model.get_sequence(self.sequence_index)
            if sequence:
                sequence.sequence = self.original_sequence
                self.model.sequence_changed.emit(self.sequence_index)


class ComplementSequenceCommand(SequenceCommand):
    """互补序列命令"""
    
    def __init__(self, sequence_index: int, model, description: str = "Complement sequence"):
        super().__init__(description)
        self.sequence_index = sequence_index
        self.model = model
        self.original_sequence = None
        
    def redo(self):
        """重做操作"""
        sequence = self.model.get_sequence(self.sequence_index)
        if sequence and sequence.sequence_type in [SequenceType.DNA, SequenceType.RNA]:
            self.original_sequence = sequence.sequence
            sequence.sequence = self._complement(sequence.sequence, sequence.sequence_type)
            self.model.sequence_changed.emit(self.sequence_index)
            
    def undo(self):
        """撤销操作"""
        if self.original_sequence:
            sequence = self.model.get_sequence(self.sequence_index)
            if sequence:
                sequence.sequence = self.original_sequence
                self.model.sequence_changed.emit(self.sequence_index)
                
    def _complement(self, sequence: str, sequence_type: SequenceType) -> str:
        """计算互补序列"""
        if sequence_type == SequenceType.DNA:
            complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        elif sequence_type == SequenceType.RNA:
            complement_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        else:
            return sequence
            
        return ''.join(complement_map.get(c.upper(), c) for c in sequence)


class ReverseComplementSequenceCommand(SequenceCommand):
    """反向互补序列命令"""
    
    def __init__(self, sequence_index: int, model, description: str = "Reverse complement sequence"):
        super().__init__(description)
        self.sequence_index = sequence_index
        self.model = model
        self.original_sequence = None
        
    def redo(self):
        """重做操作"""
        sequence = self.model.get_sequence(self.sequence_index)
        if sequence and sequence.sequence_type in [SequenceType.DNA, SequenceType.RNA]:
            self.original_sequence = sequence.sequence
            sequence.sequence = self._reverse_complement(sequence.sequence, sequence.sequence_type)
            self.model.sequence_changed.emit(self.sequence_index)
            
    def undo(self):
        """撤销操作"""
        if self.original_sequence:
            sequence = self.model.get_sequence(self.sequence_index)
            if sequence:
                sequence.sequence = self.original_sequence
                self.model.sequence_changed.emit(self.sequence_index)
                
    def _reverse_complement(self, sequence: str, sequence_type: SequenceType) -> str:
        """计算反向互补序列"""
        if sequence_type == SequenceType.DNA:
            complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        elif sequence_type == SequenceType.RNA:
            complement_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        else:
            return sequence
            
        complemented = ''.join(complement_map.get(c.upper(), c) for c in sequence)
        return complemented[::-1]


class TranslateSequenceCommand(SequenceCommand):
    """翻译序列命令"""
    
    def __init__(self, sequence_index: int, model, description: str = "Translate sequence"):
        super().__init__(description)
        self.sequence_index = sequence_index
        self.model = model
        self.original_sequence = None
        self.translated_sequence = None
        
    def redo(self):
        """重做操作"""
        sequence = self.model.get_sequence(self.sequence_index)
        if sequence and sequence.sequence_type == SequenceType.DNA:
            self.original_sequence = sequence.sequence
            self.translated_sequence = self._translate(sequence.sequence)
            sequence.sequence = self.translated_sequence
            sequence.sequence_type = SequenceType.PROTEIN
            self.model.sequence_changed.emit(self.sequence_index)
            
    def undo(self):
        """撤销操作"""
        if self.original_sequence:
            sequence = self.model.get_sequence(self.sequence_index)
            if sequence:
                sequence.sequence = self.original_sequence
                sequence.sequence_type = SequenceType.DNA
                self.model.sequence_changed.emit(self.sequence_index)
                
    def _translate(self, sequence: str) -> str:
        """翻译DNA序列为蛋白质"""
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
        seq = sequence.upper().replace('-', '')
        if len(seq) % 3 != 0:
            seq = seq[:-(len(seq) % 3)]
            
        protein = ''
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            protein += codon_table.get(codon, 'X')
            
        return protein


class BatchOperationCommand(SequenceCommand):
    """批量操作命令"""
    
    def __init__(self, operations: List[SequenceCommand], description: str = "Batch operation"):
        super().__init__(description)
        self.operations = operations
        
    def redo(self):
        """重做操作"""
        for operation in self.operations:
            operation.redo()
            
    def undo(self):
        """撤销操作"""
        for operation in reversed(self.operations):
            operation.undo()


class FindReplaceCommand(SequenceCommand):
    """查找替换命令"""
    
    def __init__(self, sequence_index: int, find_text: str, replace_text: str, 
                 model, description: str = "Find and replace"):
        super().__init__(description)
        self.sequence_index = sequence_index
        self.find_text = find_text
        self.replace_text = replace_text
        self.model = model
        self.original_sequence = None
        self.replacements = []
        
    def redo(self):
        """重做操作"""
        sequence = self.model.get_sequence(self.sequence_index)
        if sequence:
            self.original_sequence = sequence.sequence
            sequence.sequence = sequence.sequence.replace(self.find_text, self.replace_text)
            self.model.sequence_changed.emit(self.sequence_index)
            
    def undo(self):
        """撤销操作"""
        if self.original_sequence:
            sequence = self.model.get_sequence(self.sequence_index)
            if sequence:
                sequence.sequence = self.original_sequence
                self.model.sequence_changed.emit(self.sequence_index)
