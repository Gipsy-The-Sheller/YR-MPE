#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Sequence Highlighter for YR-MPE Sequence Editor
Syntax highlighting for DNA, RNA, and protein sequences
"""

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

try:
    from .sequence_models import SequenceType
except ImportError:
    from sequence_models import SequenceType


class SequenceHighlighter(QSyntaxHighlighter):
    """序列语法高亮器"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequence_type = SequenceType.UNKNOWN
        self.highlighting_rules = []
        self.setup_highlighting_rules()
        
    def set_sequence_type(self, sequence_type: SequenceType):
        """设置序列类型"""
        self.sequence_type = sequence_type
        self.setup_highlighting_rules()
        self.rehighlight()
        
    def setup_highlighting_rules(self):
        """设置高亮规则"""
        self.highlighting_rules = []
        
        if self.sequence_type == SequenceType.DNA:
            self.setup_dna_highlighting()
        elif self.sequence_type == SequenceType.RNA:
            self.setup_rna_highlighting()
        elif self.sequence_type == SequenceType.PROTEIN:
            self.setup_protein_highlighting()
        else:
            self.setup_default_highlighting()
            
    def setup_dna_highlighting(self):
        """设置DNA序列高亮规则"""
        # A - 红色
        format_a = QTextCharFormat()
        format_a.setForeground(QColor(255, 0, 0))
        format_a.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Aa]'), format_a))
        
        # T - 蓝色
        format_t = QTextCharFormat()
        format_t.setForeground(QColor(0, 0, 255))
        format_t.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Tt]'), format_t))
        
        # C - 绿色
        format_c = QTextCharFormat()
        format_c.setForeground(QColor(0, 128, 0))
        format_c.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Cc]'), format_c))
        
        # G - 橙色
        format_g = QTextCharFormat()
        format_g.setForeground(QColor(255, 165, 0))
        format_g.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Gg]'), format_g))
        
        # 模糊碱基 - 灰色
        format_ambiguous = QTextCharFormat()
        format_ambiguous.setForeground(QColor(128, 128, 128))
        format_ambiguous.setFontItalic(True)
        self.highlighting_rules.append((QRegExp(r'[NnRrYyWwSsKkMmBbDdHhVv]'), format_ambiguous))
        
        # 间隔符 - 浅灰色
        format_gap = QTextCharFormat()
        format_gap.setForeground(QColor(200, 200, 200))
        format_gap.setBackground(QColor(240, 240, 240))
        self.highlighting_rules.append((QRegExp(r'[-]'), format_gap))
        
    def setup_rna_highlighting(self):
        """设置RNA序列高亮规则"""
        # A - 红色
        format_a = QTextCharFormat()
        format_a.setForeground(QColor(255, 0, 0))
        format_a.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Aa]'), format_a))
        
        # U - 蓝色
        format_u = QTextCharFormat()
        format_u.setForeground(QColor(0, 0, 255))
        format_u.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Uu]'), format_u))
        
        # C - 绿色
        format_c = QTextCharFormat()
        format_c.setForeground(QColor(0, 128, 0))
        format_c.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Cc]'), format_c))
        
        # G - 橙色
        format_g = QTextCharFormat()
        format_g.setForeground(QColor(255, 165, 0))
        format_g.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[Gg]'), format_g))
        
        # 模糊碱基 - 灰色
        format_ambiguous = QTextCharFormat()
        format_ambiguous.setForeground(QColor(128, 128, 128))
        format_ambiguous.setFontItalic(True)
        self.highlighting_rules.append((QRegExp(r'[NnRrYyWwSsKkMmBbDdHhVv]'), format_ambiguous))
        
        # 间隔符 - 浅灰色
        format_gap = QTextCharFormat()
        format_gap.setForeground(QColor(200, 200, 200))
        format_gap.setBackground(QColor(240, 240, 240))
        self.highlighting_rules.append((QRegExp(r'[-]'), format_gap))
        
    def setup_protein_highlighting(self):
        """设置蛋白质序列高亮规则"""
        # 疏水性氨基酸 - 蓝色
        hydrophobic = ['A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W']
        format_hydrophobic = QTextCharFormat()
        format_hydrophobic.setForeground(QColor(0, 0, 255))
        format_hydrophobic.setFontWeight(QFont.Bold)
        for aa in hydrophobic:
            self.highlighting_rules.append((QRegExp(f'[{aa}{aa.lower()}]'), format_hydrophobic))
            
        # 极性氨基酸 - 绿色
        polar = ['S', 'T', 'N', 'Q']
        format_polar = QTextCharFormat()
        format_polar.setForeground(QColor(0, 128, 0))
        format_polar.setFontWeight(QFont.Bold)
        for aa in polar:
            self.highlighting_rules.append((QRegExp(f'[{aa}{aa.lower()}]'), format_polar))
            
        # 带正电荷氨基酸 - 红色
        positive = ['K', 'R', 'H']
        format_positive = QTextCharFormat()
        format_positive.setForeground(QColor(255, 0, 0))
        format_positive.setFontWeight(QFont.Bold)
        for aa in positive:
            self.highlighting_rules.append((QRegExp(f'[{aa}{aa.lower()}]'), format_positive))
            
        # 带负电荷氨基酸 - 紫色
        negative = ['D', 'E']
        format_negative = QTextCharFormat()
        format_negative.setForeground(QColor(128, 0, 128))
        format_negative.setFontWeight(QFont.Bold)
        for aa in negative:
            self.highlighting_rules.append((QRegExp(f'[{aa}{aa.lower()}]'), format_negative))
            
        # 特殊氨基酸 - 橙色
        special = ['C', 'G', 'P']
        format_special = QTextCharFormat()
        format_special.setForeground(QColor(255, 165, 0))
        format_special.setFontWeight(QFont.Bold)
        for aa in special:
            self.highlighting_rules.append((QRegExp(f'[{aa}{aa.lower()}]'), format_special))
            
        # 终止密码子 - 红色背景
        format_stop = QTextCharFormat()
        format_stop.setBackground(QColor(255, 200, 200))
        format_stop.setForeground(QColor(255, 0, 0))
        format_stop.setFontWeight(QFont.Bold)
        self.highlighting_rules.append((QRegExp(r'[*]'), format_stop))
        
        # 间隔符 - 浅灰色
        format_gap = QTextCharFormat()
        format_gap.setForeground(QColor(200, 200, 200))
        format_gap.setBackground(QColor(240, 240, 240))
        self.highlighting_rules.append((QRegExp(r'[-]'), format_gap))
        
    def setup_default_highlighting(self):
        """设置默认高亮规则"""
        # 数字 - 蓝色
        format_number = QTextCharFormat()
        format_number.setForeground(QColor(0, 0, 255))
        self.highlighting_rules.append((QRegExp(r'\d'), format_number))
        
        # 字母 - 黑色
        format_letter = QTextCharFormat()
        format_letter.setForeground(QColor(0, 0, 0))
        self.highlighting_rules.append((QRegExp(r'[A-Za-z]'), format_letter))
        
    def highlightBlock(self, text):
        """高亮文本块"""
        for pattern, format in self.highlighting_rules:
            expression = QRegExp(pattern)
            index = expression.indexIn(text)
            while index >= 0:
                length = expression.matchedLength()
                self.setFormat(index, length, format)
                index = expression.indexIn(text, index + length)


class SequenceTextEdit(QTextEdit):
    """序列文本编辑器"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.highlighter = SequenceHighlighter(self.document())
        self.sequence_type = SequenceType.UNKNOWN
        self.setup_editor()
        
    def setup_editor(self):
        """设置编辑器"""
        # 设置字体
        font = QFont("Consolas", 12)
        font.setFixedPitch(True)
        self.setFont(font)
        
        # 设置制表符宽度
        self.setTabStopWidth(40)
        
        # 设置行号
        self.setLineWrapMode(QTextEdit.NoWrap)
        
        # 设置自动完成
        self.setup_autocompletion()
        
    def setup_autocompletion(self):
        """设置自动完成"""
        # 这里可以实现自动完成功能
        pass
        
    def set_sequence_type(self, sequence_type: SequenceType):
        """设置序列类型"""
        self.sequence_type = sequence_type
        self.highlighter.set_sequence_type(sequence_type)
        
    def keyPressEvent(self, event):
        """按键事件处理"""
        # 限制输入字符
        if self.sequence_type == SequenceType.DNA:
            valid_chars = set('ATCGNMRWSYKVHDBNatcgnmrwsykvhdbn-')
        elif self.sequence_type == SequenceType.RNA:
            valid_chars = set('AUCGNMRWSYKVHDBNaucgnmrwsykvhdbn-')
        elif self.sequence_type == SequenceType.PROTEIN:
            valid_chars = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-*')
        else:
            valid_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-')
            
        if event.text() and event.text()[0] not in valid_chars:
            event.ignore()
            return
            
        super().keyPressEvent(event)
        
    def insertFromMimeData(self, source):
        """从剪贴板插入数据"""
        if source.hasText():
            text = source.text()
            # 过滤无效字符
            if self.sequence_type == SequenceType.DNA:
                valid_chars = set('ATCGNMRWSYKVHDBNatcgnmrwsykvhdbn-\n\r')
            elif self.sequence_type == SequenceType.RNA:
                valid_chars = set('AUCGNMRWSYKVHDBNaucgnmrwsykvhdbn-\n\r')
            elif self.sequence_type == SequenceType.PROTEIN:
                valid_chars = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-*\n\r')
            else:
                valid_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-\n\r')
                
            filtered_text = ''.join(c for c in text if c in valid_chars)
            if filtered_text != text:
                QMessageBox.warning(self, "Invalid Characters", 
                                  "Some characters were filtered out as they are not valid for this sequence type.")
                                  
            # 创建新的MimeData
            new_source = QMimeData()
            new_source.setText(filtered_text)
            super().insertFromMimeData(new_source)
        else:
            super().insertFromMimeData(source)
