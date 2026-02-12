#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
注释标准化规则引擎

用于将不同的基因注释名称标准化为统一的名称，
支持同义词规则和Entrez-like查询规则。
"""

import re
import json
import logging
from typing import Dict, List, Optional, Callable
from pathlib import Path

logger = logging.getLogger(__name__)


class AnnotationNormalizer:
    """
    注释标准化器
    
    支持多种规则类型：
    1. 同义词规则: alias = canonical
    2. 前缀/后缀规则: prefix* = canonical 或 *suffix = canonical
    3. 模式匹配规则: regex = canonical
    4. Entrez-like规则: [type]'gene' AND [name]'name'
    """
    
    def __init__(self):
        """初始化标准化器"""
        self.synonym_rules: Dict[str, str] = {}      # {alias: canonical}
        self.prefix_rules: List[Tuple[str, str]] = []  # [(prefix, canonical)]
        self.suffix_rules: List[Tuple[str, str]] = []  # [(suffix, canonical)]
        self.regex_rules: List[Tuple[re.Pattern, str]] = []  # [(pattern, canonical)]
        self.entrez_rules: List[str] = []             # Entrez-like规则
        self.custom_handlers: List[Callable] = []     # 自定义处理函数
        
        # 用于记录映射关系
        self.mapping_history: Dict[str, str] = {}
    
    def add_synonym_rule(self, alias: str, canonical: str):
        """
        添加同义词规则
        
        Args:
            alias: 别名
            canonical: 标准名称
        """
        self.synonym_rules[alias] = canonical
        logger.debug(f"添加同义词规则: {alias} -> {canonical}")
    
    def add_prefix_rule(self, prefix: str, canonical: str):
        """
        添加前缀规则
        
        Args:
            prefix: 前缀（会自动添加*匹配符）
            canonical: 标准名称
        """
        self.prefix_rules.append((prefix, canonical))
        logger.debug(f"添加前缀规则: {prefix}* -> {canonical}")
    
    def add_suffix_rule(self, suffix: str, canonical: str):
        """
        添加后缀规则
        
        Args:
            suffix: 后缀（会自动添加*匹配符）
            canonical: 标准名称
        """
        self.suffix_rules.append((suffix, canonical))
        logger.debug(f"添加后缀规则: *{suffix} -> {canonical}")
    
    def add_regex_rule(self, pattern: str, canonical: str):
        """
        添加正则表达式规则
        
        Args:
            pattern: 正则表达式
            canonical: 标准名称
        """
        try:
            compiled_pattern = re.compile(pattern, re.IGNORECASE)
            self.regex_rules.append((compiled_pattern, canonical))
            logger.debug(f"添加正则规则: {pattern} -> {canonical}")
        except re.error as e:
            logger.error(f"无效的正则表达式: {pattern} - {e}")
    
    def add_entrez_rule(self, rule: str):
        """
        添加Entrez-like规则
        
        Args:
            rule: 规则字符串，例如: ([type]'gene' AND [name]'cox1')
        """
        self.entrez_rules.append(rule)
        logger.debug(f"添加Entrez规则: {rule}")
    
    def add_custom_handler(self, handler: Callable[[str], Optional[str]]):
        """
        添加自定义处理函数
        
        Args:
            handler: 处理函数，接收原始名称，返回标准化名称或None
        """
        self.custom_handlers.append(handler)
    
    def normalize(self, annotation: str) -> str:
        """
        标准化注释名称
        
        应用所有规则，按优先级顺序：
        1. 精确同义词匹配
        2. 前缀规则
        3. 后缀规则
        4. 正则表达式规则
        5. Entrez-like规则
        6. 自定义处理函数
        
        Args:
            annotation: 原始注释名称
            
        Returns:
            标准化后的名称
        """
        if not annotation:
            return annotation
        
        # 1. 精确同义词匹配
        if annotation in self.synonym_rules:
            canonical = self.synonym_rules[annotation]
            self._record_mapping(annotation, canonical)
            return canonical
        
        # 2. 前缀规则
        for prefix, canonical in self.prefix_rules:
            if annotation.startswith(prefix):
                self._record_mapping(annotation, canonical)
                return canonical
        
        # 3. 后缀规则
        for suffix, canonical in self.suffix_rules:
            if annotation.endswith(suffix):
                self._record_mapping(annotation, canonical)
                return canonical
        
        # 4. 正则表达式规则
        for pattern, canonical in self.regex_rules:
            if pattern.search(annotation):
                self._record_mapping(annotation, canonical)
                return canonical
        
        # 5. Entrez-like规则
        for rule in self.entrez_rules:
            if self._match_entrez_rule(annotation, rule):
                canonical = self._extract_from_entrez_rule(annotation, rule)
                if canonical:
                    self._record_mapping(annotation, canonical)
                    return canonical
        
        # 6. 自定义处理函数
        for handler in self.custom_handlers:
            try:
                canonical = handler(annotation)
                if canonical:
                    self._record_mapping(annotation, canonical)
                    return canonical
            except Exception as e:
                logger.warning(f"自定义处理函数执行失败: {e}")
        
        # 没有匹配的规则，返回原始名称
        return annotation
    
    def _record_mapping(self, original: str, canonical: str):
        """记录映射关系"""
        self.mapping_history[original] = canonical
    
    def _match_entrez_rule(self, annotation: str, rule: str) -> bool:
        """
        匹配Entrez-like规则
        
        Args:
            annotation: 注释名称
            rule: 规则字符串
            
        Returns:
            是否匹配
        """
        # 简化的Entrez规则匹配
        # 处理NOT操作
        if " NOT " in rule:
            parts = rule.split(" NOT ")
            positive_part = parts[0]
            negative_part = parts[1]
            
            return (self._match_simple_condition(annotation, positive_part) and 
                   not self._match_simple_condition(annotation, negative_part))
        
        # 处理AND操作
        if " AND " in rule:
            parts = rule.split(" AND ")
            return all(self._match_simple_condition(annotation, part) for part in parts)
        
        # 处理OR操作
        if " OR " in rule:
            parts = rule.split(" OR ")
            return any(self._match_simple_condition(annotation, part) for part in parts)
        
        # 简单条件
        return self._match_simple_condition(annotation, rule)
    
    def _match_simple_condition(self, annotation: str, condition: str) -> bool:
        """
        匹配简单条件
        
        Args:
            annotation: 注释名称
            condition: 条件字符串
            
        Returns:
            是否匹配
        """
        # 移除[type]和[name]前缀
        condition = re.sub(r'\[(type|name)\]', '', condition).strip()
        
        # 处理引号包围的字符串
        if condition.startswith("'") and condition.endswith("'"):
            condition_value = condition[1:-1]
            return condition_value.lower() in annotation.lower()
        
        # 处理无引号的字符串
        return condition.lower() in annotation.lower()
    
    def _extract_from_entrez_rule(self, annotation: str, rule: str) -> Optional[str]:
        """
        从Entrez规则中提取标准名称
        
        Args:
            annotation: 注释名称
            rule: 规则字符串
            
        Returns:
            标准名称或None
        """
        # 简化实现：提取[name]'value'中的value
        name_match = re.search(r"\[name\]'([^']+)'", rule)
        if name_match:
            return name_match.group(1)
        
        return None
    
    def get_mapping(self) -> Dict[str, str]:
        """
        获取所有映射关系
        
        Returns:
            {original: canonical} 字典
        """
        return self.mapping_history.copy()
    
    def load_from_file(self, file_path: str):
        """
        从文件加载规则
        
        支持JSON格式:
        {
            "synonyms": {"alias1": "canonical1", ...},
            "prefixes": [["prefix1", "canonical1"], ...],
            "suffixes": [["suffix1", "canonical1"], ...],
            "regexes": [["pattern1", "canonical1"], ...],
            "entrez": ["rule1", "rule2", ...]
        }
        
        Args:
            file_path: 规则文件路径
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                rules = json.load(f)
            
            if "synonyms" in rules:
                for alias, canonical in rules["synonyms"].items():
                    self.add_synonym_rule(alias, canonical)
            
            if "prefixes" in rules:
                for prefix, canonical in rules["prefixes"]:
                    self.add_prefix_rule(prefix, canonical)
            
            if "suffixes" in rules:
                for suffix, canonical in rules["suffixes"]:
                    self.add_suffix_rule(suffix, canonical)
            
            if "regexes" in rules:
                for pattern, canonical in rules["regexes"]:
                    self.add_regex_rule(pattern, canonical)
            
            if "entrez" in rules:
                for rule in rules["entrez"]:
                    self.add_entrez_rule(rule)
            
            logger.info(f"从 {file_path} 加载了 {len(self.mapping_history)} 条规则")
            
        except Exception as e:
            logger.error(f"加载规则文件失败: {e}")
            raise
    
    def save_to_file(self, file_path: str):
        """
        保存规则到文件
        
        Args:
            file_path: 规则文件路径
        """
        try:
            rules = {
                "synonyms": self.synonym_rules,
                "prefixes": self.prefix_rules,
                "suffixes": self.suffix_rules,
                "regexes": [[pattern.pattern, canonical] for pattern, canonical in self.regex_rules],
                "entrez": self.entrez_rules
            }
            
            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(rules, f, indent=2, ensure_ascii=False)
            
            logger.info(f"规则已保存到 {file_path}")
            
        except Exception as e:
            logger.error(f"保存规则文件失败: {e}")
            raise
    
    def clear_rules(self):
        """清空所有规则"""
        self.synonym_rules.clear()
        self.prefix_rules.clear()
        self.suffix_rules.clear()
        self.regex_rules.clear()
        self.entrez_rules.clear()
        self.custom_handlers.clear()
        self.mapping_history.clear()
        logger.info("所有规则已清空")
    
    def get_rule_count(self) -> int:
        """
        获取规则总数
        
        Returns:
            规则数量
        """
        return (len(self.synonym_rules) + len(self.prefix_rules) + 
                len(self.suffix_rules) + len(self.regex_rules) + 
                len(self.entrez_rules) + len(self.custom_handlers))