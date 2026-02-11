"""
SeqDBG - 序列数据库可视化工具

用于分析多基因组的注释邻接关系，发现结构变异和注释不统一。
"""

__version__ = "1.0.0"
__author__ = "YR-MPE Team"

from .models import AnnotationNode, AnnotationGraph, GeneInfo, GeneStats
from .parser import GeneBankParser
from .normalizer import AnnotationNormalizer
from .engine import SeqDBGraphBuilder

__all__ = [
    "AnnotationNode",
    "AnnotationGraph",
    "GeneInfo",
    "GeneStats",
    "GeneBankParser",
    "AnnotationNormalizer",
    "SeqDBGraphBuilder",
]