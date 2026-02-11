#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SeqDBG解析后台线程

将GenBank解析和图构建移到后台线程，避免阻塞主线程UI。
支持进度回调显示，提供更好的用户体验。
"""

import logging
from typing import List, Optional, Dict, Callable
from pathlib import Path

from PyQt5.QtCore import QThread, pyqtSignal

from .models import AnnotationGraph, GeneStats
from .engine import SeqDBGraphBuilder

logger = logging.getLogger(__name__)


class SeqDBGParserThread(QThread):
    """
    GenBank解析后台线程
    
    在后台线程中执行GenBank文件解析和图构建操作，
    通过信号机制与主线程UI通信。
    """
    
    # 信号定义 - PyQt5 信号定义支持 Python 基本类型
    progress = pyqtSignal(str, int, int)           # 进度信号 (消息, 当前值, 总值)
    file_progress = pyqtSignal(str, int, int)       # 文件处理进度 (文件名, 当前值, 总值)
    finished = pyqtSignal(object, list)            # 完成信号 (图, 统计信息列表)
    error = pyqtSignal(str)                         # 错误信号
    status = pyqtSignal(str)                        # 状态消息
    
    def __init__(
        self,
        builder: SeqDBGraphBuilder,
        file_paths: List[str],
        extract_sequences: bool = True,
        filter_entries: Optional[List[str]] = None
    ):
        """
        初始化解析线程
        
        Args:
            builder: SeqDBGraphBuilder实例
            file_paths: GenBank文件路径列表
            extract_sequences: 是否提取序列
            filter_entries: 可选，只统计这些entries中的基因
        """
        super().__init__()
        self.builder = builder
        self.file_paths = file_paths
        self.extract_sequences = extract_sequences
        self.filter_entries = filter_entries
        self.is_cancelled = False
        
        # 进度回调函数（用于内部进度更新）
        self._progress_callback: Optional[Callable[[str, int, int], None]] = None
        
    def run(self):
        """执行解析任务"""
        try:
            self.status.emit("开始解析GenBank文件...")
            logger.info(f"开始解析 {len(self.file_paths)} 个文件")
            
            # 设置进度回调
            self._progress_callback = self._emit_file_progress
            
            # 重置构建器
            self.builder.reset()
            
            # 构建多基因组图
            graph = self.builder.build_multi_genome_graph(
                self.file_paths,
                extract_sequences=self.extract_sequences,
                progress_callback=self._progress_callback
            )
            
            if self.is_cancelled:
                self.status.emit("解析已取消")
                return
            
            # 获取统计信息
            self.status.emit("计算统计信息...")
            stats = self.builder.get_gene_stats(filter_entries=self.filter_entries)
            
            if self.is_cancelled:
                self.status.emit("解析已取消")
                return
            
            # 发送完成信号
            self.finished.emit(graph, stats)
            self.status.emit("解析完成")
            
            logger.info(f"解析完成: {len(stats)} 个基因")
            
        except Exception as e:
            logger.error(f"解析失败: {e}")
            self.error.emit(f"解析失败: {str(e)}")
    
    def _emit_file_progress(self, message: str, current: int, total: int):
        """
        发送文件处理进度信号
        
        Args:
            message: 进度消息
            current: 当前进度值
            total: 总值
        """
        self.file_progress.emit(message, current, total)
    
    def cancel(self):
        """取消解析任务"""
        self.is_cancelled = True
        self.status.emit("正在取消解析...")
        logger.info("解析任务已请求取消")


class SeqDBGIncrementalParserThread(QThread):
    """
    增量解析线程
    
    用于添加新文件而不需要完全重建图结构。
    适用于已加载部分文件后添加更多文件的场景。
    """
    
    # 信号定义 - PyQt5 信号定义支持 Python 基本类型
    progress = pyqtSignal(str, int, int)
    file_progress = pyqtSignal(str, int, int)
    finished = pyqtSignal(object, list)  # (图, 统计信息列表)
    error = pyqtSignal(str)
    status = pyqtSignal(str)
    
    def __init__(
        self,
        builder: SeqDBGraphBuilder,
        new_files: List[str],
        extract_sequences: bool = True,
        filter_entries: Optional[List[str]] = None
    ):
        """
        初始化增量解析线程
        
        Args:
            builder: SeqDBGraphBuilder实例（包含已有图）
            new_files: 新增的GenBank文件路径列表
            extract_sequences: 是否提取序列
            filter_entries: 可选，只统计这些entries中的基因
        """
        super().__init__()
        self.builder = builder
        self.new_files = new_files
        self.extract_sequences = extract_sequences
        self.filter_entries = filter_entries
        self.is_cancelled = False
        
    def run(self):
        """执行增量解析任务"""
        try:
            self.status.emit(f"开始增量解析 {len(self.new_files)} 个文件...")
            logger.info(f"开始增量解析 {len(self.new_files)} 个文件")
            
            # 加载新文件（不重置现有图）
            total_files = len(self.new_files)
            for i, file_path in enumerate(self.new_files):
                if self.is_cancelled:
                    self.status.emit("解析已取消")
                    return
                
                self.file_progress.emit(f"正在解析: {Path(file_path).name}", i + 1, total_files)
                self.builder.load_file(file_path, extract_sequences=self.extract_sequences)
            
            # 获取统计信息
            self.status.emit("计算统计信息...")
            stats = self.builder.get_gene_stats(filter_entries=self.filter_entries)
            
            if self.is_cancelled:
                self.status.emit("解析已取消")
                return
            
            # 发送完成信号
            graph = self.builder.graph
            self.finished.emit(graph, stats)
            self.status.emit("增量解析完成")
            
            logger.info(f"增量解析完成: {len(stats)} 个基因")
            
        except Exception as e:
            logger.error(f"增量解析失败: {e}")
            self.error.emit(f"增量解析失败: {str(e)}")
    
    def cancel(self):
        """取消解析任务"""
        self.is_cancelled = True
        self.status.emit("正在取消解析...")
        logger.info("增量解析任务已请求取消")