#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
力导向布局引擎

实现用于注释邻接图可视化的力导向布局算法。
"""

import math
import logging
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class NodePosition:
    """节点位置信息"""
    x: float
    y: float
    vx: float = 0.0  # 速度 X
    vy: float = 0.0  # 速度 Y


class ForceDirectedLayout:
    """
    力导向布局引擎
    
    使用力导向算法布局图节点，包括：
    - 库仑斥力（节点间排斥）
    - 胡克引力（连接节点吸引）
    - 中心引力（防止漂移）
    - 阻尼力（减少震荡）
    """
    
    def __init__(
        self,
        repulsion_strength: float = 5000.0,
        spring_strength: float = 0.001,
        spring_length: float = 100.0,
        center_strength: float = 0.01,
        damping: float = 0.85,
        max_iterations: int = 500,
        convergence_threshold: float = 0.1
    ):
        """
        初始化布局引擎
        
        Args:
            repulsion_strength: 斥力强度
            spring_strength: 弹簧引力强度
            spring_length: 理想弹簧长度
            center_strength: 中心引力强度
            damping: 阻尼系数（0-1，越小阻尼越大）
            max_iterations: 最大迭代次数
            convergence_threshold: 收敛阈值
        """
        self.repulsion_strength = repulsion_strength
        self.spring_strength = spring_strength
        self.spring_length = spring_length
        self.center_strength = center_strength
        self.damping = damping
        self.max_iterations = max_iterations
        self.convergence_threshold = convergence_threshold
        
        self.positions: Dict[str, NodePosition] = {}
        self.node_order: List[str] = []
        self.iteration = 0
        self.converged = False
        
    def initialize(
        self, 
        node_names: List[str],
        node_weights: Optional[Dict[str, float]] = None
    ):
        """
        初始化节点位置
        
        Args:
            node_names: 节点名称列表
            node_weights: 节点权重（用于调整初始布局）
        """
        self.positions = {}
        self.node_order = node_names
        self.iteration = 0
        self.converged = False
        
        n = len(node_names)
        if n == 0:
            return
        
        # 圆形初始布局
        radius = max(100, 50 * math.sqrt(n))
        
        for i, name in enumerate(node_names):
            angle = 2 * math.pi * i / n
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            self.positions[name] = NodePosition(x, y)
        
        logger.debug(f"初始化布局: {n} 个节点，半径 {radius}")
    
    def set_positions(self, positions: Dict[str, Tuple[float, float]]):
        """
        手动设置节点位置
        
        Args:
            positions: {node_name: (x, y)} 字典
        """
        for name, (x, y) in positions.items():
            self.positions[name] = NodePosition(x, y)
    
    def calculate_forces(self, edges: List[Tuple[str, str, int]]):
        """
        计算所有节点受到的力
        
        Args:
            edges: 边列表 [(source, target, count), ...]
        """
        # 重置所有节点的力（通过速度累加）
        for name in self.positions:
            self.positions[name].vx = 0
            self.positions[name].vy = 0
        
        nodes = list(self.positions.keys())
        n = len(nodes)
        
        if n < 2:
            return
        
        # 1. 计算库仑斥力（节点间排斥）
        for i in range(n):
            for j in range(i + 1, n):
                name1 = nodes[i]
                name2 = nodes[j]
                
                pos1 = self.positions[name1]
                pos2 = self.positions[name2]
                
                dx = pos1.x - pos2.x
                dy = pos1.y - pos2.y
                distance = math.sqrt(dx * dx + dy * dy)
                distance = max(0.1, distance)  # 避免除零
                
                # 斥力公式: F = k / d^2
                force = self.repulsion_strength / (distance * distance)
                
                # 斥力方向
                fx = force * dx / distance
                fy = force * dy / distance
                
                pos1.vx += fx
                pos1.vy += fy
                pos2.vx -= fx
                pos2.vy -= fy
        
        # 2. 计算胡克引力（连接节点吸引）
        for source, target, count in edges:
            if source not in self.positions or target not in self.positions:
                continue
            
            pos1 = self.positions[source]
            pos2 = self.positions[target]
            
            dx = pos1.x - pos2.x
            dy = pos1.y - pos2.y
            distance = math.sqrt(dx * dx + dy * dy)
            distance = max(0.1, distance)
            
            # 根据边的权重调整引力
            strength = self.spring_strength * (1 + count * 0.1)
            
            # 引力公式: F = k * (d - L)
            displacement = distance - self.spring_length
            force = strength * displacement
            
            # 引力方向
            fx = force * dx / distance
            fy = force * dy / distance
            
            pos1.vx -= fx
            pos1.vy -= fy
            pos2.vx += fx
            pos2.vy += fy
        
        # 3. 计算中心引力（防止漂移）
        if n > 0:
            center_x = sum(pos.x for pos in self.positions.values()) / n
            center_y = sum(pos.y for pos in self.positions.values()) / n
            
            for name, pos in self.positions.items():
                dx = pos.x - center_x
                dy = pos.y - center_y
                distance = math.sqrt(dx * dx + dy * dy)
                
                if distance > 0:
                    force = self.center_strength * distance
                    fx = force * dx / distance
                    fy = force * dy / distance
                    
                    pos.vx -= fx
                    pos.vy -= fy
    
    def update_positions(self) -> float:
        """
        根据力更新节点位置
        
        Returns:
            最大位移距离
        """
        max_displacement = 0.0
        
        for name, pos in self.positions.items():
            # 应用阻尼
            pos.vx *= self.damping
            pos.vy *= self.damping
            
            # 限制最大速度
            speed = math.sqrt(pos.vx * pos.vx + pos.vy * pos.vy)
            max_speed = 10.0
            if speed > max_speed:
                pos.vx = pos.vx / speed * max_speed
                pos.vy = pos.vy / speed * max_speed
            
            # 更新位置
            displacement = math.sqrt(pos.vx * pos.vx + pos.vy * pos.vy)
            pos.x += pos.vx
            pos.y += pos.vy
            
            max_displacement = max(max_displacement, displacement)
        
        return max_displacement
    
    def step(self, edges: List[Tuple[str, str, int]]) -> float:
        """
        执行一步布局迭代
        
        Args:
            edges: 边列表
            
        Returns:
            最大位移距离
        """
        if self.converged or self.iteration >= self.max_iterations:
            return 0.0
        
        self.calculate_forces(edges)
        max_displacement = self.update_positions()
        self.iteration += 1
        
        # 检查收敛
        if max_displacement < self.convergence_threshold and self.iteration > 100:
            self.converged = True
            logger.debug(f"布局收敛于第 {self.iteration} 次迭代")
        
        return max_displacement
    
    def layout(
        self, 
        node_names: List[str], 
        edges: List[Tuple[str, str, int]]
    ) -> Dict[str, NodePosition]:
        """
        执行完整的布局计算
        
        Args:
            node_names: 节点名称列表
            edges: 边列表
            
        Returns:
            {node_name: NodePosition} 字典
        """
        self.initialize(node_names)
        
        while not self.converged and self.iteration < self.max_iterations:
            self.step(edges)
        
        logger.info(f"布局完成: {self.iteration} 次迭代, "
                   f"{'已收敛' if self.converged else '达到最大迭代次数'}")
        
        return self.positions
    
    def get_center(self) -> Tuple[float, float]:
        """
        获取当前布局的中心点
        
        Returns:
            (center_x, center_y)
        """
        if not self.positions:
            return (0.0, 0.0)
        
        n = len(self.positions)
        center_x = sum(pos.x for pos in self.positions.values()) / n
        center_y = sum(pos.y for pos in self.positions.values()) / n
        
        return (center_x, center_y)
    
    def center_layout(self):
        """将布局中心移动到原点"""
        center_x, center_y = self.get_center()
        
        for pos in self.positions.values():
            pos.x -= center_x
            pos.y -= center_y
    
    def scale_layout(self, scale: float):
        """
        缩放布局
        
        Args:
            scale: 缩放因子
        """
        for pos in self.positions.values():
            pos.x *= scale
            pos.y *= scale
    
    def get_bounds(self) -> Tuple[float, float, float, float]:
        """
        获取布局边界
        
        Returns:
            (min_x, min_y, max_x, max_y)
        """
        if not self.positions:
            return (0.0, 0.0, 0.0, 0.0)
        
        xs = [pos.x for pos in self.positions.values()]
        ys = [pos.y for pos in self.positions.values()]
        
        return (min(xs), min(ys), max(xs), max(ys))