"""
Plugin Interface Module
定义插件基类和抽象接口
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
from pathlib import Path


class BasePlugin(ABC):
    """插件基类接口"""
    
    @classmethod
    @abstractmethod
    def get_tool_name(cls) -> str:
        """返回工具名称"""
        pass
    
    @classmethod
    @abstractmethod
    def get_supported_operations(cls) -> list:
        """返回支持的操作类型列表"""
        pass
    
    @abstractmethod
    def execute(self, 
                input_data: Dict[str, Any], 
                workspace_path: Optional[Path] = None,
                **kwargs) -> Dict[str, Any]:
        """
        执行插件操作
        
        Args:
            input_data: 输入数据字典
            workspace_path: 工作区路径（可选）
            **kwargs: 其他参数
            
        Returns:
            result: 包含结果文件路径、元数据等的字典
        """
        pass
    
    @abstractmethod
    def validate_input(self, input_data: Dict[str, Any]) -> bool:
        """验证输入数据"""
        pass