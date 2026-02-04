import json
import os
from typing import Dict, Any, Optional
from .file_operations import FileOperations

class PluginConfigManager:
    """
    插件配置管理器 - 负责读取和验证插件配置
    """
    
    def __init__(self, config_path: str = None):
        self.config_path = config_path or os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 
            "config.json"
        )
        self.file_ops = FileOperations()
        self.config_data = self._load_config()
        
    def _load_config(self) -> Dict[str, Any]:
        """加载配置文件"""
        if not os.path.exists(self.config_path):
            return {}
            
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"Failed to load config file: {e}")
            return {}
            
    def get_tool_path(self, tool_name: str) -> Optional[str]:
        """
        获取工具路径
        
        Args:
            tool_name: 工具名称
            
        Returns:
            工具完整路径，如果不存在则返回None
        """
        if not self.config_data:
            return None
            
        for tool in self.config_data:
            if tool.get("name") == tool_name:
                path = tool.get("path")
                if path:
                    # 处理相对路径
                    if path.startswith("/"):
                        path = path[1:]  # 移除开头的斜杠
                    full_path = os.path.join(
                        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                        path
                    )
                    if os.path.exists(full_path):
                        return full_path
                    else:
                        # 尝试直接使用配置中的路径
                        if os.path.exists(path):
                            return path
        return None
        
    def validate_tool_path(self, tool_name: str) -> bool:
        """验证工具路径是否存在"""
        tool_path = self.get_tool_path(tool_name)
        return tool_path is not None and os.path.exists(tool_path)
        
    def get_all_tool_paths(self) -> Dict[str, str]:
        """获取所有工具的路径映射"""
        tool_paths = {}
        if not self.config_data:
            return tool_paths
            
        for tool in self.config_data:
            tool_name = tool.get("name")
            if tool_name:
                tool_path = self.get_tool_path(tool_name)
                if tool_path:
                    tool_paths[tool_name] = tool_path
                    
        return tool_paths