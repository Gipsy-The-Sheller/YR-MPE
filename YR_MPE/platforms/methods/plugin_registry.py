"""
Plugin Registry Module
实现插件注册表和注册装饰器
"""
from typing import Dict, Type, List
from .plugin_interface import BasePlugin


class PluginRegistry:
    """插件注册表 - 单例模式"""
    
    _instance = None
    _plugins: Dict[str, Type[BasePlugin]] = {}
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    @classmethod
    def register(cls, plugin_class: Type[BasePlugin]):
        """注册插件"""
        name = plugin_class.get_tool_name()
        cls._plugins[name] = plugin_class
    
    @classmethod
    def get_plugin(cls, name: str) -> Type[BasePlugin]:
        """获取插件类"""
        return cls._plugins.get(name)
    
    @classmethod
    def get_all_plugins(cls) -> Dict[str, Type[BasePlugin]]:
        """获取所有插件"""
        return cls._plugins.copy()
    
    @classmethod
    def get_plugins_by_operation(cls, operation: str) -> List[str]:
        """根据操作类型获取插件列表"""
        result = []
        for name, plugin_class in cls._plugins.items():
            if operation in plugin_class.get_supported_operations():
                result.append(name)
        return result


def register_plugin(cls):
    """插件注册装饰器"""
    PluginRegistry.register(cls)
    return cls