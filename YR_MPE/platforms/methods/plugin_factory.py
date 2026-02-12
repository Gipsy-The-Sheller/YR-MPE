from typing import Dict, Any, Type
from .plugin_registry import PluginRegistry
from .plugin_config_manager import PluginConfigManager

class PluginFactory:
    """
    插件工厂 - 负责创建插件实例并注入必要的依赖
    """
    
    def __init__(self, registry: PluginRegistry):
        self.registry = registry
        self.config_manager = PluginConfigManager()
        
    def create_plugin(self, plugin_name: str, **kwargs) -> Any:
        """
        创建插件实例
        
        Args:
            plugin_name: 插件名称
            **kwargs: 传递给插件构造函数的参数
            
        Returns:
            插件实例
        """
        if plugin_name not in self.registry.plugins:
            raise ValueError(f"Plugin '{plugin_name}' not found in registry")
            
        plugin_class = self.registry.plugins[plugin_name]
        
        # 注入工具路径（如果插件需要）
        tool_name = self._get_tool_name_for_plugin(plugin_name)
        if tool_name:
            tool_path = self.config_manager.get_tool_path(tool_name)
            if tool_path:
                kwargs['tool_path'] = tool_path
                
        return plugin_class(**kwargs)
        
    def _get_tool_name_for_plugin(self, plugin_name: str) -> str:
        """
        根据插件名称获取对应的工具名称
        
        Args:
            plugin_name: 插件名称
            
        Returns:
            工具名称
        """
        plugin_tool_mapping = {
            'clustal_omega': 'Clustal Omega',
            'muscle5': 'Muscle5', 
            'mafft': 'MAFFT',
            'trimal': 'TrimAl',
            'gblocks': 'GBlocks',
            'model_finder': 'ModelFinder',
            'iqtree': 'IQ-TREE 3',
            'caster_site': 'CASTER-site'
        }
        
        return plugin_tool_mapping.get(plugin_name, '')