from typing import Dict, Any, Optional, Type
from .plugin_registry import PluginRegistry
from .workspace_manager import WorkspaceManager
from .plugin_adapter import PluginAdapter
import logging

class PluginManager:
    """
    插件管理器 - 负责插件的生命周期管理、配置和调用协调
    """
    
    def __init__(self, workspace_manager: WorkspaceManager):
        self.workspace_manager = workspace_manager
        self.registry = PluginRegistry()
        self.logger = logging.getLogger(__name__)
        self.active_adapters: Dict[str, PluginAdapter] = {}

    def register_plugin(self, plugin_name: str, plugin_class: type) -> None:
        """
        注册单个插件

        Args:
            plugin_name: 插件名称
            plugin_class: 插件类
        """
        self.registry._plugins[plugin_name] = plugin_class

    def register_all_plugins(self):
        """注册所有已知的插件"""
        # 序列比对插件
        self._register_alignment_plugins()
        # 修剪插件  
        self._register_trimming_plugins()
        # 模型选择插件
        self._register_model_selection_plugins()
        # 距离计算插件
        self._register_distance_plugins()
        # 系统发育插件
        self._register_phylogeny_plugins()
        # 溯祖插件
        self._register_coalescent_plugins()
        # Dataset插件
        self._register_dataset_plugins()
        
    def _register_alignment_plugins(self):
        """注册序列比对插件"""
        try:
            from YR_MPE.plugins.clustal_omega_plugin import ClustalOmegaPluginEntry
            self.register_plugin("clustal_omega", ClustalOmegaPluginEntry)
        except ImportError:
            pass
            
        try:
            from YR_MPE.plugins.muscle5_plugin import Muscle5PluginEntry  
            self.register_plugin("muscle5", Muscle5PluginEntry)
        except ImportError:
            pass
            
        try:
            from YR_MPE.plugins.mafft_plugin import MAFFTPluginEntry
            self.register_plugin("mafft", MAFFTPluginEntry)
        except ImportError:
            pass
        
    def _register_trimming_plugins(self):
        """注册修剪插件"""
        try:
            from YR_MPE.plugins.trimal_plugin import TrimAlPluginEntry
            self.register_plugin("trimal", TrimAlPluginEntry)
        except ImportError:
            pass
            
        try:
            from YR_MPE.plugins.gblocks_plugin import GBlocksPluginEntry
            self.register_plugin("gblocks", GBlocksPluginEntry)
        except ImportError:
            pass
        
    def _register_model_selection_plugins(self):
        """注册模型选择插件"""
        try:
            from YR_MPE.plugins.model_finder_plugin import ModelFinderPluginEntry
            self.register_plugin("model_finder", ModelFinderPluginEntry)
        except ImportError:
            pass
        
    def _register_distance_plugins(self):
        """注册距离计算插件"""
        try:
            from YR_MPE.plugins.ml_distance_plugin import MLDistancePluginEntry
            self.register_plugin("ml_distance", MLDistancePluginEntry)
        except ImportError:
            pass
        
    def _register_phylogeny_plugins(self):
        """注册系统发育插件"""
        try:
            from YR_MPE.plugins.iqtree_plugin import IQTreePluginEntry
            self.register_plugin("iqtree", IQTreePluginEntry)
        except ImportError:
            pass
            
        try:
            from YR_MPE.plugins.icytree import IcyTreePluginEntry
            self.register_plugin("icytree", IcyTreePluginEntry)
        except ImportError:
            pass
        
    def _register_coalescent_plugins(self):
        """注册溯祖插件"""
        try:
            from YR_MPE.plugins.caster_site_plugin import CASTERSitePluginEntry
            self.register_plugin("caster_site", CASTERSitePluginEntry)
        except ImportError:
            pass
        
    def _register_dataset_plugins(self):
        """注册Dataset插件"""
        try:
            from YR_MPE.plugins.dataset_plugin import DatasetPluginEntry
            self.register_plugin("dataset", DatasetPluginEntry)
        except ImportError:
            pass
            
    def get_available_plugins(self) -> Dict[str, type]:
        """
        获取所有可用插件
        
        Returns:
            插件名称到插件类的映射
        """
        return self.registry._plugins.copy()
        
    def create_plugin_instance(self, plugin_name: str, **kwargs) -> Any:
        """
        创建插件实例
        
        Args:
            plugin_name: 插件名称
            **kwargs: 传递给插件构造函数的参数
            
        Returns:
            插件实例
        """
        try:
            if plugin_name in self.registry._plugins:
                plugin_class = self.registry._plugins[plugin_name]
                return plugin_class(**kwargs)
            else:
                raise ValueError(f"Plugin '{plugin_name}' not found")
        except Exception as e:
            self.logger.error(f"Failed to create plugin instance {plugin_name}: {e}")
            raise
            
    def create_plugin_adapter(self, plugin_name: str, **kwargs) -> PluginAdapter:
        """
        创建插件适配器
        
        Args:
            plugin_name: 插件名称
            **kwargs: 传递给插件构造函数的参数
            
        Returns:
            插件适配器实例
        """
        plugin_instance = self.create_plugin_instance(plugin_name, **kwargs)
        adapter = PluginAdapter(plugin_instance, plugin_name)
        self.active_adapters[plugin_name] = adapter
        return adapter
        
    def execute_plugin(self, plugin_name: str, input_data: Any = None, 
                      parameters: Dict[str, Any] = None) -> PluginAdapter:
        """
        执行插件并返回适配器
        
        Args:
            plugin_name: 插件名称
            input_data: 输入数据
            parameters: 插件参数
            
        Returns:
            插件适配器实例
        """
        kwargs = {"import_from": "YR_MPEA"}
        if input_data is not None:
            kwargs["import_data"] = input_data
            
        adapter = self.create_plugin_adapter(plugin_name, **kwargs)
        return adapter