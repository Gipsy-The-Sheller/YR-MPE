import importlib
from typing import Any


class PluginFactory:
    """插件工厂 - 统一插件加载和管理"""
    
    def __init__(self):
        self._plugins = {}
    
    def load_plugin(self, plugin_path: str, class_name: str = None):
        """动态加载插件
        Args:
            plugin_path: 插件路径，例如 "YR_MPE.plugins.dataset_plugin"
            class_name: 类名，如果为None则使用模块中的默认类
        """
        try:
            module = importlib.import_module(plugin_path)
            if class_name:
                plugin_class = getattr(module, class_name)
            else:
                # 假设模块中有默认的插件类
                plugin_class = getattr(module, 'PluginEntry')
            return plugin_class()
        except (ImportError, AttributeError) as e:
            raise ImportError(f"Failed to load plugin {plugin_path}: {str(e)}")
    
    def get_sequence_viewer(self):
        """获取序列查看器插件"""
        return self.load_plugin("YR_MPE.sequence_editor", "SequenceEditorEntry")
    
    def get_dataset_manager(self):
        """获取数据集管理器"""
        return self.load_plugin("YR_MPE.platforms.methods.dataset_manager", "DatasetManager")
    
    def get_iqtree_plugin(self):
        """获取IQ-TREE插件"""
        return self.load_plugin("YR_MPE.plugins.iqtree_plugin", "IQTreePluginEntry")
    
    def get_muscle5_plugin(self):
        """获取MUSCLE5插件"""
        return self.load_plugin("YR_MPE.plugins.muscle5_plugin", "Muscle5PluginEntry")
    
    def get_mafft_plugin(self):
        """获取MAFFT插件"""
        return self.load_plugin("YR_MPE.plugins.mafft_plugin", "MAFFTPluginEntry")
    
    def get_clustal_omega_plugin(self):
        """获取Clustal Omega插件"""
        return self.load_plugin("YR_MPE.plugins.clustal_omega_plugin", "ClustalOmegaPluginEntry")
    
    def get_trimal_plugin(self):
        """获取TrimAl插件"""
        return self.load_plugin("YR_MPE.plugins.trimal_plugin", "TrimAlPluginEntry")
    
    def get_gblocks_plugin(self):
        """获取GBlocks插件"""
        return self.load_plugin("YR_MPE.plugins.gblocks_plugin", "GBlocksPluginEntry")
    
    def get_model_finder_plugin(self):
        """获取ModelFinder插件"""
        return self.load_plugin("YR_MPE.plugins.model_finder_plugin", "ModelFinderPluginEntry")
    
    def get_ml_distance_plugin(self):
        """获取ML距离计算插件"""
        return self.load_plugin("YR_MPE.plugins.ml_distance_plugin", "MLDistancePluginEntry")
    
    def get_caster_site_plugin(self):
        """获取CASTER-site插件"""
        return self.load_plugin("YR_MPE.plugins.caster_site_plugin", "CasterSitePluginEntry")
    
    def get_macse_plugin(self):
        """获取MACSE插件"""
        return self.load_plugin("YR_MPE.plugins.macse_plugin", "MACSEPluginEntry")
