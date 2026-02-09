"""
YR-MPEA 插件系统方法模块
"""

from .plugin_manager import PluginManager
from .plugin_adapter import PluginAdapter
from .plugin_executor import PluginExecutor
from .plugin_config_manager import PluginConfigManager
from .plugin_factory import PluginFactory
from .workspace_manager import WorkspaceManager
from .plugin_interface import BasePlugin
from .file_operations import FileOperations
from .ui_helpers import UIHelpers

# 数据集选择机制相关导入
from .dataset_models import (
    DatasetItem, DatasetInfo, SelectionNode,
    ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT, ITEM_TYPE_MODEL,
    ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY, ITEM_TYPE_VARIANT,
    ITEM_TYPE_COALESCENT, ITEM_TYPE_CLOCK,
    SELECTION_STATE_NONE, SELECTION_STATE_GREEN,
    SELECTION_STATE_BLUE, SELECTION_STATE_RED,
    ALL_ITEM_TYPES, ALL_SELECTION_STATES
)
from .dataset_selection_manager import DatasetSelectionManager
from .selection_engine import SelectionEngine
from .workflow_manager import WorkflowManager
from .dataset_item_button import DatasetItemButton