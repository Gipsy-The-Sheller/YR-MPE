"""
Workflow Manager Module
实现工作流的半自动调用功能
"""
from typing import Dict, List, Any, Optional, Callable
from .dataset_selection_manager import DatasetSelectionManager
from .dataset_models import (
    DatasetItem,
    ALL_ITEM_TYPES,
    SELECTION_STATE_GREEN,
    SELECTION_STATE_NONE,
    ITEM_TYPE_ALIGNMENT,
    ITEM_TYPE_MODEL,
    ITEM_TYPE_DISTANCE,
    ITEM_TYPE_PHYLOGENY,
    ITEM_TYPE_SEQUENCE
)


class WorkflowManager:
    """
    工作流管理器 - 负责根据选中的数据项自动调用相应的分析工具
    
    实现规则：
    - 只有绿色状态的数据项参与工作流调用
    - 自动验证数据项的完整性和兼容性
    - 自动构建工作流链
    """
    
    def __init__(self, selection_manager: DatasetSelectionManager):
        """
        初始化工作流管理器
        
        Args:
            selection_manager: 数据集选择管理器
        """
        self.selection_manager = selection_manager
        self.workflow_handlers: Dict[str, Callable] = {}
        self._register_default_handlers()
    
    def _register_default_handlers(self):
        """注册默认的工作流处理器"""
        # 注册各种工作流处理器
        self.workflow_handlers = {
            "align_sequence": self._handle_align_sequence,
            "compute_distance": self._handle_compute_distance,
            "build_phylogeny": self._handle_build_phylogeny,
            "find_best_model": self._handle_find_best_model,
            "trim_alignment": self._handle_trim_alignment,
        }
    
    def register_handler(self, workflow_name: str, handler: Callable):
        """
        注册自定义工作流处理器
        
        Args:
            workflow_name: 工作流名称
            handler: 处理函数
        """
        self.workflow_handlers[workflow_name] = handler
    
    def trigger_workflow(self, workflow_name: str, **kwargs) -> bool:
        """
        触发工作流
        
        Args:
            workflow_name: 工作流名称
            **kwargs: 额外参数
            
        Returns:
            是否成功触发
        """
        handler = self.workflow_handlers.get(workflow_name)
        if not handler:
            print(f"Warning: Unknown workflow '{workflow_name}'")
            return False
        
        try:
            return handler(**kwargs)
        except Exception as e:
            print(f"Error executing workflow '{workflow_name}': {str(e)}")
            return False
    
    def _get_selected_items(self, item_type: Optional[str] = None) -> List[DatasetItem]:
        """
        获取选中的数据项
        
        Args:
            item_type: 可选的数据类型过滤
            
        Returns:
            选中的数据项列表
        """
        items = self.selection_manager.get_selected_items()
        
        if item_type:
            items = [item for item in items if item.item_type == item_type]
        
        return items
    
    def _validate_workflow_input(self, required_type: str, min_count: int = 1) -> List[DatasetItem]:
        """
        验证工作流输入
        
        Args:
            required_type: 需要的数据类型
            min_count: 最少需要的数量
            
        Returns:
            验证通过的数据项列表，否则返回空列表
        """
        items = self._get_selected_items(required_type)
        
        if len(items) < min_count:
            print(f"Error: Need at least {min_count} {required_type} item(s), but got {len(items)}")
            return []
        
        return items
    
    # ========== 工作流处理器 ==========
    
    def _handle_align_sequence(self, **kwargs) -> bool:
        """
        处理序列比对工作流
        
        输入：选中的sequence数据项
        输出：alignment数据项
        """
        sequences = self._validate_workflow_input(ITEM_TYPE_SEQUENCE)
        if not sequences:
            return False
        
        # 获取比对工具参数
        aligner = kwargs.get("aligner", "mafft")
        
        print(f"Aligning {len(sequences)} sequence(s) using {aligner}...")
        
        # TODO: 实际调用比对工具
        # 这里应该调用相应的比对插件
        
        return True
    
    def _handle_compute_distance(self, **kwargs) -> bool:
        """
        处理距离计算工作流
        
        输入：选中的model数据项
        输出：distance数据项
        """
        models = self._validate_workflow_input(ITEM_TYPE_MODEL)
        if not models:
            return False
        
        print(f"Computing distances for {len(models)} model(s)...")
        
        # TODO: 实际调用距离计算工具
        
        return True
    
    def _handle_build_phylogeny(self, **kwargs) -> bool:
        """
        处理系统发育树构建工作流
        
        输入：选中的alignment或model+distance数据项
        输出：phylogeny数据项
        """
        # 优先使用alignment
        alignments = self._validate_workflow_input(ITEM_TYPE_ALIGNMENT)
        
        if alignments:
            print(f"Building phylogeny from {len(alignments)} alignment(s)...")
            # TODO: 实际调用系统发育树构建工具
            return True
        
        # 如果没有alignment，尝试使用model+distance
        models = self._validate_workflow_input(ITEM_TYPE_MODEL)
        distances = self._validate_workflow_input(ITEM_TYPE_DISTANCE)
        
        if models and distances:
            print(f"Building phylogeny from {len(models)} model(s) and {len(distances)} distance(s)...")
            # TODO: 实际调用系统发育树构建工具
            return True
        
        print("Error: Need either alignment or model+distance to build phylogeny")
        return False
    
    def _handle_find_best_model(self, **kwargs) -> bool:
        """
        处理模型选择工作流
        
        输入：选中的alignment数据项
        输出：model数据项
        """
        alignments = self._validate_workflow_input(ITEM_TYPE_ALIGNMENT)
        if not alignments:
            return False
        
        print(f"Finding best model for {len(alignments)} alignment(s)...")
        
        # TODO: 实际调用模型选择工具
        
        return True
    
    def _handle_trim_alignment(self, **kwargs) -> bool:
        """
        处理序列修剪工作流
        
        输入：选中的alignment数据项
        输出：alignment数据项（修剪后的）
        """
        alignments = self._validate_workflow_input(ITEM_TYPE_ALIGNMENT)
        if not alignments:
            return False
        
        trimmer = kwargs.get("trimmer", "trimal")
        
        print(f"Trimming {len(alignments)} alignment(s) using {trimmer}...")
        
        # TODO: 实际调用序列修剪工具
        
        return True
    
    # ========== 辅助方法 ==========
    
    def get_available_workflows(self) -> List[str]:
        """
        获取所有可用的工作流
        
        Returns:
            工作流名称列表
        """
        return list(self.workflow_handlers.keys())
    
    def can_execute_workflow(self, workflow_name: str) -> bool:
        """
        检查是否可以执行指定的工作流
        
        Args:
            workflow_name: 工作流名称
            
        Returns:
            是否可以执行
        """
        if workflow_name not in self.workflow_handlers:
            return False
        
        # 根据工作流类型检查是否有足够的选中项
        # TODO: 实现更详细的检查
        
        return True
    
    def get_workflow_requirements(self, workflow_name: str) -> Dict[str, Any]:
        """
        获取工作流的执行要求
        
        Args:
            workflow_name: 工作流名称
            
        Returns:
            包含要求信息的字典
        """
        requirements = {
            "align_sequence": {
                "required_type": ITEM_TYPE_SEQUENCE,
                "min_count": 1,
                "description": "Select sequence(s) to align"
            },
            "compute_distance": {
                "required_type": ITEM_TYPE_MODEL,
                "min_count": 1,
                "description": "Select model(s) to compute distances"
            },
            "build_phylogeny": {
                "required_type": ITEM_TYPE_ALIGNMENT,
                "min_count": 1,
                "alternative": {"required_types": [ITEM_TYPE_MODEL, ITEM_TYPE_DISTANCE], "min_count": [1, 1]},
                "description": "Select alignment(s) or model(s)+distance(s) to build phylogeny"
            },
            "find_best_model": {
                "required_type": ITEM_TYPE_ALIGNMENT,
                "min_count": 1,
                "description": "Select alignment(s) to find best model"
            },
            "trim_alignment": {
                "required_type": ITEM_TYPE_ALIGNMENT,
                "min_count": 1,
                "description": "Select alignment(s) to trim"
            }
        }
        
        return requirements.get(workflow_name, {})