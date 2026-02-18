"""
选择策略引擎
实现复杂的选择规则和策略
"""
from typing import Dict, List, Set, Optional, Tuple
from collections import deque

from .dataset_models import (
    DatasetItem, DatasetInfo, SelectionNode,
    ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT, ITEM_TYPE_MODEL,
    ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY, ITEM_TYPE_VARIANT,
    ITEM_TYPE_COALESCENT, ITEM_TYPE_CLOCK, ITEM_TYPE_CHAIN,
    SELECTION_STATE_NONE, SELECTION_STATE_GREEN,
    SELECTION_STATE_BLUE, SELECTION_STATE_RED,
    DATA_DEPENDENCIES, ITEM_TYPE_TO_GROUP, OWNERSHIP_GROUPS
)


class SelectionEngine:
    """选择策略引擎 - 实现复杂的选择规则"""
    
    def __init__(self, selection_manager):
        """
        初始化选择引擎
        
        Args:
            selection_manager: DatasetSelectionManager实例
        """
        self.manager = selection_manager
    
    # ========== 主要选择方法 ==========
    
    def select_item(self, item_id: str, append: bool = False) -> bool:
        """
        选中数据项（基于所有权 UUID 机制）
        
        规则：
        1. 清除同类型的不同 UUID 的激活项（确保同一类型只有一个激活）
        2. 清除不同类型且不同 UUID 的激活项
        3. 保留不同类型但相同 UUID 的激活项
        4. 将当前项设置为绿色
        5. 高亮同 UUID 的其他未激活项为蓝色
        
        Args:
            item_id: 数据项ID
            append: 是否追加选择（暂不使用，保留接口兼容）
        
        Returns:
            是否成功选择
        """
        item = self.manager.get_item(item_id)
        if not item:
            return False
        
        # 获取当前项的类型分组
        current_group = ITEM_TYPE_TO_GROUP.get(item.item_type)
        
        # 获取所有当前激活（绿色）的数据项
        green_items = self.manager.get_items_by_state(SELECTION_STATE_GREEN)
        
        # 遍历所有绿色项，根据规则处理
        for green_item in green_items:
            # 如果是当前项本身，跳过
            if green_item.id == item.id:
                continue
            
            # 检查是否是同一类型分组
            green_item_group = ITEM_TYPE_TO_GROUP.get(green_item.item_type)
            
            if green_item_group == current_group:
                # 同一类型分组：清除激活（确保同一类型只有一个激活）
                self._set_item_state(green_item, SELECTION_STATE_NONE, "")
            elif green_item.ownership_uuid != item.ownership_uuid:
                # 不同类型且不同 UUID：清除激活
                self._set_item_state(green_item, SELECTION_STATE_NONE, "")
            # else: 不同类型但相同 UUID：保留激活状态
        
        # 将当前项设置为绿色
        self._set_item_state(item, SELECTION_STATE_GREEN, "Direct selection")
        
        # 基于所有权 UUID 的高亮机制
        if item.ownership_uuid:
            self._apply_ownership_highlighting(item)
        
        # 选择祖先节点
        self._select_ancestors(item)
        
        # 执行冲突检测
        self._check_conflicts(item)
        
        return True
    
    def _apply_ownership_highlighting(self, item: DatasetItem):
        """
        应用基于所有权 UUID 的高亮机制
        
        规则：
        1. 高亮同 UUID 的其他未激活项目为蓝色（上下文高亮）
        2. 不改变已激活（绿色）的项目
        """
        # 遍历所有数据项
        for other_item in self.manager.items.values():
            if other_item.id == item.id:
                continue
            
            # 如果有相同的 UUID，设置为蓝色（上下文高亮）
            # 但不要改变已经激活的项目
            if (other_item.ownership_uuid and 
                other_item.ownership_uuid == item.ownership_uuid and
                other_item.selection_state != SELECTION_STATE_GREEN):
                self._set_item_state(other_item, SELECTION_STATE_BLUE,
                                   f"Same ownership as {item.get_name()}")
    
    def _perform_auto_selection(self, item: DatasetItem):
        """
        执行自动关联选择（规则2）
        
        当选择某个数据项时，根据其类型自动选择相关的数据项
        """
        # 选择所有祖先节点
        self._select_ancestors(item)
        
        # 选择相关项（根据数据类型）
        self._select_related_items(item)
        
        # 选择所属数据集
        if item.dataset_id:
            self._select_dataset_context(item.dataset_id)
    
    def _select_ancestors(self, item: DatasetItem):
        """
        选择所有祖先节点
        
        从当前项向上遍历，将所有父项设置为绿色
        """
        current_id = item.parent_id
        
        while current_id:
            parent_item = self.manager.get_item(current_id)
            if not parent_item:
                break
            
            self._set_item_state(parent_item, SELECTION_STATE_GREEN, 
                               f"Ancestor of {item.get_name()}")
            
            current_id = parent_item.parent_id
    
    def _select_related_items(self, item: DatasetItem):
        """
        选择相关项（规则2）
        
        根据数据类型选择相关的数据项：
        - alignment: 选择model和variant
        - model: 选择distance和phylogeny
        - distance: 选择phylogeny
        - phylogeny: 选择clock
        """
        # 定义每种数据类型的相关类型
        related_types_mapping = {
            ITEM_TYPE_ALIGNMENT: [ITEM_TYPE_MODEL, ITEM_TYPE_VARIANT],
            ITEM_TYPE_MODEL: [ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY],
            ITEM_TYPE_DISTANCE: [ITEM_TYPE_PHYLOGENY],
            ITEM_TYPE_PHYLOGENY: [ITEM_TYPE_CLOCK],
            ITEM_TYPE_SEQUENCE: [],  # sequence没有相关项
            ITEM_TYPE_VARIANT: [],   # variant是终端项
            ITEM_TYPE_COALESCENT: [], # coalescent是终端项
            ITEM_TYPE_CLOCK: []      # clock是终端项
        }
        
        related_types = related_types_mapping.get(item.item_type, [])
        
        for related_type in related_types:
            self._select_latest_related(item, related_type)
    
    def _select_latest_related(self, item: DatasetItem, related_type: str):
        """
        选择指定类型的最新相关项（规则3）
        
        对于同类型的多个子项，只有最新的为green，其他为blue
        """
        # 获取同数据集、相同parent_id的指定类型的数据项
        related_items = [
            ri for ri in self.manager.items.values()
            if (ri.dataset_id == item.dataset_id and 
                ri.item_type == related_type and 
                ri.parent_id == item.id)
        ]
        
        if not related_items:
            return
        
        # 按创建时间排序，最新的为绿色，其他为蓝色（规则3）
        related_items.sort(key=lambda x: x.created_at, reverse=True)
        
        # 选中最新的
        latest = related_items[0]
        self._set_item_state(latest, SELECTION_STATE_GREEN,
                           f"Latest {related_type} of {item.get_name()}")
        
        # 其他设置为蓝色
        for ri in related_items[1:]:
            self._set_item_state(ri, SELECTION_STATE_BLUE,
                               f"Related {related_type} of {item.get_name()}")
    
    def _select_dataset_context(self, dataset_id: str):
        """
        选择数据集上下文（规则4）
        
        将同一数据集的其他不相关数据项设置为blue
        """
        if dataset_id not in self.manager.datasets:
            return
        
        # 获取数据集的所有数据项
        all_items = self.manager.get_items_by_dataset(dataset_id)
        
        # 获取所有绿色和蓝色的数据项
        selected_items = [
            item for item in all_items 
            if item.selection_state in [SELECTION_STATE_GREEN, SELECTION_STATE_BLUE]
        ]
        
        # 获取未选中的数据项
        unselected_items = [
            item for item in all_items 
            if item.selection_state == SELECTION_STATE_NONE
        ]
        
        # 将未选中的数据项设置为蓝色（规则4：上下文高亮）
        for item in unselected_items:
            self._set_item_state(item, SELECTION_STATE_BLUE,
                               f"Context of dataset {dataset_id}")
    
    def _check_conflicts(self, item: DatasetItem):
        """
        检查冲突（规则6）
        
        检查是否存在冲突：
        - 属于不同数据集
        - 数据不兼容
        - 依赖关系不完整
        """
        conflicts = []
        
        # 检查所有绿色数据项是否属于同一数据集
        green_items = self.manager.get_items_by_state(SELECTION_STATE_GREEN)
        dataset_ids = set(i.dataset_id for i in green_items)
        
        if len(dataset_ids) > 1:
            # 存在跨数据集选择
            conflicts.append("Items belong to different datasets")
            
            # 将非当前数据集的数据项设置为红色
            for green_item in green_items:
                if green_item.dataset_id != item.dataset_id:
                    self._set_item_state(green_item, SELECTION_STATE_RED,
                                       f"Conflict: Different dataset than {item.get_name()}")
        
        # 检查数据兼容性
        for green_item in green_items:
            if green_item.id != item.id and not item.is_compatible_with(green_item):
                conflicts.append(f"Incompatible data: {item.get_name()} and {green_item.get_name()}")
                self._set_item_state(green_item, SELECTION_STATE_RED,
                                   f"Conflict: Incompatible with {item.get_name()}")
        
        # 检查依赖关系完整性
        dependency_conflicts = self._check_dependency_completeness(green_items)
        if dependency_conflicts:
            conflicts.extend(dependency_conflicts)
    
    def _check_dependency_completeness(self, items: List[DatasetItem]) -> List[str]:
        """
        检查依赖关系完整性
        
        Args:
            items: 要检查的数据项列表
        
        Returns:
            冲突错误列表
        """
        conflicts = []
        
        for item in items:
            # 获取该类型所需的依赖
            required_types = DATA_DEPENDENCIES.get(item.item_type, [])
            
            # 对于某些操作，依赖是可选的
            # 例如：phylogeny可以基于alignment，也可以基于model+distance
            if item.item_type == ITEM_TYPE_PHYLOGENY:
                # 检查是否有alignment或(model+distance)
                has_alignment = any(
                    i.item_type == ITEM_TYPE_ALIGNMENT and i.parent_id == item.parent_id
                    for i in items
                )
                has_model = any(
                    i.item_type == ITEM_TYPE_MODEL and i.parent_id == item.parent_id
                    for i in items
                )
                has_distance = any(
                    i.item_type == ITEM_TYPE_DISTANCE and i.parent_id == item.parent_id
                    for i in items
                )
                
                if not (has_alignment or (has_model and has_distance)):
                    conflicts.append(
                        f"Phylogeny {item.get_name()} requires alignment or model+distance"
                    )
        
        return conflicts
    
    def clear_all_selections(self):
        """清除所有选择"""
        for item in self.manager.items.values():
            self._set_item_state(item, SELECTION_STATE_NONE, "")
    
    def _set_item_state(self, item: DatasetItem, state: str, reason: str):
        """
        设置数据项的选择状态
        
        Args:
            item: 数据项
            state: 选择状态
            reason: 选择原因
        """
        item.selection_state = state
        item.selection_reason = reason
        
        # 更新选中项集合
        if state == SELECTION_STATE_GREEN:
            self.manager.selected_items.add(item.id)
        elif item.id in self.manager.selected_items:
            self.manager.selected_items.remove(item.id)
    
    # ========== 特殊选择方法 ==========
    
    def select_dataset(self, dataset_id: str):
        """
        选择整个数据集
        
        将数据集的所有数据项设置为蓝色（规则4）
        """
        if dataset_id not in self.manager.datasets:
            return
        
        # 清除所有选择
        self.clear_all_selections()
        
        # 将数据集的所有数据项设置为蓝色
        items = self.manager.get_items_by_dataset(dataset_id)
        for item in items:
            self._set_item_state(item, SELECTION_STATE_BLUE,
                               f"Dataset: {self.manager.datasets[dataset_id].name}")
    
    def select_coalescent(self, coalescent_id: str):
        """
        选择溯祖分析结果
        
        特殊规则：当选择Coalescent时，自动选择所有参与溯祖分析的Phylogeny及其祖先项
        
        Args:
            coalescent_id: Coalescent数据项ID
        """
        coalescent_item = self.manager.get_item(coalescent_id)
        if not coalescent_item or coalescent_item.item_type != ITEM_TYPE_COALESCENT:
            return
        
        # 清除所有选择
        self.clear_all_selections()
        
        # 选择Coalescent
        self._set_item_state(coalescent_item, SELECTION_STATE_GREEN, "Coalescent selection")
        
        # 选择所有相关的Phylogeny
        # 假设coalescent的parent_id指向第一个phylogeny
        # 实际实现需要根据具体的数据结构调整
        phylogeny_id = coalescent_item.parent_id
        if phylogeny_id:
            phylogeny_item = self.manager.get_item(phylogeny_id)
            if phylogeny_item and phylogeny_item.item_type == ITEM_TYPE_PHYLOGENY:
                self._set_item_state(phylogeny_item, SELECTION_STATE_GREEN,
                                   f"Phylogeny for {coalescent_item.get_name()}")
                
                # 选择Phylogeny的祖先
                self._select_ancestors(phylogeny_item)
                
                # 如果有多个Phylogeny参与Coalescent，需要额外处理
                # TODO: 根据实际数据结构调整
    
    def select_by_type(self, dataset_id: str, item_type: str):
        """
        选择指定数据集的指定类型的所有数据项
        
        Args:
            dataset_id: 数据集ID
            item_type: 数据项类型
        """
        if dataset_id not in self.manager.datasets:
            return
        
        # 清除所有选择
        self.clear_all_selections()
        
        # 获取指定类型的所有数据项
        items = self.manager.get_items_by_type(item_type)
        dataset_items = [item for item in items if item.dataset_id == dataset_id]
        
        # 将所有数据项设置为绿色
        for item in dataset_items:
            self._set_item_state(item, SELECTION_STATE_GREEN,
                               f"Selected by type: {item_type}")
    
    # ========== 验证和检查方法 ==========
    
    def validate_selection(self) -> Tuple[bool, List[str]]:
        """
        验证当前选择的有效性
        
        Returns:
            (是否有效, 错误列表)
        """
        green_items = self.manager.get_items_by_state(SELECTION_STATE_GREEN)
        
        if not green_items:
            return False, ["No items selected"]
        
        errors = []
        
        # 检查是否都属于同一数据集
        dataset_ids = set(item.dataset_id for item in green_items)
        if len(dataset_ids) > 1:
            errors.append("Items belong to different datasets")
        
        # 检查数据有效性
        for item in green_items:
            if not item.is_valid:
                errors.extend(item.validation_errors)
        
        # 检查依赖关系
        dependency_errors = self._check_dependency_completeness(green_items)
        errors.extend(dependency_errors)
        
        return len(errors) == 0, errors
    
    def get_selection_summary(self) -> Dict[str, any]:
        """
        获取当前选择的摘要信息
        
        Returns:
            选择摘要字典
        """
        green_items = self.manager.get_items_by_state(SELECTION_STATE_GREEN)
        blue_items = self.manager.get_items_by_state(SELECTION_STATE_BLUE)
        red_items = self.manager.get_items_by_state(SELECTION_STATE_RED)
        
        return {
            "green_count": len(green_items),
            "blue_count": len(blue_items),
            "red_count": len(red_items),
            "by_type": {
                item_type: len([i for i in green_items if i.item_type == item_type])
                for item_type in [ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT, ITEM_TYPE_MODEL,
                                 ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY, ITEM_TYPE_VARIANT,
                                 ITEM_TYPE_COALESCENT, ITEM_TYPE_CLOCK]
            },
            "datasets": list(set(item.dataset_id for item in green_items))
        }