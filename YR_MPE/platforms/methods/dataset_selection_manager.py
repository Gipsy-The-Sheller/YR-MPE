"""
数据集选择管理器
实现数据所有权与选择机制的核心管理逻辑
"""
import os
import json
from datetime import datetime
from typing import Dict, List, Any, Optional, Set, Tuple
from collections import deque

from .dataset_models import (
    DatasetItem, DatasetInfo, SelectionNode,
    ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT, ITEM_TYPE_MODEL,
    ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY, ITEM_TYPE_TREESET,
    ITEM_TYPE_VARIANT, ITEM_TYPE_COALESCENT, ITEM_TYPE_CLOCK,
    ALL_ITEM_TYPES, ALL_SELECTION_STATES,
    SELECTION_STATE_NONE, SELECTION_STATE_GREEN,
    SELECTION_STATE_BLUE, SELECTION_STATE_RED,
    DATA_DEPENDENCIES
)


class DatasetSelectionManager:
    """数据集选择管理器 - 负责管理数据集、数据项和选择状态"""
    
    def __init__(self, disable_auto_save: bool = False):
        # 数据存储
        self.datasets: Dict[str, DatasetInfo] = {}  # 数据集ID -> DatasetInfo
        self.items: Dict[str, DatasetItem] = {}     # 数据项ID -> DatasetItem
        self.selection_trees: Dict[str, SelectionNode] = {}  # 数据集ID -> 选择树根节点
        
        # 当前选中的数据项
        self.selected_items: Set[str] = set()
        
        # 管理器状态
        self.workspace_path: Optional[str] = None
        self.version = "2.0"
        self.disable_auto_save = disable_auto_save  # 禁用自动保存
        
        # 初始化选择引擎（延迟初始化）
        self._selection_engine = None
        
    # ========== 数据集管理 ==========
    
    def create_dataset(self, name: str, description: str = "", is_multigene: bool = False) -> str:
        """创建新数据集"""
        dataset = DatasetInfo()
        dataset.name = name
        dataset.description = description
        dataset.is_multigene = is_multigene
        
        self.datasets[dataset.id] = dataset
        self._save_state()
        
        return dataset.id
    
    def get_dataset(self, dataset_id: str) -> Optional[DatasetInfo]:
        """获取数据集"""
        return self.datasets.get(dataset_id)
    
    def delete_dataset(self, dataset_id: str) -> bool:
        """删除数据集及其所有数据项"""
        if dataset_id not in self.datasets:
            return False
        
        # 删除数据集的所有数据项
        dataset = self.datasets[dataset_id]
        for item_id in dataset.items:
            if item_id in self.items:
                del self.items[item_id]
        
        # 删除选择树
        if dataset_id in self.selection_trees:
            del self.selection_trees[dataset_id]
        
        # 删除数据集
        del self.datasets[dataset_id]
        
        # 从选中项中移除
        self.selected_items = {item_id for item_id in self.selected_items 
                              if item_id in self.items}
        
        self._save_state()
        return True
    
    def get_all_datasets(self) -> List[DatasetInfo]:
        """获取所有数据集"""
        return list(self.datasets.values())
    
    # ========== 数据项管理 ==========
    
    def add_item(self, item: DatasetItem, dataset_id: str) -> bool:
        """添加数据项到数据集"""
        if dataset_id not in self.datasets:
            return False
        
        # 验证数据项
        if not item.validate():
            return False
        
        # 设置数据集ID
        item.dataset_id = dataset_id
        
        # 添加到存储
        self.items[item.id] = item
        
        # 添加到数据集
        self.datasets[dataset_id].add_item(item.id)
        
        # 更新数据集的统计信息
        dataset = self.datasets[dataset_id]
        if item.item_type == "sequence":
            dataset.loci_count += 1
            if item.sequence_count > dataset.taxa_count:
                dataset.taxa_count = item.sequence_count
        
        # 重建选择树
        self._build_selection_tree(dataset_id)
        
        self._save_state()
        return True
    
    def get_item(self, item_id: str) -> Optional[DatasetItem]:
        """获取数据项"""
        return self.items.get(item_id)
    
    def delete_item(self, item_id: str) -> bool:
        """删除数据项"""
        if item_id not in self.items:
            return False
        
        item = self.items[item_id]
        dataset_id = item.dataset_id
        
        # 从数据集中移除
        if dataset_id in self.datasets:
            self.datasets[dataset_id].remove_item(item_id)
        
        # 删除数据项
        del self.items[item_id]
        
        # 从选中项中移除
        if item_id in self.selected_items:
            self.selected_items.remove(item_id)
        
        # 重建选择树
        if dataset_id:
            self._build_selection_tree(dataset_id)
        
        self._save_state()
        return True
    
    def get_items_by_type(self, item_type: str) -> List[DatasetItem]:
        """获取指定类型的所有数据项"""
        return [item for item in self.items.values() if item.item_type == item_type]
    
    def get_items_by_dataset(self, dataset_id: str) -> List[DatasetItem]:
        """获取指定数据集的所有数据项"""
        if dataset_id not in self.datasets:
            return []
        
        dataset = self.datasets[dataset_id]
        return [self.items[item_id] for item_id in dataset.items if item_id in self.items]
    
    def get_items_by_state(self, selection_state: str) -> List[DatasetItem]:
        """获取指定选择状态的所有数据项"""
        return [item for item in self.items.values() if item.selection_state == selection_state]
    
    # ========== 选择树管理 ==========
    
    def _build_selection_tree(self, dataset_id: str) -> SelectionNode:
        """构建选择树"""
        if dataset_id not in self.datasets:
            return None
        
        # 获取数据集的所有数据项
        dataset_items = self.get_items_by_dataset(dataset_id)
        
        if not dataset_items:
            return None
        
        # 创建节点映射
        nodes: Dict[str, SelectionNode] = {}
        root_nodes: List[SelectionNode] = []
        
        # 为每个数据项创建节点
        for item in dataset_items:
            nodes[item.id] = SelectionNode(item.id)
        
        # 构建父子关系
        for item in dataset_items:
            node = nodes[item.id]
            
            if item.parent_id is None or item.parent_id not in nodes:
                # 根节点
                root_nodes.append(node)
            else:
                # 子节点
                parent_node = nodes[item.parent_id]
                parent_node.add_child(node)
        
        # 如果只有一个根节点，使用它作为选择树的根
        if len(root_nodes) == 1:
            self.selection_trees[dataset_id] = root_nodes[0]
        else:
            # 如果有多个根节点，创建一个虚拟根节点
            virtual_root = SelectionNode(f"virtual_{dataset_id}")
            for root in root_nodes:
                virtual_root.add_child(root)
            self.selection_trees[dataset_id] = virtual_root
        
        return self.selection_trees[dataset_id]
    
    def get_selection_tree(self, dataset_id: str) -> Optional[SelectionNode]:
        """获取选择树"""
        if dataset_id not in self.selection_trees:
            self._build_selection_tree(dataset_id)
        return self.selection_trees.get(dataset_id)
    
    # ========== 选择机制 ==========
    
    def select_item(self, item_id: str, append: bool = False) -> bool:
        """选中数据项（直接选择）"""
        if item_id not in self.items:
            return False
        
        item = self.items[item_id]
        
        # 使用 selection_engine 进行选择（包含所有权 UUID 机制）
        return self.selection_engine.select_item(item_id, append)
    
    def _perform_auto_selection(self, item: DatasetItem):
        """执行自动关联选择"""
        # 选择所有祖先节点
        self._select_ancestors(item)
        
        # 选择相关项（根据数据类型）
        self._select_related_items(item)
        
        # 更新UI通知（TODO: 需要实现信号机制）
        pass
    
    def _select_ancestors(self, item: DatasetItem):
        """选择所有祖先节点"""
        current_id = item.parent_id
        
        while current_id and current_id in self.items:
            parent_item = self.items[current_id]
            parent_item.selection_state = SELECTION_STATE_GREEN
            parent_item.selection_reason = f"Ancestor of {item.id}"
            self.selected_items.add(current_id)
            
            current_id = parent_item.parent_id
    
    def _select_related_items(self, item: DatasetItem):
        """选择相关项"""
        # 根据数据类型选择相关项
        related_types = []
        
        if item.item_type == ITEM_TYPE_ALIGNMENT:
            related_types = [ITEM_TYPE_MODEL, ITEM_TYPE_VARIANT]
        elif item.item_type == ITEM_TYPE_MODEL:
            related_types = [ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY]
        elif item.item_type == ITEM_TYPE_DISTANCE:
            related_types = [ITEM_TYPE_PHYLOGENY]
        elif item.item_type == ITEM_TYPE_PHYLOGENY:
            related_types = [ITEM_TYPE_CLOCK]
        
        # 获取同数据集的相同parent_id的相关项
        for related_type in related_types:
            related_items = [
                ri for ri in self.items.values()
                if (ri.dataset_id == item.dataset_id and 
                    ri.item_type == related_type and 
                    ri.parent_id == item.id)
            ]
            
            if related_items:
                # 按创建时间排序，最新的为绿色，其他为蓝色
                related_items.sort(key=lambda x: x.created_at, reverse=True)
                
                # 选中最新的
                latest = related_items[0]
                latest.selection_state = SELECTION_STATE_GREEN
                latest.selection_reason = f"Latest {related_type} of {item.id}"
                self.selected_items.add(latest.id)
                
                # 其他设置为蓝色
                for ri in related_items[1:]:
                    ri.selection_state = SELECTION_STATE_BLUE
                    ri.selection_reason = f"Related {related_type} of {item.id}"
                    self.selected_items.discard(ri.id)  # 不包含在选中项中
    
    def clear_all_selections(self):
        """清除所有选择"""
        for item in self.items.values():
            item.selection_state = SELECTION_STATE_NONE
            item.selection_reason = ""
        
        self.selected_items.clear()
    
    def get_selected_items(self) -> List[DatasetItem]:
        """获取所有选中的数据项（绿色状态）"""
        return [item for item in self.items.values() 
                if item.selection_state == SELECTION_STATE_GREEN]
    
    # ========== 兼容性方法 ==========
    
    def add_alignment(self, alignment_data: Any, dataset_id: str = None) -> str:
        """添加比对结果（兼容旧接口）"""
        # 创建新的DatasetItem
        item = DatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
        
        # 根据数据类型设置属性
        if isinstance(alignment_data, dict):
            sequences = alignment_data.get('data', alignment_data.get('sequences', []))
        else:
            sequences = alignment_data
        
        if sequences and len(sequences) > 0:
            item.sequences = sequences
            item.sequence_count = len(sequences)
            # 安全地获取序列长度
            if hasattr(sequences[0], 'seq'):
                try:
                    item.length = len(str(sequences[0].seq))
                except:
                    item.length = 0
            else:
                item.length = 0
            item.is_aligned = True
        
        item.loci_name = f"Alignment_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        item.name = item.loci_name
        
        # 添加到数据集
        if dataset_id:
            self.add_item(item, dataset_id)
        else:
            # 如果没有指定数据集，创建一个临时数据集
            if not self.datasets:
                dataset_id = self.create_dataset("Default Dataset")
            self.add_item(item, dataset_id)
        
        return item.id
    
    def add_model(self, model_data: Any, parent_id: str = None, dataset_id: str = None) -> str:
        """添加模型结果（兼容旧接口）"""
        item = DatasetItem(item_type=ITEM_TYPE_MODEL)
        item.parent_id = parent_id
        
        # 存储模型数据
        if isinstance(model_data, dict):
            item.data = model_data
        else:
            item.data = {'content': model_data}
        
        item.loci_name = f"Model_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        item.name = item.loci_name
        
        # 添加到数据集
        if not dataset_id and parent_id:
            parent_item = self.items.get(parent_id)
            if parent_item:
                dataset_id = parent_item.dataset_id
        
        if dataset_id:
            self.add_item(item, dataset_id)
        
        return item.id
    
    def add_distance(self, distance_data: Any, parent_id: str = None, dataset_id: str = None) -> str:
        """添加距离矩阵（兼容旧接口）"""
        item = DatasetItem(item_type=ITEM_TYPE_DISTANCE)
        item.parent_id = parent_id
        
        # 存储距离数据
        if isinstance(distance_data, dict):
            item.data = distance_data
        else:
            item.data = {'content': distance_data}
        
        item.loci_name = f"Distance_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        item.name = item.loci_name
        
        # 添加到数据集
        if not dataset_id and parent_id:
            parent_item = self.items.get(parent_id)
            if parent_item:
                dataset_id = parent_item.dataset_id
        
        if dataset_id:
            self.add_item(item, dataset_id)
        
        return item.id
    
    def add_phylogeny(self, phylogeny_data: Any, parent_id: str = None, dataset_id: str = None) -> str:
        """添加系统发育树（兼容旧接口）"""
        item = DatasetItem(item_type=ITEM_TYPE_PHYLOGENY)
        item.parent_id = parent_id
        
        # 存储树数据
        if isinstance(phylogeny_data, dict):
            item.data = phylogeny_data
        else:
            item.data = {'content': phylogeny_data}
        
        item.loci_name = f"Phylogeny_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        item.name = item.loci_name
        
        # 添加到数据集
        if not dataset_id and parent_id:
            parent_item = self.items.get(parent_id)
            if parent_item:
                dataset_id = parent_item.dataset_id
        
        if dataset_id:
            self.add_item(item, dataset_id)
        
        return item.id
    
    # ========== 持久化 ==========
    
    def set_workspace_path(self, path: str):
        """设置工作区路径"""
        self.workspace_path = path
        
        # 加载保存的状态
        state_file = os.path.join(path, "dataset_selection_state.json")
        if os.path.exists(state_file):
            self._load_state()
    
    def _save_state(self):
        """保存状态到文件"""
        if self.disable_auto_save:
            return  # 禁用自动保存

        if not self.workspace_path:
            return

        state_file = os.path.join(self.workspace_path, "dataset_selection_state.json")

        state = {
            "version": self.version,
            "datasets": {dataset_id: dataset.to_dict()
                        for dataset_id, dataset in self.datasets.items()},
            "items": {item_id: item.to_dict()
                     for item_id, item in self.items.items()},
            "selected_items": list(self.selected_items)
        }

        with open(state_file, 'w', encoding='utf-8') as f:
            json.dump(state, f, indent=2, ensure_ascii=False)
    
    def _load_state(self):
        """从文件加载状态"""
        if not self.workspace_path:
            return

        state_file = os.path.join(self.workspace_path, "dataset_selection_state.json")
        if not os.path.exists(state_file):
            return

        try:
            
            with open(state_file, 'r', encoding='utf-8') as f:
                state = json.load(f)

            # 加载数据集
            for dataset_id, dataset_data in state.get("datasets", {}).items():
                self.datasets[dataset_id] = DatasetInfo.from_dict(dataset_data)
                if 'dataset_items' in dataset_data.get('settings', {}):
                    saved_items = dataset_data['settings']['dataset_items']

            # 加载数据项
            for item_id, item_data in state.get("items", {}).items():
                self.items[item_id] = DatasetItem.from_dict(item_data)


            # 如果JSON中没有items，尝试从dataset.settings中恢复
            if len(self.items) == 0:
                for dataset_id, dataset in self.datasets.items():
                    if 'dataset_items' in dataset.settings:
                        saved_items = dataset.settings['dataset_items']
                        from .dataset_models import ITEM_TYPE_ALIGNMENT
                        for item_data in saved_items:
                            # 创建新的DatasetItem
                            new_item = DatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
                            new_item.dataset_id = dataset_id
                            new_item.loci_name = item_data.get('loci_name', '')
                            new_item.file_path = item_data.get('file_path', '')
                            new_item.length = item_data.get('length', 0)
                            new_item.sequence_count = item_data.get('sequence_count', 0)
                            new_item.is_aligned = item_data.get('is_aligned', False)
                            
                            # 设置选择状态
                            if item_data.get('selected', False):
                                new_item.selection_state = SELECTION_STATE_GREEN
                                self.selected_items.add(new_item.id)
                            else:
                                new_item.selection_state = SELECTION_STATE_NONE
                            
                            # 添加到items
                            self.items[new_item.id] = new_item
                            
                            # 添加到dataset.items
                            dataset.add_item(new_item.id)
                            

            # 加载选中项（如果JSON中有）
            if "selected_items" in state:
                loaded_selected = state.get("selected_items", [])
                self.selected_items = set(loaded_selected)

            # 重建选择树
            for dataset_id in self.datasets.keys():
                self._build_selection_tree(dataset_id)

        except Exception as e:
            import traceback
            traceback.print_exc()
    
    def export_to_dict(self) -> Dict[str, Any]:
        """导出为字典格式"""
        return {
            "version": self.version,
            "datasets": {dataset_id: dataset.to_dict() 
                        for dataset_id, dataset in self.datasets.items()},
            "items": {item_id: item.to_dict() 
                     for item_id, item in self.items.items()},
            "selected_items": list(self.selected_items)
        }
    
    def import_from_dict(self, data: Dict[str, Any]):
        """从字典格式导入"""
        self.version = data.get("version", "2.0")
        
        # 导入数据集
        for dataset_id, dataset_data in data.get("datasets", {}).items():
            self.datasets[dataset_id] = DatasetInfo.from_dict(dataset_data)
        
        # 导入数据项
        for item_id, item_data in data.get("items", {}).items():
            self.items[item_id] = DatasetItem.from_dict(item_data)
        
        # 导入选中项
        self.selected_items = set(data.get("selected_items", []))
        
        # 重建选择树
        for dataset_id in self.datasets.keys():
            self._build_selection_tree(dataset_id)
    
    @property
    def selection_engine(self):
        """获取选择引擎实例（延迟初始化）"""
        if self._selection_engine is None:
            from .selection_engine import SelectionEngine
            self._selection_engine = SelectionEngine(self)
        return self._selection_engine