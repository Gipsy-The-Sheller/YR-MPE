"""
数据集模型模块
实现数据所有权与选择机制的核心数据模型
"""
import uuid
from datetime import datetime
from typing import List, Dict, Any, Optional, Set
from Bio.SeqRecord import SeqRecord


# 数据类型常量
ITEM_TYPE_SEQUENCE = "sequence"
ITEM_TYPE_ALIGNMENT = "alignment"
ITEM_TYPE_MODEL = "model"
ITEM_TYPE_DISTANCE = "distance"
ITEM_TYPE_PHYLOGENY = "phylogeny"
ITEM_TYPE_VARIANT = "variant"
ITEM_TYPE_COALESCENT = "coalescent"
ITEM_TYPE_CLOCK = "clock"
ITEM_TYPE_CHAIN = "chain"

ALL_ITEM_TYPES = [
    ITEM_TYPE_SEQUENCE,
    ITEM_TYPE_ALIGNMENT,
    ITEM_TYPE_MODEL,
    ITEM_TYPE_DISTANCE,
    ITEM_TYPE_PHYLOGENY,
    ITEM_TYPE_VARIANT,
    ITEM_TYPE_COALESCENT,
    ITEM_TYPE_CLOCK,
    ITEM_TYPE_CHAIN
]

# 选择状态常量
SELECTION_STATE_NONE = "none"
SELECTION_STATE_GREEN = "green"
SELECTION_STATE_BLUE = "blue"
SELECTION_STATE_RED = "red"

ALL_SELECTION_STATES = [
    SELECTION_STATE_NONE,
    SELECTION_STATE_GREEN,
    SELECTION_STATE_BLUE,
    SELECTION_STATE_RED
]

# 数据依赖关系
DATA_DEPENDENCIES = {
    ITEM_TYPE_SEQUENCE: [],
    ITEM_TYPE_ALIGNMENT: [ITEM_TYPE_SEQUENCE],
    ITEM_TYPE_MODEL: [ITEM_TYPE_ALIGNMENT],
    ITEM_TYPE_DISTANCE: [ITEM_TYPE_MODEL],
    ITEM_TYPE_PHYLOGENY: [ITEM_TYPE_ALIGNMENT, ITEM_TYPE_MODEL, ITEM_TYPE_DISTANCE],
    ITEM_TYPE_VARIANT: [ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT],
    ITEM_TYPE_COALESCENT: [ITEM_TYPE_PHYLOGENY],
    ITEM_TYPE_CLOCK: [ITEM_TYPE_PHYLOGENY],
    ITEM_TYPE_CHAIN: [ITEM_TYPE_ALIGNMENT]
}

# 所有权分组（同一组内只能同时激活一个项目）
# 根据用户要求：
# - Dataset / Alignment / Sequence 算一种
# - distance 算一种
# - model 算一种
# - tree (phylogeny) 算一种
# - trace plot (chain) 算一种
OWNERSHIP_GROUPS = {
    "sequence_group": [ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT],
    "distance_group": [ITEM_TYPE_DISTANCE],
    "model_group": [ITEM_TYPE_MODEL],
    "phylogeny_group": [ITEM_TYPE_PHYLOGENY],
    "chain_group": [ITEM_TYPE_CHAIN],
    "variant_group": [ITEM_TYPE_VARIANT],
    "coalescent_group": [ITEM_TYPE_COALESCENT],
    "clock_group": [ITEM_TYPE_CLOCK]
}

# 反向映射：item_type -> group_name
ITEM_TYPE_TO_GROUP = {}
for group_name, item_types in OWNERSHIP_GROUPS.items():
    for item_type in item_types:
        ITEM_TYPE_TO_GROUP[item_type] = group_name


class DatasetItem:
    """数据集项模型 - 扩展版本，支持数据所有权和选择机制"""
    
    def __init__(self, item_type: str = ITEM_TYPE_SEQUENCE):
        # 基础属性
        self.id: str = str(uuid.uuid4())  # 唯一标识符
        self.dataset_id: str = ""          # 所属数据集ID
        self.loci_name: str = ""           # 位点名称/文件名
        self.file_path: str = ""           # 原始文件路径
        self.item_type: str = item_type    # 数据项类型
        self.parent_id: Optional[str] = None  # 父项ID
        self.created_at: datetime = datetime.now()
        self.modified_at: datetime = datetime.now()
        
        # 序列数据（如果是sequence或alignment类型）
        self.sequences: List[SeqRecord] = []
        self.is_aligned: bool = False
        
        # 元数据
        self.length: int = 0                  # 序列/比对长度
        self.sequence_count: int = 0          # 序列数量
        self.format: str = ""                 # 文件格式：'fasta', 'phylip', 'nexus', 'newick'
        
        # 特定类型的数据
        self.data: Dict[str, Any] = {}        # 存储特定类型的数据
        self.metadata: Dict[str, Any] = {}    # 额外的元数据
        
        # 选择状态
        self.selection_state: str = SELECTION_STATE_NONE
        self.selection_reason: str = ""       # 选择原因（用于调试和日志）
        
        # 所有权 UUID（用于高亮机制）
        self.ownership_uuid: str = ""         # 数据所有权标识，同 UUID 的数据会一起高亮
        
        # 验证状态
        self.is_valid: bool = True
        self.validation_errors: List[str] = []
        
        # 兼容性：保留旧属性
        self.selected = False
        self.name = ""
        
    def __str__(self):
        return f"DatasetItem(id={self.id[:8]}..., type={self.item_type}, name={self.loci_name})"
        
    def __repr__(self):
        return self.__str__()
    
    def update_modified_time(self):
        """更新修改时间"""
        self.modified_at = datetime.now()
        
    def set_name(self, name: str):
        """设置名称"""
        self.name = name
        self.loci_name = name
        self.update_modified_time()
        
    def get_name(self) -> str:
        """获取名称"""
        return self.name or self.loci_name
        
    def get_dependencies(self) -> List[str]:
        """获取该数据项的依赖数据项ID列表"""
        # 这个方法需要在DatasetManager中实现，因为它需要访问所有数据项
        # 这里只返回空列表作为占位符
        return []
        
    def get_dependents(self) -> List[str]:
        """获取依赖于该数据项的其他数据项ID列表"""
        # 这个方法需要在DatasetManager中实现，因为它需要访问所有数据项
        # 这里只返回空列表作为占位符
        return []
        
    def is_compatible_with(self, other: 'DatasetItem') -> bool:
        """检查与另一个数据项是否兼容"""
        # 检查是否属于同一数据集
        if self.dataset_id != other.dataset_id:
            return False
        
        # 检查数据是否匹配（例如：序列数量一致）
        if self.sequence_count > 0 and other.sequence_count > 0:
            if self.sequence_count != other.sequence_count:
                return False
        
        return True
    
    def validate(self) -> bool:
        """验证数据项的有效性"""
        self.validation_errors = []
        
        # 基本验证
        if not self.id:
            self.validation_errors.append("Missing ID")
        
        if not self.item_type or self.item_type not in ALL_ITEM_TYPES:
            self.validation_errors.append(f"Invalid item type: {self.item_type}")
        
        if not self.loci_name:
            self.validation_errors.append("Missing loci name")
        
        # 类型特定验证
        if self.item_type in [ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT]:
            if not self.sequences:
                self.validation_errors.append(f"{self.item_type} has no sequences")
            elif self.sequence_count != len(self.sequences):
                self.validation_errors.append(f"Sequence count mismatch: declared={self.sequence_count}, actual={len(self.sequences)}")
        
        if self.item_type == ITEM_TYPE_ALIGNMENT and not self.is_aligned:
            self.validation_errors.append("Alignment marked as not aligned")
        
        self.is_valid = len(self.validation_errors) == 0
        return self.is_valid
        
    def to_dict(self) -> Dict[str, Any]:
        """转换为字典格式（用于序列化）"""
        return {
            "id": self.id,
            "dataset_id": self.dataset_id,
            "loci_name": self.loci_name,
            "file_path": self.file_path,
            "item_type": self.item_type,
            "parent_id": self.parent_id,
            "created_at": self.created_at.isoformat(),
            "modified_at": self.modified_at.isoformat(),
            "is_aligned": self.is_aligned,
            "length": self.length,
            "sequence_count": self.sequence_count,
            "format": self.format,
            "data": self.data,
            "metadata": self.metadata,
            "selection_state": self.selection_state,
            "selection_reason": self.selection_reason,
            "ownership_uuid": self.ownership_uuid,
            "is_valid": self.is_valid,
            "validation_errors": self.validation_errors
        }
        
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'DatasetItem':
        """从字典格式创建实例（用于反序列化）"""
        item = cls(item_type=data.get("item_type", ITEM_TYPE_SEQUENCE))
        item.id = data.get("id", str(uuid.uuid4()))
        item.dataset_id = data.get("dataset_id", "")
        item.loci_name = data.get("loci_name", "")
        item.file_path = data.get("file_path", "")
        item.parent_id = data.get("parent_id")
        
        if "created_at" in data:
            item.created_at = datetime.fromisoformat(data["created_at"])
        if "modified_at" in data:
            item.modified_at = datetime.fromisoformat(data["modified_at"])
        
        item.is_aligned = data.get("is_aligned", False)
        item.length = data.get("length", 0)
        item.sequence_count = data.get("sequence_count", 0)
        item.format = data.get("format", "")
        item.data = data.get("data", {})
        item.metadata = data.get("metadata", {})
        item.selection_state = data.get("selection_state", SELECTION_STATE_NONE)
        item.selection_reason = data.get("selection_reason", "")
        item.ownership_uuid = data.get("ownership_uuid", "")
        item.is_valid = data.get("is_valid", True)
        item.validation_errors = data.get("validation_errors", [])
        
        # 兼容性
        item.name = item.loci_name
        
        return item


class DatasetInfo:
    """数据集信息"""
    
    def __init__(self):
        self.id: str = str(uuid.uuid4())
        self.name: str = ""
        self.description: str = ""
        self.created_at: datetime = datetime.now()
        self.modified_at: datetime = datetime.now()
        self.is_multigene: bool = False
        self.items: List[str] = []  # 包含的数据项ID列表
        self.settings: Dict[str, Any] = {}
        
        # 多基因数据集特有属性
        self.partition_scheme: Optional[str] = None
        self.loci_count: int = 0
        self.taxa_count: int = 0
        
        # 选择状态（用于高亮显示）
        self.selection_state: str = SELECTION_STATE_NONE
        
    def update_modified_time(self):
        """更新修改时间"""
        self.modified_at = datetime.now()
        
    def add_item(self, item_id: str):
        """添加数据项ID"""
        if item_id not in self.items:
            self.items.append(item_id)
            self.update_modified_time()
            
    def remove_item(self, item_id: str):
        """移除数据项ID"""
        if item_id in self.items:
            self.items.remove(item_id)
            self.update_modified_time()
            
    def get_items_by_type(self, all_items: Dict[str, DatasetItem], item_type: str) -> List[DatasetItem]:
        """获取指定类型的所有数据项"""
        return [all_items[item_id] for item_id in self.items 
                if item_id in all_items and all_items[item_id].item_type == item_type]
        
    def to_dict(self) -> Dict[str, Any]:
        """转换为字典格式"""
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "created_at": self.created_at.isoformat(),
            "modified_at": self.modified_at.isoformat(),
            "is_multigene": self.is_multigene,
            "items": self.items,
            "settings": self.settings,
            "partition_scheme": self.partition_scheme,
            "loci_count": self.loci_count,
            "taxa_count": self.taxa_count
        }
        
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'DatasetInfo':
        """从字典格式创建实例"""
        info = cls()
        info.id = data.get("id", str(uuid.uuid4()))
        info.name = data.get("name", "")
        info.description = data.get("description", "")
        
        if "created_at" in data:
            info.created_at = datetime.fromisoformat(data["created_at"])
        if "modified_at" in data:
            info.modified_at = datetime.fromisoformat(data["modified_at"])
        
        info.is_multigene = data.get("is_multigene", False)
        info.items = data.get("items", [])
        info.settings = data.get("settings", {})
        info.partition_scheme = data.get("partition_scheme")
        info.loci_count = data.get("loci_count", 0)
        info.taxa_count = data.get("taxa_count", 0)
        
        return info


class SelectionNode:
    """选择树节点"""
    
    def __init__(self, item_id: str):
        self.item_id: str = item_id
        self.children: List['SelectionNode'] = []
        self.parent: Optional['SelectionNode'] = None
        self.state: str = SELECTION_STATE_NONE
        
    def add_child(self, child: 'SelectionNode'):
        """添加子节点"""
        child.parent = self
        self.children.append(child)
        
    def remove_child(self, child: 'SelectionNode'):
        """移除子节点"""
        if child in self.children:
            child.parent = None
            self.children.remove(child)
        
    def get_leaf_nodes(self) -> List['SelectionNode']:
        """获取所有叶子节点"""
        leaves = []
        if not self.children:
            leaves.append(self)
        else:
            for child in self.children:
                leaves.extend(child.get_leaf_nodes())
        return leaves
        
    def get_all_nodes(self) -> List['SelectionNode']:
        """获取所有节点（包括自己和所有后代）"""
        nodes = [self]
        for child in self.children:
            nodes.extend(child.get_all_nodes())
        return nodes
        
    def get_latest_child_by_type(self, item_type: str, all_items: Dict[str, DatasetItem]) -> Optional['SelectionNode']:
        """获取指定类型的最新子节点"""
        # 筛选出指定类型的子节点
        typed_children = [child for child in self.children 
                         if child.item_id in all_items and all_items[child.item_id].item_type == item_type]
        
        if not typed_children:
            return None
        
        # 按创建时间排序，返回最新的
        typed_children.sort(key=lambda node: all_items[node.item_id].created_at, reverse=True)
        return typed_children[0]
        
    def get_depth(self) -> int:
        """获取节点深度（根节点为0）"""
        depth = 0
        current = self.parent
        while current:
            depth += 1
            current = current.parent
        return depth
        
    def __str__(self):
        return f"SelectionNode(item_id={self.item_id[:8]}..., state={self.state}, children={len(self.children)})"


class ChainItem(DatasetItem):
    """MCMC链文件数据项"""
    
    def __init__(self, file_paths: list = None, run_number: int = 1, chain_count: int = 1, tool: str = ""):
        super().__init__(item_type=ITEM_TYPE_CHAIN)
        
        # MCMC链特定属性
        self.file_paths = file_paths if file_paths else []  # 链文件路径列表（多个.p文件）
        self.run_number = run_number                        # 运行编号（run1, run2...）
        self.chain_count = chain_count                      # 链数量
        self.tool = tool                                    # 生成工具（mrbayes/phylobayes）
        
        # 自动设置loci_name
        self.loci_name = f"{tool}_run{run_number}"
        
        # 存储链的统计信息
        self.generation_count = 0           # 代数
        self.burnin = 0                     # burn-in值
        self.mean_likelihood = 0.0          # 平均似然值
        self.convergence_check = False      # 是否收敛
        
    def __str__(self):
        return f"ChainItem(run={self.run_number}, chains={self.chain_count}, tool={self.tool}, files={len(self.file_paths)})"
    
    def __repr__(self):
        return self.__str__()