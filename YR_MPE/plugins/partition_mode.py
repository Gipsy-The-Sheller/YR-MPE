# partition_mode.py
#
# Copyright (c) 2026 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# YR-MPE program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
分区模式共享组件
提供分区定义、模式转换和文件导入导出功能
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
from enum import Enum
import os
import re
import tempfile


class PartitionMode(Enum):
    """分区模式枚举"""
    EL = "EL"  # Edge-linked Proportional (-p)
    TL = "TL"  # Edge-linked Equal (-q)
    EUL = "EUL"  # Edge-unlinked (-Q)
    TUL = "TUL"  # Separate Tree Inference (-S)
    
    def get_iqtree_parameter(self) -> str:
        """获取IQ-TREE对应的命令行参数"""
        mapping = {
            PartitionMode.EL: "-p",
            PartitionMode.TL: "-q",
            PartitionMode.EUL: "-Q",
            PartitionMode.TUL: "-S"
        }
        return mapping[self]
    
    def get_description(self) -> str:
        """获取模式描述"""
        descriptions = {
            PartitionMode.EL: "Edge-linked Proportional - 共享拓扑和共享分支长度（按比例缩放）",
            PartitionMode.TL: "Edge-linked Equal - 共享拓扑和相等分支长度",
            PartitionMode.EUL: "Edge-unlinked - 共享拓扑和独立分支长度",
            PartitionMode.TUL: "Separate Tree Inference - 独立拓扑和独立分支长度"
        }
        return descriptions[self]


@dataclass
class PartitionDefinition:
    """分区定义数据类"""
    name: str  # 分区名称
    file_path: str  # 文件路径（如果从单个文件导入不同位点范围，则为空字符串）
    seq_type: str  # 序列类型 (DNA, AA, CODON, MORPH)
    model_range: str  # 模型范围 (如 "1-1000" 或 "gene1.fas:1-1000")
    selected_model: Optional[str] = None  # 选定模型
    
    def __post_init__(self):
        """初始化后处理"""
        # 确保seq_type是大写的
        self.seq_type = self.seq_type.upper()
        
        # 处理文件路径为空的情况
        if self.file_path is None:
            self.file_path = ""
    
    def get_display_range(self) -> str:
        """获取显示用的范围字符串"""
        if self.file_path:
            # 如果有文件路径，显示 "file:range"
            filename = os.path.basename(self.file_path)
            return f"{filename}:{self.model_range}"
        else:
            # 否则只显示范围
            return self.model_range
    
    def validate(self) -> Tuple[bool, str]:
        """验证分区定义的有效性"""
        # 检查名称
        if not self.name or not self.name.strip():
            return False, "分区名称不能为空"
        
        # 检查序列类型
        valid_seq_types = ["DNA", "AA", "CODON", "MORPH", "STANDARD"]
        if self.seq_type not in valid_seq_types:
            return False, f"无效的序列类型: {self.seq_type}，必须是 {', '.join(valid_seq_types)} 之一"
        
        # 检查模型范围
        if not self.model_range or not self.model_range.strip():
            return False, "模型范围不能为空"
        
        # 检查模型范围格式
        if ":" in self.model_range:
            # 跨文件格式: "file:range"
            parts = self.model_range.split(":", 1)
            if len(parts) != 2:
                return False, f"无效的模型范围格式: {self.model_range}"
            range_part = parts[1]
        else:
            range_part = self.model_range
        
        # 验证范围格式
        range_pattern = r'^(\d+-\d+|\d+)$'
        if not re.match(range_pattern, range_part):
            return False, f"无效的范围格式: {range_part}，应该是 '1-1000' 或 '1000' 格式"
        
        # 检查文件路径（如果提供了）
        if self.file_path and not os.path.exists(self.file_path):
            return False, f"文件不存在: {self.file_path}"
        
        return True, ""


class PartitionModeConverter:
    """分区模式转换器"""
    
    # MrBayes支持的分区模式
    MRBAYES_SUPPORTED = [PartitionMode.EL, PartitionMode.EUL]
    
    @classmethod
    def convert_to_mrbayes(cls, source_mode: PartitionMode) -> Tuple[PartitionMode, List[str]]:
        """
        将分区模式转换为MrBayes兼容的格式
        
        Args:
            source_mode: 源分区模式
            
        Returns:
            tuple: (目标模式, 警告信息列表)
        """
        warnings = []
        
        if source_mode not in PartitionMode:
            return PartitionMode.EL, [f"⚠️ 未知分区模式: {source_mode}"]
        
        if source_mode in cls.MRBAYES_SUPPORTED:
            # 完全支持，无需转换
            return source_mode, warnings
        else:
            # 需要转换
            target_mode = PartitionMode.EL
            
            if source_mode == PartitionMode.TUL:
                warnings.append("⚠️ 转换警告: MrBayes不支持拓扑解链（TUL）模式")
                warnings.append("⚠️ 转换说明: 已转换为边解链（EUL）模式")
                warnings.append("⚠️ 影响范围: 所有分区将共享相同的拓扑结构")
                warnings.append("⚠️ 建议: 如果确实需要独立拓扑，请分别分析每个分区")
                warnings.append("⚠️ 替代方案: 使用IQ-TREE的-S模式进行独立树推断")
                target_mode = PartitionMode.EUL
            elif source_mode == PartitionMode.TL:
                warnings.append("⚠️ 转换说明: TL模式转换为EL模式")
                warnings.append("⚠️ 影响: 分支长度将从相等变为按比例缩放")
                warnings.append("⚠️ 结果: 不同分区可能有不同的总体进化速率")
            
            return target_mode, warnings
    
    @classmethod
    def is_supported(cls, mode: PartitionMode, tool: str) -> bool:
        """
        检查工具是否支持指定的分区模式
        
        Args:
            mode: 分区模式
            tool: 工具名称 (iqtree/mrbayes/modelfinder)
            
        Returns:
            bool: 是否支持
        """
        if tool == 'mrbayes':
            return mode in cls.MRBAYES_SUPPORTED
        elif tool in ['iqtree', 'modelfinder']:
            return mode in PartitionMode
        return False
    
    @classmethod
    def get_supported_modes(cls, tool: str) -> List[PartitionMode]:
        """获取工具支持的所有分区模式"""
        if tool == 'mrbayes':
            return cls.MRBAYES_SUPPORTED
        elif tool in ['iqtree', 'modelfinder']:
            return list(PartitionMode)
        return []


class PartitionFileIO:
    """分区文件导入导出处理器"""
    
    @staticmethod
    def parse_nexus_partition(file_path: str) -> List[PartitionDefinition]:
        """解析NEXUS格式分区文件（支持完整格式）"""
        partitions = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 检查是否是NEXUS格式
            if not content.lower().startswith('#nexus'):
                raise ValueError("文件不是NEXUS格式")
            
            # 查找sets块
            sets_pattern = r'begin\s+sets\s*;(.*?)end\s*;'
            sets_match = re.search(sets_pattern, content, re.IGNORECASE | re.DOTALL)
            
            if not sets_match:
                raise ValueError("未找到sets块")
            
            sets_content = sets_match.group(1)
            
            # 解析charset定义
            charset_pattern = r'charset\s+(\w+)\s*=\s*([^;]+);'
            charset_matches = re.findall(charset_pattern, sets_content, re.IGNORECASE)
            
            # 解析charpartition命令（如果存在）
            charpartition_pattern = r'charpartition\s+\w+\s*=\s*([^;]+);'
            charpartition_match = re.search(charpartition_pattern, sets_content, re.IGNORECASE)
            
            partition_models = {}
            if charpartition_match:
                # 解析每个分区的模型定义
                model_defs = charpartition_match.group(1)
                for model_def in model_defs.split(','):
                    model_def = model_def.strip()
                    if ':' in model_def:
                        model, partition_name = model_def.split(':', 1)
                        partition_models[partition_name.strip()] = model.strip()
            
            # 获取文件所在目录（用于相对路径解析）
            file_dir = os.path.dirname(file_path)
            
            for name, range_def in charset_matches:
                name = name.strip()
                range_def = range_def.strip()
                
                # 检测是否包含文件路径（格式: filename:range）
                file_path = ""
                seq_type = "DNA"  # 默认DNA类型
                if ':' in range_def and not re.match(r'^\d', range_def):
                    # 可能是文件路径格式
                    parts = range_def.split(':', 1)
                    if len(parts) == 2:
                        potential_file = parts[0].strip()
                        range_def = parts[1].strip()
                        
                        # 检查文件是否存在（支持相对路径）
                        if os.path.exists(potential_file):
                            file_path = os.path.abspath(potential_file)
                        elif os.path.exists(os.path.join(file_dir, potential_file)):
                            file_path = os.path.abspath(os.path.join(file_dir, potential_file))
                        else:
                            # 如果文件不存在，仍然保留路径，让用户处理
                            file_path = potential_file
                
                # 获取该分区的模型（如果定义了）
                selected_model = partition_models.get(name, None)
                
                # 尝试从charpartition推断序列类型
                if selected_model:
                    if selected_model.startswith("LG") or selected_model.startswith("WAG") or selected_model.startswith("JTT"):
                        seq_type = "AA"
                    elif selected_model.startswith("DNA") or selected_model.startswith("GTR"):
                        seq_type = "DNA"
                
                partition = PartitionDefinition(
                    name=name,
                    file_path=file_path,
                    seq_type=seq_type,
                    model_range=range_def,
                    selected_model=selected_model
                )
                partitions.append(partition)
            
            if not partitions:
                raise ValueError("未找到有效的分区定义")
            
            return partitions
            
        except Exception as e:
            raise ValueError(f"解析NEXUS分区文件失败: {str(e)}")
    
    @staticmethod
    def parse_raxml_partition(file_path: str) -> List[PartitionDefinition]:
        """解析RAxML格式分区文件"""
        partitions = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    # 格式: DATATYPE, NAME = RANGE
                    parts = line.split(',', 2)
                    if len(parts) == 3:
                        seq_type = parts[0].strip()
                        name_range = parts[1].split('=', 1)
                        if len(name_range) == 2:
                            name = name_range[0].strip()
                            range_def = name_range[1].strip()
                            
                            partition = PartitionDefinition(
                                name=name,
                                file_path="",
                                seq_type=seq_type.upper(),
                                model_range=range_def
                            )
                            partitions.append(partition)
            
            if not partitions:
                raise ValueError("未找到有效的分区定义")
            
            return partitions
            
        except Exception as e:
            raise ValueError(f"解析RAxML分区文件失败: {str(e)}")

    @staticmethod
    def parse_nexus_full(file_path: str) -> Tuple[Dict, List[PartitionDefinition], Dict]:
        """
        完整解析 NEXUS 文件（序列+分区）

        Args:
            file_path: NEXUS 文件路径

        Returns:
            tuple: (sequences: Dict[str, str],
                    partitions: List[PartitionDefinition],
                    metadata: Dict)

        Raises:
            FileNotFoundError: 文件不存在
            ValueError: 文件格式错误或解析失败
        """
        import sys
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from YR_MPE.platforms.methods.import_partitioned_nexus import (
            parse_charset_range, import_partitioned_nexus
        )

        # 使用现有的 import_partitioned_nexus 函数获取数据
        try:
            # 直接调用 import_partitioned_nexus 获取数据
            # 注意：这需要处理 DatasetItem 的转换
            from YR_MPE.platforms.methods.dataset_manager import DatasetItem
            dataset_items, partition_scheme, summary = import_partitioned_nexus(file_path)

            # 提取序列数据
            sequences = {}
            if dataset_items:
                for item in dataset_items:
                    # 这里需要处理多个位点的情况
                    # 由于 import_partitioned_nexus 返回的是拆分后的位点，
                    # 我们需要重新组合它们
                    pass

            # 重新解析原始文件以获取完整的序列数据
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # 解析 data 块
            data_pattern = r'begin\s+data\s*;(.*?)end\s*;'
            data_match = re.search(data_pattern, content, re.IGNORECASE | re.DOTALL)

            if not data_match:
                raise ValueError("未找到 data 块")

            data_content = data_match.group(1)

            # 提取 datatype
            datatype_pattern = r'datatype=(\w+)'
            datatype_match = re.search(datatype_pattern, data_content, re.IGNORECASE)
            datatype = datatype_match.group(1).upper() if datatype_match else 'DNA'

            # 提取 matrix
            matrix_pattern = r'matrix\s*(.*?)end\s*;'
            matrix_match = re.search(matrix_pattern, data_content, re.IGNORECASE | re.DOTALL)

            if not matrix_match:
                raise ValueError("未找到 matrix 定义")

            matrix_content = matrix_match.group(1)

            # 解析序列
            sequences = {}
            for line in matrix_content.split('\n'):
                line = line.strip()
                if not line or line.startswith(';'):
                    continue

                parts = line.split(maxsplit=1)
                if len(parts) >= 2:
                    seq_name = parts[0]
                    seq_content = parts[1].replace(' ', '')
                    sequences[seq_name] = seq_content

            # 解析分区定义
            partitions = PartitionFileIO.parse_nexus_partition(file_path)

            # 元数据
            metadata = {
                'source_file': file_path,
                'ntax': len(sequences),
                'nchar': len(list(sequences.values())[0]) if sequences else 0,
                'datatype': datatype,
                'partition_count': len(partitions),
                'partition_scheme': partition_scheme
            }

            return sequences, partitions, metadata

        except Exception as e:
            raise ValueError(f"解析完整 NEXUS 文件失败: {str(e)}")

    @staticmethod
    def export_nexus_partition(partitions: List[PartitionDefinition], file_path: str):
        """导出NEXUS格式分区文件"""
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write("#nexus\n")
                f.write("begin sets;\n")
                
                # 写入charset定义
                for partition in partitions:
                    if partition.file_path:
                        f.write(f"    charset {partition.name} = {partition.file_path}:{partition.model_range};\n")
                    else:
                        f.write(f"    charset {partition.name} = {partition.model_range};\n")
                
                # 写入charpartition（如果有选定模型）
                models = [p for p in partitions if p.selected_model]
                if models:
                    f.write("    charpartition models = ")
                    model_defs = []
                    for partition in partitions:
                        if partition.selected_model:
                            model_defs.append(f"{partition.selected_model}:{partition.name}")
                        else:
                            model_defs.append(f"DNA+G:{partition.name}")
                    f.write(", ".join(model_defs))
                    f.write(";\n")
                
                f.write("end;\n")
                
        except Exception as e:
            raise IOError(f"导出NEXUS分区文件失败: {str(e)}")
    
    @staticmethod
    def export_raxml_partition(partitions: List[PartitionDefinition], file_path: str):
        """导出RAxML格式分区文件"""
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                for partition in partitions:
                    f.write(f"{partition.seq_type}, {partition.name} = {partition.model_range}\n")
                    
        except Exception as e:
            raise IOError(f"导出RAxML分区文件失败: {str(e)}")
    
    @staticmethod
    def validate_partition_file(file_path: str) -> Tuple[bool, str, List[PartitionDefinition]]:
        """
        验证分区文件并返回分区定义
        
        Args:
            file_path: 分区文件路径
            
        Returns:
            tuple: (是否有效, 错误消息, 分区定义列表)
        """
        try:
            # 尝试解析为NEXUS格式
            try:
                partitions = PartitionFileIO.parse_nexus_partition(file_path)
                # 验证每个分区
                for partition in partitions:
                    valid, error = partition.validate()
                    if not valid:
                        return False, f"分区 '{partition.name}' 验证失败: {error}", []
                return True, "", partitions
            except:
                pass
            
            # 尝试解析为RAxML格式
            try:
                partitions = PartitionFileIO.parse_raxml_partition(file_path)
                # 验证每个分区
                for partition in partitions:
                    valid, error = partition.validate()
                    if not valid:
                        return False, f"分区 '{partition.name}' 验证失败: {error}", []
                return True, "", partitions
            except:
                pass
            
            return False, "无法识别的分区文件格式，请使用NEXUS或RAxML格式", []
            
        except Exception as e:
            return False, f"验证分区文件失败: {str(e)}", []


# 便捷函数
def create_partition_from_dict(data: Dict) -> PartitionDefinition:
    """从字典创建分区定义"""
    return PartitionDefinition(
        name=data.get('name', ''),
        file_path=data.get('file_path', ''),
        seq_type=data.get('seq_type', 'DNA'),
        model_range=data.get('model_range', ''),
        selected_model=data.get('selected_model', None)
    )


def partition_to_dict(partition: PartitionDefinition) -> Dict:
    """将分区定义转换为字典"""
    return {
        'name': partition.name,
        'file_path': partition.file_path,
        'seq_type': partition.seq_type,
        'model_range': partition.model_range,
        'selected_model': partition.selected_model
    }


# ==================== MrBayes 分区相关类 ====================

@dataclass
class MrBayesPartitionDefinition:
    """MrBayes分区定义"""
    name: str              # 分区名称（如subset1）
    range: str             # 位点范围（如1-675\3 3784-4893\3）
    seq_type: str          # 序列类型：DNA 或 Protein
    
    # DNA模型参数
    nst: int = 6           # 替换模型数（1-6），仅DNA使用
    
    # Protein模型参数
    aamodel: str = "mixed" # 氨基酸替代模型，仅Protein使用
    
    # 通用参数
    rates: str = "gamma"           # 速率异质性
    ngammacat: int = 4             # Gamma类别数
    
    def __post_init__(self):
        """初始化后处理"""
        self.seq_type = self.seq_type.upper()
        
    def get_model_display(self) -> str:
        """获取模型显示名称"""
        if self.seq_type == "DNA":
            nst_names = {1: "JC69", 2: "HKY85", 6: "GTR"}
            return f"{nst_names.get(self.nst, f'nst={self.nst}')}"
        else:
            return self.aamodel
    
    def validate(self) -> Tuple[bool, str]:
        """验证分区定义的有效性"""
        if not self.name or not self.name.strip():
            return False, "分区名称不能为空"
        
        if self.seq_type not in ["DNA", "PROTEIN"]:
            return False, f"无效的序列类型: {self.seq_type}，必须是 DNA 或 PROTEIN"
        
        if self.seq_type == "DNA":
            if self.nst not in [1, 2, 6]:
                return False, f"无效的nst值: {self.nst}，DNA必须是1、2或6"
        else:
            valid_models = ["BLOSUM62", "BLOSUM", "WAG", "LG", "GTR", "JONES", "MTREV", "POISSON", "MIXED"]
            if self.aamodel.upper() not in valid_models:
                return False, f"无效的氨基酸模型: {self.aamodel}"
        
        if not self.range or not self.range.strip():
            return False, "位点范围不能为空"
        
        return True, ""


class MrBayesModelConverter:
    """MrBayes模型转换器"""
    
    # DNA nst映射（简化版本）
    DNA_NST_MODELS = {
        1: "JC69",
        2: "HKY85", 
        6: "GTR"
    }
    
    # ModelFinder核苷酸模型到MrBayes的完整映射
    # 映射格式: "Model": {"nst": int, "statefreq": str}
    # statefreq可选值: "fixed(equal)", "fixed(empirical)", "dirichlet"
    MODELFINDER_TO_MRBAYES_DNA = {
        # 等速率模型 (nst=1)
        "JC": {"nst": 1, "statefreq": "fixed(equal)"},
        "JC69": {"nst": 1, "statefreq": "fixed(equal)"},
        "F81": {"nst": 1, "statefreq": "dirichlet"},
        
        # 双参数模型 (nst=2)
        "K80": {"nst": 2, "statefreq": "fixed(equal)"},
        "K2P": {"nst": 2, "statefreq": "fixed(equal)"},
        "TNe": {"nst": 2, "statefreq": "fixed(equal)"},
        "K81": {"nst": 2, "statefreq": "fixed(equal)"},
        "K3P": {"nst": 2, "statefreq": "fixed(equal)"},
        "TPM2": {"nst": 2, "statefreq": "fixed(equal)"},
        "TPM3": {"nst": 2, "statefreq": "fixed(equal)"},
        
        # 三参数模型 (nst=2, 但有频率)
        "K81u": {"nst": 2, "statefreq": "dirichlet"},
        "K3Pu": {"nst": 2, "statefreq": "dirichlet"},
        "TPM2u": {"nst": 2, "statefreq": "dirichlet"},
        "TPM3u": {"nst": 2, "statefreq": "dirichlet"},
        
        # 四参数模型 (nst=6)
        "HKY": {"nst": 6, "statefreq": "dirichlet"},
        "HKY85": {"nst": 6, "statefreq": "dirichlet"},
        "TN": {"nst": 6, "statefreq": "dirichlet"},
        "TN93": {"nst": 6, "statefreq": "dirichlet"},
        "TIM": {"nst": 6, "statefreq": "dirichlet"},
        "TIM2": {"nst": 6, "statefreq": "dirichlet"},
        "TIM3": {"nst": 6, "statefreq": "dirichlet"},
        "TVM": {"nst": 6, "statefreq": "dirichlet"},
        
        # 等频率的4参数模型 (nst=6, equal frequencies)
        "TIMe": {"nst": 6, "statefreq": "fixed(equal)"},
        "TIM2e": {"nst": 6, "statefreq": "fixed(equal)"},
        "TIM3e": {"nst": 6, "statefreq": "fixed(equal)"},
        "TVMe": {"nst": 6, "statefreq": "fixed(equal)"},
        "SYM": {"nst": 6, "statefreq": "fixed(equal)"},
        
        # 完整模型 (nst=6, dirichlet)
        "GTR": {"nst": 6, "statefreq": "dirichlet"},
    }
    
    # Protein模型映射
    PROTEIN_MODELS = [
        "BLOSUM62",
        "BLOSUM",
        "WAG",
        "LG",
        "GTR",
        "JONES",
        "MTREV",
        "POISSON",
        "MIXED"
    ]
    
    # ModelFinder蛋白质模型到MrBayes的映射
    MODELFINDER_TO_MRBAYES_PROTEIN = {
        "BLOSUM62": "blosum",
        "BLOSUM": "blosum",
        "WAG": "wag",
        "LG": "lg",
        "JTT": "jones",
        "JONES": "jones",
        "GTR": "gtr",
        "MTREV": "mtrev",
        "POISSON": "poisson",
    }
    
    # 速率异质性选项
    RATES_OPTIONS = ["equal", "gamma", "invgamma", "propinv", "lnorm", "adgamma"]
    
    @classmethod
    def parse_modelfinder_model(cls, model_string: str) -> Tuple[str, Dict[str, str]]:
        """
        解析ModelFinder模型字符串，返回基础模型和附加参数
        
        Args:
            model_string: ModelFinder模型字符串（如 "GTR+F+I+G4"）
            
        Returns:
            tuple: (基础模型名称, 附加参数字典)
        """
        # 分割模型字符串
        parts = model_string.split("+")
        base_model = parts[0].upper()
        
        # 解析附加参数
        params = {}
        for part in parts[1:]:
            part_upper = part.upper()
            if part_upper == "I":
                params["invariant"] = "true"
            elif part_upper == "F":
                params["freq"] = "empirical"
            elif part_upper == "FO":
                params["freq"] = "ml"
            elif part_upper == "FQ":
                params["freq"] = "equal"
            elif part_upper.startswith("G"):
                try:
                    params["gamma_categories"] = int(part_upper[1:])
                except (ValueError, IndexError):
                    params["gamma_categories"] = 4
            elif part_upper == "R":
                params["freerate"] = "true"
            elif part_upper == "ASC":
                params["ascertainment"] = "true"
        
        return base_model, params
    
    @classmethod
    def convert_dna_model_to_mrbayes(cls, model_string: str) -> Dict[str, any]:
        """
        将ModelFinder DNA模型转换为MrBayes格式
        
        Args:
            model_string: ModelFinder模型字符串（如 "GTR+F+I+G4"）
            
        Returns:
            dict: {"nst": int, "statefreq": str, "rates": str, "ngammacat": int, "use_propinv": bool}
                  如果模型无法转换，返回None
        """
        # 解析模型字符串
        base_model, params = cls.parse_modelfinder_model(model_string)
        
        # 查找基础模型的映射
        if base_model not in cls.MODELFINDER_TO_MRBAYES_DNA:
            return None
        
        mrbayes_params = cls.MODELFINDER_TO_MRBAYES_DNA[base_model].copy()
        
        # 处理频率参数
        # +F: 使用经验频率（fixed(empirical)），覆盖原来的 dirichlet
        # +FO: 估计频率（estimated(dirichlet)，保持不变）
        # +FQ: 等频率（fixed(equal)，保持不变）
        if params.get("freq") == "empirical":  # +F
            mrbayes_params["statefreq"] = "fixed(empirical)"
        # 对于 +FO 和 +FQ，statefreq 保持原样
        
        # 处理速率异质性
        # +G → Gamma (+G)
        # +I → PropInv (+I)
        # +G+I → Invgamma (+G+I)
        rates = "equal"
        ngammacat = 4
        
        if "gamma_categories" in params:  # +G
            rates = "gamma"
            ngammacat = params["gamma_categories"]
        
        if params.get("invariant") == "true":  # +I
            if rates == "gamma":
                # 同时有 +G 和 +I，使用 Invgamma
                rates = "invgamma"
            else:
                # 只有 +I，使用 PropInv
                rates = "propinv"
        
        mrbayes_params.update({
            "rates": rates,
            "ngammacat": ngammacat
        })
        
        return mrbayes_params
    
    @classmethod
    def convert_protein_model_to_mrbayes(cls, model_string: str) -> Tuple[str, List[str]]:
        """
        将ModelFinder蛋白质模型转换为MrBayes格式
        
        Args:
            model_string: ModelFinder模型字符串
            
        Returns:
            tuple: (MrBayes模型名称, 警告信息列表)
                  如果模型无法直接匹配，返回 (None, 警告信息)
        """
        # 解析模型字符串
        base_model, _ = cls.parse_modelfinder_model(model_string)
        
        # 尝试匹配MrBayes支持的模型
        if base_model in cls.MODELFINDER_TO_MRBAYES_PROTEIN:
            return cls.MODELFINDER_TO_MRBAYES_PROTEIN[base_model], []
        
        # 无法匹配，返回警告
        warnings = [
            f"⚠️ 模型警告: MrBayes不支持 '{base_model}' 蛋白质模型",
            f"⚠️ MrBayes支持的蛋白质模型: {', '.join(cls.MODELFINDER_TO_MRBAYES_PROTEIN.keys())}",
            "⚠️ 请选择以下选项之一:"
        ]
        
        return None, warnings
    
    @classmethod
    def convert_model_to_mrbayes(cls, model_string: str, seq_type: str = "DNA") -> Tuple[Dict, List[str]]:
        """
        将ModelFinder模型转换为MrBayes格式
        
        Args:
            model_string: ModelFinder模型字符串
            seq_type: 序列类型 ("DNA" 或 "PROTEIN")
            
        Returns:
            tuple: (MrBayes参数字典或模型名称, 警告信息列表)
                  对于DNA: 返回包含nst、statefreq等参数的字典
                  对于PROTEIN: 返回模型名称字符串
        """
        warnings = []
        
        if seq_type.upper() == "DNA":
            mrbayes_params = cls.convert_dna_model_to_mrbayes(model_string)
            if mrbayes_params is None:
                warnings.append(f"⚠️ 模型警告: 无法将 '{model_string}' 转换为MrBayes格式")
                warnings.append("⚠️ 将使用默认GTR模型")
                mrbayes_params = cls.MODELFINDER_TO_MRBAYES_DNA["GTR"].copy()
                mrbayes_params.update({"rates": "gamma", "ngammacat": 4, "use_propinv": False})
            return mrbayes_params, warnings
        
        elif seq_type.upper() == "PROTEIN":
            mrbayes_model, model_warnings = cls.convert_protein_model_to_mrbayes(model_string)
            if model_warnings:
                warnings.extend(model_warnings)
            return mrbayes_model, warnings
        
        else:
            warnings.append(f"⚠️ 错误: 不支持的序列类型 '{seq_type}'")
            return None, warnings
    
    @classmethod
    def generate_partition_commands(cls, partitions: List[MrBayesPartitionDefinition], 
                                      mode: PartitionMode = PartitionMode.EL) -> List[str]:
        """
        生成MrBayes分区命令
        
        Args:
            partitions: 分区定义列表
            mode: 分区模式（EL=EUL, TL=EL, TUL=EUL, EUL=EUL）
            
        Returns:
            命令列表
        """
        commands = []
        
        if not partitions:
            return commands
        
        # 生成charset定义
        for i, partition in enumerate(partitions, 1):
            commands.append(f"charset {partition.name} = {partition.range};")
        
        # 生成partition定义
        partition_names = ", ".join([p.name for p in partitions])
        commands.append(f"partition Names = {len(partitions)}:{partition_names};")
        commands.append("set partition=Names;")
        
        # 根据序列类型分组处理
        dna_partitions = [p for p in partitions if p.seq_type == "DNA"]
        protein_partitions = [p for p in partitions if p.seq_type == "PROTEIN"]
        
        # 处理DNA分区
        if dna_partitions:
            dna_indices = [i+1 for i, p in enumerate(partitions) if p.seq_type == "DNA"]
            dna_indices_str = ",".join([f"({i})" for i in dna_indices])
            commands.append(f"lset applyto=({dna_indices_str}) nucmodel=4by4;")
        
        # 处理Protein分区
        if protein_partitions:
            prot_indices = [i+1 for i, p in enumerate(partitions) if p.seq_type == "PROTEIN"]
            prot_indices_str = ",".join([f"({i})" for i in prot_indices])
            commands.append(f"lset applyto=({prot_indices_str}) nucmodel=protein;")
        
        # 为每个分区设置具体模型参数
        for i, partition in enumerate(partitions, 1):
            if partition.seq_type == "DNA":
                # DNA: 设置nst
                commands.append(f"lset applyto=({i}) nst={partition.nst};")
            else:
                # Protein: 设置aamodel
                if partition.aamodel.upper() == "MIXED":
                    commands.append(f"prset applyto=({i}) aamodelpr=mixed;")
                else:
                    commands.append(f"prset applyto=({i}) aamodelpr=fixed({partition.aamodel.lower()});")
            
            # 设置速率异质性
            if partition.rates == "equal":
                commands.append(f"lset applyto=({i}) rates=equal;")
            elif partition.rates in ["gamma", "invgamma", "lnorm", "adgamma"]:
                commands.append(f"lset applyto=({i}) rates={partition.rates} ngammacat={partition.ngammacat};")
            elif partition.rates == "propinv":
                commands.append(f"lset applyto=({i}) rates=propinv;")
        
        # Edge-unlinked模式：添加unlink命令
        if mode == PartitionMode.EUL or mode == PartitionMode.TUL:
            commands.append("unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all) tratio=(all);")
            commands.append("unlink brlens=(all);")
        
        return commands
    
    @classmethod
    def parse_mrbayes_block(cls, mrbayes_block: str) -> Tuple[List[MrBayesPartitionDefinition], PartitionMode]:
        """
        解析MrBayes数据块中的分区定义
        
        Args:
            mrbayes_block: MrBayes数据块内容
            
        Returns:
            (分区定义列表, 分区模式)
        """
        partitions = []
        mode = PartitionMode.EL  # 默认Edge-linked
        
        lines = mrbayes_block.split('\n')
        
        # 解析charset定义
        charset_pattern = r'charset\s+(\w+)\s*=\s*([^;]+);'
        charsets = {}
        for line in lines:
            match = re.match(charset_pattern, line.strip(), re.IGNORECASE)
            if match:
                name = match.group(1)
                range_def = match.group(2).strip()
                charsets[name] = range_def
        
        # 解析partition定义
        partition_pattern = r'partition\s+(\w+)\s*=\s*(\d+):([^;]+);'
        partition_match = None
        for line in lines:
            match = re.match(partition_pattern, line.strip(), re.IGNORECASE)
            if match:
                partition_match = match
                break
        
        if partition_match:
            partition_names = [name.strip() for name in partition_match.group(3).split(',')]
            for name in partition_names:
                if name in charsets:
                    partitions.append(MrBayesPartitionDefinition(
                        name=name,
                        range=charsets[name],
                        seq_type="DNA",  # 默认DNA，需要后续推断
                        nst=6,
                        rates="gamma",
                        ngammacat=4
                    ))
        
        # 检测是否为Edge-unlinked模式
        if "unlink" in mrbayes_block.lower():
            mode = PartitionMode.EUL
        
        return partitions, mode


def mrbayes_partition_to_dict(partition: MrBayesPartitionDefinition) -> Dict:
    """将MrBayes分区定义转换为字典"""
    return {
        'name': partition.name,
        'range': partition.range,
        'seq_type': partition.seq_type,
        'nst': partition.nst,
        'aamodel': partition.aamodel,
        'rates': partition.rates,
        'ngammacat': partition.ngammacat
    }


def mrbayes_partition_from_dict(data: Dict) -> MrBayesPartitionDefinition:
    """从字典创建MrBayes分区定义"""
    return MrBayesPartitionDefinition(
        name=data.get('name', ''),
        range=data.get('range', ''),
        seq_type=data.get('seq_type', 'DNA'),
        nst=data.get('nst', 6),
        aamodel=data.get('aamodel', 'mixed'),
        rates=data.get('rates', 'gamma'),
        ngammacat=data.get('ngammacat', 4)
    )


# ==================== Dataset 到 ModelFinder 分区转换 ====================

def dataset_to_model_finder_partitions(dataset_info, dataset_items: List) -> Tuple[List[PartitionDefinition], PartitionMode, Dict]:
    """
    将 dataset 转换为 ModelFinder 分区定义
    
    Args:
        dataset_info: DatasetInfo 对象，包含 partition_scheme 和 settings
        dataset_items: List[DatasetItem]，包含序列数据
    
    Returns:
        tuple: (分区定义列表, 分区模式, 元数据字典)
    
    Raises:
        ValueError: 当 dataset_info 或 partition_scheme 无效时
    """
    import sys
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    
    from YR_MPE.platforms.methods.dataset_models import DatasetItem
    
    # 验证输入
    if not dataset_info:
        raise ValueError("dataset_info 不能为空")
    
    if not dataset_items:
        raise ValueError("dataset_items 不能为空")
    
    # 1. 解析 partition_scheme（NEXUS 格式）
    partition_scheme = dataset_info.partition_scheme
    if not partition_scheme:
        # 如果没有 partition_scheme，创建单个分区
        if dataset_items:
            first_item = dataset_items[0]
            partitions = [PartitionDefinition(
                name=first_item.loci_name or "full",
                file_path=first_item.file_path or "",
                seq_type=_infer_seq_type(first_item),
                model_range=f"1-{first_item.length}"
            )]
        else:
            raise ValueError("既没有 partition_scheme 也没有 dataset_items")
    else:
        # 解析 NEXUS 格式的 partition_scheme
        partitions = _parse_partition_scheme_to_definitions(partition_scheme, dataset_items)
    
    # 2. 从 settings 中读取 linkage 设置
    settings = dataset_info.settings if hasattr(dataset_info, 'settings') else {}
    topo_linked = settings.get('topo_linked', False)
    edge_linked = settings.get('edge_linked', False)
    
    # 3. 将 linkage 设置映射到 PartitionMode
    # 映射规则（来自 model_finder_plugin.py:185-198）：
    # - topo_linked=False → TUL (Separate Tree)
    # - topo_linked=True, edge_linked=False → EUL (Edge-unlinked)
    # - topo_linked=True, edge_linked=True → TL (Edge-linked Equal)
    if not topo_linked:
        partition_mode = PartitionMode.TUL  # Separate Tree
    elif edge_linked:
        partition_mode = PartitionMode.TL   # Edge-linked Equal
    else:
        partition_mode = PartitionMode.EUL  # Edge-unlinked
    
    # 4. 为每个分区添加序列信息
    partitions = _enrich_partitions_with_sequences(partitions, dataset_items)
    
    # 5. 构建元数据
    metadata = {
        'dataset_id': dataset_info.id if hasattr(dataset_info, 'id') else '',
        'dataset_name': dataset_info.name if hasattr(dataset_info, 'name') else '',
        'partition_count': len(partitions),
        'topo_linked': topo_linked,
        'edge_linked': edge_linked,
        'partition_mode': partition_mode.value,
        'total_taxa': dataset_items[0].sequence_count if dataset_items else 0
    }
    
    return partitions, partition_mode, metadata


def _parse_partition_scheme_to_definitions(partition_scheme: str, dataset_items: List) -> List[PartitionDefinition]:
    """
    解析 NEXUS 格式的 partition_scheme 并转换为 PartitionDefinition 列表
    
    Args:
        partition_scheme: NEXUS 格式的分区方案字符串
        dataset_items: DatasetItem 列表，用于获取序列类型等信息
    
    Returns:
        PartitionDefinition 列表
    """
    partitions = []
    
    # 解析 charset 定义
    charset_pattern = r'charset\s+(\S+)\s*=\s*([^;]+);'
    charset_matches = re.findall(charset_pattern, partition_scheme, re.IGNORECASE)
    
    # 创建名称到 dataset_item 的映射
    item_map = {item.loci_name: item for item in dataset_items if hasattr(item, 'loci_name')}
    
    for name, range_def in charset_matches:
        name = name.strip()
        range_def = range_def.strip()
        
        # 跳过外部文件引用（暂时不支持）
        if ':' in range_def and not re.match(r'^\d', range_def):
            continue
        
        # 获取对应的 dataset_item 以推断序列类型
        seq_type = "DNA"  # 默认
        if name in item_map:
            seq_type = _infer_seq_type(item_map[name])
        
        # 创建分区定义
        partition = PartitionDefinition(
            name=name,
            file_path="",  # 从同一个文件导入
            seq_type=seq_type,
            model_range=range_def,
            selected_model=None
        )
        
        partitions.append(partition)
    
    # 如果没有找到分区，创建单个分区
    if not partitions and dataset_items:
        first_item = dataset_items[0]
        partitions.append(PartitionDefinition(
            name=first_item.loci_name or "full",
            file_path=first_item.file_path or "",
            seq_type=_infer_seq_type(first_item),
            model_range=f"1-{first_item.length}"
        ))
    
    return partitions


def _enrich_partitions_with_sequences(partitions: List[PartitionDefinition], dataset_items: List) -> List[PartitionDefinition]:
    """
    为分区定义添加序列信息
    
    Args:
        partitions: PartitionDefinition 列表
        dataset_items: DatasetItem 列表
    
    Returns:
        更新后的 PartitionDefinition 列表
    """
    # 创建名称到 dataset_item 的映射
    item_map = {item.loci_name: item for item in dataset_items if hasattr(item, 'loci_name')}
    
    for partition in partitions:
        if partition.name in item_map:
            item = item_map[partition.name]
            # 更新文件路径
            if hasattr(item, 'file_path') and item.file_path:
                partition.file_path = item.file_path
            
            # 如果有选定的模型，添加到分区定义
            if hasattr(item, 'data') and 'selected_model' in item.data:
                partition.selected_model = item.data['selected_model']
    
    return partitions


def _infer_seq_type(dataset_item) -> str:
    """
    从 DatasetItem 推断序列类型
    
    Args:
        dataset_item: DatasetItem 对象
    
    Returns:
        序列类型 (DNA, AA, CODON, MORPH, STANDARD)
    """
    # 检查 metadata
    if hasattr(dataset_item, 'metadata') and 'seq_type' in dataset_item.metadata:
        seq_type = dataset_item.metadata['seq_type'].upper()
        if seq_type in ['DNA', 'AA', 'CODON', 'MORPH', 'STANDARD']:
            return seq_type
    
    # 检查 data 字段
    if hasattr(dataset_item, 'data') and 'seq_type' in dataset_item.data:
        seq_type = dataset_item.data['seq_type'].upper()
        if seq_type in ['DNA', 'AA', 'CODON', 'MORPH', 'STANDARD']:
            return seq_type
    
    # 从序列内容推断
    if hasattr(dataset_item, 'sequences') and dataset_item.sequences:
        first_seq = str(dataset_item.sequences[0].seq).upper()
        
        # 检查是否是氨基酸序列
        aa_chars = set('ACDEFGHIKLMNPQRSTVWY*')
        dna_chars = set('ACGTN-')
        
        # 统计氨基酸字符出现次数
        aa_count = sum(1 for c in first_seq if c in aa_chars)
        # 统计 DNA 字符出现次数
        dna_count = sum(1 for c in first_seq if c in dna_chars)
        
        # 如果氨基酸字符占主导，判定为蛋白质
        total_valid = aa_count + dna_count
        if total_valid > 0 and aa_count / total_valid > 0.5:
            return "AA"
    
    # 默认返回 DNA
    return "DNA"


def get_dataset_partition_sequences(dataset_info, dataset_items: List) -> Dict[str, List]:
    """
    获取每个分区的序列数据
    
    Args:
        dataset_info: DatasetInfo 对象
        dataset_items: List[DatasetItem]
    
    Returns:
        字典: {分区名称: [SeqRecord1, SeqRecord2, ...]}
    """
    sequences_dict = {}
    
    # 创建名称到 dataset_item 的映射
    item_map = {item.loci_name: item for item in dataset_items if hasattr(item, 'loci_name')}
    
    # 解析 partition_scheme 获取分区名称
    partition_scheme = dataset_info.partition_scheme if hasattr(dataset_info, 'partition_scheme') else ""
    if partition_scheme:
        charset_pattern = r'charset\s+(\S+)\s*='
        partition_names = re.findall(charset_pattern, partition_scheme, re.IGNORECASE)
    else:
        partition_names = list(item_map.keys())
    
    # 为每个分区提取序列
    for name in partition_names:
        if name in item_map:
            item = item_map[name]
            if hasattr(item, 'sequences') and item.sequences:
                sequences_dict[name] = item.sequences
    
    return sequences_dict