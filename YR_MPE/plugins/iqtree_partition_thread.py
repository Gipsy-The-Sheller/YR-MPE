# iqtree_partition_thread.py
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
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
IQ-TREE 分区模式系统发育推断线程类
处理分区模式下的IQ-TREE系统发育树构建
"""

import os
import tempfile
from typing import List, Optional, Dict
from Bio import SeqIO

from ..templates.base_process_thread import BaseProcessThread
from .partition_mode import PartitionDefinition, PartitionMode


class IQTreePartitionThread(BaseProcessThread):
    """IQ-TREE 分区系统发育推断线程类"""
    
    def __init__(
        self,
        tool_path: str,
        input_file: str,
        partitions: List[PartitionDefinition],
        mode: PartitionMode,
        model_params: Dict,
        bootstrap_enabled: bool = True,
        bootstrap_replicates: int = 1000,
        sampling_type: str = "SITE",
        output_prefix: Optional[str] = None
    ):
        """
        初始化IQ-TREE分区线程
        
        Args:
            tool_path: IQ-TREE可执行文件路径
            input_file: 输入比对文件路径
            partitions: 分区定义列表
            mode: 分区模式
            model_params: 模型参数字典
            bootstrap_enabled: 是否启用Bootstrap
            bootstrap_replicates: Bootstrap重复次数
            sampling_type: 抽样类型（SITE/GENE/GENESITE）
            output_prefix: 输出文件前缀
        """
        super().__init__(tool_path, [input_file], [], [])
        self.input_file = input_file
        self.partitions = partitions
        self.mode = mode
        self.model_params = model_params
        self.bootstrap_enabled = bootstrap_enabled
        self.bootstrap_replicates = bootstrap_replicates
        self.sampling_type = sampling_type
        self.output_prefix = output_prefix or os.path.splitext(os.path.basename(input_file))[0]
    
    def get_tool_name(self):
        """返回工具名称"""
        return "IQ-TREE Partition Phylogeny"
    
    def execute_commands(self):
        """执行IQ-TREE分区命令"""
        try:
            output_files = []
            
            self.progress.emit("Preparing partition analysis...")
            self.console_output.emit("Starting IQ-TREE partition phylogenetic inference", "info")
            
            # 创建临时分区文件
            partition_file = self.create_temp_partition_file()
            
            # 构建IQ-TREE命令
            cmd = self.build_iqtree_command(partition_file)
            
            self.console_output.emit(f"Executing: {' '.join(cmd)}", "info")
            
            # 执行命令
            result = self.execute_command(cmd)
            
            if result.returncode != 0:
                self.error.emit(f"IQ-TREE execution failed: {result.stderr}")
                return
            
            # 收集输出文件
            output_files = self.collect_output_files()
            
            if not output_files:
                self.error.emit("No output files generated")
                return
            
            self.progress.emit("Phylogenetic inference completed")
            self.console_output.emit(f"Successfully generated {len(output_files)} output file(s)", "info")
            
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"IQ-TREE partition analysis exception: {str(e)}")
    
    def create_temp_partition_file(self) -> str:
        """创建临时分区定义文件"""
        temp_file = self.create_temp_file(suffix='.nex')
        
        with open(temp_file, 'w') as f:
            f.write("#nexus\n")
            f.write("begin sets;\n")
            
            for partition in self.partitions:
                # 写入字符集定义
                if partition.file_path:
                    # 如果有单独的文件，使用文件路径格式
                    f.write(f"    charset {partition.name} = {partition.file_path}: {partition.model_range};\n")
                else:
                    # 否则使用位点范围格式
                    f.write(f"    charset {partition.name} = {partition.model_range};\n")
            
            f.write("end;\n")
        
        return temp_file
    
    def build_iqtree_command(self, partition_file: str) -> List[str]:
        """构建IQ-TREE分区命令"""
        cmd = [self.tool_path, "-s", self.input_file]
        
        # 添加分区模式参数
        mode_flag = self.get_mode_flag()
        cmd.extend([mode_flag, partition_file])
        
        # 添加模型参数
        model = self.build_model_string()
        cmd.extend(["-m", model])
        
        # 添加Bootstrap参数
        if self.bootstrap_enabled:
            cmd.extend(["-bb", str(self.bootstrap_replicates)])
        
        # 添加aBayes
        cmd.append("-abayes")
        
        # 添加线程数
        threads = self.model_params.get('threads', 1)
        if threads > 1:
            cmd.extend(["-nt", str(threads)])
        else:
            cmd.append("-nt AUTO")
        
        # 添加输出前缀
        cmd.extend(["--prefix", self.output_prefix])
        
        # 添加redo参数（重新计算）
        cmd.append("-redo")
        
        return cmd
    
    def get_mode_flag(self) -> str:
        """获取分区模式标志"""
        mode_map = {
            PartitionMode.EL: "-p",   # Edge-linked proportional
            PartitionMode.EUL: "-Q",  # Edge-unlinked
            PartitionMode.TUL: "-S"   # Separate tree (topo unlinked)
        }
        return mode_map.get(self.mode, "-p")
    
    def build_model_string(self) -> str:
        """构建模型字符串"""
        # 如果分区中已指定模型，使用MFP+MERGE进行模型选择和合并
        # 否则使用用户指定的模型
        
        # 检查是否有分区指定了模型
        has_models = any(p.selected_model for p in self.partitions)
        
        if has_models:
            # 如果有指定模型，使用MFP+MERGE
            return "MFP+MERGE"
        else:
            # 否则使用用户指定的模型
            seq_type = self.model_params.get('seq_type', 'AUTO')
            
            if seq_type == "AUTO":
                return "MFP+MERGE"
            else:
                # 构建模型字符串
                model = self.model_params.get('model', 'GTR+G')
                return model
    
    def collect_output_files(self) -> List[str]:
        """收集输出文件"""
        output_files = []
        
        # 查找.treefile文件
        treefile = f"{self.output_prefix}.treefile"
        if os.path.exists(treefile):
            output_files.append(treefile)
        
        # 查找.iqtree文件
        iqtree_file = f"{self.output_prefix}.iqtree"
        if os.path.exists(iqtree_file):
            output_files.append(iqtree_file)
        
        # 查找.log文件
        log_file = f"{self.output_prefix}.log"
        if os.path.exists(log_file):
            output_files.append(log_file)
        
        # 查找.partition文件（如果有）
        partition_file = f"{self.output_prefix}.partition"
        if os.path.exists(partition_file):
            output_files.append(partition_file)
        
        return output_files