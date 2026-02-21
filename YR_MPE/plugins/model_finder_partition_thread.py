# model_finder_partition_thread.py
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
ModelFinder分区模式线程类
处理分区模型查找的命令构建和执行
"""

import os
import tempfile
from typing import List, Dict, Optional
from ..templates.base_process_thread import BaseProcessThread
from .partition_mode import PartitionDefinition, PartitionMode, PartitionFileIO


class ModelFinderPartitionThread(BaseProcessThread):
    """ModelFinder分区模型查找线程类"""
    
    def __init__(self, tool_path: str, input_files: List[str], 
                 partitions: List[PartitionDefinition], 
                 partition_mode: PartitionMode,
                 rcluster: bool = False,
                 rcluster_percent: Optional[int] = None,
                 extra_params: Optional[List[str]] = None):
        """
        初始化ModelFinder分区线程
        
        Args:
            tool_path: IQ-TREE可执行文件路径
            input_files: 输入文件列表
            partitions: 分区定义列表
            partition_mode: 分区模式
            rcluster: 是否使用分层聚类加速
            rcluster_percent: 分层聚类百分比
            extra_params: 额外的命令行参数
        """
        super().__init__(tool_path, input_files, extra_params or [])
        self.partitions = partitions
        self.partition_mode = partition_mode
        self.rcluster = rcluster
        self.rcluster_percent = rcluster_percent
    
    def get_tool_name(self) -> str:
        """返回工具名称"""
        return "IQ-TREE ModelFinder with Partition Mode"
    
    def create_temp_partition_file(self) -> str:
        """创建临时分区定义文件"""
        try:
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.nex', delete=False)
            
            # 写入NEXUS格式的分区定义
            PartitionFileIO.export_nexus_partition(self.partitions, temp_file.name)
            
            temp_file.close()
            return temp_file.name
            
        except Exception as e:
            self.error.emit(f"Failed to create partition file: {str(e)}")
            raise
    
    def build_iqtree_command(self, input_file: str, partition_file: str) -> List[str]:
        """
        构建IQ-TREE分区命令
        
        Args:
            input_file: 输入文件路径
            partition_file: 分区文件路径
            
        Returns:
            IQ-TREE命令列表
        """
        cmd = [self.tool_path]
        
        # 添加输入文件
        cmd.extend(["-s", input_file])
        
        # 添加分区模式参数
        mode_param = self.partition_mode.get_iqtree_parameter()
        cmd.extend([mode_param, partition_file])
        
        # 添加模型测试参数 (使用MFP+MERGE，IQ-TREE 3推荐)
        cmd.extend(["-m", "MFP+MERGE"])
        
        # 添加性能优化参数
        cmd.extend(["-T", "AUTO"])  # 自动检测CPU核心数
        
        # 可选：添加分层聚类加速（适用于大型数据集）
        if self.rcluster and self.rcluster_percent:
            cmd.extend(["--rcluster", str(self.rcluster_percent)])
        
        # 添加输出前缀
        input_base = os.path.splitext(os.path.basename(input_file))[0]
        cmd.extend(["--prefix", f"{input_base}_partitioned"])
        
        # 添加额外参数
        if self.parameters:
            cmd.extend(self.parameters)
        
        return cmd
    
    def execute_commands(self):
        """执行ModelFinder分区命令"""
        try:
            output_files = []
            html_files = []
            
            # 创建临时分区文件
            partition_file = self.create_temp_partition_file()
            
            # 处理每个输入文件
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                
                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(
                    f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", 
                    "info"
                )
                
                # 构建命令
                cmd = self.build_iqtree_command(input_file, partition_file)
                
                # 显示命令
                self.console_output.emit(f"Command: {' '.join(cmd)}", "info")
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(
                        f"ModelFinder partition execution failed for file {i+1}: {result.stderr}"
                    )
                    return
                
                # 查找生成的输出文件
                input_base = os.path.splitext(os.path.basename(input_file))[0]
                output_prefix = f"{input_base}_partitioned"
                
                # 查找.iqtree文件
                iqtree_file = f"{output_prefix}.iqtree"
                if os.path.exists(iqtree_file):
                    output_files.append(iqtree_file)
                else:
                    self.console_output.emit(
                        f"Warning: Could not find .iqtree file for {input_file}", 
                        "warning"
                    )
                
                # 查找.log文件
                log_file = f"{output_prefix}.log"
                if os.path.exists(log_file):
                    output_files.append(log_file)
                
                # 查找.best_model.nex文件（如果存在）
                best_model_file = f"{output_prefix}.best_model.nex"
                if os.path.exists(best_model_file):
                    output_files.append(best_model_file)
            
            # 删除临时分区文件
            try:
                if os.path.exists(partition_file):
                    os.remove(partition_file)
            except:
                pass
            
            self.progress.emit("ModelFinder partition analysis completed")
            self.finished.emit(output_files, html_files)
            
        except Exception as e:
            self.error.emit(f"ModelFinder partition analysis exception: {str(e)}")


class PartitionResultParser:
    """分区结果解析器"""
    
    @staticmethod
    def parse_partition_results(iqtree_file: str) -> Dict:
        """
        解析分区模型查找结果（IQ-TREE 3格式）
        
        Args:
            iqtree_file: IQ-TREE输出文件路径
            
        Returns:
            包含解析结果的字典
        """
        results = {
            'overall_model': '',
            'partition_models': {},
            'log_likelihood': 0.0,
            'aic': 0.0,
            'aicc': 0.0,
            'bic': 0.0,
            'num_partitions': 0
        }
        
        try:
            with open(iqtree_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # 解析整体模型（IQ-TREE 3格式）
            overall_patterns = [
                r"Best-fit partition scheme: (.+)",
                r"Best-fit model for full alignment: (.+)",
                r"Best partition scheme found: (.+)"
            ]
            for pattern in overall_patterns:
                match = re.search(pattern, content)
                if match:
                    results['overall_model'] = match.group(1).strip()
                    break
            
            # 解析各分区模型
            partition_patterns = [
                r"Partition (\d+): (.+) \((.+)\) --> (.+)",
                r"Partition (\d+): (.+) --> (.+)",
                r"Best model for (.+): (.+)"
            ]
            
            import re
            
            for pattern in partition_patterns:
                matches = re.findall(pattern, content)
                if matches:
                    for match in matches:
                        if len(match) == 4:
                            partition_num, partition_name, partition_desc, model = match
                            partition_name = partition_name.strip()
                            results['partition_models'][partition_name] = {
                                'model': model.strip(),
                                'description': partition_desc.strip()
                            }
                        elif len(match) == 3:
                            partition_num, partition_name, model = match
                            partition_name = partition_name.strip()
                            results['partition_models'][partition_name] = {
                                'model': model.strip(),
                                'description': ''
                            }
                        elif len(match) == 2:
                            partition_name, model = match
                            partition_name = partition_name.strip()
                            results['partition_models'][partition_name] = {
                                'model': model.strip(),
                                'description': ''
                            }
                    break
            
            # 更新分区数量
            results['num_partitions'] = len(results['partition_models'])
            
            # 解析统计信息（IQ-TREE 3格式）
            stats_patterns = [
                r"Log-likelihood:\s+([+-]?[\d\.]+).*?AICc:\s+([+-]?[\d\.]+).*?BIC:\s+([+-]?[\d\.]+)",
                r"Log-likelihood:\s+([+-]?[\d\.]+).*?AIC:\s+([+-]?[\d\.]+).*?AICc:\s+([+-]?[\d\.]+).*?BIC:\s+([+-]?[\d\.]+)"
            ]
            for pattern in stats_patterns:
                match = re.search(pattern, content, re.DOTALL)
                if match:
                    if len(match.groups()) == 3:
                        results['log_likelihood'] = float(match.group(1))
                        results['aicc'] = float(match.group(2))
                        results['bic'] = float(match.group(3))
                    elif len(match.groups()) == 4:
                        results['log_likelihood'] = float(match.group(1))
                        results['aic'] = float(match.group(2))
                        results['aicc'] = float(match.group(3))
                        results['bic'] = float(match.group(4))
                    break
            
        except Exception as e:
            raise ValueError(f"Failed to parse IQ-TREE file: {str(e)}")
        
        return results