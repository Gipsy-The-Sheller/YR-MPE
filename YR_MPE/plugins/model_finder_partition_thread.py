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
import re
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
                 global_params: Optional[Dict] = None):
        """
        初始化ModelFinder分区线程
        
        Args:
            tool_path: IQ-TREE可执行文件路径
            input_files: 输入文件列表
            partitions: 分区定义列表
            partition_mode: 分区模式
            rcluster: 是否使用分层聚类加速
            rcluster_percent: 分层聚类百分比
            global_params: 全局参数字典，包含 seq_type, criterion, threads
        """
        super().__init__(tool_path, input_files, [])
        self.partitions = partitions
        self.partition_mode = partition_mode
        self.rcluster = rcluster
        self.rcluster_percent = rcluster_percent
        self.global_params = global_params or {}
    
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
        
        # 添加序列类型参数（从全局参数获取）
        seq_type = self.global_params.get('seq_type', 'AUTO')
        if seq_type != "AUTO":
            seq_type_code = {"dna": "DNA", "prot": "AA"}[seq_type.lower()]
            cmd.extend(["-st", seq_type_code])
        
        # 添加准则参数（从全局参数获取）
        criterion = self.global_params.get('criterion', 'BIC')
        if criterion == "BIC":
            pass  # BIC是默认值，不需要额外参数
        elif criterion == "AICc":
            cmd.extend(["-AICc"])
        elif criterion == "AIC":
            cmd.extend(["-AIC"])
        
        # 添加模型集参数
        cmd.extend(["-mset", "ALL"])  # 测试所有模型
        
        # 添加线程数参数（从全局参数获取）
        threads = self.global_params.get('threads', 1)
        if threads > 1:
            cmd.extend(["-nt", str(threads)])
        
        # 可选：添加分层聚类加速（适用于大型数据集）
        if self.rcluster and self.rcluster_percent:
            cmd.extend(["--rcluster", str(self.rcluster_percent)])
        
        # 添加输出前缀
        input_base = os.path.splitext(os.path.basename(input_file))[0]
        cmd.extend(["--prefix", f"{input_base}_partitioned"])
        
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
                
                # 查找.best_scheme.nex文件（如果存在，包含分区合并后的最终方案）
                best_scheme_file = f"{output_prefix}.best_scheme.nex"
                if os.path.exists(best_scheme_file):
                    output_files.append(best_scheme_file)
                
                # 查找.best_model.nex文件（如果存在，旧版本IQ-TREE使用）
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
            
            # 解析统计信息（IQ-TREE 3格式）- 尝试两种模式
            try:
                stats_patterns = [
                    r"Log-likelihood:\s+([+-]?[\d\.]+).*?AICc:\s+([+-]?[\d\.]+).*?BIC:\s+([+-]?[\d\.]+)",
                    r"Log-likelihood:\s+([+-]?[\d\.]+).*?AIC:\s+([+-]?[\d\.]+).*?AICc:\s+([+-]?[\d\.]+).*?BIC:\s+([+-]?[\d\.]+)"
                ]
                
                # 先尝试第一种模式
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
                pass  # 统计信息解析失败不影响主要结果
            
            # 解析 "TREE USED FOR ModelFinder" 部分的详细信息（优先级更高）
            try:
                # Log-likelihood of the tree: -1519.5139 (s.e. 51.6783)
                tree_logl_pattern = r"Log-likelihood of the tree:\s+([+-]?[\d\.]+)\s+\(s\.e\.\s+([+-]?[\d\.]+)\)"
                tree_logl_match = re.search(tree_logl_pattern, content)
                if tree_logl_match:
                    results['log_likelihood'] = float(tree_logl_match.group(1))
                    results['tree_log_likelihood'] = float(tree_logl_match.group(1))
                    results['tree_log_likelihood_se'] = float(tree_logl_match.group(2))
                
                # Unconstrained log-likelihood (without tree): -2518.4649
                unconstrained_logl_pattern = r"Unconstrained log-likelihood \(without tree\):\s+([+-]?[\d\.]+)"
                unconstrained_logl_match = re.search(unconstrained_logl_pattern, content)
                if unconstrained_logl_match:
                    results['unconstrained_log_likelihood'] = float(unconstrained_logl_match.group(1))
                
                # Number of free parameters (#branches + #model parameters): 65
                free_params_pattern = r"Number of free parameters.*?:\s+(\d+)"
                free_params_match = re.search(free_params_pattern, content)
                if free_params_match:
                    results['free_parameters'] = int(free_params_match.group(1))
                
                # Akaike information criterion (AIC) score: 3169.0278
                aic_score_pattern = r"Akaike information criterion \(AIC\) score:\s+([+-]?[\d\.]+)"
                aic_score_match = re.search(aic_score_pattern, content)
                if aic_score_match:
                    results['aic_score'] = float(aic_score_match.group(1))
                    results['aic'] = float(aic_score_match.group(1))
                
                # Corrected Akaike information criterion (AICc) score: 3182.4340
                aicc_score_pattern = r"Corrected Akaike information criterion \(AICc\) score:\s+([+-]?[\d\.]+)"
                aicc_score_match = re.search(aicc_score_pattern, content)
                if aicc_score_match:
                    results['aicc_score'] = float(aicc_score_match.group(1))
                    results['aicc'] = float(aicc_score_match.group(1))
                
                # Bayesian information criterion (BIC) score: 3465.4028
                bic_score_pattern = r"Bayesian information criterion \(BIC\) score:\s+([+-]?[\d\.]+)"
                bic_score_match = re.search(bic_score_pattern, content)
                if bic_score_match:
                    results['bic_score'] = float(bic_score_match.group(1))
                    results['bic'] = float(bic_score_match.group(1))
            except Exception as e:
                pass  # 统计信息解析失败不影响主要结果
            
            return results
            
        except Exception as e:
            import traceback
            error_msg = f"Failed to parse IQ-TREE file: {str(e)}\nTraceback: {traceback.format_exc()}"
            raise ValueError(error_msg)
            
            # Total tree length (sum of branch lengths): 0.2410
            tree_length_pattern = r"Total tree length.*?:\s+([+-]?[\d\.]+)"
            tree_length_match = re.search(tree_length_pattern, content)
            if tree_length_match:
                results['tree_length'] = float(tree_length_match.group(1))
            
            # Sum of internal branch lengths: 0.0263 (10.8978% of tree length)
            internal_branch_length_pattern = r"Sum of internal branch lengths:\s+([+-]?[\d\.]+)\s+\(([\d\.]+)%\s+of tree length\)"
            internal_branch_match = re.search(internal_branch_length_pattern, content)
            if internal_branch_match:
                results['internal_branch_length'] = float(internal_branch_match.group(1))
                results['internal_branch_length_percent'] = float(internal_branch_match.group(2))
            
            # 解析每个分区的 Speed 和完整参数
            # 格式: ID  Model           Speed  Parameters
            #        1  TIM2e+G4       1.0000  TIM2e{2.2081,4.37413,11.5914}+FQ+G4{0.275195}
            substitution_process_pattern = r"SUBSTITUTION PROCESS\s*-+\s*Edge-linked.*?\s+ID\s+Model\s+Speed\s+Parameters\s*(.*?)(?=\n[A-Z]|\Z)"
            substitution_process_match = re.search(substitution_process_pattern, content, re.DOTALL)
            if substitution_process_match:
                partition_details = substitution_process_match.group(1).strip()
                # 每行格式: ID  Model  Speed  Parameters
                for line in partition_details.split('\n'):
                    line = line.strip()
                    if not line or not line[0].isdigit():
                        continue
                    parts = line.split()
                    if len(parts) >= 4:
                        # 假设格式: ID Model Speed Parameters
                        # 参数部分可能包含空格，所以需要特殊处理
                        id_num = parts[0]
                        model_name = parts[1]
                        speed = parts[2]
                        # Parameters 是剩余部分
                        parameters = ' '.join(parts[3:])
                        
                        # 找到对应的分区并添加 Speed 和 Parameters
                        for partition_name, model_info in results['partition_models'].items():
                            if model_info.get('model', '') == model_name:
                                model_info['speed'] = float(speed)
                                model_info['parameters'] = parameters
                                break
            
        except Exception as e:
            raise ValueError(f"Failed to parse IQ-TREE file: {str(e)}")
    
    @staticmethod
    def parse_best_scheme_nex(best_scheme_file: str) -> Dict:
        """
        从 .best_scheme.nex 文件解析最终的分区方案（推荐使用，因为它包含分区合并后的结果）
        
        Args:
            best_scheme_file: .best_scheme.nex 文件路径
            
        Returns:
            包含解析结果的字典
        """
        results = {
            'overall_model': '',
            'partition_models': {},
            'charset_ranges': {},  # 新增：存储每个 charset 的位点范围
            'log_likelihood': 0.0,
            'aic': 0.0,
            'aicc': 0.0,
            'bic': 0.0,
            'num_partitions': 0
        }
        
        try:
            with open(best_scheme_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # 检查是否是NEXUS格式
            if not content.lower().startswith('#nexus'):
                raise ValueError("File is not in NEXUS format")
            
            # 步骤1：解析 charset 定义（获取分区名称和位点范围）
            # 格式: charset partition_name = 1-100  200-300;
            charset_pattern = r'charset\s+([^\s=]+)\s*=\s*([^;]+);'
            charset_matches = re.findall(charset_pattern, content, re.IGNORECASE)
            
            # 存储每个 charset 的名称和位点范围列表
            charset_data = {}
            for charset_name, range_str in charset_matches:
                charset_name = charset_name.strip()
                # 解析位点范围（可能有多个，用空格分隔）
                ranges = [r.strip() for r in range_str.split()]
                charset_data[charset_name] = ranges
                results['charset_ranges'][charset_name] = ranges
            
            # 步骤2：解析 charpartition 命令（获取模型分配）
            # 格式: charpartition mymodels = MODEL: partition_name, MODEL: partition_name;
            charpartition_pattern = r'charpartition\s+\w+\s*=\s*([^;]+);'
            charpartition_match = re.search(charpartition_pattern, content, re.IGNORECASE | re.DOTALL)
            
            if charpartition_match:
                partition_defs = charpartition_match.group(1)
                
                # 逐行解析分区定义
                lines = partition_defs.strip().split('\n')
                
                for line in lines:
                    line = line.strip()
                    
                    # 跳过空行和注释行
                    if not line or line.startswith('#'):
                        continue
                    
                    # 移除尾部的逗号
                    line = line.rstrip(',')
                    
                    # 使用正则匹配: MODEL: partition_name1, partition_name2
                    partition_match = re.match(r'^(.+?)\s*:\s*(.+)$', line)
                    
                    if partition_match:
                        model_str = partition_match.group(1).strip()
                        partition_names_str = partition_match.group(2).strip()
                        
                        # 分割多个分区名称（可能有逗号）
                        partition_names = [name.strip() for name in partition_names_str.split(',')]
                        
                        # 为每个分区分配模型
                        for partition_name in partition_names:
                            if partition_name:  # 确保名称不为空
                                results['partition_models'][partition_name] = {
                                    'model': model_str,
                                    'description': ''
                                }
                
                # 更新分区数量
                results['num_partitions'] = len(results['partition_models'])
                
                # 构建整体模型描述
                if results['partition_models']:
                    models = [pm['model'] for pm in results['partition_models'].values()]
                    # 如果所有分区使用相同模型，显示一个；否则显示数量
                    unique_models = list(set(models))
                    if len(unique_models) == 1:
                        results['overall_model'] = unique_models[0]
                    else:
                        results['overall_model'] = f"{len(unique_models)} different models"
            
            # 从 .iqtree 文件读取统计信息（如果存在）
            iqtree_file = best_scheme_file.replace('.best_scheme.nex', '.iqtree')
            if os.path.exists(iqtree_file):
                try:
                    with open(iqtree_file, 'r', encoding='utf-8', errors='ignore') as f:
                        iqtree_content = f.read()
                    
                    # 解析统计信息
                    stats_patterns = [
                        r"Log-likelihood:\s+([+-]?[\d\.]+).*?AICc:\s+([+-]?[\d\.]+).*?BIC:\s+([+-]?[\d\.]+)",
                        r"Log-likelihood:\s+([+-]?[\d\.]+).*?AIC:\s+([+-]?[\d\.]+).*?AICc:\s+([+-]?[\d\.]+).*?BIC:\s+([+-]?[\d\.]+)"
                    ]
                    for pattern in stats_patterns:
                        match = re.search(pattern, iqtree_content, re.DOTALL)
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
                except:
                    pass  # 统计信息解析失败不影响主要结果
            
        except Exception as e:
            raise ValueError(f"Failed to parse .best_scheme.nex file: {str(e)}")
        
        return results
        
        return results