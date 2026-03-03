from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, 
                             QPushButton, QLineEdit, QTextEdit, 
                             QLabel, QComboBox, QCheckBox, QRadioButton, 
                             QSpinBox, QDoubleSpinBox, QScrollArea, 
                             QFrame, QTextEdit, QToolButton, QDialog,
                             QGroupBox, QSizePolicy, QFormLayout, QGridLayout,
                             QFileDialog, QMessageBox, QApplication,
                             QListWidget, QTableWidget, QTableWidgetItem,
                             QHeaderView, QAbstractItemView)
from PyQt5.QtCore import Qt, pyqtSignal, QRegExp
from PyQt5.QtGui import QFont, QTextCursor, QSyntaxHighlighter, QTextCharFormat, QColor, QPalette
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from ..platforms.methods.dataset_models import ChainItem
from .partition_mode import PartitionMode, MrBayesPartitionDefinition, MrBayesModelConverter
from .mrbayes_partition_ui import MrBayesPartitionDialog
import os
import tempfile
import json
from Bio import SeqIO
import re

class MPIBeagleSettingsDialog(QDialog):
    """MPI & BEAGLE Settings Dialog"""
    
    def __init__(self, parent=None, use_mpi=True, 
                 use_beagle=True, beagle_device='CPU', 
                 beagle_precision='Double', beagle_scaling='Dynamic'):
        super().__init__(parent)
        self.setWindowTitle("MPI & BEAGLE Settings")
        self.setMinimumSize(400, 300)
        
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # MPI Settings
        mpi_group = QGroupBox("MPI Settings")
        mpi_layout = QVBoxLayout()
        mpi_group.setLayout(mpi_layout)
        
        self.use_mpi = QCheckBox("Use MPI")
        self.use_mpi.setChecked(use_mpi)
        mpi_layout.addWidget(self.use_mpi)
        
        # 添加说明文本
        mpi_info = QLabel("MPI will use (Runs × Chains) processors")
        mpi_info.setStyleSheet("color: #666; font-size: 11px;")
        mpi_layout.addWidget(mpi_info)
        
        layout.addWidget(mpi_group)
        
        # BEAGLE Settings
        beagle_group = QGroupBox("BEAGLE Settings")
        beagle_layout = QVBoxLayout()
        beagle_group.setLayout(beagle_layout)
        
        self.use_beagle = QCheckBox("Use BEAGLE")
        self.use_beagle.setChecked(use_beagle)
        self.use_beagle.stateChanged.connect(self.on_beagle_toggled)
        beagle_layout.addWidget(self.use_beagle)
        
        beagle_form_layout = QFormLayout()
        
        self.beagle_device_combo = QComboBox()
        self.beagle_device_combo.addItems(["CPU", "GPU"])
        self.beagle_device_combo.setCurrentText(beagle_device)
        beagle_form_layout.addRow("Device:", self.beagle_device_combo)
        
        self.beagle_precision_combo = QComboBox()
        self.beagle_precision_combo.addItems(["Double", "Single"])
        self.beagle_precision_combo.setCurrentText(beagle_precision)
        beagle_form_layout.addRow("Precision:", self.beagle_precision_combo)
        
        self.beagle_scaling_combo = QComboBox()
        self.beagle_scaling_combo.addItems(["Dynamic", "Always"])
        self.beagle_scaling_combo.setCurrentText(beagle_scaling)
        beagle_form_layout.addRow("Scaling:", self.beagle_scaling_combo)
        
        beagle_layout.addLayout(beagle_form_layout)
        
        layout.addWidget(beagle_group)
        
        # Buttons
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        button_layout.addWidget(ok_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
        
        # Initial state
        self.on_beagle_toggled()
    
    def on_beagle_toggled(self):
        """BEAGLE复选框状态改变时处理"""
        enabled = self.use_beagle.isChecked()
        self.beagle_device_combo.setEnabled(enabled)
        self.beagle_precision_combo.setEnabled(enabled)
        self.beagle_scaling_combo.setEnabled(enabled)
    
    def get_settings(self):
        """获取设置值"""
        return {
            'use_mpi': self.use_mpi.isChecked(),
            'use_beagle': self.use_beagle.isChecked(),
            'beagle_device': self.beagle_device_combo.currentText(),
            'beagle_precision': self.beagle_precision_combo.currentText(),
            'beagle_scaling': self.beagle_scaling_combo.currentText()
        }

class MrBayesThread(BaseProcessThread):
    """MrBayes系统发育推断线程类"""
    
    def __init__(self, tool_path, mpirun_path, input_files, parameters, imported_files=None, run_data_block_checked=True, 
                 use_partition_mode=False, partition_definitions=None, partition_mode=None, workdir=None, use_mpi=False):
        super().__init__(tool_path, input_files, parameters, imported_files, workdir=workdir)
        self.mpirun_path = mpirun_path
        self.run_data_block_checked = run_data_block_checked
        self.use_mpi = use_mpi
        # 分区模式相关
        self.use_partition_mode = use_partition_mode
        self.partition_definitions = partition_definitions or []
        self.partition_mode = partition_mode
    
    def get_tool_name(self):
        """返回工具名称"""
        return "MrBayes-MPI-BEAGLE Phylogeny"
    
    def execute_commands(self):
        """执行MrBayes命令"""
        try:
            output_files = []
            
            # 分别处理每个输入文件
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                    
                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", "info")
                
                # 获取参数
                params = self.parameters.copy()
                
                # 获取MCMC参数
                if len(params) > 0 and isinstance(params[0], dict):
                    param_dict = params[0]
                    nchains = param_dict.get('nchains', 4)
                    nruns = param_dict.get('nruns', 2)
                else:
                    nchains = 4
                    nruns = 2
                
                # 计算MPI进程数
                # MPI进程数 = nchains * nruns（每个进程处理一条链）
                total_threads = nchains * nruns
                
                # 生成MrBayes NEXUS文件
                nexus_file = self._generate_mrbayes_script(input_file, params)
                
                # 构建命令
                if self.use_mpi and total_threads > 1 and self.mpirun_path:
                    # 使用MPI并行版本
                    cmd = [
                        self.mpirun_path,
                        "-np", str(total_threads),
                        self.tool_path
                    ]
                else:
                    # 使用单线程版本
                    cmd = [self.tool_path]
                
                # 添加输入文件
                cmd.append(nexus_file)
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"MrBayes execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的输出文件（在nexus_file所在目录）
                nexus_dir = os.path.dirname(nexus_file)
                nexus_name = os.path.splitext(os.path.basename(nexus_file))[0]
                
                # 查找.con.tre文件
                con_tre_files = [f for f in os.listdir(nexus_dir) if f.startswith(nexus_name) and f.endswith('.con.tre')]
                for f in con_tre_files:
                    output_files.append(os.path.join(nexus_dir, f))
                
                # 查找.p文件（chain文件）
                chain_files = [f for f in os.listdir(nexus_dir) if f.startswith(nexus_name) and '.run' in f and f.endswith('.p')]
                for f in chain_files:
                    output_files.append(os.path.join(nexus_dir, f))
                
                # 查找.t文件
                t_files = [f for f in os.listdir(nexus_dir) if f.startswith(nexus_name) and '.run' in f and f.endswith('.t')]
                for f in t_files:
                    output_files.append(os.path.join(nexus_dir, f))
                        
            self.progress.emit("Phylogenetic inference completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Phylogenetic inference exception: {str(e)}")
    
    def _generate_mrbayes_script(self, input_file, params):
        """生成MrBayes命令脚本（完整的NEXUS文件）"""
        
        # 检查输入文件格式
        file_ext = os.path.splitext(input_file)[1].lower()
        if file_ext in ['.nex', '.nexus']:
            # 检查NEXUS文件是否包含MrBayes块
            has_mrbayes_block = self._check_nexus_for_mrbayes_block(input_file)
            
            if has_mrbayes_block and self.run_data_block_checked:
                # 如果文件包含MrBayes块且checkbox被选中，直接使用原文件
                return input_file
            elif has_mrbayes_block and not self.run_data_block_checked:
                # 如果文件包含MrBayes块但checkbox未被选中，移除原有的MrBayes块
                nexus_file = self._remove_mrbayes_block_from_nexus(input_file)
                # 在这个文件后面添加我们的MrBayes命令块
                return self._add_mrbayes_block_to_nexus(nexus_file, params)
            else:
                # 文件不包含MrBayes块，添加我们的MrBayes命令块
                return self._add_mrbayes_block_to_nexus(input_file, params)
        else:
            # 非NEXUS格式，转换为NEXUS格式
            nexus_file = self._convert_to_nexus(input_file)
            # 添加MrBayes命令块
            return self._add_mrbayes_block_to_nexus(nexus_file, params)
    
    def _add_mrbayes_block_to_nexus(self, nexus_file, params):
        """向NEXUS文件添加MrBayes命令块"""
        try:
            # 读取原NEXUS文件内容
            with open(nexus_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 生成MrBayes命令块
            mrbayes_commands = []
            mrbayes_commands.append("begin mrbayes;")
            
            # 检查是否启用分区模式
            if self.use_partition_mode and self.partition_definitions:
                # 分区模式
                partition_commands = MrBayesModelConverter.generate_partition_commands(
                    self.partition_definitions, 
                    self.partition_mode
                )
                mrbayes_commands.extend(partition_commands)
            else:
                # 单一模型模式
                # 设置模型参数
                mrbayes_commands.extend(self._get_model_commands(params))
            
            # 设置分子钟定年参数
            mrbayes_commands.extend(self._get_clock_commands(params))
            
            # 设置MCMC参数
            mrbayes_commands.extend(self._get_mcmc_commands(params))
            
            # 运行MCMC
            mrbayes_commands.append("mcmc;")
            
            # 总结结果
            mrbayes_commands.extend(self._get_summary_commands(params))
            
            # 结束MrBayes块
            mrbayes_commands.append("end;")
            
            # 添加退出命令
            mrbayes_commands.append("quit;")
            
            # 将MrBayes命令块添加到NEXUS文件末尾
            mrbayes_block = "\n".join(mrbayes_commands)
            new_content = content + "\n" + mrbayes_block
            
            # 创建临时文件保存新内容
            output_file = self.create_temp_file(suffix='.nex')
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(new_content)
            
            return output_file
        except Exception as e:
            self.error.emit(f"Error adding MrBayes block to NEXUS file: {str(e)}")
            return nexus_file
    
    def _check_nexus_for_mrbayes_block(self, nexus_file):
        """检查NEXUS文件是否包含MrBayes块"""
        try:
            with open(nexus_file, 'r') as f:
                content = f.read().lower()
                # 检查是否包含MrBayes块
                return bool(content.find('begin mrbayes;'))
        except Exception:
            self.error.emit("Error reading NEXUS file")
            return False
    
    def _remove_mrbayes_block_from_nexus(self, nexus_file):
        """从NEXUS文件中移除MrBayes块，保留数据块"""
        try:
            with open(nexus_file, 'r') as f:
                content = f.read()
            
            # 分割为lines，逐行处理
            lines = content.splitlines()
            
            # 标记是否在MrBayes块内
            inside_mrbayes_block = False
            processed_lines = []
            
            for line in lines:
                stripped_line = line.strip().lower()
                
                # 检查MrBayes块开始
                if stripped_line.startswith('begin mrbayes;') or stripped_line.startswith('begin mrbayes'):
                    inside_mrbayes_block = True
                    continue  # 跳过MrBayes块开始行
                
                # 检查MrBayes块结束
                if inside_mrbayes_block and stripped_line.startswith('end;'):
                    inside_mrbayes_block = False
                    continue  # 跳过MrBayes块结束行
                
                # 如果在MrBayes块内，跳过该行
                if inside_mrbayes_block:
                    continue
                
                # 否则保留该行
                processed_lines.append(line)
            
            # 创建临时文件保存处理后的内容
            temp_nexus_file = self.create_temp_file(suffix='.nex')
            with open(temp_nexus_file, 'w') as f:
                f.write('\n'.join(processed_lines))
            
            return temp_nexus_file
        except Exception:
            # 如果处理失败，返回原文件
            return nexus_file
    
    def _convert_to_nexus(self, input_file):
        """将输入文件转换为NEXUS格式（interleave格式）"""
        # 确定输入文件格式
        file_ext = os.path.splitext(input_file)[1].lower()
        if file_ext in ['.fas', '.fasta', '.fa', '.fna']:
            input_format = 'fasta'
        elif file_ext in ['.phy', '.phylip']:
            input_format = 'phylip'
        elif file_ext in ['.nex', '.nexus']:
            # 如果已经是NEXUS格式，检查是否是interleave格式
            # 如果不是interleave格式，需要重新转换
            with open(input_file, 'r', encoding='utf-8') as f:
                content = f.read()
                # 检查是否是interleave格式（检查matrix块的格式）
                if 'interleave' in content.lower() or 'interleaved' in content.lower():
                    return input_file
            # 如果不是interleave格式，继续转换
            input_format = 'nexus'
        else:
            # 默认尝试fasta格式
            input_format = 'fasta'
        
        # 生成临时NEXUS文件
        output_file = self.create_temp_file(suffix='.nex')
        
        try:
            # 使用BioPython进行格式转换
            sequences = list(SeqIO.parse(input_file, input_format))
            
            # 清理分类名称中的特殊字符
            for seq in sequences:
                seq.id = self._clean_taxon_name(seq.id)
                if seq.name:
                    seq.name = self._clean_taxon_name(seq.name)
                if hasattr(seq, 'description') and seq.description:
                    seq.description = self._clean_taxon_name(seq.description)
            
            # 为序列添加分子类型（默认DNA）
            for seq in sequences:
                seq.annotations['molecule_type'] = 'DNA'
            
            # 使用'nexus'格式写入
            SeqIO.write(sequences, output_file, 'nexus')
            
            # 手动转换为interleave格式
            self._convert_to_interleave_format(output_file)
            
            return output_file
        except Exception as e:
            self.error.emit(f"Error converting file to NEXUS format: {str(e)}")
            # 如果转换失败，返回原始文件
            return input_file
    
    def _clean_taxon_name(self, name):
        """清理分类名称中的特殊字符，替换为下划线"""
        if not name:
            return name
        
        # 将空格和标点符号替换为下划线
        import string
        # 允许的字符：字母、数字、下划线
        allowed_chars = set(string.ascii_letters + string.digits + '_')
        
        cleaned = []
        for char in name:
            if char in allowed_chars:
                cleaned.append(char)
            else:
                cleaned.append('_')  # 将其他字符替换为下划线
        
        # 避免连续的下划线
        result = ''.join(cleaned)
        while '__' in result:
            result = result.replace('__', '_')
        
        # 避免以非字母开头
        if result and result[0] not in string.ascii_letters:
            result = 'T_' + result
        
        return result
    
    def _convert_to_interleave_format(self, nexus_file):
        """将sequential NEXUS文件转换为interleave格式"""
        try:
            with open(nexus_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 检查format行中是否已经有interleave关键字
            if 'interleave' not in content.lower():
                # 如果没有，在format行中添加interleave关键字
                content = content.replace('format datatype=', 'format datatype=')
                # 在format行的末尾添加interleave关键字（在分号之前）
                content = re.sub(r'(format\s+datatype=\w+\s+missing=\?\s+gap=-\s*)(;)', r'\1interleave \2', content)
            
            with open(nexus_file, 'w', encoding='utf-8') as f:
                f.write(content)
        except Exception as e:
            print(f"Warning: Failed to convert to interleave format: {e}")
    
    def _get_clock_commands(self, params):
        """获取分子钟设置命令"""
        commands = []
        
        # 处理参数，可能是列表或字典
        if isinstance(params, dict):
            param_dict = params
        elif isinstance(params, list) and len(params) > 0 and isinstance(params[0], dict):
            param_dict = params[0]
        else:
            param_dict = {}
        
        enable_dating = param_dict.get('enable_dating', False)
        if not enable_dating:
            return commands
        
        # 启用分子钟
        commands.append("prset brlenspr=clock:uniform;")
        
        # 设置松弛分子钟类型
        clock_model = param_dict.get('clock_model', 'strict')
        clock_params = param_dict.get('clock_params', {})
        
        if clock_model == 'tk02':
            # TK02参数
            tk02varpr = clock_params.get('tk02varpr', 'exponential')
            tk02varpr_mean = clock_params.get('tk02varpr_mean', 1.00)
            if tk02varpr == 'fixed':
                commands.append(f"prset clockvarpr=tk02 tk02varpr=fixed({tk02varpr_mean});")
            elif tk02varpr == 'exponential':
                commands.append(f"prset clockvarpr=tk02 tk02varpr=exponential({tk02varpr_mean});")
            else:  # uniform
                commands.append(f"prset clockvarpr=tk02 tk02varpr=uniform({tk02varpr_mean});")
        
        elif clock_model == 'igr':
            # IGR参数
            igrvarpr = clock_params.get('igrvarpr', 'exponential')
            igrvarpr_mean = clock_params.get('igrvarpr_mean', 10.00)
            if igrvarpr == 'fixed':
                commands.append(f"prset clockvarpr=igr igrvarpr=fixed({igrvarpr_mean});")
            elif igrvarpr == 'exponential':
                commands.append(f"prset clockvarpr=igr igrvarpr=exponential({igrvarpr_mean});")
            else:  # uniform
                commands.append(f"prset clockvarpr=igr igrvarpr=uniform({igrvarpr_mean});")
        
        elif clock_model == 'cpp':
            # CPP参数
            cppratepr = clock_params.get('cppratepr', 'exponential')
            cppratepr_mean = clock_params.get('cppratepr_mean', 0.10)
            cppmultdevpr_value = clock_params.get('cppmultdevpr_value', 0.40)
            
            if cppratepr == 'fixed':
                commands.append(f"prset clockvarpr=cpp cppratepr=fixed({cppratepr_mean});")
            else:  # exponential
                commands.append(f"prset clockvarpr=cpp cppratepr=exponential({cppratepr_mean});")
            
            commands.append(f"prset cppmultdevpr=fixed({cppmultdevpr_value});")
        
        # strict 不需要额外参数
        
        # 设置校准节点
        calibrations = param_dict.get('calibrations', [])
        if calibrations:
            for i, calib in enumerate(calibrations):
                node_name = f"calibration_{i+1}"
                
                # 检查是否使用约束
                use_constraint = calib.get('use_constraint', True)
                
                # 定义分类群集合（仅在use_constraint为True时）
                if use_constraint:
                    taxa = calib.get('taxa', [])
                    if taxa:
                        # MrBayes的constraint命令：taxa之间用逗号分隔，不能有空格
                        taxa_str = ",".join(taxa)
                        commands.append(f"constraint {node_name} = {taxa_str};")
                    
                    # 设置校准先验
                    prior_type = calib.get('prior_type', 'fixed')
                    calib_params = calib.get('params', {})
                    
                    if prior_type == 'fixed':
                        value = calib_params.get('fixed_value', 0)
                        commands.append(f"calibrate {node_name} = fixed({value});")
                    elif prior_type == 'uniform':
                        min_val = calib_params.get('uniform_min', 0)
                        max_val = calib_params.get('uniform_max', 0)
                        commands.append(f"calibrate {node_name} = unif({min_val},{max_val});")
                    elif prior_type == 'offset_exp':
                        offset = calib_params.get('offset_exp_offset', 0)
                        rate = calib_params.get('offset_exp_rate', 1.0)
                        commands.append(f"calibrate {node_name} = offsetexp({offset},{rate});")
                    elif prior_type == 'offset_gamma':
                        offset = calib_params.get('offset_gamma_offset', 0)
                        alpha = calib_params.get('offset_gamma_alpha', 1.0)
                        beta = calib_params.get('offset_gamma_beta', 1.0)
                        commands.append(f"calibrate {node_name} = offsetgamma({offset},{alpha},{beta});")
                    elif prior_type == 'offset_lognorm':
                        offset = calib_params.get('offset_lognorm_offset', 0)
                        mean = calib_params.get('offset_lognorm_mean', 0.0)
                        std = calib_params.get('offset_lognorm_std', 1.0)
                        commands.append(f"calibrate {node_name} = offsetlognormal({offset},{mean},{std});")
                    elif prior_type == 'truncated_normal':
                        offset = calib_params.get('trunc_norm_offset', 0)
                        mean = calib_params.get('trunc_norm_mean', 0.0)
                        std = calib_params.get('trunc_norm_std', 1.0)
                        commands.append(f"calibrate {node_name} = truncatednormal({offset},{mean},{std});")
                else:
                    # 不使用约束，只定义节点名称（不生成constraint命令）
                    # 但仍需要设置校准先验
                    # 在MrBayes中，如果不定义constraint，需要使用其他方式来指定校准节点
                    # 这种情况下，我们仍然可以设置校准，但需要使用不同的方法
                    pass
            
            # 如果有任何校准节点，启用校准
            # 检查是否有使用约束的校准节点
            has_constraint_calibrations = any(calib.get('use_constraint', True) for calib in calibrations)
            if has_constraint_calibrations:
                commands.append("prset nodeagepr=calibrated;")
        
        # 设置外群
        outgroup = param_dict.get('outgroup', '')
        if outgroup:
            # 外群可能是单个或多个分类群，用逗号分隔
            # 先清理每个外群名称
            outgroup_taxa = [self._clean_taxon_name(t.strip()) for t in outgroup.split(',') if t.strip()]
            if outgroup_taxa:
                # 定义ingroup constraint（单数）
                outgroup_str = ",".join(outgroup_taxa)
                commands.append(f"constraint ingroup = {outgroup_str};")
                # 应用constraint（复数，引用constraint集合）
                commands.append("prset topologypr=constraints(ingroups);")
        
        return commands
    
    def _get_model_commands(self, params):
        """获取模型设置命令"""
        commands = []
        
        # 数据类型 - params现在直接是字典
        if isinstance(params, dict):
            param_dict = params
        elif isinstance(params, list) and len(params) > 0 and isinstance(params[0], dict):
            param_dict = params[0]
        else:
            param_dict = {}

        datatype = param_dict.get('datatype', 'DNA')
        
        if datatype == 'DNA':
            # DNA模型设置
            nst = param_dict.get('nst', 6)
            statefreq = param_dict.get('statefreq', 'estimated(dirichlet)')
            
            # 设置nucmodel和nst
            commands.append(f"lset nucmodel=4by4 nst={nst};")
            
            # 设置状态频率
            if statefreq == 'fixed(equal)':
                commands.append("prset statefreqpr=fixed(equal);")
            elif statefreq == 'fixed(empirical)':
                commands.append("prset statefreqpr=fixed(empirical);")
            else:
                commands.append("prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);")
        else:
            # 蛋白质模型设置
            aamodel = param_dict.get('aamodel', 'mixed')
            statefreq = param_dict.get('statefreq', 'estimated(dirichlet)')
            
            # 设置nucmodel为protein
            commands.append("lset nucmodel=protein;")
            
            # 设置氨基酸模型
            if aamodel == 'mixed':
                commands.append("prset aamodelpr=mixed;")
            else:
                commands.append(f"prset aamodelpr=fixed({aamodel.lower()});")
            
            # 设置状态频率
            if statefreq == 'fixed(equal)':
                commands.append("prset statefreqpr=fixed(equal);")
            elif statefreq == 'fixed(empirical)':
                commands.append("prset statefreqpr=fixed(empirical);")
            else:
                commands.append("prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);")
        
        # 设置速率异质性
        rates = param_dict.get('rates', 'equal')
        ngammacat = param_dict.get('ngammacat', 4)
        
        if rates == 'equal':
            commands.append("lset rates=equal;")
        elif rates == 'gamma':
            commands.append(f"lset rates=gamma ngammacat={ngammacat};")
        elif rates == 'invgamma':
            commands.append(f"lset rates=invgamma ngammacat={ngammacat};")
        elif rates == 'propinv':
            commands.append("lset rates=propinv;")
        elif rates == 'lnorm':
            commands.append(f"lset rates=lnorm nlnormcat={ngammacat};")
        elif rates == 'adgamma':
            commands.append(f"lset rates=adgamma ngammacat={ngammacat};")
        
        # BEAGLE设置
        use_beagle = param_dict.get('use_beagle', True)
        if use_beagle:
            beagle_device = param_dict.get('beagle_device', 'cpu')
            beagle_precision = param_dict.get('beagle_precision', 'double')
            beagle_scaling = param_dict.get('beagle_scaling', 'dynamic')
            
            commands.append(f"set usebeagle=yes beagledevice={beagle_device} beagleprecision={beagle_precision} beaglescaling={beagle_scaling};")
        
        return commands
    
    def _get_mcmc_commands(self, params):
        """获取MCMC设置命令"""
        commands = []
        
        # 处理参数，可能是列表或字典
        if isinstance(params, dict):
            param_dict = params
        elif isinstance(params, list) and len(params) > 0 and isinstance(params[0], dict):
            param_dict = params[0]
        else:
            param_dict = {}

        ngen = param_dict.get('ngen', 1000000)
        samplefreq = param_dict.get('samplefreq', 1000)
        nchains = param_dict.get('nchains', 4)
        nruns = param_dict.get('nruns', 2)
        
        # 设置MCMC参数
        commands.append(f"mcmcp ngen={ngen} samplefreq={samplefreq} nchains={nchains} nruns={nruns} printfreq={samplefreq} savebrlens=yes checkpoint=yes checkfreq=5000;")
        
        return commands
    
    def _get_summary_commands(self, params):
        """获取总结命令"""
        commands = []
        
        # 处理参数，可能是列表或字典
        if isinstance(params, dict):
            param_dict = params
        elif isinstance(params, list) and len(params) > 0 and isinstance(params[0], dict):
            param_dict = params[0]
        else:
            param_dict = {}

        # Burn-in设置
        burnin_as_fraction = param_dict.get('burnin_as_fraction', True)
        if burnin_as_fraction:
            burnin_frac = param_dict.get('burnin_frac', 0.25)
            commands.append(f"sump relburnin=yes burninfrac={burnin_frac};")
        else:
            burnin_states = param_dict.get('burnin_states', 1000)
            commands.append(f"sump burnin={burnin_states};")
        
        # Consensus树设置
        contype = param_dict.get('contype', 'majorityrule')
        conformat = param_dict.get('conformat', 'figtree')
        
        if burnin_as_fraction:
            burnin_frac = param_dict.get('burnin_frac', 0.25)
            commands.append(f"sumt relburnin=yes burninfrac={burnin_frac} contype={contype} conformat={conformat};")
        else:
            burnin_states = param_dict.get('burnin_states', 1000)
            commands.append(f"sumt burnin={burnin_states} contype={contype} conformat={conformat};")
        
        return commands


class MrBayesPlugin(BasePlugin):
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_phylogeny_result_signal = pyqtSignal(dict)  # 导出系统发育树结果信号
    export_chain_result_signal = pyqtSignal(object)  # 导出MCMC链文件信号
    
    def __init__(self, import_from=None, import_data=None, workdir=None, imported_model=None, model_conversion_result=None, seq_type="DNA"):
        # 存储导入的模型信息 - 必须在 super().__init__() 之前，因为父类会调用 init_ui()
        self.imported_model = imported_model
        self.model_conversion_result = model_conversion_result
        self.imported_seq_type = seq_type
        
        super().__init__(import_from, import_data, workdir=workdir)
        
        # 初始化变量
        if not hasattr(self, 'imported_files'):
            self.imported_files = []
        if not hasattr(self, 'file_tags'):
            self.file_tags = []
        
        # 分区模式相关变量
        # 注意：use_partition_mode在UI中作为QCheckBox创建，这里不初始化
        self.partition_definitions = []
        self.partition_mode = PartitionMode.EL
        
        # 处理Dataset Manager导入的数据
        if import_from == "DATASET_MANAGER" and import_data is not None:
            self.handle_import_data(import_data)
    
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "MrBayes-MPI-BEAGLE Phylogeny"
        self.tool_name = "MrBayes-MPI-BEAGLE"
        self.citation = [
            """Ronquist, F., Teslenko, M., van der Mark, P., Ayres, D. L., Darling, A., Höhna, S., Larget, B., Liu, L., Suchard, M. A., & Huelsenbeck, J. P. (2012). MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. Systematic biology, 61(3), 539–542. https://doi.org/10.1093/sysbio/sys029""",
            """Daniel L Ayres, Michael P Cummings, Guy Baele, Aaron E Darling, Paul O Lewis, David L Swofford, John P Huelsenbeck, Philippe Lemey, Andrew Rambaut, Marc A Suchard, BEAGLE 3: Improved Performance, Scaling, and Usability for a High-Performance Computing Library for Statistical Phylogenetics, Systematic Biology, Volume 68, Issue 6, November 2019, Pages 1052–1061, https://doi.org/10.1093/sysbio/syz020"""
        ]
        self.input_types = {"PHYLIP": ["phy"], "NEXUS": ["nex", "nexus"], "Chain File": [""]}
        self.output_types = {"Chain File": [".chain"], "Tree File": [".con.tre"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)

        # File input
        input_group = QGroupBox("Input")
        input_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select sequence files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)

        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # QCheckBox
        # If input NEXUS file contains a MrBayes Data Block, run it directly?
        self.run_data_block = QCheckBox("Run MrBayes Data Block directly if present")
        self.run_data_block.setChecked(True)
        input_layout.addRow("", self.run_data_block)

        # Phylogenetic parameters
        phy_parameters_group = QGroupBox("Phylogenetic parameters")
        phy_parameters_layout = QFormLayout()
        phy_parameters_group.setLayout(phy_parameters_layout)
        layout.addWidget(phy_parameters_group)

        # Format of MrBayes data block
        # begin mrbayes;
        # ...
        # end;

        # Rate param numbers
        self.rate_params_layout = QHBoxLayout()

        self.datatype_combo = QComboBox()
        self.datatype_combo.addItems(["DNA", "Protein"]) # default: DNA
        # self.datatype_combo.currentIndexChanged.connect(self.on_datatype_changed)
        self.datatype_combo.setCurrentIndex(0)

        type_label = QLabel("Type:")
        type_label.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)

        self.rate_params_layout.addWidget(type_label)
        self.rate_params_layout.addWidget(self.datatype_combo)

        # for DNA data
        # lset nst = <dna_rate_num>, statefreq = <fixed(equal) / fixed(empirical)> // <none - estimated>;
        self.dnadata_widget = QWidget()
        self.dnadata_widget.setContentsMargins(0, 0, 0, 0)
        self.dnadata_layout = QHBoxLayout()
        self.dnadata_layout.setContentsMargins(0, 0, 0, 0)
        self.dnadata_widget.setLayout(self.dnadata_layout)

        self.dna_rate_num_spinbox = QSpinBox()
        self.dna_rate_num_spinbox.setRange(1, 6)
        self.dna_rate_num_spinbox.setValue(6)

        self.state_freq_pr_combo_dna = QComboBox()
        self.state_freq_pr_combo_dna.addItems(["estimated(dirichlet)", "fixed(equal)", "fixed(empirical)"])
        self.state_freq_pr_combo_dna.setCurrentIndex(0)

        self.dnadata_layout.addWidget(QLabel("Subst. num:"))
        self.dnadata_layout.addWidget(self.dna_rate_num_spinbox)
        self.dnadata_layout.addWidget(QLabel("State freq.:"))
        self.dnadata_layout.addWidget(self.state_freq_pr_combo_dna)

        self.rate_params_layout.addWidget(self.dnadata_widget)

        # for Protein data
        # prset aamodelpr = fixed(...) / mixed;
        self.prodata_widget = QWidget()
        self.prodata_widget.setContentsMargins(0, 0, 0, 0)
        self.prodata_widget.setVisible(False)
        self.prodata_layout = QHBoxLayout()
        self.prodata_layout.setContentsMargins(0, 0, 0, 0)
        self.prodata_widget.setLayout(self.prodata_layout)
        
        # self.prodata_layout.setVisible(False)
        self.prot_model_combo = QComboBox()
        self.prot_model_combo.addItems(["Blosum62", "Blosum", "Wag", "Lg", "gtr", "jones", "mtrev", "Poisson", "mixed"])

        self.prodata_layout.addWidget(QLabel("Model"))
        self.prodata_layout.addWidget(self.prot_model_combo)
        self.prodata_layout.addWidget(QLabel("State freq.:"))
        

        self.state_freq_pr_combo_prot = QComboBox()
        self.state_freq_pr_combo_prot.addItems(["estimated(dirichlet)", "fixed(equal)", "fixed(empirical)"])
        self.state_freq_pr_combo_prot.setCurrentIndex(0)

        self.prodata_layout.addWidget(self.state_freq_pr_combo_prot)

        self.rate_params_layout.addWidget(self.prodata_widget)

        phy_parameters_layout.addRow("Model Settings:", self.rate_params_layout)
        
        # Rate Heterogenity
        # lset rates = euqal/gamma/invgamma/propinv/lnorm/adgamma Ngammacat=<gamma categories>;
        rate_hetero_layout = QHBoxLayout()
        self.rate_hetero_combo = QComboBox()
        self.rate_hetero_combo.addItems(["Equal", "Gamma (+G)", "InvGamma (+G+I)", "PropInv (+I)", "Lognormal", "Adgamma"])

        self.gamma_categories_spinbox = QSpinBox()
        self.gamma_categories_spinbox.setRange(1, 100)
        self.gamma_categories_spinbox.setValue(4)

        # rate_hetero_layout.addWidget(QLabel("Rate Var.:"))
        rate_hetero_layout.addWidget(self.rate_hetero_combo)
        rate_hetero_layout.addWidget(QLabel("Gamma Categories:"))
        rate_hetero_layout.addWidget(self.gamma_categories_spinbox)

        phy_parameters_layout.addRow("Rate Heterogenity:", rate_hetero_layout)


        
        # MPI & BEAGLE settings button
        mpi_beagle_layout = QHBoxLayout()
        self.mpi_beagle_btn = QPushButton("MPI & BEAGLE Settings...")
        self.mpi_beagle_btn.clicked.connect(self.open_mpi_beagle_dialog)
        mpi_beagle_layout.addWidget(self.mpi_beagle_btn)
        mpi_beagle_layout.addStretch()
        layout.addLayout(mpi_beagle_layout)

        # Initialize MPI & BEAGLE settings with default values
        self.use_mpi = True
        self.use_beagle = True
        self.beagle_device = 'CPU'
        self.beagle_precision = 'Double'
        self.beagle_scaling = 'Dynamic'

        # Partition mode settings
        partition_layout = QHBoxLayout()
        self.use_partition_mode = QCheckBox("Enable Partition Mode")
        self.use_partition_mode.setChecked(False)
        self.use_partition_mode.stateChanged.connect(self.on_partition_mode_toggled)
        
        self.partition_config_btn = QPushButton("Configure Partitions...")
        self.partition_config_btn.setEnabled(False)
        self.partition_config_btn.clicked.connect(self.open_partition_config_dialog)
        
        partition_layout.addWidget(self.use_partition_mode)
        partition_layout.addWidget(self.partition_config_btn)
        partition_layout.addStretch()
        layout.addLayout(partition_layout)
        
        # 分区状态标签（初始隐藏）
        self.partition_status_label = QLabel("No partitions defined")
        self.partition_status_label.setStyleSheet("color: #6c757d; font-style: italic;")
        self.partition_status_label.setVisible(False)
        layout.addWidget(self.partition_status_label)
        
        # 分子钟定年设置
        dating_layout = QHBoxLayout()
        self.enable_dating = QCheckBox("Enable Phylogenetic Dating")
        self.enable_dating.setChecked(False)
        self.enable_dating.stateChanged.connect(self.on_dating_toggled)
        
        self.dating_settings_btn = QPushButton("Dating Settings...")
        self.dating_settings_btn.setEnabled(False)
        self.dating_settings_btn.clicked.connect(self.open_dating_settings_dialog)
        
        dating_layout.addWidget(self.enable_dating)
        dating_layout.addWidget(self.dating_settings_btn)
        dating_layout.addStretch()
        layout.addLayout(dating_layout)
        
        # 分子钟定年状态标签
        self.dating_status_label = QLabel("Dating disabled")
        self.dating_status_label.setStyleSheet("color: #6c757d; font-style: italic;")
        self.dating_status_label.setVisible(False)
        layout.addWidget(self.dating_status_label)
        
        # 初始化分子钟定年参数
        self.dating_enabled = False
        self.clock_model = "strict"
        self.clock_params = {}
        self.calibrations = []
        self.outgroup = []  # 改为列表，支持多个外群

        mcmc_params_group = QGroupBox("MCMC settings")
        mcmc_params_layout = QFormLayout()
        mcmc_params_group.setLayout(mcmc_params_layout)

        layout.addWidget(mcmc_params_group)

        # generation; sampling frequency; run num; chain num;
        # mcmcp ngen=* samplefreq=* printfreq=<samplefreq> nchains=* nruns=* savebrlens=yes checkpoint=yes checkfreq=5000;
        self.generation_spinbox = QLineEdit()
        self.generation_spinbox.setText("1000000")
        self.generation_spinbox.setMinimumWidth(100)

        self.sampling_frequency_spinbox = QLineEdit()
        self.sampling_frequency_spinbox.setText("1000")
        self.sampling_frequency_spinbox.setMinimumWidth(100)

        self.run_num_spinbox = QSpinBox()
        self.run_num_spinbox.setMinimum(1)
        self.run_num_spinbox.setValue(2)

        self.chain_num_spinbox = QSpinBox()
        self.chain_num_spinbox.setMinimum(1)
        self.chain_num_spinbox.setValue(4) # MC^3

        mcmcp_widget = QWidget()
        mcmcp_widget.setContentsMargins(0, 0, 0, 0)
        mcmcp_layout = QVBoxLayout()
        mcmcp_layout.setContentsMargins(0, 0, 0, 0)
        mcmcp_layout.setSpacing(5)
        mcmcp_widget.setLayout(mcmcp_layout)

        # 第一行：Generations和Sampling Freq.
        row1_widget = QWidget()
        row1_widget.setContentsMargins(0, 0, 0, 0)
        row1_layout = QHBoxLayout()
        row1_layout.setContentsMargins(0, 0, 0, 0)
        row1_widget.setLayout(row1_layout)

        row1_layout.addWidget(QLabel("Generations:"))
        row1_layout.addWidget(self.generation_spinbox)
        row1_layout.addWidget(QLabel("Sampling Freq.:"))
        row1_layout.addWidget(self.sampling_frequency_spinbox)
        # row1_layout.addStretch()

        # 第二行：Runs和Chains
        row2_widget = QWidget()
        row2_widget.setContentsMargins(0, 0, 0, 0)
        row2_layout = QHBoxLayout()
        row2_layout.setContentsMargins(0, 0, 0, 0)
        row2_widget.setLayout(row2_layout)

        row2_layout.addWidget(QLabel("Runs:"))
        row2_layout.addWidget(self.run_num_spinbox)
        row2_layout.addWidget(QLabel("Chains:"))
        row2_layout.addWidget(self.chain_num_spinbox)
        # row2_layout.addStretch()

        mcmcp_layout.addWidget(row1_widget)
        mcmcp_layout.addWidget(row2_widget)

        mcmc_params_layout.addRow("MCMC:", mcmcp_widget)

        # summary and consensus
        # sump <relburnin=yes burninfrac=*> / <burnin=*>;
        # sumt <relburnin=yes burninfrac=*> / <burnin=*> conformat=<Simple/FigTree> contype=*;

        sum_widget = QWidget()
        sum_layout = QHBoxLayout()
        sum_widget.setContentsMargins(0, 0, 0, 0)
        sum_layout.setContentsMargins(0, 0, 0, 0)
        sum_widget.setLayout(sum_layout)

        self.contree_type = QComboBox()
        self.contree_type.addItems(["Majority-rule", "Halfcompat", "Allcompat"])

        self.contree_format = QComboBox()
        self.contree_format.addItems(["FigTree", "Simple"])

        sum_layout.addWidget(QLabel("Method:"))
        sum_layout.addWidget(self.contree_type)
        sum_layout.addWidget(QLabel("Output format:"))
        sum_layout.addWidget(self.contree_format)

        mcmc_params_layout.addRow("Consensus:", sum_widget)

        burnin_widget = QWidget()
        burnin_layout = QHBoxLayout()
        burnin_widget.setContentsMargins(0, 0, 0, 0)
        burnin_layout.setContentsMargins(0, 0, 0, 0)
        burnin_widget.setLayout(burnin_layout)

        self.burnin_as_fraction = QRadioButton("Fraction:")
        self.burnin_as_fraction.setChecked(True)

        self.burnin_fraction = QDoubleSpinBox()
        self.burnin_fraction.setRange(0, 1)
        self.burnin_fraction.setSingleStep(0.01)
        self.burnin_fraction.setValue(0.25)

        self.burnin_as_states = QRadioButton("States:")
        self.burnin_as_states.setChecked(False)

        self.burnin_states = QSpinBox()
        self.burnin_states.setRange(0, 2147483647)
        self.burnin_states.setValue(1000)

        burnin_layout.addWidget(self.burnin_as_fraction)
        burnin_layout.addWidget(self.burnin_fraction)
        burnin_layout.addWidget(self.burnin_as_states)
        burnin_layout.addWidget(self.burnin_states)

        mcmc_params_layout.addRow("Burn-in:", burnin_widget)



        mb_data_block_widget = QWidget()
        mb_data_block_widget.setContentsMargins(0,0,0,0)
        mb_data_block_layout = QHBoxLayout()
        mb_data_block_layout.setContentsMargins(0,0,0,0)
        mb_data_block_widget.setLayout(mb_data_block_layout)

        show_mb_data_block = QPushButton("Show MrBayes data block") # use QDialog to showcase MrBayes data block (Font: Consolas; same style as Console)
        show_mb_data_block.clicked.connect(self.show_mb_data_block)
        copy_mb_data_block = QPushButton("Copy MrBayes data block")
        copy_mb_data_block.clicked.connect(self.copy_mb_data_block)

        mb_data_block_layout.addWidget(show_mb_data_block)
        mb_data_block_layout.addWidget(copy_mb_data_block)

        layout.addWidget(mb_data_block_widget)
        
        # 连接数据类型切换
        self.datatype_combo.currentIndexChanged.connect(self.on_datatype_changed)
        
        # 处理导入的数据
        if self.import_file:
            self.file_path_edit.setText(self.import_file)
            input_group.setVisible(False)   
        elif hasattr(self, 'imported_files') and self.imported_files:
            # 显示导入的文件
            for file_path in self.imported_files:
                self.add_file_tag(file_path)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
        
        # 应用导入的模型参数
        self._apply_imported_model()
    
    def _apply_imported_model(self):
        """应用导入的模型参数到UI"""
        if not self.imported_model:
            return
        
        from .partition_mode import MrBayesModelConverter
        
        # 解析模型
        base_model, params = MrBayesModelConverter.parse_modelfinder_model(self.imported_model)
        
        if self.imported_seq_type.upper() == "DNA":
            # DNA模型
            if self.model_conversion_result and isinstance(self.model_conversion_result, dict):
                # 应用转换后的参数
                mrbayes_params = self.model_conversion_result
                
                # 设置数据类型为DNA
                self.datatype_combo.setCurrentIndex(0)
                
                # 设置nst
                if "nst" in mrbayes_params:
                    self.dna_rate_num_spinbox.setValue(mrbayes_params["nst"])
                
                # 设置statefreq
                if "statefreq" in mrbayes_params:
                    statefreq = mrbayes_params["statefreq"]
                    if statefreq == "fixed(equal)":
                        self.state_freq_pr_combo_dna.setCurrentIndex(1)
                    elif statefreq == "fixed(empirical)":
                        self.state_freq_pr_combo_dna.setCurrentIndex(2)
                    else:  # dirichlet
                        self.state_freq_pr_combo_dna.setCurrentIndex(0)
                
                # 设置速率异质性
                if "rates" in mrbayes_params:
                    rates = mrbayes_params["rates"]
                    if rates == "equal":
                        self.rate_hetero_combo.setCurrentIndex(0)
                    elif rates == "gamma":
                        self.rate_hetero_combo.setCurrentIndex(1)
                    elif rates == "invgamma":
                        self.rate_hetero_combo.setCurrentIndex(2)
                    elif rates == "propinv":
                        self.rate_hetero_combo.setCurrentIndex(3)
                    elif rates == "lnorm":
                        self.rate_hetero_combo.setCurrentIndex(4)
                    elif rates == "adgamma":
                        self.rate_hetero_combo.setCurrentIndex(5)
                
                # 设置gamma类别数
                if "ngammacat" in mrbayes_params:
                    self.gamma_categories_spinbox.setValue(mrbayes_params["ngammacat"])
        
        elif self.imported_seq_type.upper() == "PROTEIN":
            # 蛋白质模型
            if self.model_conversion_result:
                # 设置数据类型为Protein
                self.datatype_combo.setCurrentIndex(1)
                
                # 设置氨基酸模型
                if self.model_conversion_result == "mixed":
                    self.prot_model_combo.setCurrentText("mixed")
                else:
                    # 查找模型在组合框中的位置
                    model_index = self.prot_model_combo.findText(self.model_conversion_result, Qt.MatchFixedString)
                    if model_index >= 0:
                        self.prot_model_combo.setCurrentIndex(model_index)
    
    def handle_import_data(self, import_data):
        """处理从Dataset Manager导入的数据"""
        if not isinstance(import_data, dict):
            return
        
        dataset_items = import_data.get('dataset_items', [])
        dataset_config = import_data.get('dataset_config', {})
        
        # 筛选出 alignment 类型的 items
        from ..platforms.methods.dataset_models import ITEM_TYPE_ALIGNMENT
        selected_items = [item for item in dataset_items if item.item_type == ITEM_TYPE_ALIGNMENT]
        
        if len(selected_items) == 0:
            # 没有alignment，不导入
            return
        
        if len(selected_items) == 1:
            # 单一alignment，直接使用
            self.import_file = selected_items[0].file_path
            self.imported_files = [selected_items[0].file_path]
        else:
            # 多个alignment，合并成超级矩阵
            try:
                # 检查所有alignment是否都已对齐
                unaligned_items = [item for item in selected_items if not item.is_aligned]
                if unaligned_items:
                    unaligned_names = [item.loci_name for item in unaligned_items]
                    warning_msg = f"The following partitions are not aligned:\n"
                    warning_msg += "\n".join(f"  - {name}" for name in unaligned_names)
                    warning_msg += "\n\nPlease align these partitions first or select only aligned partitions."
                    self.add_console_message(f"Warning: {warning_msg}", "warning")
                
                # 获取dataset settings
                topo_linked = dataset_config.get('topo_linked', False)
                edge_linked = dataset_config.get('edge_linked', False)
                
                # 映射到MrBayes分区模式
                # MrBayes只支持EL和EUL，不支持TUL（拓扑解链）
                if not topo_linked:
                    # TUL转换为EUL，给出警告
                    self.add_console_message("⚠️ 转换警告: MrBayes不支持拓扑解链（TUL）模式", "warning")
                    self.add_console_message("⚠️ 转换说明: 已转换为边解链（EUL）模式", "warning")
                    self.partition_mode = PartitionMode.EUL
                elif edge_linked:
                    # TL转换为EL，给出警告
                    self.add_console_message("⚠️ 转换说明: TL模式转换为EL模式", "warning")
                    self.partition_mode = PartitionMode.EL
                else:
                    self.partition_mode = PartitionMode.EUL
                
                # 合并所有选中的 partition 为 supermatrix
                supermatrix_sequences = {}
                partition_definitions = []
                current_pos = 1
                
                # 获取所有序列名称
                if selected_items:
                    all_seq_names = [seq.id for seq in selected_items[0].sequences]
                else:
                    all_seq_names = []
                
                # 合并序列
                for seq_name in all_seq_names:
                    supermatrix_seq = ""
                    for item in selected_items:
                        # 找到对应的序列
                        seq = next((s for s in item.sequences if s.id == seq_name), None)
                        if seq:
                            supermatrix_seq += str(seq.seq)
                        else:
                            # 如果某个 partition 缺少该序列，用 ? 填充
                            supermatrix_seq += "?" * item.length
                    supermatrix_sequences[seq_name] = supermatrix_seq
                
                # 创建 supermatrix 临时文件
                temp_file = self.create_temp_file(suffix='.phy')
                
                # 写入PHYLIP格式（MrBayes推荐格式）
                with open(temp_file, 'w') as f:
                    # 第一行：序列数量和序列长度
                    num_sequences = len(supermatrix_sequences)
                    seq_length = len(list(supermatrix_sequences.values())[0]) if supermatrix_sequences else 0
                    f.write(f"{num_sequences} {seq_length}\n")
                    
                    # 写入序列数据
                    for seq_name, seq_content in supermatrix_sequences.items():
                        # 序列名称限制为10个字符（PHYLIP格式）
                        short_name = seq_name[:10].ljust(10)
                        f.write(f"{short_name} {seq_content}\n")
                
                self.import_file = temp_file
                self.imported_files = [temp_file]
                
                # 计算 partition 坐标并创建 MrBayes 分区定义
                for item in selected_items:
                    end_pos = current_pos + item.length - 1
                    
                    # 检测序列类型
                    seq_type = self._detect_sequence_type(item)
                    
                    # 创建MrBayes分区定义
                    partition_def = MrBayesPartitionDefinition(
                        name=item.loci_name,
                        range=f"{current_pos}-{end_pos}",
                        seq_type=seq_type,
                        nst=6,  # 默认GTR
                        aamodel="mixed",  # 默认mixed
                        rates="gamma",  # 默认gamma
                        ngammacat=4  # 默认4个gamma类别
                    )
                    partition_definitions.append(partition_def)
                    current_pos = end_pos + 1
                
                self.partition_definitions = partition_definitions
                
                # 启用分区模式
                if hasattr(self, 'use_partition_mode'):
                    self.use_partition_mode.setChecked(True)
                    self.partition_config_btn.setEnabled(True)
                
                # 更新分区状态显示
                self.update_partition_status()
                
                # 添加控制台消息
                self.add_console_message(f"Dataset imported: {len(selected_items)} partitions", "info")
                self.add_console_message(f"Partition mode: {self.partition_mode.value}", "info")
                self.add_console_message(f"Supermatrix created: {len(supermatrix_sequences)} taxa", "info")
                
            except Exception as e:
                self.add_console_message(f"Failed to import dataset: {str(e)}", "error")
                import traceback
                self.add_console_message(traceback.format_exc(), "error")
                self.import_file = None
                self.imported_files = []
                self.partition_definitions = []
        
        # 更新UI显示导入的文件
        if hasattr(self, 'file_path_edit') and self.file_path_edit and self.import_file:
            self.file_path_edit.setText(self.import_file)
    
    def _detect_sequence_type(self, dataset_item) -> str:
        """从DatasetItem检测序列类型"""
        # 检查metadata
        if hasattr(dataset_item, 'metadata') and 'seq_type' in dataset_item.metadata:
            seq_type = dataset_item.metadata['seq_type'].upper()
            if seq_type in ['DNA', 'AA', 'PROTEIN', 'CODON']:
                if seq_type == 'AA' or seq_type == 'PROTEIN':
                    return 'PROTEIN'
                return 'DNA'
        
        # 检查data字段
        if hasattr(dataset_item, 'data') and 'seq_type' in dataset_item.data:
            seq_type = dataset_item.data['seq_type'].upper()
            if seq_type in ['DNA', 'AA', 'PROTEIN', 'CODON']:
                if seq_type == 'AA' or seq_type == 'PROTEIN':
                    return 'PROTEIN'
                return 'DNA'
        
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
                return 'PROTEIN'
        
        # 默认返回 DNA
        return 'DNA'
    
    def update_partition_status(self):
        """更新分区状态显示"""
        if not hasattr(self, 'partition_status_label'):
            return
        
        if not self.partition_definitions:
            self.partition_status_label.setText("No partitions defined")
            self.partition_status_label.setStyleSheet("color: #6c757d; font-style: italic;")
        else:
            status_text = f"{len(self.partition_definitions)} partition(s) configured"
            status_text += f" [{self.partition_mode.value} mode]"
            self.partition_status_label.setText(status_text)
            self.partition_status_label.setStyleSheet("color: #28a745; font-weight: bold;")
    
    def apply_partition_model_config(self, partition_model_config):
        """
        应用分区模型配置（用模型信息覆盖分区定义）
        
        Args:
            partition_model_config: dict {
                'partition_mode': str,
                'partitions': list,
                'model_count': int
            }
        """
        from .partition_mode import MrBayesModelConverter
        
        # 更新分区模式
        mode_map = {
            'EL': PartitionMode.EL,
            'TL': PartitionMode.EL,  # MrBayes不支持TL，转换为EL
            'EUL': PartitionMode.EUL,
            'TUL': PartitionMode.EUL  # MrBayes不支持TUL，转换为EUL
        }
        input_mode = partition_model_config.get('partition_mode', 'EL')
        converted_mode = mode_map.get(input_mode, PartitionMode.EL)
        
        if converted_mode != input_mode:
            self.add_console_message(f"⚠️ 分区模式转换: {input_mode} -> {converted_mode.value} (MrBayes不支持{input_mode})", "warning")
        
        self.partition_mode = converted_mode
        
        # 用模型信息更新分区定义
        model_partitions = partition_model_config.get('partitions', [])
        if model_partitions and self.partition_definitions:
            # 确保分区数量匹配
            if len(model_partitions) == len(self.partition_definitions):
                for i, model_part in enumerate(model_partitions):
                    # 支持两种字段名：'model' 和 'best_model'
                    model_code = model_part.get('model') or model_part.get('best_model')
                    if model_code:
                        partition = self.partition_definitions[i]
                        
                        # 使用MrBayesModelConverter转换模型
                        if partition.seq_type == "DNA":
                            mrbayes_params, warnings = MrBayesModelConverter.convert_model_to_mrbayes(
                                model_code, "DNA"
                            )
                            if warnings:
                                for warning in warnings:
                                    self.add_console_message(f"⚠️ 分区{i+1}警告: {warning}", "warning")
                            
                            if mrbayes_params:
                                partition.nst = mrbayes_params.get('nst', 6)
                                partition.rates = mrbayes_params.get('rates', 'gamma')
                                partition.ngammacat = mrbayes_params.get('ngammacat', 4)
                            
                            self.add_console_message(f"  Partition {i+1}: {partition.name} -> {model_code}", "info")
                        else:  # Protein
                            mrbayes_model, warnings = MrBayesModelConverter.convert_model_to_mrbayes(
                                model_code, "PROTEIN"
                            )
                            if warnings:
                                for warning in warnings:
                                    self.add_console_message(f"⚠️ 分区{i+1}警告: {warning}", "warning")
                            
                            if mrbayes_model:
                                partition.aamodel = mrbayes_model
                            
                            self.add_console_message(f"  Partition {i+1}: {partition.name} -> {model_code} ({mrbayes_model})", "info")
                    
                    if model_part.get('seq_type'):
                        # 更新序列类型
                        seq_type = model_part['seq_type'].upper()
                        if seq_type in ['DNA', 'PROTEIN', 'AA']:
                            if seq_type == 'AA':
                                seq_type = 'PROTEIN'
                            self.partition_definitions[i].seq_type = seq_type
            else:
                # 分区数量不匹配，警告用户
                self.add_console_message(
                    f"Warning: Model has {len(model_partitions)} partitions but dataset has {len(self.partition_definitions)} partitions. Partition count mismatch!",
                    "warning"
                )
        
        # 更新分区状态显示
        self.update_partition_status()

    # ============ MrBayesPlugin 业务逻辑方法 ============
    
    def on_datatype_changed(self, index):
        """数据类型改变时的处理"""
        if index == 0:  # DNA
            self.dnadata_widget.setVisible(True)
            self.prodata_widget.setVisible(False)
        else:  # Protein
            self.dnadata_widget.setVisible(False)
            self.prodata_widget.setVisible(True)
    
    def open_mpi_beagle_dialog(self):
        """打开MPI & BEAGLE设置对话框"""
        dialog = MPIBeagleSettingsDialog(
            self,
            use_mpi=self.use_mpi,
            use_beagle=self.use_beagle,
            beagle_device=self.beagle_device,
            beagle_precision=self.beagle_precision,
            beagle_scaling=self.beagle_scaling
        )
        
        if dialog.exec_() == QDialog.Accepted:
            settings = dialog.get_settings()
            self.use_mpi = settings['use_mpi']
            self.use_beagle = settings['use_beagle']
            self.beagle_device = settings['beagle_device']
            self.beagle_precision = settings['beagle_precision']
            self.beagle_scaling = settings['beagle_scaling']
    
    def on_partition_mode_toggled(self):
        """分区模式复选框状态改变时处理"""
        enabled = self.use_partition_mode.isChecked()
        self.partition_config_btn.setEnabled(enabled)
    
    def open_partition_config_dialog(self):
        """打开分区配置对话框"""
        dialog = MrBayesPartitionDialog(
            self,
            partitions=self.partition_definitions,
            mode=self.partition_mode
        )
        
        if dialog.exec_() == QDialog.Accepted:
            self.partition_definitions = dialog.get_partitions()
            self.partition_mode = dialog.get_mode()
            self.add_console_message(
                f"Partition mode enabled with {len(self.partition_definitions)} partition(s), "
                f"mode: {self.partition_mode.value}", "info"
            )
    
    def on_dating_toggled(self):
        """分子钟定年复选框状态改变时处理"""
        enabled = self.enable_dating.isChecked()
        self.dating_settings_btn.setEnabled(enabled)
        self.dating_enabled = enabled
        
        if not enabled:
            self.dating_status_label.setVisible(False)
        else:
            self.dating_status_label.setVisible(True)
            self.update_dating_status()
    
    def open_dating_settings_dialog(self):
        """打开分子钟定年设置对话框"""
        # 从输入文件中提取taxa列表
        taxa_list = self._extract_taxa_from_input()
        
        # 创建对话框
        dialog = MolecularDatingSettingsDialog(
            self,
            clock_model=self.clock_model,
            clock_params=self.clock_params.copy(),
            calibrations=self.calibrations.copy(),
            outgroup=self.outgroup,
            taxa_list=taxa_list
        )
        
        if dialog.exec_() == QDialog.Accepted:
            # 获取设置
            self.clock_model = dialog.get_clock_model()
            self.clock_params = dialog.get_clock_params()
            self.calibrations = dialog.get_calibrations()
            self.outgroup = dialog.get_outgroup()
            
            # 更新状态显示
            self.update_dating_status()
            
            # 添加控制台消息
            calib_count = len(self.calibrations)
            self.add_console_message(
                f"Phylogenetic dating settings updated: clock={self.clock_model}, "
                f"calibrations={calib_count}", "info"
            )
    
    def _extract_taxa_from_input(self):
        """从输入文件中提取taxa列表"""
        taxa_list = []
        
        try:
            # 优先从导入的文件中提取
            files_to_check = []
            if self.imported_files:
                files_to_check = self.imported_files
            elif self.import_file:
                files_to_check = [self.import_file]
            
            for file_path in files_to_check:
                if not os.path.exists(file_path):
                    continue
                
                # 根据文件类型解析
                file_ext = os.path.splitext(file_path)[1].lower()
                
                if file_ext in ['.phy', '.phylip']:
                    # PHYLIP格式
                    taxa_list = self._extract_taxa_from_phylip(file_path)
                elif file_ext in ['.fas', '.fasta', '.fa', '.fna']:
                    # FASTA格式
                    taxa_list = self._extract_taxa_from_fasta(file_path)
                elif file_ext in ['.nex', '.nexus']:
                    # NEXUS格式
                    taxa_list = self._extract_taxa_from_nexus(file_path)
                
                # 如果成功提取，就使用第一个文件的taxa
                if taxa_list:
                    break
                    
        except Exception as e:
            print(f"Error extracting taxa from input: {e}")
        
        return taxa_list
    
    def _extract_taxa_from_phylip(self, file_path):
        """从PHYLIP文件中提取taxa列表"""
        taxa_list = []
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                # 跳过第一行（序列数量和长度）
                for line in lines[1:]:
                    line = line.strip()
                    if not line:
                        continue
                    # PHYLIP格式：名称和序列用空格分隔
                    # 提取完整的名称（第一个空格之前的部分）
                    parts = line.split()
                    if parts:
                        taxon = parts[0]
                        # 清理名称，确保与MrBayes使用的名称一致
                        taxon = self._clean_taxon_name(taxon)
                        if taxon:
                            taxa_list.append(taxon)
        except Exception as e:
            print(f"Error parsing PHYLIP file: {e}")
        return taxa_list
    
    def _extract_taxa_from_fasta(self, file_path):
        """从FASTA文件中提取taxa列表"""
        taxa_list = []
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        taxon = line[1:].split()[0]  # 去掉>并取第一个词
                        # 清理名称，确保与MrBayes使用的名称一致
                        taxon = self._clean_taxon_name(taxon)
                        if taxon:
                            taxa_list.append(taxon)
        except Exception as e:
            print(f"Error parsing FASTA file: {e}")
        return taxa_list
    
    def _extract_taxa_from_nexus(self, file_path):
        """从NEXUS文件中提取taxa列表"""
        taxa_list = []
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
                # 查找matrix块
                matrix_match = re.search(r'matrix\s*(.*?)end\s*;', content, re.IGNORECASE | re.DOTALL)
                if matrix_match:
                    matrix_content = matrix_match.group(1)
                    for line in matrix_content.split('\n'):
                        line = line.strip()
                        if not line or line.startswith(';'):
                            continue
                        # 提取分类名称
                        parts = line.split()
                        if parts:
                            taxon = parts[0]
                            if taxon:
                                taxa_list.append(taxon)
        except Exception as e:
            print(f"Error parsing NEXUS file: {e}")
        return taxa_list
    
    def update_dating_status(self):
        """更新分子钟定年状态显示"""
        if not hasattr(self, 'dating_status_label'):
            return
        
        if not self.dating_enabled:
            self.dating_status_label.setText("Dating disabled")
            self.dating_status_label.setStyleSheet("color: #6c757d; font-style: italic;")
        else:
            status_text = f"Clock: {self.clock_model.upper()}"
            calib_count = len(self.calibrations)
            status_text += f" | Calibrations: {calib_count}"
            if self.outgroup:
                status_text += f" | Outgroup: {self.outgroup}"
            
            self.dating_status_label.setText(status_text)
            self.dating_status_label.setStyleSheet("color: #28a745; font-weight: bold;")
    
    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.phy *.phylip *.nex *.nexus);;All files (*)"
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select alignment files", "", file_filter)
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def add_file_tag(self, file_path):
        """添加文件标签"""
        self.imported_files.append(file_path)
        
        # Create file tag widget
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.Box)
        tag_widget.setStyleSheet("""
            QFrame {
                background-color: #e9ecef;
                border-radius: 15px;
                margin: 2px;
            }
        """)
        
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)
        
        # Get display name (handle duplicate names)
        display_name = self.get_display_name(file_path)
        
        # File name label
        name_label = QLabel(display_name)
        name_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(name_label)
        
        # Close button
        close_btn = QPushButton("×")
        close_btn.setFixedSize(20, 20)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: transparent;
                border: none;
                color: #6c757d;
                font-weight: bold;
                font-size: 14px;
            }
            QPushButton:hover {
                color: #dc3545;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_file_tag(tag_widget, file_path))
        tag_layout.addWidget(close_btn)
        
        # Add to container
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append(tag_widget)
        
        # Show container
        self.file_tags_container.setVisible(True)
    
    def remove_file_tag(self, tag_widget, file_path):
        """移除文件标签"""
        # Remove from lists
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        if tag_widget in self.file_tags:
            self.file_tags.remove(tag_widget)
        
        # Remove widget
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Update file path display
        if len(self.imported_files) == 0:
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
        elif len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        else:
            self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def get_display_name(self, file_path):
        """获取文件显示名称（处理重名情况）"""
        base_name = os.path.basename(file_path)
        count = sum(1 for f in self.imported_files if os.path.basename(f) == base_name)
        if count > 1:
            name, ext = os.path.splitext(base_name)
            dir_name = os.path.basename(os.path.dirname(file_path))
            return f"{name} ({dir_name}){ext}"
        return base_name
    
    def prepare_input_files(self):
        """准备输入文件"""
        try:
            input_files = []
            
            # 如果有导入的文件，使用它们
            if self.imported_files:
                for file_path in self.imported_files:
                    if os.path.exists(file_path):
                        input_files.append(file_path)
                    else:
                        QMessageBox.warning(self, "Warning", f"File does not exist: {file_path}")
                return input_files
            elif self.import_file:
                return [self.import_file]
            
            # 否则，从文本输入创建临时文件（如果有的话）
            if hasattr(self, 'sequence_text') and self.sequence_text.toPlainText().strip():
                sequence_text = self.sequence_text.toPlainText().strip()
                temp_file = self.create_temp_file(suffix='.phy')
                with open(temp_file, 'w') as f:
                    f.write(sequence_text)
                return [temp_file]
            
            return None
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
    
    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 数据类型
        datatype = self.datatype_combo.currentText()
        params.append({'datatype': datatype})
        
        # 模型参数
        if datatype == "DNA":
            # DNA模型
            nst = self.dna_rate_num_spinbox.value()
            statefreq = self.state_freq_pr_combo_dna.currentText()
            params[-1].update({
                'nst': nst,
                'statefreq': statefreq
            })
        else:
            # 蛋白质模型
            aamodel = self.prot_model_combo.currentText().lower()
            statefreq = self.state_freq_pr_combo_prot.currentText()
            params[-1].update({
                'aamodel': aamodel,
                'statefreq': statefreq
            })
        
        # 速率异质性
        rates_mapping = {
            "Equal": "equal",
            "Gamma (+G)": "gamma",
            "InvGamma (+G+I)": "invgamma",
            "PropInv (+I)": "propinv",
            "Lognormal": "lnorm",
            "Adgamma": "adgamma"
        }
        rates = rates_mapping.get(self.rate_hetero_combo.currentText(), "equal")
        ngammacat = self.gamma_categories_spinbox.value()
        params[-1].update({
            'rates': rates,
            'ngammacat': ngammacat
        })
        
        # BEAGLE设置
        params[-1].update({
            'use_beagle': self.use_beagle,
            'beagle_device': self.beagle_device.lower(),
            'beagle_precision': self.beagle_precision.lower(),
            'beagle_scaling': self.beagle_scaling.lower()
        })
        
        # MCMC参数
        nchains = self.chain_num_spinbox.value()
        params[-1].update({
            'ngen': int(self.generation_spinbox.text()),
            'samplefreq': int(self.sampling_frequency_spinbox.text()),
            'nchains': nchains,
            'nruns': self.run_num_spinbox.value()
        })
        
        # Consensus参数
        contype_mapping = {
            "Majority-rule": "majorityrule",
            "Halfcompat": "halfcompat",
            "Allcompat": "allcompat"
        }
        contype = contype_mapping.get(self.contree_type.currentText(), "majorityrule")
        conformat = self.contree_format.currentText().lower()
        
        # Burn-in参数
        burnin_as_fraction = self.burnin_as_fraction.isChecked()
        if burnin_as_fraction:
            burnin_frac = self.burnin_fraction.value()
            params[-1].update({
                'burnin_as_fraction': True,
                'burnin_frac': burnin_frac,
                'contype': contype,
                'conformat': conformat
            })
        else:
            burnin_states = self.burnin_states.value()
            params[-1].update({
                'burnin_as_fraction': False,
                'burnin_states': burnin_states,
                'contype': contype,
                'conformat': conformat
            })
        
        # 分子钟定年参数
        params[-1].update({
            'enable_dating': self.dating_enabled,
            'clock_model': self.clock_model,
            'clock_params': self.clock_params,
            'calibrations': self.calibrations,
            'outgroup': self.outgroup
        })
        
        return params
    
    def run_analysis(self):
        """运行分析"""
        # 检查输入
        if not self.imported_files and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide alignment files!")
            return
            
        # 检查MrBayes可执行文件是否存在
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "MrBayes executable file not found! Please check config.json.")
            return
            
        # 检查MPI可执行文件是否存在（如果使用MPI）
        mpirun_path = None
        if self.use_mpi:
            mpirun_path = self.get_mpirun_path()
            if not mpirun_path or not os.path.exists(mpirun_path):
                QMessageBox.critical(self, "Error", "MPI executable file not found! Please check config.json for 'MPIRun' entry.")
                return
            
        # 添加控制台消息
        self.add_console_message("Starting phylogenetic inference with MrBayes...", "info")
        
        # 准备输入文件
        input_files = self.prepare_input_files()
        if not input_files:
            return
            
        # 更新UI状态
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 未知进度
        
        # 获取参数
        params = self.get_parameters()
        
        # 在单独的线程中运行MrBayes
        self.analysis_thread = MrBayesThread(
            self.tool_path, mpirun_path, input_files, params, self.imported_files,
            self.run_data_block.isChecked(),
            use_partition_mode=self.use_partition_mode.isChecked(),
            partition_definitions=self.partition_definitions,
            partition_mode=self.partition_mode,
            workdir=self.workdir,
            use_mpi=self.use_mpi
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()
    
    def get_mpirun_path(self):
        """获取MPI运行程序路径"""
        config_path = os.path.join(self.plugin_path, 'config.json')
        
        if not os.path.exists(config_path):
            return None
            
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                config = json.load(f)
            
            # 查找MPI运行程序路径
            for tool in config:
                if tool.get('name', '').lower() in ['mpirun', 'mpiexec', 'mpi']:
                    return os.path.join(self.plugin_path, './'+tool['path'])
                    
        except Exception as e:
            print(f"Error reading MPI path from config: {e}")
            
        return None
    
    def analysis_finished(self, output_files, html_files):
        """分析完成处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 保存输出文件
        self.current_output_files = output_files
        
        # 检测并发送MCMC链文件
        self._detect_and_emit_chain_files(output_files)
        
        # 检测并发送系统发育树文件
        self._detect_and_emit_phylogeny_files(output_files)
        
        # 显示结果
        self.display_results(output_files)
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        # 添加控制台消息
        self.add_console_message(f"Phylogenetic inference completed successfully! Found {len(output_files)} result file(s)", "info")
        
        QMessageBox.information(self, "Completed", "Phylogenetic inference completed!")
    
    def _detect_and_emit_phylogeny_files(self, output_files):
        """检测系统发育树文件并发送信号"""
        # 查找.con.tre文件（共识树文件）
        con_tre_files = [f for f in output_files if f.endswith('.con.tre')]
        
        if con_tre_files:
            tree_file = con_tre_files[0]
            # 创建系统发育树数据字典
            phylogeny_data = {
                'file_path': tree_file,
                'file_type': 'newick',
                'tool': 'mrbayes',
                'tree_type': 'consensus'
            }
            # 发送信号
            self.export_phylogeny_result_signal.emit(phylogeny_data)
            self.add_console_message(f"Detected consensus tree: {os.path.basename(tree_file)}", "info")
    
    def _detect_and_emit_chain_files(self, output_files):
        """检测MCMC链文件并发送信号"""
        # 查找所有.p文件（MrBayes的chain文件）
        chain_files = [f for f in output_files if f.endswith('.p')]
        
        if not chain_files:
            return
        
        # 统计run数量
        run_numbers = set()
        for chain_file in chain_files:
            filename = os.path.basename(chain_file)
            if '.run' in filename:
                try:
                    run_part = filename.split('.run')[1].split('.')[0]
                    run_number = int(run_part)
                    run_numbers.add(run_number)
                except (ValueError, IndexError):
                    pass
        
        run_count = len(run_numbers)
        chain_count = len(chain_files)
        
        # 创建一个ChainItem，包含所有链文件
        chain_item = ChainItem(
            file_paths=chain_files,
            run_number=1,  # 合并所有run，使用1作为统一编号
            chain_count=chain_count,
            tool="mrbayes"
        )
        
        # 发送信号
        self.export_chain_result_signal.emit(chain_item)
        self.add_console_message(f"Detected MCMC chain files: {run_count} run(s), {chain_count} chain file(s) total", "info")
    
    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        self.add_console_message(f"Analysis error: {error_message}", "error")
        QMessageBox.critical(self, "Error", f"Analysis failed: {error_message}")
    
    def stop_analysis(self):
        """停止分析"""
        if self.analysis_thread and self.analysis_thread.isRunning():
            # 停止线程
            self.analysis_thread.stop()
            self.add_console_message("Analysis stopped by user", "info")
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "Analysis has been stopped.")
    
    def get_mb_data_block(self):
        """生成MrBayes数据块内容"""
        try:
            params_list = self.get_parameters()  # 获取参数列表
            # 第一个参数是MPI并行数，第二个元素才是参数字典
            if len(params_list) > 1 and isinstance(params_list[1], dict):
                params = params_list[1]
            elif len(params_list) > 0 and isinstance(params_list[0], dict):
                params = params_list[0]
            else:
                # 如果参数获取失败，使用默认值
                self.add_console_message("Using default parameters for data block", "warning")
                params = self.get_default_parameters()
        except (IndexError, TypeError) as e:
            # 如果参数获取失败，使用默认值
            self.add_console_message(f"Error getting parameters, using defaults: {e}", "warning")
            params = self.get_default_parameters()
        
        commands = []
        commands.append("begin mrbayes;")
        
        # 检查是否启用分区模式
        if self.use_partition_mode and self.partition_definitions:
            # 分区模式
            partition_commands = MrBayesModelConverter.generate_partition_commands(
                self.partition_definitions, 
                self.partition_mode
            )
            commands.extend(partition_commands)
        else:
            # 单一模型模式
            # 模型设置
            if params.get('datatype', 'DNA') == 'DNA':
                nst = params.get('nst', 6)
                statefreq = params.get('statefreq', 'estimated(dirichlet)')
                commands.append(f"lset nucmodel=4by4 nst={nst};")
                
                if statefreq == 'fixed(equal)':
                    commands.append("prset statefreqpr=fixed(equal);")
                elif statefreq == 'fixed(empirical)':
                    commands.append("prset statefreqpr=fixed(empirical);")
                else:
                    commands.append("prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);")
            else:
                aamodel = params.get('aamodel', 'mixed')
                statefreq = params.get('statefreq', 'estimated(dirichlet)')
                commands.append("lset nucmodel=protein;")
                
                if aamodel == 'mixed':
                    commands.append("prset aamodelpr=mixed;")
                else:
                    commands.append(f"prset aamodelpr=fixed({aamodel});")
                
                if statefreq == 'fixed(equal)':
                    commands.append("prset statefreqpr=fixed(equal);")
                elif statefreq == 'fixed(empirical)':
                    commands.append("prset statefreqpr=fixed(empirical);")
                else:
                    commands.append("prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0);")
            
            # 速率异质性
            rates = params.get('rates', 'equal')
            ngammacat = params.get('ngammacat', 4)
            
            if rates == 'equal':
                commands.append("lset rates=equal;")
            elif rates == 'gamma':
                commands.append(f"lset rates=gamma ngammacat={ngammacat};")
            elif rates == 'invgamma':
                commands.append(f"lset rates=invgamma ngammacat={ngammacat};")
            elif rates == 'propinv':
                commands.append("lset rates=propinv;")
            elif rates == 'lnorm':
                commands.append(f"lset rates=lnorm nlnormcat={ngammacat};")
            elif rates == 'adgamma':
                commands.append(f"lset rates=adgamma ngammacat={ngammacat};")
        
        # BEAGLE设置
        if params.get('use_beagle', True):
            commands.append(f"set usebeagle=yes beagledevice={params.get('beagle_device', 'cpu')} beagleprecision={params.get('beagle_precision', 'double')} beaglescaling={params.get('beagle_scaling', 'dynamic')};")
        
        # MCMC设置
        commands.append(f"mcmcp ngen={params.get('ngen', 1000000)} samplefreq={params.get('samplefreq', 1000)} nchains={params.get('nchains', 4)} nruns={params.get('nruns', 2)} printfreq={params.get('samplefreq', 1000)} savebrlens=yes checkpoint=yes checkfreq=5000;")
        
        # 运行MCMC
        commands.append("mcmc;")
        
        # 总结结果
        if params.get('burnin_as_fraction', True):
            burnin_frac = params.get('burnin_frac', 0.25)
            commands.append(f"sump relburnin=yes burninfrac={burnin_frac};")
            commands.append(f"sumt relburnin=yes burninfrac={burnin_frac} contype={params.get('contype', 'majorityrule')} conformat={params.get('conformat', 'figtree')};")
        else:
            burnin_states = params.get('burnin_states', 1000)
            commands.append(f"sump burnin={burnin_states};")
            commands.append(f"sumt burnin={burnin_states} contype={params.get('contype', 'majorityrule')} conformat={params.get('conformat', 'figtree')};")
        
        commands.append("end;")
        
        return "\n".join(commands)
    
    def get_default_parameters(self):
        """获取默认参数字典"""
        return {
            'datatype': 'DNA',
            'nst': 6,
            'statefreq': 'estimated(dirichlet)',
            'aamodel': 'mixed',
            'rates': 'gamma',
            'ngammacat': 4,
            'use_beagle': True,
            'beagle_device': 'cpu',
            'beagle_precision': 'double',
            'beagle_scaling': 'dynamic',
            'ngen': 1000000,
            'samplefreq': 1000,
            'nchains': 4,
            'nruns': 2,
            'burnin_as_fraction': True,
            'burnin_frac': 0.25,
            'burnin_states': 1000,
            'contype': 'majorityrule',
            'conformat': 'figtree'
        }

    def show_mb_data_block(self):
        """显示MrBayes数据块"""
        dialog = QDialog(self)
        dialog.setWindowTitle("MrBayes Data Block")
        dialog.setMinimumSize(600, 400)
        
        layout = QVBoxLayout()
        
        text_edit = QTextEdit()
        text_edit.setFont(QFont("Consolas", 10))
        text_edit.setStyleSheet("""
            QTextEdit {
                background-color: #272822;
                color: #f8f8f2;
                border: 1px solid #3c3c3c;
            }
        """)
        
        try:
            data_block = self.get_mb_data_block()
            text_edit.setPlainText(data_block)
            
            # 添加语法高亮
            highlighter = MrBayesHighlighter(text_edit.document())
            
        except Exception as e:
            error_msg = f"Error generating MrBayes data block: {str(e)}"
            self.add_console_message(error_msg, "error")
            text_edit.setPlainText(error_msg)
        
        layout.addWidget(text_edit)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)
        
        dialog.setLayout(layout)
        dialog.exec_()

    def copy_mb_data_block(self):
        """复制MrBayes数据块到剪贴板"""
        try:
            data_block = self.get_mb_data_block()
            QApplication.clipboard().setText(data_block)
            self.add_console_message("MrBayes data block copied to clipboard", "info")
            QMessageBox.information(self, "Success", "MrBayes data block copied to clipboard!")
        except Exception as e:
            error_msg = f"Error copying MrBayes data block: {str(e)}"
            self.add_console_message(error_msg, "error")
            QMessageBox.critical(self, "Error", error_msg)
    
    def display_results(self, output_files):
        """显示分析结果，使用IcyTree插件显示系统发育树"""
        if not output_files:
            self.output_info.setText("No output files generated")
            return

        # 查找.con.tre文件（共识树文件）
        con_tre_files = [f for f in output_files if f.endswith('.con.tre')]
        
        if con_tre_files:
            tree_file = con_tre_files[0]
            try:
                # 读取树文件内容
                with open(tree_file, 'r') as f:
                    tree_content = f.read().strip()
                
                # 确保树内容不为空
                if not tree_content:
                    self.output_info.setText("Tree file is empty")
                    return
                
                # 导入IcyTree插件
                from ..icytree import IcyTreePlugin
                import os
                
                # 创建IcyTree插件实例
                plugin_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '')
                icytree_plugin = IcyTreePlugin(plugin_path=plugin_path)
                
                # 设置Newick字符串并显示
                icytree_plugin.set_newick_string(tree_content)
                
                # 在输出标签页中显示IcyTree
                output_layout = self.output_tab.layout()
                if output_layout:
                    # 隐藏web_view（如果存在）
                    if hasattr(self, 'web_view'):
                        self.web_view.setVisible(False)
                    
                    # 查找并移除之前的IcyTree插件（如果有）
                    for i in reversed(range(output_layout.count())):
                        widget = output_layout.itemAt(i).widget()
                        if widget and widget != self.output_info and widget != self.report_combo.parentWidget():
                            widget.setParent(None)
                    
                    # 添加IcyTree插件到输出标签页
                    output_layout.addWidget(icytree_plugin)
                
                self.output_info.setText(f"Consensus tree visualization ready: {os.path.basename(tree_file)}")
                
            except ImportError:
                # 如果无法导入IcyTree插件，显示错误信息
                self.output_info.setText("Error: IcyTree plugin not available")
                
            except Exception as e:
                error_msg = f"Error processing tree file: {str(e)}"
                self.output_info.setText(error_msg)
                self.add_console_message(error_msg, "error")
        else:
            # 没有找到共识树文件，显示信息
            self.output_info.setText(f"No consensus tree (.con.tre) found. Generated {len(output_files)} file(s).")

class MolecularDatingSettingsDialog(QDialog):
    """分子钟定年设置对话框"""
    
    def __init__(self, parent=None, clock_model="strict", clock_params=None, calibrations=None, outgroup="", taxa_list=None):
        super().__init__(parent)
        self.setWindowTitle("Molecular Dating Settings")
        self.setMinimumSize(700, 600)
        
        self.clock_model = clock_model
        self.clock_params = clock_params or {}
        self.calibrations = calibrations or []
        self.outgroup = outgroup if isinstance(outgroup, list) else [outgroup] if outgroup else []
        self.taxa_list = taxa_list or []
        
        self.init_ui()
        self.load_data()
    
    def init_ui(self):
        """初始化UI"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # 分子钟模型选择
        clock_group = QGroupBox("Clock Model")
        clock_layout = QFormLayout()
        clock_group.setLayout(clock_layout)
        
        self.clock_model_combo = QComboBox()
        self.clock_model_combo.addItems([
            "Strict Clock",
            "Independent Gamma Rates (IGR)",
            "Compound Poisson Process (CPP)",
            "Thorne-Kishino 2002 (TK02)"
        ])
        self.clock_model_combo.currentIndexChanged.connect(self.on_clock_model_changed)
        clock_layout.addRow("Clock Model:", self.clock_model_combo)
        
        # IGR参数 - 默认隐藏
        self.igr_widget = QWidget()
        self.igr_layout = QHBoxLayout()
        self.igr_layout.setContentsMargins(0, 0, 0, 0)
        self.igr_widget.setLayout(self.igr_layout)
        
        self.igrvarpr_combo = QComboBox()
        self.igrvarpr_combo.addItems(["Fixed", "Exponential", "Uniform"])
        self.igrvarpr_combo.setMinimumWidth(100)
        self.igr_layout.addWidget(QLabel("Distribution:"))
        self.igr_layout.addWidget(self.igrvarpr_combo)
        
        self.igrvarpr_mean_spinbox = QDoubleSpinBox()
        self.igrvarpr_mean_spinbox.setRange(0.001, 1000.0)
        self.igrvarpr_mean_spinbox.setValue(10.00)
        self.igrvarpr_mean_spinbox.setSingleStep(0.1)
        self.igrvarpr_mean_spinbox.setDecimals(2)
        self.igr_layout.addWidget(QLabel("Mean:"))
        self.igr_layout.addWidget(self.igrvarpr_mean_spinbox)
        
        self.igr_layout.addStretch()
        clock_layout.addRow("", self.igr_widget)
        
        # TK02参数 - 默认隐藏
        self.tk02_widget = QWidget()
        self.tk02_layout = QHBoxLayout()
        self.tk02_layout.setContentsMargins(0, 0, 0, 0)
        self.tk02_widget.setLayout(self.tk02_layout)
        
        self.tk02varpr_combo = QComboBox()
        self.tk02varpr_combo.addItems(["Fixed", "Exponential", "Uniform"])
        self.tk02varpr_combo.setMinimumWidth(100)
        self.tk02_layout.addWidget(QLabel("Distribution:"))
        self.tk02_layout.addWidget(self.tk02varpr_combo)
        
        self.tk02varpr_mean_spinbox = QDoubleSpinBox()
        self.tk02varpr_mean_spinbox.setRange(0.001, 1000.0)
        self.tk02varpr_mean_spinbox.setValue(1.00)
        self.tk02varpr_mean_spinbox.setSingleStep(0.1)
        self.tk02varpr_mean_spinbox.setDecimals(2)
        self.tk02_layout.addWidget(QLabel("Mean:"))
        self.tk02_layout.addWidget(self.tk02varpr_mean_spinbox)
        
        self.tk02_layout.addStretch()
        clock_layout.addRow("", self.tk02_widget)
        
        # CPP参数 - 默认隐藏
        self.cpp_widget = QWidget()
        self.cpp_layout = QVBoxLayout()
        self.cpp_layout.setContentsMargins(0, 0, 0, 0)
        self.cpp_widget.setLayout(self.cpp_layout)
        
        cpp_row1 = QHBoxLayout()
        cpp_row1.addWidget(QLabel("Rate Distribution:"))
        self.cppratepr_combo = QComboBox()
        self.cppratepr_combo.addItems(["Fixed", "Exponential"])
        self.cppratepr_combo.setMinimumWidth(100)
        cpp_row1.addWidget(self.cppratepr_combo)
        cpp_row1.addWidget(QLabel("Rate Mean:"))
        self.cppratepr_mean_spinbox = QDoubleSpinBox()
        self.cppratepr_mean_spinbox.setRange(0.001, 1000.0)
        self.cppratepr_mean_spinbox.setValue(0.10)
        self.cppratepr_mean_spinbox.setSingleStep(0.01)
        self.cppratepr_mean_spinbox.setDecimals(2)
        cpp_row1.addWidget(self.cppratepr_mean_spinbox)
        cpp_row1.addStretch()
        self.cpp_layout.addLayout(cpp_row1)
        
        cpp_row2 = QHBoxLayout()
        cpp_row2.addWidget(QLabel("Multiplier Value:"))
        self.cppmultdevpr_value_spinbox = QDoubleSpinBox()
        self.cppmultdevpr_value_spinbox.setRange(0.001, 1000.0)
        self.cppmultdevpr_value_spinbox.setValue(0.40)
        self.cppmultdevpr_value_spinbox.setSingleStep(0.01)
        self.cppmultdevpr_value_spinbox.setDecimals(2)
        cpp_row2.addWidget(self.cppmultdevpr_value_spinbox)
        cpp_row2.addStretch()
        self.cpp_layout.addLayout(cpp_row2)
        
        clock_layout.addRow("", self.cpp_widget)
        
        layout.addWidget(clock_group)
        
        # 外群设置 - 多选列表
        outgroup_group = QGroupBox("Outgroup Settings")
        outgroup_layout = QVBoxLayout()
        outgroup_group.setLayout(outgroup_layout)
        
        outgroup_info_layout = QHBoxLayout()
        outgroup_info_layout.addWidget(QLabel("Outgroup(s):"))
        self.outgroup_count_label = QLabel("(0 selected)")
        self.outgroup_count_label.setStyleSheet("color: gray;")
        outgroup_info_layout.addWidget(self.outgroup_count_label)
        outgroup_info_layout.addStretch()
        outgroup_layout.addLayout(outgroup_info_layout)
        
        # 创建带复选框的列表
        self.outgroup_list_widget = QTableWidget()
        self.outgroup_list_widget.setColumnCount(2)
        self.outgroup_list_widget.setHorizontalHeaderLabels(["", "Taxon"])
        self.outgroup_list_widget.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.outgroup_list_widget.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.outgroup_list_widget.verticalHeader().setVisible(False)
        self.outgroup_list_widget.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.outgroup_list_widget.setMaximumHeight(200)
        outgroup_layout.addWidget(self.outgroup_list_widget)
        
        # 初始化外群列表
        self.init_outgroup_list()
        
        layout.addWidget(outgroup_group)
        
        # 校准节点设置 - 表格形式
        calibration_group = QGroupBox("Calibration Nodes")
        calibration_layout = QVBoxLayout()
        calibration_group.setLayout(calibration_layout)
        
        # 按钮行
        calib_btn_layout = QHBoxLayout()
        self.add_calibration_btn = QPushButton("+ Add")
        self.add_calibration_btn.clicked.connect(self.add_calibration)
        self.edit_calibration_btn = QPushButton("Edit")
        self.edit_calibration_btn.clicked.connect(self.edit_calibration)
        self.remove_calibration_btn = QPushButton("- Remove")
        self.remove_calibration_btn.clicked.connect(self.remove_calibration)
        calib_btn_layout.addWidget(self.add_calibration_btn)
        calib_btn_layout.addWidget(self.edit_calibration_btn)
        calib_btn_layout.addWidget(self.remove_calibration_btn)
        calib_btn_layout.addStretch()
        calibration_layout.addLayout(calib_btn_layout)
        
        # Calibration表格
        self.calibration_table = QTableWidget()
        self.calibration_table.setColumnCount(5)
        self.calibration_table.setHorizontalHeaderLabels([
            "Constraint", "Node Name", "Taxa", "Prior Type", "Parameters"
        ])
        self.calibration_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.calibration_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        calibration_layout.addWidget(self.calibration_table)
        
        self.update_calibration_table()
        
        layout.addWidget(calibration_group)
        
        # 按钮
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        button_layout.addWidget(ok_button)
        
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(cancel_button)
        
        layout.addLayout(button_layout)
        
        # 初始状态
        self.on_clock_model_changed()
    
    def init_outgroup_list(self):
        """初始化外群列表"""
        self.outgroup_list_widget.setRowCount(len(self.taxa_list))
        
        for i, taxon in enumerate(self.taxa_list):
            # 复选框
            checkbox_item = QTableWidgetItem()
            checkbox_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            if taxon in self.outgroup:
                checkbox_item.setCheckState(Qt.Checked)
            else:
                checkbox_item.setCheckState(Qt.Unchecked)
            self.outgroup_list_widget.setItem(i, 0, checkbox_item)
            
            # 分类名称
            name_item = QTableWidgetItem(taxon)
            self.outgroup_list_widget.setItem(i, 1, name_item)
        
        # 连接信号
        self.outgroup_list_widget.itemChanged.connect(self.on_outgroup_changed)
        
        self.update_outgroup_count()
    
    def on_outgroup_changed(self, item):
        """外群复选框改变时处理"""
        if item.column() == 0:
            self.update_outgroup_count()
    
    def update_outgroup_count(self):
        """更新外群计数"""
        count = 0
        for i in range(self.outgroup_list_widget.rowCount()):
            item = self.outgroup_list_widget.item(i, 0)
            if item and item.checkState() == Qt.Checked:
                count += 1
        
        self.outgroup_count_label.setText(f"({count} selected)")
    
    def on_clock_model_changed(self):
        """分子钟模型改变时更新UI"""
        model = self.clock_model_combo.currentText()
        
        # 隐藏所有参数
        self.igr_widget.setVisible(False)
        self.tk02_widget.setVisible(False)
        self.cpp_widget.setVisible(False)
        
        # 根据选择的模型显示相应的参数
        if "IGR" in model:
            self.igr_widget.setVisible(True)
        elif "TK02" in model:
            self.tk02_widget.setVisible(True)
        elif "CPP" in model:
            self.cpp_widget.setVisible(True)
        # Strict Clock不需要额外参数
    
    def load_data(self):
        """加载数据"""
        # 设置分子钟模型
        model_map = {
            "strict": "Strict Clock",
            "igr": "Independent Gamma Rates (IGR)",
            "cpp": "Compound Poisson Process (CPP)",
            "tk02": "Thorne-Kishino 2002 (TK02)"
        }
        self.clock_model_combo.setCurrentText(model_map.get(self.clock_model, "Strict Clock"))
        
        # 加载参数
        if self.clock_model == "igr":
            self.igrvarpr_combo.setCurrentText(self.clock_params.get('igrvarpr', 'Exponential'))
            self.igrvarpr_mean_spinbox.setValue(self.clock_params.get('igrvarpr_mean', 10.00))
        elif self.clock_model == "tk02":
            self.tk02varpr_combo.setCurrentText(self.clock_params.get('tk02varpr', 'Exponential'))
            self.tk02varpr_mean_spinbox.setValue(self.clock_params.get('tk02varpr_mean', 1.00))
        elif self.clock_model == "cpp":
            self.cppratepr_combo.setCurrentText(self.clock_params.get('cppratepr', 'Exponential'))
            self.cppratepr_mean_spinbox.setValue(self.clock_params.get('cppratepr_mean', 0.10))
            self.cppmultdevpr_value_spinbox.setValue(self.clock_params.get('cppmultdevpr_value', 0.40))
        
        # 更新校准节点表格
        self.update_calibration_table()
    
    def update_calibration_table(self):
        """更新校准节点表格"""
        self.calibration_table.setRowCount(len(self.calibrations))
        
        for row, calib in enumerate(self.calibrations):
            # 复选框（constraint）
            checkbox_item = QTableWidgetItem()
            checkbox_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            checkbox_item.setCheckState(Qt.Checked if calib.get('use_constraint', False) else Qt.Unchecked)
            self.calibration_table.setItem(row, 0, checkbox_item)
            
            # 节点名称
            name_item = QTableWidgetItem(calib.get('name', ''))
            self.calibration_table.setItem(row, 1, name_item)
            
            # 分类群
            taxa_str = ", ".join(calib.get('taxa', []))
            taxa_item = QTableWidgetItem(taxa_str)
            self.calibration_table.setItem(row, 2, taxa_item)
            
            # 先验类型
            prior_type_item = QTableWidgetItem(calib.get('prior_type', ''))
            self.calibration_table.setItem(row, 3, prior_type_item)
            
            # 参数
            params = calib.get('params', {})
            param_str = self._format_params(calib.get('prior_type', ''), params)
            param_item = QTableWidgetItem(param_str)
            self.calibration_table.setItem(row, 4, param_item)
    
    def _format_params(self, prior_type, params):
        """格式化参数显示"""
        if prior_type == 'fixed':
            return f"fixed({params.get('fixed_value', 0)})"
        elif prior_type == 'uniform':
            return f"unif({params.get('uniform_min', 0)}, {params.get('uniform_max', 0)})"
        elif prior_type == 'offset_exp':
            return f"offsetexp({params.get('offset_exp_offset', 0)}, {params.get('offset_exp_rate', 1.0)})"
        elif prior_type == 'offset_gamma':
            return f"offsetgamma({params.get('offset_gamma_offset', 0)}, {params.get('offset_gamma_alpha', 1.0)}, {params.get('offset_gamma_beta', 1.0)})"
        elif prior_type == 'offset_lognorm':
            return f"offsetlognormal({params.get('offset_lognorm_offset', 0)}, {params.get('offset_lognorm_mean', 0.0)}, {params.get('offset_lognorm_std', 1.0)})"
        elif prior_type == 'truncated_normal':
            return f"truncatednormal({params.get('trunc_norm_offset', 0)}, {params.get('trunc_norm_mean', 0.0)}, {params.get('trunc_norm_std', 1.0)})"
        return ""
    
    def add_calibration(self):
        """添加校准节点"""
        dialog = CalibrationEditDialog(self, self.taxa_list)
        if dialog.exec_() == QDialog.Accepted:
            calib = dialog.get_calibration()
            if calib:
                self.calibrations.append(calib)
                self.update_calibration_table()
    
    def edit_calibration(self):
        """编辑校准节点"""
        current_row = self.calibration_table.currentRow()
        if current_row < 0:
            QMessageBox.warning(self, "Warning", "Please select a calibration node to edit!")
            return
        
        calib = self.calibrations[current_row]
        dialog = CalibrationEditDialog(self, self.taxa_list, calib)
        if dialog.exec_() == QDialog.Accepted:
            new_calib = dialog.get_calibration()
            if new_calib:
                self.calibrations[current_row] = new_calib
                self.update_calibration_table()
    
    def remove_calibration(self):
        """移除校准节点"""
        current_row = self.calibration_table.currentRow()
        if current_row < 0:
            QMessageBox.warning(self, "Warning", "Please select a calibration node to remove!")
            return
        
        calib = self.calibrations[current_row]
        reply = QMessageBox.question(
            self, "Confirm Removal",
            f"Are you sure you want to remove calibration node '{calib.get('name', '')}'?",
            QMessageBox.Yes | QMessageBox.No
        )
        
        if reply == QMessageBox.Yes:
            del self.calibrations[current_row]
            self.update_calibration_table()
    
    def get_clock_model(self) -> str:
        """获取分子钟模型"""
        model_text = self.clock_model_combo.currentText()
        model_map = {
            "Strict Clock": "strict",
            "Independent Gamma Rates (IGR)": "igr",
            "Compound Poisson Process (CPP)": "cpp",
            "Thorne-Kishino 2002 (TK02)": "tk02"
        }
        return model_map.get(model_text, "strict")
    
    def get_clock_params(self) -> dict:
        """获取分子钟参数"""
        model = self.get_clock_model()
        params = {}
        
        if model == "igr":
            params = {
                'igrvarpr': self.igrvarpr_combo.currentText().lower(),
                'igrvarpr_mean': self.igrvarpr_mean_spinbox.value()
            }
        elif model == "tk02":
            params = {
                'tk02varpr': self.tk02varpr_combo.currentText().lower(),
                'tk02varpr_mean': self.tk02varpr_mean_spinbox.value()
            }
        elif model == "cpp":
            params = {
                'cppratepr': self.cppratepr_combo.currentText().lower(),
                'cppratepr_mean': self.cppratepr_mean_spinbox.value(),
                'cppmultdevpr_value': self.cppmultdevpr_value_spinbox.value()
            }
        
        return params
    
    def get_calibrations(self) -> list:
        """获取校准节点列表"""
        result = []
        for i in range(self.calibration_table.rowCount()):
            checkbox_item = self.calibration_table.item(i, 0)
            use_constraint = checkbox_item.checkState() == Qt.Checked
            
            calib = self.calibrations[i].copy()
            calib['use_constraint'] = use_constraint
            result.append(calib)
        
        return result
    
    def get_outgroup(self) -> str:
        """获取外群"""
        selected_taxa = []
        for i in range(self.outgroup_list_widget.rowCount()):
            item = self.outgroup_list_widget.item(i, 0)
            if item and item.checkState() == Qt.Checked:
                taxon_item = self.outgroup_list_widget.item(i, 1)
                if taxon_item:
                    selected_taxa.append(taxon_item.text())
        
        return ", ".join(selected_taxa)


class CalibrationEditDialog(QDialog):
    """校准节点编辑对话框"""
    
    def __init__(self, parent=None, taxa_list=None, calibration=None):
        super().__init__(parent)
        self.setWindowTitle("Edit Calibration Node")
        self.setMinimumSize(500, 450)
        
        self.taxa_list = taxa_list or []
        self.calibration = calibration or {}
        
        self.init_ui()
        self.load_data()
    
    def init_ui(self):
        """初始化UI"""
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        form_layout = QFormLayout()
        
        # 节点名称
        self.name_edit = QLineEdit()
        form_layout.addRow("Node Name:", self.name_edit)
        
        # 分类群选择
        self.taxa_list_widget = QListWidget()
        self.taxa_list_widget.setSelectionMode(QListWidget.MultiSelection)
        self.taxa_list_widget.setMaximumHeight(150)
        if self.taxa_list:
            self.taxa_list_widget.addItems(self.taxa_list)
        form_layout.addRow("Taxa:", self.taxa_list_widget)
        
        # 校准先验类型
        self.prior_type_combo = QComboBox()
        self.prior_type_combo.addItems([
            "Fixed Value",
            "Uniform Distribution",
            "Offset Exponential",
            "Offset Gamma",
            "Offset Lognormal",
            "Truncated Normal"
        ])
        self.prior_type_combo.currentIndexChanged.connect(self.on_prior_type_changed)
        form_layout.addRow("Prior Type:", self.prior_type_combo)
        
        # 参数设置区域
        self.params_widget = QWidget()
        self.params_layout = QFormLayout()
        self.params_widget.setLayout(self.params_layout)
        
        # Fixed Value参数
        self.fixed_value_spinbox = QDoubleSpinBox()
        self.fixed_value_spinbox.setRange(0, 1000000)
        self.fixed_value_spinbox.setDecimals(6)
        fixed_form_layout = QFormLayout()
        fixed_form_layout.addRow("Value:", self.fixed_value_spinbox)
        self.fixed_widget = QWidget()
        self.fixed_widget.setLayout(fixed_form_layout)
        
        # Uniform参数
        self.uniform_min_spinbox = QDoubleSpinBox()
        self.uniform_min_spinbox.setRange(0, 1000000)
        self.uniform_max_spinbox = QDoubleSpinBox()
        self.uniform_max_spinbox.setRange(0, 1000000)
        uniform_form_layout = QFormLayout()
        uniform_form_layout.addRow("Min:", self.uniform_min_spinbox)
        uniform_form_layout.addRow("Max:", self.uniform_max_spinbox)
        self.uniform_widget = QWidget()
        self.uniform_widget.setLayout(uniform_form_layout)
        
        # Offset Exponential参数
        self.offset_exp_offset_spinbox = QDoubleSpinBox()
        self.offset_exp_offset_spinbox.setRange(0, 1000000)
        self.offset_exp_rate_spinbox = QDoubleSpinBox()
        self.offset_exp_rate_spinbox.setRange(0.0001, 1000)
        self.offset_exp_rate_spinbox.setValue(1.0)
        offset_exp_form_layout = QFormLayout()
        offset_exp_form_layout.addRow("Offset:", self.offset_exp_offset_spinbox)
        offset_exp_form_layout.addRow("Rate (λ):", self.offset_exp_rate_spinbox)
        self.offset_exp_widget = QWidget()
        self.offset_exp_widget.setLayout(offset_exp_form_layout)
        
        # Offset Gamma参数
        self.offset_gamma_offset_spinbox = QDoubleSpinBox()
        self.offset_gamma_offset_spinbox.setRange(0, 1000000)
        self.offset_gamma_alpha_spinbox = QDoubleSpinBox()
        self.offset_gamma_alpha_spinbox.setRange(0.01, 1000)
        self.offset_gamma_alpha_spinbox.setValue(1.0)
        self.offset_gamma_beta_spinbox = QDoubleSpinBox()
        self.offset_gamma_beta_spinbox.setRange(0.0001, 1000)
        self.offset_gamma_beta_spinbox.setValue(1.0)
        offset_gamma_form_layout = QFormLayout()
        offset_gamma_form_layout.addRow("Offset:", self.offset_gamma_offset_spinbox)
        offset_gamma_form_layout.addRow("Alpha (α):", self.offset_gamma_alpha_spinbox)
        offset_gamma_form_layout.addRow("Beta (β):", self.offset_gamma_beta_spinbox)
        self.offset_gamma_widget = QWidget()
        self.offset_gamma_widget.setLayout(offset_gamma_form_layout)
        
        # Offset Lognormal参数
        self.offset_lognorm_offset_spinbox = QDoubleSpinBox()
        self.offset_lognorm_offset_spinbox.setRange(0, 1000000)
        self.offset_lognorm_mean_spinbox = QDoubleSpinBox()
        self.offset_lognorm_mean_spinbox.setRange(-100, 100)
        self.offset_lognorm_mean_spinbox.setValue(0.0)
        self.offset_lognorm_std_spinbox = QDoubleSpinBox()
        self.offset_lognorm_std_spinbox.setRange(0.01, 100)
        self.offset_lognorm_std_spinbox.setValue(1.0)
        offset_lognorm_form_layout = QFormLayout()
        offset_lognorm_form_layout.addRow("Offset:", self.offset_lognorm_offset_spinbox)
        offset_lognorm_form_layout.addRow("Mean (μ):", self.offset_lognorm_mean_spinbox)
        offset_lognorm_form_layout.addRow("Std (σ):", self.offset_lognorm_std_spinbox)
        self.offset_lognorm_widget = QWidget()
        self.offset_lognorm_widget.setLayout(offset_lognorm_form_layout)
        
        # Truncated Normal参数
        self.trunc_norm_offset_spinbox = QDoubleSpinBox()
        self.trunc_norm_offset_spinbox.setRange(0, 1000000)
        self.trunc_norm_mean_spinbox = QDoubleSpinBox()
        self.trunc_norm_mean_spinbox.setRange(-100, 100)
        self.trunc_norm_mean_spinbox.setValue(0.0)
        self.trunc_norm_std_spinbox = QDoubleSpinBox()
        self.trunc_norm_std_spinbox.setRange(0.01, 100)
        self.trunc_norm_std_spinbox.setValue(1.0)
        trunc_norm_form_layout = QFormLayout()
        trunc_norm_form_layout.addRow("Offset:", self.trunc_norm_offset_spinbox)
        trunc_norm_form_layout.addRow("Mean (μ):", self.trunc_norm_mean_spinbox)
        trunc_norm_form_layout.addRow("Std (σ):", self.trunc_norm_std_spinbox)
        self.trunc_norm_widget = QWidget()
        self.trunc_norm_widget.setLayout(trunc_norm_form_layout)
        
        # 将所有参数widget添加到params_layout（而不是params_form_layout）
        self.params_layout.addRow("", self.fixed_widget)
        self.params_layout.addRow("", self.uniform_widget)
        self.params_layout.addRow("", self.offset_exp_widget)
        self.params_layout.addRow("", self.offset_gamma_widget)
        self.params_layout.addRow("", self.offset_lognorm_widget)
        self.params_layout.addRow("", self.trunc_norm_widget)
        
        form_layout.addRow("Parameters:", self.params_widget)
        
        layout.addLayout(form_layout)
        
        # 按钮
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(ok_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)
        
        # 初始状态
        self.on_prior_type_changed()
    
    def on_prior_type_changed(self):
        """先验类型改变时更新UI"""
        prior_type = self.prior_type_combo.currentText()
        
        # 隐藏所有参数widget
        self.fixed_widget.setVisible(False)
        self.uniform_widget.setVisible(False)
        self.offset_exp_widget.setVisible(False)
        self.offset_gamma_widget.setVisible(False)
        self.offset_lognorm_widget.setVisible(False)
        self.trunc_norm_widget.setVisible(False)
        
        # 根据选择的先验类型显示相应的参数
        if prior_type == "Fixed Value":
            self.fixed_widget.setVisible(True)
        elif prior_type == "Uniform Distribution":
            self.uniform_widget.setVisible(True)
        elif prior_type == "Offset Exponential":
            self.offset_exp_widget.setVisible(True)
        elif prior_type == "Offset Gamma":
            self.offset_gamma_widget.setVisible(True)
        elif prior_type == "Offset Lognormal":
            self.offset_lognorm_widget.setVisible(True)
        elif prior_type == "Truncated Normal":
            self.trunc_norm_widget.setVisible(True)
    
    def load_data(self):
        """加载数据"""
        self.name_edit.setText(self.calibration.get('name', ''))
        
        # 加载分类群选择
        selected_taxa = self.calibration.get('taxa', [])
        for i in range(self.taxa_list_widget.count()):
            item = self.taxa_list_widget.item(i)
            if item.text() in selected_taxa:
                item.setSelected(True)
        
        # 加载先验类型
        prior_type = self.calibration.get('prior_type', '')
        if prior_type == 'fixed':
            self.prior_type_combo.setCurrentText("Fixed Value")
        elif prior_type == 'uniform':
            self.prior_type_combo.setCurrentText("Uniform Distribution")
        elif prior_type == 'offset_exp':
            self.prior_type_combo.setCurrentText("Offset Exponential")
        elif prior_type == 'offset_gamma':
            self.prior_type_combo.setCurrentText("Offset Gamma")
        elif prior_type == 'offset_lognorm':
            self.prior_type_combo.setCurrentText("Offset Lognormal")
        elif prior_type == 'truncated_normal':
            self.prior_type_combo.setCurrentText("Truncated Normal")
        
        # 加载参数
        params = self.calibration.get('params', {})
        if 'fixed_value' in params:
            self.fixed_value_spinbox.setValue(params['fixed_value'])
        if 'uniform_min' in params:
            self.uniform_min_spinbox.setValue(params['uniform_min'])
        if 'uniform_max' in params:
            self.uniform_max_spinbox.setValue(params['uniform_max'])
        if 'offset_exp_offset' in params:
            self.offset_exp_offset_spinbox.setValue(params['offset_exp_offset'])
        if 'offset_exp_rate' in params:
            self.offset_exp_rate_spinbox.setValue(params['offset_exp_rate'])
        if 'offset_gamma_offset' in params:
            self.offset_gamma_offset_spinbox.setValue(params['offset_gamma_offset'])
        if 'offset_gamma_alpha' in params:
            self.offset_gamma_alpha_spinbox.setValue(params['offset_gamma_alpha'])
        if 'offset_gamma_beta' in params:
            self.offset_gamma_beta_spinbox.setValue(params['offset_gamma_beta'])
        if 'offset_lognorm_offset' in params:
            self.offset_lognorm_offset_spinbox.setValue(params['offset_lognorm_offset'])
        if 'offset_lognorm_mean' in params:
            self.offset_lognorm_mean_spinbox.setValue(params['offset_lognorm_mean'])
        if 'offset_lognorm_std' in params:
            self.offset_lognorm_std_spinbox.setValue(params['offset_lognorm_std'])
        if 'trunc_norm_offset' in params:
            self.trunc_norm_offset_spinbox.setValue(params['trunc_norm_offset'])
        if 'trunc_norm_mean' in params:
            self.trunc_norm_mean_spinbox.setValue(params['trunc_norm_mean'])
        if 'trunc_norm_std' in params:
            self.trunc_norm_std_spinbox.setValue(params['trunc_norm_std'])
    
    def get_calibration(self):
        """获取校准节点数据"""
        name = self.name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Warning", "Node name cannot be empty!")
            return None
        
        # 获取选中的分类群
        selected_items = self.taxa_list_widget.selectedItems()
        taxa = [item.text() for item in selected_items]
        
        if not taxa:
            QMessageBox.warning(self, "Warning", "Please select at least one taxon!")
            return None
        
        # 获取先验类型
        prior_type_map = {
            "Fixed Value": "fixed",
            "Uniform Distribution": "uniform",
            "Offset Exponential": "offset_exp",
            "Offset Gamma": "offset_gamma",
            "Offset Lognormal": "offset_lognorm",
            "Truncated Normal": "truncated_normal"
        }
        prior_type = prior_type_map.get(self.prior_type_combo.currentText(), "fixed")
        
        # 获取参数
        params = {}
        if prior_type == "fixed":
            params['fixed_value'] = self.fixed_value_spinbox.value()
        elif prior_type == "uniform":
            params['uniform_min'] = self.uniform_min_spinbox.value()
            params['uniform_max'] = self.uniform_max_spinbox.value()
        elif prior_type == "offset_exp":
            params['offset_exp_offset'] = self.offset_exp_offset_spinbox.value()
            params['offset_exp_rate'] = self.offset_exp_rate_spinbox.value()
        elif prior_type == "offset_gamma":
            params['offset_gamma_offset'] = self.offset_gamma_offset_spinbox.value()
            params['offset_gamma_alpha'] = self.offset_gamma_alpha_spinbox.value()
            params['offset_gamma_beta'] = self.offset_gamma_beta_spinbox.value()
        elif prior_type == "offset_lognorm":
            params['offset_lognorm_offset'] = self.offset_lognorm_offset_spinbox.value()
            params['offset_lognorm_mean'] = self.offset_lognorm_mean_spinbox.value()
            params['offset_lognorm_std'] = self.offset_lognorm_std_spinbox.value()
        elif prior_type == "truncated_normal":
            params['trunc_norm_offset'] = self.trunc_norm_offset_spinbox.value()
            params['trunc_norm_mean'] = self.trunc_norm_mean_spinbox.value()
            params['trunc_norm_std'] = self.trunc_norm_std_spinbox.value()
        
        return {
            'name': name,
            'taxa': taxa,
            'prior_type': prior_type,
            'params': params,
            'use_constraint': False  # 默认不设置拓扑约束
        }


class MrBayesHighlighter(QSyntaxHighlighter):
    """MrBayes语法高亮器，使用Monokai主题"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # 定义格式
        self.begin_end_format = QTextCharFormat()
        self.begin_end_format.setForeground(QColor("#66D9EF"))  # 青色
        self.begin_end_format.setFontWeight(QFont.Bold)
        
        self.command_format = QTextCharFormat()
        self.command_format.setForeground(QColor("#F92772"))     # 红橙色
        self.command_format.setFontWeight(QFont.Bold)
        
        self.setting_format = QTextCharFormat()
        self.setting_format.setForeground(QColor("#E6DB74"))     # 黄色
        
        self.value_format = QTextCharFormat()
        self.value_format.setForeground(QColor("#AE81FF"))       # 紫色
        
        # 编译正则表达式
        self.begin_end_regex = QRegExp(r"^\s*(begin mrbayes|end)\s*;?\s*$")
        self.command_regex = QRegExp(r"\b(lset|mcmcp|sump|sumt|prset|set|mcmc|charset|partition)\b")
        self.settings_regex = QRegExp(r"([a-zA-Z]+)=([^\s;]+)")
        self.settings_regex_spacing = QRegExp(r"([a-zA-Z]+) = ([^\s;]+)")
    
    def highlightBlock(self, text):
        # 高亮 begin/end 行
        pos = self.begin_end_regex.indexIn(text)
        if pos >= 0:
            self.setFormat(0, len(text), self.begin_end_format)
            return  # 整行已经高亮，不需要再处理
        
        # 高亮命令
        pos = 0
        index = self.command_regex.indexIn(text, pos)
        while index >= 0:
            length = self.command_regex.matchedLength()
            self.setFormat(index, length, self.command_format)
            pos = index + length
            index = self.command_regex.indexIn(text, pos)
        
        # 高亮设置项和值
        pos = 0
        index = self.settings_regex.indexIn(text, pos)
        while index >= 0:
            matched_str = self.settings_regex.cap(0)
            setting_part = self.settings_regex.cap(1)
            value_part = self.settings_regex.cap(2)
            
            # 高亮设置项（等号前的部分）
            setting_pos = text.find(setting_part, index)
            self.setFormat(setting_pos, len(setting_part), self.setting_format)
            
            # 高亮值（等号后的部分）
            value_pos = text.find(value_part, setting_pos + len(setting_part))
            self.setFormat(value_pos, len(value_part), self.value_format)
            
            pos = index + len(matched_str)
            index = self.settings_regex.indexIn(text, pos)

        index = self.settings_regex_spacing.indexIn(text, pos)
        while index >= 0:
            matched_str = self.settings_regex_spacing.cap(0)
            setting_part = self.settings_regex_spacing.cap(1)
            value_part = self.settings_regex_spacing.cap(2)
            
            # 高亮设置项（等号前的部分）
            setting_pos = text.find(setting_part, index)
            self.setFormat(setting_pos, len(setting_part), self.setting_format)
            
            # 高亮值（等号后的部分）
            value_pos = text.find(value_part, setting_pos + len(setting_part))
            self.setFormat(value_pos, len(value_part), self.value_format)
            
            pos = index + len(matched_str)
            index = self.settings_regex_spacing.indexIn(text, pos)

    def copy_mb_data_block(self):
        """复制MrBayes数据块到剪贴板"""
        try:
            data_block = self.get_mb_data_block()
            QApplication.clipboard().setText(data_block)
            self.add_console_message("MrBayes data block copied to clipboard", "info")
            QMessageBox.information(self, "Success", "MrBayes data block copied to clipboard!")
        except Exception as e:
            error_msg = f"Error copying MrBayes data block: {str(e)}"
            self.add_console_message(error_msg, "error")
            QMessageBox.critical(self, "Error", error_msg)
    
    def display_results(self, output_files):
        """显示分析结果，使用IcyTree插件显示系统发育树"""
        if not output_files:
            self.output_info.setText("No output files generated")
            return

        # 查找.con.tre文件（共识树文件）
        con_tre_files = [f for f in output_files if f.endswith('.con.tre')]
        
        if con_tre_files:
            tree_file = con_tre_files[0]
            try:
                # 读取树文件内容
                with open(tree_file, 'r') as f:
                    tree_content = f.read().strip()
                
                # 确保树内容不为空
                if not tree_content:
                    self.output_info.setText("Tree file is empty")
                    return
                
                # 导入IcyTree插件
                from ..icytree import IcyTreePlugin
                import os
                
                # 创建IcyTree插件实例
                plugin_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '')
                icytree_plugin = IcyTreePlugin(plugin_path=plugin_path)
                
                # 设置Newick字符串并显示
                icytree_plugin.set_newick_string(tree_content)
                
                # 在输出标签页中显示IcyTree
                output_layout = self.output_tab.layout()
                if output_layout:
                    # 清除现有部件
                    for i in reversed(range(output_layout.count())):
                        widget = output_layout.itemAt(i).widget()
                        if widget and widget != self.output_info:
                            widget.setParent(None)
                
                # 添加IcyTree插件到输出标签页
                output_layout.addWidget(icytree_plugin)
                
                self.output_info.setText(f"Consensus tree visualization ready: {os.path.basename(tree_file)}")
                
            except ImportError:
                # 如果无法导入IcyTree插件，回退到显示文件列表
                self.output_info.setText("IcyTree plugin not available, showing file list instead")
                
                html_content = f"""
                <!DOCTYPE html>
                <html>
                <head>
                    <title>MrBayes Result</title>
                    <style>
                        body {{ font-family: Arial, sans-serif; margin: 20px; }}
                        h1 {{ color: #2c3e50; }}
                        .info {{ background-color: #e8f4f8; padding: 10px; border-radius: 5px; }}
                    </style>
                </head>
                <body>
                    <h1>MrBayes Analysis Result</h1>
                    <div class="info">
                        <p><strong>Consensus Tree File:</strong> {os.path.basename(tree_file)}</p>
                        <p><strong>Total Output Files:</strong> {len(output_files)}</p>
                    </div>
                    <h2>Tree Content (Newick Format)</h2>
                    <pre>{tree_content}</pre>
                    <h2>All Output Files</h2>
                    <ul>
                """
                
                for f in output_files:
                    html_content += f'<li>{os.path.basename(f)}</li>'
                
                html_content += """
                    </ul>
                </body>
                </html>
                """
                
                html_file = self.create_temp_file(suffix='.html')
                with open(html_file, 'w', encoding='utf-8') as f:
                    f.write(html_content)
                
                self.reports = [html_file]
                self.update_report_combo()
                self.show_current_report()
                
            except Exception as e:
                error_msg = f"Error processing tree file: {str(e)}"
                self.output_info.setText(error_msg)
                self.add_console_message(error_msg, "error")
                
                # 创建错误信息的HTML页面
                html_content = f"""
                <!DOCTYPE html>
                <html>
                <head>
                    <title>Error Display</title>
                    <style>
                        body {{ 
                            font-family: Arial, sans-serif; 
                            margin: 20px; 
                            background-color: #ffe6e6; 
                            color: #d00;
                        }}
                    </style>
                </head>
                <body>
                    <h1>Error Processing Tree File</h1>
                    <p>{error_msg}</p>
                </body>
                </html>
                """
                
                html_file = self.create_temp_file(suffix='.html')
                with open(html_file, 'w', encoding='utf-8') as f:
                    f.write(html_content)
                
                self.reports = [html_file]
                self.update_report_combo()
                self.show_current_report()
        else:
            self.output_info.setText(f"No consensus tree found. Generated {len(output_files)} file(s).")
            
            # 显示所有文件
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>MrBayes Result</title>
                <style>
                    body {{ font-family: Arial, sans-serif; margin: 20px; }}
                    h1 {{ color: #2c3e50; }}
                    .info {{ background-color: #e8f4f8; padding: 10px; border-radius: 5px; }}
                </style>
            </head>
            <body>
                <h1>MrBayes Analysis Result</h1>
                <div class="info">
                    <p><strong>Total Output Files:</strong> {len(output_files)}</p>
                </div>
                <h2>Output Files</h2>
                <ul>
            """
            
            for f in output_files:
                html_content += f'<li>{os.path.basename(f)}</li>'
            
            html_content += """
                </ul>
            </body>
            </html>
            """
            
            html_file = self.create_temp_file(suffix='.html')
            with open(html_file, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            self.reports = [html_file]
            self.update_report_combo()
            self.show_current_report()

class MrBayesPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.plugin_path = plugin_path
        # self.config = config_loader()
    
    def run(self, import_from=None, import_data=None, workdir=None, imported_model=None, model_conversion_result=None, seq_type="DNA"):
        return MrBayesPlugin(import_from, import_data, workdir=workdir, imported_model=imported_model, model_conversion_result=model_conversion_result, seq_type=seq_type)
