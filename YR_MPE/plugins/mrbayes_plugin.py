from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, 
                             QPushButton, QLineEdit, QTextEdit, 
                             QLabel, QComboBox, QCheckBox, QRadioButton, 
                             QSpinBox, QDoubleSpinBox, QScrollArea, 
                             QFrame, QTextEdit, QToolButton, QDialog,
                             QGroupBox, QSizePolicy, QFormLayout, QGridLayout,
                             QFileDialog, QMessageBox, QApplication)
from PyQt5.QtCore import Qt, pyqtSignal, QRegExp
from PyQt5.QtGui import QFont, QTextCursor, QSyntaxHighlighter, QTextCharFormat, QColor, QPalette
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import os
import tempfile
import json
from Bio import SeqIO
import re

class MrBayesThread(BaseProcessThread):
    """MrBayes系统发育推断线程类"""
    
    def __init__(self, tool_path, mpirun_path, input_files, parameters, imported_files=None, run_data_block_checked=True):
        super().__init__(tool_path, input_files, parameters, imported_files)
        self.mpirun_path = mpirun_path
        self.run_data_block_checked = run_data_block_checked
    
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
                
                # 获取MPI并行数
                mpi_parallel = params.pop(0)  # 第一个参数是并行数
                
                # 构建MrBayes命令脚本
                mb_script = self._generate_mrbayes_script(input_file, params)
                
                # 创建临时脚本文件
                script_file = self.create_temp_file(suffix='.nexus')
                with open(script_file, 'w') as f:
                    f.write(mb_script)
                
                # 构建命令
                if mpi_parallel > 1:
                    # 使用MPI并行版本
                    cmd = [
                        self.mpirun_path,
                        "-np", str(mpi_parallel),
                        self.tool_path
                    ]
                else:
                    # 使用单线程版本
                    cmd = [self.tool_path]
                
                # 添加输入文件
                cmd.append(script_file)
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"MrBayes execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的输出文件
                input_dir = os.path.dirname(input_file)
                input_name = os.path.splitext(os.path.basename(input_file))[0]
                
                # 查找.con.tre文件
                con_tre_files = [f for f in os.listdir(input_dir) if f.startswith(input_name) and f.endswith('.con.tre')]
                for f in con_tre_files:
                    output_files.append(os.path.join(input_dir, f))
                
                # 查找.chain文件
                chain_files = [f for f in os.listdir(input_dir) if f.startswith(input_name) and '.run' in f and f.endswith('.p')]
                for f in chain_files:
                    output_files.append(os.path.join(input_dir, f))
                
                # 查找.t文件
                t_files = [f for f in os.listdir(input_dir) if f.startswith(input_name) and '.run' in f and f.endswith('.t')]
                for f in t_files:
                    output_files.append(os.path.join(input_dir, f))
                        
            self.progress.emit("Phylogenetic inference completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Phylogenetic inference exception: {str(e)}")
    
    def _generate_mrbayes_script(self, input_file, params):
        """生成MrBayes命令脚本"""
        script = []
        
        # 检查输入文件格式
        file_ext = os.path.splitext(input_file)[1].lower()
        if file_ext in ['.nex', '.nexus']:
            # 检查NEXUS文件是否包含MrBayes块
            has_mrbayes_block = self._check_nexus_for_mrbayes_block(input_file)
            
            if has_mrbayes_block and self.run_data_block_checked:
                # 如果文件包含MrBayes块且checkbox被选中，直接运行
                script.append(f"execute {input_file};")
                # 结束，因为我们不需要添加额外的MrBayes命令
                script.append("quit;")
                return "\n".join(script)
            elif has_mrbayes_block and not self.run_data_block_checked:
                # 如果文件包含MrBayes块但checkbox未被选中，移除原有的MrBayes块
                input_file = self._remove_mrbayes_block_from_nexus(input_file)
                script.append(f"execute {input_file};")
            else:
                # 文件不包含MrBayes块，使用我们的配置
                script.append(f"execute {input_file};")
        else:
            # 非NEXUS格式，需要转换
            nexus_file = self._convert_to_nexus(input_file)
            script.append(f"execute {nexus_file};")
        
        # 添加MrBayes命令块
        script.append("begin mrbayes;")
        
        # 设置模型参数
        script.extend(self._get_model_commands(params))
        
        # 设置MCMC参数
        script.extend(self._get_mcmc_commands(params))
        
        # 运行MCMC
        script.append("mcmc;")
        
        # 总结结果
        script.extend(self._get_summary_commands(params))
        
        # 结束MrBayes块
        script.append("end;")
        
        # 添加退出命令
        script.append("quit;")
        
        return "\n".join(script)
    
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
        """将输入文件转换为NEXUS格式"""
        # 确定输入文件格式
        file_ext = os.path.splitext(input_file)[1].lower()
        if file_ext in ['.fas', '.fasta', '.fa', '.fna']:
            input_format = 'fasta'
        elif file_ext in ['.phy', '.phylip']:
            input_format = 'phylip'
        elif file_ext in ['.nex', '.nexus']:
            # 如果已经是NEXUS格式，直接返回
            return input_file
        else:
            # 默认尝试fasta格式
            input_format = 'fasta'
        
        # 生成临时NEXUS文件
        output_file = self.create_temp_file(suffix='.nex')
        
        try:
            # 使用BioPython进行格式转换
            sequences = list(SeqIO.parse(input_file, input_format))
            SeqIO.write(sequences, output_file, 'nexus')
            return output_file
        except Exception as e:
            self.error.emit(f"Error converting file to NEXUS format: {str(e)}")
            # 如果转换失败，返回原始文件
            return input_file
    
    def _get_model_commands(self, params):
        """获取模型设置命令"""
        commands = []
        
        # 数据类型 - params可能是列表或字典，如果是列表则取第一个元素之后的字典
        if isinstance(params, list) and len(params) > 0:
            # 如果是列表，取第一个参数（MPI并行数）之后的参数字典
            if len(params) > 1 and isinstance(params[1], dict):
                param_dict = params[1]
            elif len(params) > 0 and isinstance(params[0], dict):
                param_dict = params[0]
            else:
                # 如果参数格式不对，使用默认值
                param_dict = {}
        else:
            param_dict = params if isinstance(params, dict) else {}

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
        if isinstance(params, list) and len(params) > 0:
            # 如果是列表，取第一个参数（MPI并行数）之后的参数字典
            if len(params) > 1 and isinstance(params[1], dict):
                param_dict = params[1]
            elif len(params) > 0 and isinstance(params[0], dict):
                param_dict = params[0]
            else:
                # 如果参数格式不对，使用默认值
                param_dict = {}
        else:
            param_dict = params if isinstance(params, dict) else {}

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
        if isinstance(params, list) and len(params) > 0:
            # 如果是列表，取第一个参数（MPI并行数）之后的参数字典
            if len(params) > 1 and isinstance(params[1], dict):
                param_dict = params[1]
            elif len(params) > 0 and isinstance(params[0], dict):
                param_dict = params[0]
            else:
                # 如果参数格式不对，使用默认值
                param_dict = {}
        else:
            param_dict = params if isinstance(params, dict) else {}

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
    def __init__(self, import_from=None, import_data=None):
        super().__init__(import_from, import_data)
        
        # 初始化变量
        if not hasattr(self, 'imported_files'):
            self.imported_files = []
        if not hasattr(self, 'file_tags'):
            self.file_tags = []
    
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

        self.rate_params_layout.addWidget(QLabel("Type:"))
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

        self.prodata_layout.addWidget(QLabel("Protein model"))
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


        
        # MPI & BEAGLE settings
        mpi_beagle_group = QGroupBox("MPI && BEAGLE settings")
        mpi_beagle_layout = QFormLayout()
        mpi_beagle_group.setLayout(mpi_beagle_layout)
        layout.addWidget(mpi_beagle_group)

        # MPI parallels
        mpi_params_layout = QHBoxLayout()
        self.use_mpi = QCheckBox("Use MPI")
        self.use_mpi.setChecked(True)
        # Number of Parallels
        self.n_mpi_parallel = QSpinBox()
        n_cpu_cores = os.cpu_count()
        # Set the maximum value to the number of CPU cores
        self.n_mpi_parallel.setMaximum(n_cpu_cores)
        # default value: max(n_cpu_cores / 4, 1)
        self.n_mpi_parallel.setValue(max(int(n_cpu_cores / 4), 1))

        mpi_params_layout.addWidget(self.use_mpi)
        mpi_params_layout.addWidget(QLabel("Number of Parallels:"))
        mpi_params_layout.addWidget(self.n_mpi_parallel)
        mpi_beagle_layout.addRow("MPI settings:", mpi_params_layout) # <mpirun [UNIX] / mpiexec [WINDOWS]> -np <number>

        self.beagle_params_layout = QGridLayout()
        # use beagle?

        # set usebeagle=yes/no beagledevice=cpu/gpu beagleprecision=double/single beaglescaling=dynamic/always;
        self.use_beagle = QCheckBox("Use BEAGLE") 
        self.use_beagle.setChecked(True)

        self.beagle_device_label = QLabel("Device:")
        self.beagle_device_label.setAlignment(Qt.AlignRight)
        self.beagle_device_combo = QComboBox()
        self.beagle_device_combo.addItems(["CPU", "GPU"])
        self.beagle_device_combo.setCurrentText("CPU")

        self.beagle_precision_label = QLabel("Precision:")
        self.beagle_precision_label.setAlignment(Qt.AlignRight)
        self.beagle_precision_combo = QComboBox()
        self.beagle_precision_combo.addItems(["Double", "Single"])
        self.beagle_precision_combo.setCurrentText("Double")

        self.beagle_scaling_label = QLabel("Scaling:")
        self.beagle_scaling_label.setAlignment(Qt.AlignRight)
        self.beagle_scaling_combo = QComboBox()
        self.beagle_scaling_combo.addItems(["Dynamic", "Always"])
        self.beagle_scaling_combo.setCurrentText("Dynamic")

        self.beagle_params_layout.addWidget(self.use_beagle, 0, 0)
        self.beagle_params_layout.addWidget(self.beagle_device_label, 0, 2)
        self.beagle_params_layout.addWidget(self.beagle_device_combo, 0, 3)
        self.beagle_params_layout.addWidget(self.beagle_precision_label, 1, 0)
        self.beagle_params_layout.addWidget(self.beagle_precision_combo, 1, 1)
        self.beagle_params_layout.addWidget(self.beagle_scaling_label, 1, 2)
        self.beagle_params_layout.addWidget(self.beagle_scaling_combo, 1, 3)

        mpi_beagle_layout.addRow("BEAGLE settings", self.beagle_params_layout)

        mcmc_params_group = QGroupBox("MCMC settings")
        mcmc_params_layout = QFormLayout()
        mcmc_params_group.setLayout(mcmc_params_layout)

        layout.addWidget(mcmc_params_group)

        # generation; sampling frequency; run num; chain num;
        # mcmcp ngen=* samplefreq=* printfreq=<samplefreq> nchains=* nruns=* savebrlens=yes checkpoint=yes checkfreq=5000;
        self.generation_spinbox = QSpinBox()
        self.generation_spinbox.setMinimum(1)
        self.generation_spinbox.setValue(1000000)
        self.generation_spinbox.setMaximum(2147483647)

        self.sampling_frequency_spinbox = QSpinBox()
        self.sampling_frequency_spinbox.setMinimum(1)
        self.sampling_frequency_spinbox.setValue(1000)
        self.sampling_frequency_spinbox.setMaximum(2147483647)

        self.run_num_spinbox = QSpinBox()
        self.run_num_spinbox.setMinimum(1)
        self.run_num_spinbox.setValue(2)

        self.chain_num_spinbox = QSpinBox()
        self.chain_num_spinbox.setMinimum(1)
        self.chain_num_spinbox.setValue(4) # MC^3

        mcmcp_widget = QWidget()
        mcmcp_layout = QHBoxLayout()
        mcmcp_widget.setContentsMargins(0, 0, 0, 0)
        mcmcp_layout.setContentsMargins(0, 0, 0, 0)
        mcmcp_widget.setLayout(mcmcp_layout)

        mcmcp_layout.addWidget(QLabel("Generations:"))
        mcmcp_layout.addWidget(self.generation_spinbox)
        mcmcp_layout.addWidget(QLabel("Sampling Freq.:"))
        mcmcp_layout.addWidget(self.sampling_frequency_spinbox)
        mcmcp_layout.addWidget(QLabel("Runs:"))
        mcmcp_layout.addWidget(self.run_num_spinbox)
        mcmcp_layout.addWidget(QLabel("Chains:"))
        mcmcp_layout.addWidget(self.chain_num_spinbox)

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

    # ============ MrBayesPlugin 业务逻辑方法 ============
    
    def on_datatype_changed(self, index):
        """数据类型改变时的处理"""
        if index == 0:  # DNA
            self.dnadata_widget.setVisible(True)
            self.prodata_widget.setVisible(False)
        else:  # Protein
            self.dnadata_widget.setVisible(False)
            self.prodata_widget.setVisible(True)
    
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
        
        # 添加MPI并行数作为第一个参数
        if self.use_mpi.isChecked():
            params.append(self.n_mpi_parallel.value())
        else:
            params.append(1)
        
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
            'use_beagle': self.use_beagle.isChecked(),
            'beagle_device': self.beagle_device_combo.currentText().lower(),
            'beagle_precision': self.beagle_precision_combo.currentText().lower(),
            'beagle_scaling': self.beagle_scaling_combo.currentText().lower()
        })
        
        # MCMC参数
        params[-1].update({
            'ngen': self.generation_spinbox.value(),
            'samplefreq': self.sampling_frequency_spinbox.value(),
            'nchains': self.chain_num_spinbox.value(),
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
        if self.use_mpi.isChecked():
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
        
        # 在单独的线程中运行MrBayes
        self.analysis_thread = MrBayesThread(
            self.tool_path, mpirun_path, input_files, self.get_parameters(), self.imported_files,
            self.run_data_block.isChecked()
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
        
        # 显示结果
        self.display_results(output_files)
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        # 添加控制台消息
        self.add_console_message(f"Phylogenetic inference completed successfully! Found {len(output_files)} result file(s)", "info")
        
        QMessageBox.information(self, "Completed", "Phylogenetic inference completed!")
    
    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        self.add_console_message(f"Analysis error: {error_message}", "error")
        QMessageBox.critical(self, "Error", f"Analysis failed: {error_message}")
    
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
        self.command_regex = QRegExp(r"\b(lset|mcmcp|sump|sumt|prset|set|mcmc)\b")
        self.settings_regex = QRegExp(r"([a-zA-Z]+)=([^\s;]+)")
    
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
        """显示分析结果"""
        if not output_files:
            self.output_info.setText("No output files generated")
            return
        
        # 查找.con.tre文件
        con_tre_files = [f for f in output_files if f.endswith('.con.tre')]
        
        if con_tre_files:
            # 读取并显示树文件
            try:
                with open(con_tre_files[0], 'r') as f:
                    tree_content = f.read()
                self.output_info.setText(f"Consensus tree: {os.path.basename(con_tre_files[0])}")
                
                # 生成简单的HTML预览
                html_content = f"""
                <!DOCTYPE html>
                <html>
                <head>
                    <title>MrBayes Result</title>
                    <style>
                        body {{ font-family: Arial, sans-serif; margin: 20px; }}
                        h1 {{ color: #2c3e50; }}
                        .info {{ background-color: #e8f4f8; padding: 10px; border-radius: 5px; }}
                        pre {{ background-color: #f8f9fa; padding: 10px; border-radius: 5px; overflow-x: auto; }}
                    </style>
                </head>
                <body>
                    <h1>MrBayes Analysis Result</h1>
                    <div class="info">
                        <p><strong>Consensus Tree File:</strong> {os.path.basename(con_tre_files[0])}</p>
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
                
                # 创建临时HTML文件
                html_file = self.create_temp_file(suffix='.html')
                with open(html_file, 'w', encoding='utf-8') as f:
                    f.write(html_content)
                
                self.reports = [html_file]
                self.update_report_combo()
                self.show_current_report()
                
            except Exception as e:
                self.output_info.setText(f"Error reading tree file: {str(e)}")
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
    
    def run(self, import_from=None, import_data=None):
        return MrBayesPlugin(import_from, import_data)
