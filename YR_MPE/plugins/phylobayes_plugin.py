from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QGroupBox, 
                             QFormLayout, QLineEdit, QPushButton, 
                             QFrame, QComboBox, QSpinBox,
                             QRadioButton, QLabel, QCheckBox, QTextEdit, 
                             QTabWidget, QApplication, QFileDialog, 
                             QMessageBox, QDoubleSpinBox, QSizePolicy)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
import tempfile
import os
import re
from typing import List, Optional
import subprocess
from Bio import SeqIO

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread

import logging


class PhyloBayesThread(BaseProcessThread):
    """PhyloBayes-MPI系统发育推断线程类"""
    
    def __init__(self, tool_path, mpirun_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
        self.mpirun_path = mpirun_path
    
    def get_tool_name(self):
        """返回工具名称"""
        return "PhyloBayes-MPI Phylogeny"
        
    def execute_commands(self):
        """执行PhyloBayes-MPI系统发育推断命令"""
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
                chain_name = os.path.join(os.path.abspath(input_file), os.path.basename(input_file) + "_chain")
                
                # 构建MPI命令
                cmd = [
                    self.mpirun_path,
                    "-np", str(mpi_parallel),
                    self.tool_path,
                    "-d", input_file,
                    *params,  # 其他参数
                    chain_name  # 链名称作为最后一个参数
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"PhyloBayes-MPI execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的链文件
                # 根据链名称模式查找输出文件
                import glob
                chain_pattern = f"{chain_name}*"
                chain_files = glob.glob(os.path.join(os.path.dirname(input_file), chain_pattern))
                
                if chain_files:
                    output_files.extend(chain_files)
                else:
                    self.console_output.emit(f"Warning: Could not find chain files for {input_file}", "warning")
                        
            self.progress.emit("Phylogenetic inference completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Phylogenetic inference exception: {str(e)}")


class PhyloBayesPlugin(BasePlugin):
    """PhyloBayes-MPI系统发育推断插件类"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_model_result_signal = pyqtSignal(dict)  # 导出模型结果信号
    export_phylogeny_result_signal = pyqtSignal(dict)  # 导出系统发育树结果信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化PhyloBayes插件"""
        super().__init__(import_from, import_data)
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "PhyloBayes-MPI Phylogeny"
        self.tool_name = "PhyloBayes-MPI"
        self.citation = [
            """Lartillot N., Philippe H. A Bayesian Mixture Model for Across-Site Heterogeneities in the Amino-Acid Replacement Process. Molecular Biology and Evolution 2004 21(6): 1095-1109.""",
            """Lartillot N., Philippe H. Computing Bayes factors using thermodynamic integration. Systematic Biology 2006 55:195-207.""",
            """Lartillot N., Brinkmann H., Philippe H. Suppression of long-branch attraction artefacts in the animal phylogeny using a site-heterogeneous model. BMC Evolutionary Biology 2007 Feb 8;7 Suppl 1:S4.""",
            """Nicolas Lartillot, Thomas Lepage and Samuel Blanquart. PhyloBayes 3: a Bayesian software package for phylogenetic reconstruction and molecular dating. Bioinformatics 2009 25(17): 2286–2288."""
        ]
        self.input_types = {"PHYLIP": ["phy"], "Chain File": [""]}
        self.output_types = {"Chain File": [""], "Tree File": [".con.tre"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
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
        
        # QRadioButton
        # input is a chain file / Phylip alignment file?
        self.is_alignment_file = QRadioButton("Phylip Alignment File")
        self.is_chain_file = QRadioButton("Chain File")
        self.is_alignment_file.setChecked(True)
        input_types = QHBoxLayout()
        input_types.addWidget(self.is_alignment_file)
        input_types.addWidget(self.is_chain_file)

        input_layout.addRow("Input type:", input_types)


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
            
        
        params_group = QGroupBox("Phylogenetic Parameters")
        params_layout = QFormLayout()
        params_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # Site heterogeneity model
        self.use_CAT_model = QCheckBox("Apply")
        self.use_CAT_model.setChecked(True)
        self.use_inf_CAT = QRadioButton("Infinite mixture (CAT)") # -cat
        self.use_inf_CAT.setChecked(True)
        self.use_finite_CAT = QRadioButton("Finite mixture (nCAT)") # -ncat <mixture number>
        self.mix_number_spinbox = QSpinBox()
        self.mix_number_spinbox.setValue(10)
        site_heterogeneity_layout = QHBoxLayout()
        site_heterogeneity_layout.addWidget(self.use_CAT_model) 
        site_heterogeneity_layout.addWidget(self.use_inf_CAT)
        site_heterogeneity_layout.addWidget(self.use_finite_CAT)
        site_heterogeneity_layout.addWidget(self.mix_number_spinbox)

        params_layout.addRow("Site Heterogeneity Model:", site_heterogeneity_layout)

        # Rate heterogeneity model
        self.use_Gamma_model = QCheckBox("Apply")
        self.use_Gamma_model.setChecked(True)
        gamma_categories_label = QLabel("Gamma Categories:") # -dgam <gamma categories>
        self.gamma_categories_spinbox = QSpinBox()
        self.gamma_categories_spinbox.setValue(4)
        gamma_model_layout = QHBoxLayout()
        gamma_model_layout.addWidget(self.use_Gamma_model)
        gamma_model_layout.addWidget(gamma_categories_label)
        gamma_model_layout.addWidget(self.gamma_categories_spinbox)
        params_layout.addRow("Rate Heterogeneity Model:", gamma_model_layout)

        # Substitution model
        self.subst_model_combo = QComboBox()
        self.subst_model_combo.addItem("GTR") # -gtr
        self.subst_model_combo.addItem("Poisson") # -poisson
        self.subst_model_combo.addItem("JTT (Protein only)") # -jtt
        self.subst_model_combo.addItem("WAG (Protein only)") # -wag
        self.subst_model_combo.addItem("LG (Protein only)") # -lg
        self.subst_model_combo.setCurrentText("GTR")
        params_layout.addRow("Substitution Model:", self.subst_model_combo)


        # run group
        run_group = QGroupBox("MCMC & Run Settings")
        run_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        run_layout = QFormLayout()
        run_group.setLayout(run_layout)
        layout.addWidget(run_group)

        # Number of Parallels
        self.n_mpi_parallel = QSpinBox()
        n_cpu_cores = os.cpu_count()
        # Set the maximum value to the number of CPU cores
        self.n_mpi_parallel.setMaximum(n_cpu_cores)
        # default value: max(n_cpu_cores / 4, 1)
        self.n_mpi_parallel.setValue(max(int(n_cpu_cores / 4), 1))
        run_layout.addRow("Number of MPI Parallels:", self.n_mpi_parallel) # <mpirun [UNIX] / mpiexec [WINDOWS]> -np <number>

        # Number of MCMC Generations [-x <freq> **<chainlength>**]
        self.n_mcmc_generations = QSpinBox()
        self.n_mcmc_generations.setMinimum(-1)  # Allow -1 for infinite
        self.n_mcmc_generations.setMaximum(1000000000)  # Large number as upper limit
        self.n_mcmc_generations.setValue(10000)  # Changed default from -1 to a reasonable value
        run_layout.addRow("MCMC Generations:", self.n_mcmc_generations)

        # Sampling Frequency [-x **<freq>** <chainlength>]
        self.sampling_frequency = QSpinBox()
        self.sampling_frequency.setValue(1)
        run_layout.addRow("Sampling Frequency:", self.sampling_frequency)

        # layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets

    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.phy *.phylip);;All files (*)"
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
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
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
    
    def on_text_changed(self):
        """文本输入变化处理"""
        if self.sequence_text.toPlainText().strip():
            # 如果有文本输入，清空文件选择
            self.imported_files.clear()
            for tag in self.file_tags:
                self.file_tags_layout.removeWidget(tag)
                tag.deleteLater()
            self.file_tags.clear()
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
        elif not self.imported_files:
            # 如果既没有文本也没有文件，恢复初始状态
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)

    def setup_output_tab(self):
        """设置输出标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 输出预览
        self.output_preview = QTextEdit()
        self.output_preview.setReadOnly(True)
        self.output_preview.setFont(QFont("Courier", 10))
        layout.addWidget(QLabel("Phylogenetic Tree Preview (Newick Format):"))
        layout.addWidget(self.output_preview)

    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)

    def prepare_input_files(self):
        """Prepare input files from multiple sources"""
        try:
            input_files = []
            
            # If there are imported files, use them individually
            if self.imported_files:
                for file_path in self.imported_files:
                    if os.path.exists(file_path):
                        input_files.append(file_path)
                    else:
                        QMessageBox.warning(self, "Warning", f"File does not exist: {file_path}")
                return input_files
            elif self.import_file:
                return [self.import_file]
            
            # Otherwise, use text input to create a temporary file
            sequence_text = self.sequence_text.toPlainText().strip()
            if not sequence_text and not self.import_file:
                QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
                return None
                
            # Create temporary file in PHYLIP format
            temp_file = self.create_temp_file(suffix='.phy')
            with open(temp_file, 'w') as f:
                f.write(sequence_text)
            return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None

    def get_parameters(self):
        """获取命令行参数"""
        params = []

        # 添加MPI并行数作为第一个参数
        params.append(str(self.n_mpi_parallel.value()))

        # Site heterogeneity model
        if self.use_CAT_model.isChecked():
            if self.use_inf_CAT.isChecked():
                params.append("-cat")  # Infinite CAT model
            else:
                params.extend(["-ncat", str(self.mix_number_spinbox.value())])  # Finite CAT model

        # Rate heterogeneity model
        if self.use_Gamma_model.isChecked():
            params.extend(["-dgam", str(self.gamma_categories_spinbox.value())])

        # Substitution model
        subst_model = self.subst_model_combo.currentText()
        if subst_model == "GTR":
            params.append("-gtr")
        elif subst_model == "Poisson":
            params.append("-poisson")
        elif subst_model == "JTT (Protein only)":
            params.append("-jtt")
        elif subst_model == "WAG (Protein only)":
            params.append("-wag")
        elif subst_model == "LG (Protein only)":
            params.append("-lg")

        # MCMC settings
        if self.n_mcmc_generations.value() > 0:
            params.extend(["-x", str(self.sampling_frequency.value()), str(self.n_mcmc_generations.value())])
        elif self.n_mcmc_generations.value() == -1:
            # For infinite generations, we still need to set sampling frequency
            params.extend(["-x", str(self.sampling_frequency.value()), "1000000"])  # Use large number instead of infinity

        return params

    def run_analysis(self):
        """运行分析"""
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return
            
        # 检查PhyloBayes可执行文件是否存在
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "PhyloBayes executable file not found! Please check config.json.")
            return
            
        # 检查MPI可执行文件是否存在
        mpirun_path = self.get_mpirun_path()
        if not mpirun_path or not os.path.exists(mpirun_path):
            QMessageBox.critical(self, "Error", "MPI executable file not found! Please check config.json for 'MPIRun' entry.")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting phylogenetic inference with PhyloBayes-MPI...", "info")
        
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
        
        # 在单独的线程中运行PhyloBayes
        self.analysis_thread = PhyloBayesThread(
            self.tool_path, mpirun_path, input_files, self.get_parameters(), self.imported_files
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()

    def get_mpirun_path(self):
        """获取MPI运行程序路径"""
        import json
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
            logging.error(f"Error reading MPI path from config: {e}")
            
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
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from == "YR_MPEA":
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        
        QMessageBox.information(self, "Completed", "Phylogenetic inference completed!")

    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 添加控制台消息
        self.add_console_message(f"Phylogenetic inference failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"Phylogenetic inference failed: {error_message}")

    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "Phylogenetic inference has been aborted.")

    def display_results(self, output_files):
        """显示结果"""
        if not output_files:
            return
            
        # 查找.tree文件并显示
        treefile = None
        for file in output_files:
            if file.endswith('.tre') or file.endswith('.con.tre'):
                treefile = file
                break
                
        if treefile:
            try:
                with open(treefile, 'r', encoding='utf-8') as f:
                    content = f.read()
                    self.output_preview.setPlainText(content)
            except Exception as e:
                self.add_console_message(f"Error displaying results: {str(e)}", "error")

    def import_to_platform(self):
        """将结果导入到当前平台"""
        if not hasattr(self, 'current_output_files') or not self.current_output_files:
            QMessageBox.warning(self, "Warning", "No phylogenetic results to import.")
            return
            
        try:
            # 读取系统发育树文件内容
            phylogenies = []
            for output_file in self.current_output_files:
                if output_file.endswith('.tre') or output_file.endswith('.con.tre'):
                    with open(output_file, 'r', encoding='utf-8') as f:
                        content = f.read()
                        phylogenies.append({
                            'filename': os.path.basename(output_file),
                            'content': content
                        })
            
            if not phylogenies:
                QMessageBox.warning(self, "Warning", "No phylogenetic trees found in results.")
                return
                
            # 发送信号将数据导入到平台
            self.import_phylogenies_to_platform(phylogenies)
            
            # 显示成功消息
            QMessageBox.information(self, "Success", f"Successfully imported {len(phylogenies)} phylogenetic tree(s) to the platform.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import phylogenetic trees: {str(e)}")

    def import_phylogenies_to_platform(self, phylogenies):
        """将系统发育树导入到平台的工作区"""
        # 发送信号将数据导入到平台
        self.export_phylogeny_result_signal.emit({"type": "phylogeny", "data": phylogenies})

    def handle_import_data(self, import_data):
        """处理从YR-MPEA导入的数据"""
        if isinstance(import_data, list):
            # 创建临时文件来存储导入的序列数据
            temp_file = self.create_temp_file(suffix='.phy')
            with open(temp_file, 'w') as f:
                # 转换为PHYLIP格式
                f.write(f"{len(import_data)} {len(import_data[0].seq) if import_data else 0}\n")
                for seq in import_data:
                    f.write(f"{seq.id:<10} {seq.seq}\n")
            self.temp_files.append(temp_file)
            self.import_file = temp_file
            self.imported_files = [temp_file]
            
            # 更新UI显示导入的文件
            if hasattr(self, 'file_path_edit') and self.file_path_edit:
                self.file_path_edit.setText(temp_file)
        else:
            self.import_file = None
            self.imported_files = []


# 插件入口点
class PhyloBayesPluginEntry:
    """PhyloBayes插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return PhyloBayesPlugin(import_from=import_from, import_data=import_data)