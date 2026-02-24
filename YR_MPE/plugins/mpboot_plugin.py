# mpboot_plugin.py
#
# Copyright (c) 2025 Zhi-Jie Xu
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

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QTextEdit,
                             QTabWidget, QToolButton, QApplication, QFrame, QDoubleSpinBox)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
import tempfile
import os
from typing import List, Optional, Dict

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class MPBootThread(BaseProcessThread):
    """MPBoot系统发育推断线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "MPBoot Phylogeny"
        
    def execute_commands(self):
        """执行MPBoot系统发育推断命令"""
        try:
            output_files = []
            
            # 分别处理每个输入文件
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                    
                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", "info")
                
                # 构建命令
                cmd = [
                    self.tool_path,
                    "-s", input_file,
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"MPBoot execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的.treefile文件
                treefile = input_file + ".treefile"
                if os.path.exists(treefile):
                    output_files.append(treefile)
                else:
                    self.console_output.emit(f"Warning: Could not find .treefile for {input_file}", "warning")
                        
                # 查找生成的.mpboot文件
                mpboot_file = input_file + ".mpboot"
                if os.path.exists(mpboot_file):
                    output_files.append(mpboot_file)
                
                # 查找生成的.log文件
                log_file = input_file + ".log"
                if os.path.exists(log_file):
                    output_files.append(log_file)
                
                # 查找生成的.contree文件
                contree_file = input_file + ".contree"
                if os.path.exists(contree_file):
                    output_files.append(contree_file)
            
            self.progress.emit("Phylogenetic inference completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Phylogenetic inference exception: {str(e)}")


class MPBootPlugin(BasePlugin):
    """MPBoot系统发育推断插件类"""
    
    # 定义信号
    export_phylogeny_result_signal = pyqtSignal(dict)  # 导出系统发育树结果信号
    
    def __init__(self, import_from=None, import_data=None, **kwargs):
        """初始化MPBoot插件"""
        super().__init__(import_from, import_data, **kwargs)
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "MPBoot Phylogeny"
        self.tool_name = "MPBoot"
        self.citation = [
            """D.T. Hoang, L.S. Vinh, T. Flouri, A. Stamatakis, A. von Haeseler, and B.Q. Minh (2018) MPBoot: fast phylogenetic maximum parsimony tree inference and bootstrap approximation. BMC Evol. Biol., 18, 11. https://doi.org/10.1186/s12862-018-1131-3"""
        ]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "PHYLIP": ["phy"]}
        self.output_types = {"Newick Tree": ".treefile", "MPBoot Report": ".mpboot"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    
    def setup_input_tab(self):
        """设置输入和参数标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # 输入组
        input_group = QGroupBox("Input")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # 文件输入
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select alignment files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)
        
        # 文件标签容器
        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # 文本输入
        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste alignment in FASTA or NBRF/PIR format...")
        self.sequence_text.setMaximumHeight(200)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        input_layout.addRow("Sequence text:", self.sequence_text)
        
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
            
            # 隐藏文本输入框
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
        
        # 基础参数组
        basic_params_group = QGroupBox("Basic Parameters")
        basic_params_layout = QFormLayout()
        basic_params_group.setLayout(basic_params_layout)
        layout.addWidget(basic_params_group)
        
        # 序列类型
        self.seq_type_combo = QComboBox()
        self.seq_type_combo.addItems(["AUTO", "BIN", "DNA", "AA", "CODON", "MORPH"])
        basic_params_layout.addRow("Sequence Type:", self.seq_type_combo)
        
        # 线程数
        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 5)
        self.threads_spinbox.setValue(1)
        self.threads_spinbox.setSpecialValueText("AUTO")
        basic_params_layout.addRow("Threads:", self.threads_spinbox)
        
        # 树搜索参数组
        tree_search_group = QGroupBox("Tree Search Parameters")
        tree_search_layout = QFormLayout()
        tree_search_group.setLayout(tree_search_layout)
        layout.addWidget(tree_search_group)
        
        # 初始简约树数量
        self.initial_tree_spinbox = QSpinBox()
        self.initial_tree_spinbox.setRange(1, 10000)
        self.initial_tree_spinbox.setValue(100)
        tree_search_layout.addRow("Initial Trees:", self.initial_tree_spinbox)
        
        # 最佳简约树数量
        self.best_tree_spinbox = QSpinBox()
        self.best_tree_spinbox.setRange(1, 1000)
        self.best_tree_spinbox.setValue(20)
        tree_search_layout.addRow("Best Trees:", self.best_tree_spinbox)
        
        # 候选树集合大小
        self.can_tree_spinbox = QSpinBox()
        self.can_tree_spinbox.setRange(1, 100)
        self.can_tree_spinbox.setValue(5)
        tree_search_layout.addRow("Candidate Trees:", self.can_tree_spinbox)
        
        # 迭代模式
        iter_layout = QHBoxLayout()
        self.iters_combo = QComboBox()
        self.iters_combo.addItems(["auto", "fixed"])
        self.iters_combo.currentTextChanged.connect(self.on_iters_changed)
        self.fixed_iters_value_spinbox = QSpinBox()
        self.fixed_iters_value_spinbox.setRange(1, 10000)
        self.fixed_iters_value_spinbox.setValue(100)
        self.fixed_iters_value_spinbox.setEnabled(False)
        iter_layout.addWidget(self.iters_combo)
        iter_layout.addWidget(self.fixed_iters_value_spinbox)
        tree_search_layout.addRow("Iterations:", iter_layout)
        
        # 树搜索方法
        self.method_combo = QComboBox()
        self.method_combo.addItems(["auto", "IQP", "IQPNNI (old)"])
        tree_search_layout.addRow("Method:", self.method_combo)
        
        # 扰动强度
        self.perturbation_doublespinbox = QDoubleSpinBox()
        self.perturbation_doublespinbox.setRange(0.0, 1.0)
        self.perturbation_doublespinbox.setSingleStep(0.1)
        self.perturbation_doublespinbox.setValue(0.5)
        tree_search_layout.addRow("Perturbation Strength:", self.perturbation_doublespinbox)
        
        # MPBoot 参数组
        mpboot_group = QGroupBox("MPBoot Parameters")
        mpboot_layout = QFormLayout()
        mpboot_group.setLayout(mpboot_layout)
        layout.addWidget(mpboot_group)
        
        # Ratchet 参数
        ratchet_layout = QHBoxLayout()
        self.ratchet_iter_spinbox = QSpinBox()
        self.ratchet_iter_spinbox.setRange(0, 100)
        self.ratchet_iter_spinbox.setValue(1)
        self.ratchet_iter_spinbox.setSpecialValueText("Default")
        mpboot_layout.addRow("Ratchet Iterations:", self.ratchet_iter_spinbox)
        
        self.ratchet_wgt_spinbox = QSpinBox()
        self.ratchet_wgt_spinbox.setRange(0, 100)
        self.ratchet_wgt_spinbox.setValue(1)
        mpboot_layout.addRow("Ratchet Weight:", self.ratchet_wgt_spinbox)
        
        self.ratchet_percent_spinbox = QSpinBox()
        self.ratchet_percent_spinbox.setRange(0, 100)
        self.ratchet_percent_spinbox.setValue(50)
        mpboot_layout.addRow("Ratchet Percent:", self.ratchet_percent_spinbox)
        
        self.ratchet_off_checkbox = QCheckBox("Turn off Ratchet")
        mpboot_layout.addRow("", self.ratchet_off_checkbox)
        
        # SPR 半径
        self.spr_rad_spinbox = QSpinBox()
        self.spr_rad_spinbox.setRange(1, 100)
        self.spr_rad_spinbox.setValue(3)
        mpboot_layout.addRow("SPR Radius:", self.spr_rad_spinbox)
        
        # 候选树截断
        self.cand_cutoff_spinbox = QSpinBox()
        self.cand_cutoff_spinbox.setRange(1, 100)
        self.cand_cutoff_spinbox.setValue(10)
        mpboot_layout.addRow("Candidate Cutoff (%):", self.cand_cutoff_spinbox)
        
        # 其他 MPBoot 选项
        self.opt_btree_off_checkbox = QCheckBox("Turn off tree refinement")
        mpboot_layout.addRow("", self.opt_btree_off_checkbox)
        
        self.nni_pars_checkbox = QCheckBox("Use NNI instead of SPR")
        mpboot_layout.addRow("", self.nni_pars_checkbox)
        
        # Bootstrap 参数组
        bootstrap_group = QGroupBox("Bootstrap Parameters")
        bootstrap_layout = QFormLayout()
        bootstrap_group.setLayout(bootstrap_layout)
        layout.addWidget(bootstrap_group)
        
        # UFBoot 参数
        ufboot_layout = QHBoxLayout()
        self.ufboot_checkbox = QCheckBox("Enable Bootstrap")
        self.ufboot_checkbox.setChecked(True)
        self.ufboot_spinbox = QSpinBox()
        self.ufboot_spinbox.setRange(100, 100000)
        self.ufboot_spinbox.setValue(1000)
        self.ufboot_spinbox.setEnabled(True)
        self.ufboot_checkbox.stateChanged.connect(
            lambda state: self.ufboot_spinbox.setEnabled(state == Qt.Checked)
        )
        ufboot_layout.addWidget(self.ufboot_checkbox)
        ufboot_layout.addWidget(self.ufboot_spinbox)
        bootstrap_layout.addRow("Bootstrap:", ufboot_layout)
        
        # 最大迭代次数
        self.nm_spinbox = QSpinBox()
        self.nm_spinbox.setRange(100, 10000)
        self.nm_spinbox.setValue(1000)
        bootstrap_layout.addRow("Max Iterations:", self.nm_spinbox)
        
        # 停止规则迭代次数
        self.nstep_spinbox = QSpinBox()
        self.nstep_spinbox.setRange(10, 1000)
        self.nstep_spinbox.setValue(100)
        bootstrap_layout.addRow("Stop Iterations:", self.nstep_spinbox)
        
        # 最小相关系数
        self.bcor_doublespinbox = QDoubleSpinBox()
        self.bcor_doublespinbox.setRange(0.0, 1.0)
        self.bcor_doublespinbox.setSingleStep(0.01)
        self.bcor_doublespinbox.setValue(0.99)
        bootstrap_layout.addRow("Min Correlation:", self.bcor_doublespinbox)
        
        # RELL epsilon
        self.beps_doublespinbox = QDoubleSpinBox()
        self.beps_doublespinbox.setRange(0.0, 1.0)
        self.beps_doublespinbox.setSingleStep(0.1)
        self.beps_doublespinbox.setValue(0.5)
        bootstrap_layout.addRow("RELL Epsilon:", self.beps_doublespinbox)
        
        # Consensus 参数
        consensus_group = QGroupBox("Consensus Parameters")
        consensus_layout = QFormLayout()
        consensus_group.setLayout(consensus_layout)
        layout.addWidget(consensus_group)
        
        # 最小分裂支持阈值
        self.t_threshold_doublespinbox = QDoubleSpinBox()
        self.t_threshold_doublespinbox.setRange(0.0, 1.0)
        self.t_threshold_doublespinbox.setSingleStep(0.1)
        self.t_threshold_doublespinbox.setValue(0.0)
        consensus_layout.addRow("Threshold:", self.t_threshold_doublespinbox)
        
        # Burnin
        self.bi_spinbox = QSpinBox()
        self.bi_spinbox.setRange(0, 10000)
        self.bi_spinbox.setValue(0)
        consensus_layout.addRow("Burnin:", self.bi_spinbox)
        
        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets
    
    def on_iters_changed(self, text):
        """处理迭代模式变化"""
        self.fixed_iters_value_spinbox.setEnabled(text == "fixed")
    
    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.fas *.fasta *.fa *.fna *.phy *.nexus);;All files (*)"
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
            
            # 隐藏文本输入框
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
    
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
                
            # Create temporary file
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                f.write(sequence_text)
            return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
    
    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 序列类型
        seq_type = self.seq_type_combo.currentText()
        if seq_type != "AUTO":
            params.extend(["-st", seq_type])
        
        # 树搜索参数
        params.extend(["-numpars", str(self.initial_tree_spinbox.value())])
        params.extend(["-toppars", str(self.best_tree_spinbox.value())])
        params.extend(["-numcand", str(self.can_tree_spinbox.value())])
        
        # 迭代模式
        if self.iters_combo.currentText() == "fixed":
            params.extend(["-n", str(self.fixed_iters_value_spinbox.value())])
        
        # 扰动强度
        params.extend(["-pers", str(self.perturbation_doublespinbox.value())])
        
        # 树搜索方法
        method = self.method_combo.currentText()
        if method == "IQP":
            params.extend(["-iqp"])
        elif method == "IQPNNI (old)":
            params.extend(["-iqpnni"])
        
        # MPBoot 参数
        if self.ratchet_iter_spinbox.value() > 0:
            params.extend(["-ratchet_iter", str(self.ratchet_iter_spinbox.value())])
        
        if self.ratchet_wgt_spinbox.value() > 0:
            params.extend(["-ratchet_wgt", str(self.ratchet_wgt_spinbox.value())])
        
        params.extend(["-ratchet_percent", str(self.ratchet_percent_spinbox.value())])
        
        if self.ratchet_off_checkbox.isChecked():
            params.extend(["-ratchet_off"])
        
        params.extend(["-spr_rad", str(self.spr_rad_spinbox.value())])
        params.extend(["-cand_cutoff", str(self.cand_cutoff_spinbox.value())])
        
        if self.opt_btree_off_checkbox.isChecked():
            params.extend(["-opt_btree_off"])
        
        if self.nni_pars_checkbox.isChecked():
            params.extend(["-nni_pars"])
        
        # Bootstrap 参数
        if self.ufboot_checkbox.isChecked():
            params.extend(["-bb", str(self.ufboot_spinbox.value())])
            params.extend(["-nm", str(self.nm_spinbox.value())])
            params.extend(["-nstep", str(self.nstep_spinbox.value())])
            params.extend(["-bcor", str(self.bcor_doublespinbox.value())])
            params.extend(["-beps", str(self.beps_doublespinbox.value())])
        
        # Consensus 参数
        threshold = self.t_threshold_doublespinbox.value()
        if threshold > 0.0:
            params.extend(["-t", str(threshold)])
        
        if self.bi_spinbox.value() > 0:
            params.extend(["-bi", str(self.bi_spinbox.value())])
        
        # 计算一致性树
        params.extend(["-con"])
        
        # 线程数
        threads = self.threads_spinbox.value()
        if threads > 1:
            params.extend(["-nt", str(threads)])
        
        return params
    
    def run_analysis(self):
        """运行分析"""
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "MPBoot executable file not found!")
            return
        
        # 添加控制台消息
        self.add_console_message("Starting MPBoot phylogenetic inference...", "info")
        
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
        
        # 创建并启动线程
        self.analysis_thread = MPBootThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
        )
        
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()
    
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
        """显示结果，使用IcyTree展示系统发育树"""
        if not output_files:
            QMessageBox.information(self, "error", "No output files generated")
            return
            
        # 查找.treefile文件并显示
        treefile = None
        for file in output_files:
            if file.endswith('.treefile'):
                treefile = file
                break
                
        if treefile:
            try:
                # 读取树文件内容
                with open(treefile, 'r', encoding='utf-8') as f:
                    tree_content = f.read().strip()
                
                # 确保树内容不为空
                if not tree_content:
                    QMessageBox.information(self, "error", "Tree file is empty")
                    return
                
                # 导入IcyTree插件
                from .icytree import IcyTreePlugin
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
                
                QMessageBox.information(self, "success", f"Phylogenetic tree visualization ready: {os.path.basename(treefile)}")
                
            except ImportError:
                # 如果无法导入IcyTree插件，显示错误信息
                QMessageBox.information(self, "error", "Error: IcyTree plugin not available")
                
            except Exception as e:
                error_msg = f"Error processing tree file: {str(e)}"
                QMessageBox.information(self, "error", error_msg)
                self.add_console_message(error_msg, "error")
        else:
            # 没有找到树文件，显示信息
            QMessageBox.information(self, "error", f"No treefile found. Generated {len(output_files)} file(s).")
    
    def import_to_platform(self):
        """将结果导入到当前平台"""
        if not hasattr(self, 'current_output_files') or not self.current_output_files:
            QMessageBox.warning(self, "Warning", "No phylogenetic results to import.")
            return

        try:
            # 读取系统发育树文件内容
            phylogenies = []
            for output_file in self.current_output_files:
                if output_file.endswith('.treefile'):
                    with open(output_file, 'r', encoding='utf-8') as f:
                        content = f.read()
                        phylogenies.append({
                            'filename': os.path.basename(output_file),
                            'content': content,
                            'file_path': output_file  # 添加文件路径
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
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                for seq in import_data:
                    f.write(f">{seq.id}\n{seq.seq}\n")
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
class MPBootPluginEntry:
    """MPBoot插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None, **kwargs):
        return MPBootPlugin(import_from=import_from, import_data=import_data, **kwargs)