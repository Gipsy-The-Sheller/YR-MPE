# iqtree_plugin.py
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
                             QTabWidget, QToolButton, QApplication, QFrame, QDoubleSpinBox, QDialog)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QTextCursor, QIcon
import tempfile
import os
import re
from pathlib import Path
from typing import List, Optional, Dict
from Bio import SeqIO

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from .components.alignment_common import collect_partition_sequences
from .components.model_parameter_parser import parse_iqtree_file
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess

# 导入分区模式相关模块
from .partition_mode import PartitionDefinition, PartitionMode
from .iqtree_partition_ui import IQTreePartitionDialog
from .iqtree_partition_thread import IQTreePartitionThread


class IQTreeThread(BaseProcessThread):
    """IQ-TREE phylogenetic inference worker."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters,
                 imported_files=None, workdir=None):
        super().__init__(
            tool_path, input_files, parameters, imported_files, workdir=workdir
        )
    
    def get_tool_name(self):
        """返回工具名称"""
        return "IQ-TREE Phylogeny"
        
    def execute_commands(self):
        """执行IQ-TREE系统发育推断命令"""
        try:
            output_files = []
            html_files = []
            
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
                result = self.execute_command(cmd, cwd=str(Path(input_file).parent))
                detail = (result.stderr or "").strip()
                if result.returncode != 0:
                    raise RuntimeError(
                        "IQ-TREE failed for %s%s"
                        % (os.path.basename(input_file),
                           ": " + detail if detail else "")
                    )
                
                # 查找生成的.treefile文件
                treefile = input_file + ".treefile"
                if os.path.exists(treefile):
                    output_files.append(treefile)
                else:
                    # 尝试其他可能的命名方式
                    alternate_treefile = input_file + ".phy.treefile"
                    if os.path.exists(alternate_treefile):
                        output_files.append(alternate_treefile)
                    else:
                        self.console_output.emit(f"Warning: Could not find .treefile for {input_file}", "warning")
                        
                # 查找生成的.iqtree文件
                iqtree_file = input_file + ".iqtree"
                if os.path.exists(iqtree_file):
                    output_files.append(iqtree_file)
                
                # 查找生成的.log文件
                log_file = input_file + ".log"
                if os.path.exists(log_file):
                    output_files.append(log_file)
            
            self.progress.emit("Phylogenetic inference completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])
            
        except Exception as exc:
            self.error.emit("Phylogenetic inference exception: %s" % exc)


class IQTreePlugin(BasePlugin):
    """IQ-TREE系统发育推断插件类"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_model_result_signal = pyqtSignal(dict)  # 导出模型结果信号
    export_phylogeny_result_signal = pyqtSignal(dict)  # 导出系统发育树结果信号
    
    def __init__(self, import_from=None, import_data=None, **kwargs):
        """初始化IQ-TREE插件"""
        super().__init__(import_from, import_data, **kwargs)
        
        # 初始化分区模式相关变量
        self.partition_mode_enabled = False
        self.partition_definitions: List[PartitionDefinition] = []
        self.partition_mode = PartitionMode.EL
        self.partition_bootstrap_enabled = True
        self.partition_bootstrap_replicates = 1000
        self.partition_rcluster_percent = None
        self.partition_dialog = None
        
        if import_data is not None and (
            not isinstance(import_data, list) or import_from != "YR_MPEA"
        ):
            self.handle_import_data(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "IQ-TREE Phylogeny"
        self.tool_name = "IQ-TREE 3"
        self.citation = []
        self.citation_keys = ["wong2026iqtree"]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "PHYLIP": ["phy"]}
        self.output_types = {"Newick Tree": ".treefile", "IQ-TREE Report": ".iqtree"}
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
        
        # 分区模式复选框
        self.partition_mode_checkbox = QCheckBox("Enable Partition Mode")
        self.partition_mode_checkbox.setChecked(False)
        self.partition_mode_checkbox.stateChanged.connect(self.on_partition_mode_toggled)
        input_layout.addRow("", self.partition_mode_checkbox)
        
        # 分区模式配置按钮（初始隐藏）
        partition_config_layout = QHBoxLayout()
        self.partition_config_btn = QPushButton("Configure Partitions")
        self.partition_config_btn.clicked.connect(self.open_partition_dialog)
        self.partition_config_btn.setVisible(False)
        partition_config_layout.addWidget(self.partition_config_btn)
        partition_config_layout.addStretch()
        input_layout.addRow("", partition_config_layout)
        
        # 分区状态标签（初始隐藏）
        self.partition_status_label = QLabel("No partitions defined")
        self.partition_status_label.setStyleSheet("color: #6c757d; font-style: italic;")
        self.partition_status_label.setVisible(False)
        input_layout.addRow("", self.partition_status_label)
        
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
        self.seq_type_combo.addItems(["AUTO", "DNA", "AA"])
        basic_params_layout.addRow("Sequence Type:", self.seq_type_combo)
        
        # 线程数
        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 5)
        self.threads_spinbox.setValue(1)
        self.threads_spinbox.setSpecialValueText("AUTO")
        basic_params_layout.addRow("Threads:", self.threads_spinbox)
        
        # 创建水平布局用于并列显示模型参数组和分支支持分析组
        horizontal_layout = QHBoxLayout()
        
        # 模型参数组
        model_group = QGroupBox("Substitution Model Options")
        model_layout = QFormLayout()
        model_group.setLayout(model_layout)
        
        # 模型选择
        self.model_combo = QComboBox()
        # DNA模型
        dna_models = [
            ["JC69", "JC"], ["F81"], ["K2P", "K80"], ["HKY85", "HKY"], ["TNe"], ["TN93", "TN"], ["K3P", "K81"], ["K81u", "K3Pu"], ["TPM2"], ["TPM2u"], 
            ["TPM3"], ["TPM3u"], ["TIM"], ["TIMe"], ["TIM2"], ["TIM2e"], ["TIM3"], ["TIM3e"], 
            ["TVM"], ["TVMe"], ["SYM"], ["GTR"]
        ]
        # 蛋白质模型
        protein_models = [
            "Blosum62", "cpREV", "Dayhoff", "DCMut", "EAL", "ELM", "FLAVI", "FLU", "GTR20", "HIVb", "HIVw", "JTT",
            "JTTDCMut", "LG", "mtART", "mtMAM", "mtREV", "mtZOA", "mtMet", "mtVer", "mtInv", 
            "NQ.bird", "NQ.insect", "NQ.mammal", "NQ.pfam", "NQ.plant", "NQ.yeast",
            "Poisson", "PMB", 
            "Q.bird", "Q.insect", "Q.mammal", "Q.pfam", "Q.plant", "Q.yeast",
            "rtREV", "VT", "WAG"
        ]
        
        # 创建模型映射字典和显示列表
        self.model_map = {}  # 用于存储别名到主名称的映射
        model_display_items = ['auto']

        # # 添加"DNA"分界
        # NUC_Item = QStandardItemModel("Nucleotides (DNA/RNA)")
        # NUC_Item.setFlags(Qt.NoItemFlags)
        # # 设置无法选中
        # # self.model_combo.setItemData(self.model_combo.count() - 1, Qt.ItemIsSelectable, False)
        # self.model_combo.addItem(Sepe)
        
        # 处理DNA模型
        for model_entry in dna_models:
            if isinstance(model_entry, list) and len(model_entry) == 2:
                main_model, alias = model_entry[0], model_entry[1]
                model_display_items.append(f"{main_model} ({alias})")  # 添加带别名的显示项
            elif isinstance(model_entry, list) and len(model_entry) == 1:
                model_display_items.append(model_entry[0])
            else:
                model_display_items.append(model_entry)
        
        # # 添加"AA"分界
        # model_display_items.append("Amino Acids (Protein)")
        # self.model_combo.setItemData(self.model_combo.count() - 1, Qt.ItemIsSelectable, False)

        # 处理蛋白质模型（暂时没有别名）
        for model in protein_models:
            model_display_items.append(model)
        
        # 添加所有模型到组合框
        self.model_combo.addItems(model_display_items)
        self.model_combo.currentTextChanged.connect(self.on_model_changed)
        model_layout.addRow("Substitution Model:", self.model_combo)
        
        # Rate heterogeneity参数（合并Gamma和FreeRate）
        rate_layout = QHBoxLayout()
        self.rate_combo = QComboBox()
        self.rate_combo.addItems(["Equal", "Gamma Distribution (+G)", "FreeRate Model (+R)"])
        self.rate_combo.currentTextChanged.connect(self.on_rate_combo_changed)
        rate_layout.addWidget(self.rate_combo)
        
        self.rate_categories_spinbox = QSpinBox()
        self.rate_categories_spinbox.setRange(1, 10)
        self.rate_categories_spinbox.setValue(4)
        self.rate_categories_spinbox.setEnabled(False)
        rate_layout.addWidget(self.rate_categories_spinbox)
        model_layout.addRow("Rate Heterogeneity:", rate_layout)
        
        # Invariable Sites参数
        self.invar_checkbox = QCheckBox("+I")
        model_layout.addRow("Invariable Sites:", self.invar_checkbox)

        self.ascertain_bias_checkbox = QCheckBox("+ASC")
        model_layout.addRow("Aecertain Bias:", self.ascertain_bias_checkbox)

        # # Empirical 参数
        # self.empirical_checkbox = QCheckBox("+F")
        # model_layout.addRow("Empirical Frequencies:", self.empirical_checkbox)

        # state freq selection
        self.state_freq_combo = QComboBox()
        self.state_freq_combo.addItems(["Estimated", "Empirical (+F)", "ML-optimized (+FO)", "Equal (+FQ)"])
        model_layout.addRow("State Freq.:", self.state_freq_combo)
        
        # 分支支持分析组
        bootstrap_group = QGroupBox("Branch Support Analysis")
        bootstrap_layout = QFormLayout()
        bootstrap_group.setLayout(bootstrap_layout)
        
        # UFBoot参数
        ufboot_layout = QHBoxLayout()
        self.ufboot_checkbox = QCheckBox("UFBoot")
        self.ufboot_checkbox.setChecked(True)
        self.ufboot_spinbox = QSpinBox()
        self.ufboot_spinbox.setRange(1000, 100000)
        self.ufboot_spinbox.setSingleStep(100)
        self.ufboot_spinbox.setValue(10000)
        self.ufboot_spinbox.setEnabled(True)
        self.ufboot_checkbox.stateChanged.connect(
            lambda state: self.ufboot_spinbox.setEnabled(state == Qt.Checked)
        )
        ufboot_layout.addWidget(self.ufboot_checkbox)
        ufboot_layout.addWidget(self.ufboot_spinbox)
        bootstrap_layout.addRow("Bootstrap:", ufboot_layout)
        
        # aBayes参数
        self.abayes_checkbox = QCheckBox("Approximate Bayes Test")
        bootstrap_layout.addRow("aBayes:", self.abayes_checkbox)
        
        # SH-aLRT参数
        sh_alrt_layout = QHBoxLayout()
        self.sh_alrt_checkbox = QCheckBox("SH-aLRT")
        self.sh_alrt_spinbox = QSpinBox()
        self.sh_alrt_spinbox.setRange(100, 10000)
        self.sh_alrt_spinbox.setValue(1000)
        self.sh_alrt_spinbox.setEnabled(False)
        self.sh_alrt_checkbox.stateChanged.connect(
            lambda state: self.sh_alrt_spinbox.setEnabled(state == Qt.Checked)
        )
        sh_alrt_layout.addWidget(self.sh_alrt_checkbox)
        sh_alrt_layout.addWidget(self.sh_alrt_spinbox)
        bootstrap_layout.addRow("LRT Test:", sh_alrt_layout)
        
        # 将两个组添加到水平布局中
        horizontal_layout.addWidget(model_group)
        horizontal_layout.addWidget(bootstrap_group)
        
        # 将水平布局添加到主布局中
        layout.addLayout(horizontal_layout)
        
        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets
    
    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.fas *.fasta *.fa *.fna *.phy);;All files (*)"
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

    def on_model_changed(self):
        if self.model_combo.currentText() == "auto":
            # disable +I +G +F parameters
            self.rate_combo.setEnabled(False)
            self.invar_checkbox.setEnabled(False)
            self.state_freq_combo.setEnabled(False)
            self.rate_combo.setCurrentText("Equal")
            self.invar_checkbox.setChecked(False)
            self.state_freq_combo.setCurrentText('Estimated')
        else:
            # enable +I +G +F parameters
            self.rate_combo.setEnabled(True)
            self.invar_checkbox.setEnabled(True)
            self.state_freq_combo.setEnabled(True)
    
    def on_rate_combo_changed(self, rate_option):
        """处理 Rate Heterogeneity 组合框变化"""
        # Both +G and +R take a category count: -m GTR+G8 / -m GTR+R8.
        self.rate_categories_spinbox.setEnabled(
            rate_option in ("Gamma Distribution (+G)", "FreeRate Model (+R)")
        )
    
    def setup_output_tab(self):
        """设置输出标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # Add output info label
        self.output_info = QLabel("No tree available yet.")
        self.output_info.setAlignment(Qt.AlignCenter)
        self.output_info.setStyleSheet("color: #6c757d; padding: 20px;")
        layout.addWidget(self.output_info)
        
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
    
    def _model_choices(self):
        choices = {}
        for index in range(self.model_combo.count()):
            text = self.model_combo.itemText(index)
            if text == "auto":
                choices[text] = "MFP"
            elif " (" in text and text.endswith(")"):
                choices[text] = text.split(" (")[0]
            else:
                choices[text] = text
        return choices

    def control_bindings(self):
        return [
            {"id": "input", "widget": "file_path_edit",
             "label": "Input alignment", "kind": "path",
             "cli": ["-s", "{value}"]},
            {"id": "seq_type", "widget": "seq_type_combo",
             "label": "Sequence type", "kind": "choice",
             "choices": {"AUTO": "", "DNA": "DNA", "AA": "AA"},
             "cli": ["-st", "{value}"]},
            {"id": "threads", "widget": "threads_spinbox",
             "label": "Threads", "kind": "int",
             "cli": ["-nt", "{value}"],
             "omit_values": [1]},
            {"id": "model", "widget": "model_combo",
             "label": "Substitution model", "kind": "choice",
             "choices": self._model_choices(),
             "cli": ["-m", "{value}"],
             "notes": "Base model only; rate, invariant-site, ascertainment "
                      "and state-frequency suffixes are appended to -m."},
            {"id": "rate_heterogeneity", "widget": "rate_combo",
             "label": "Rate heterogeneity", "kind": "choice",
             "choices": {
                 "Equal": [],
                 "Gamma Distribution (+G)": ["+G"],
                 "FreeRate Model (+R)": ["+R"],
             },
             "cli": ["{value}"],
             "notes": "Equal adds nothing; +G and +R both append the "
                      "category count (for example +G8 or +R8)."},
            {"id": "rate_categories", "widget": "rate_categories_spinbox",
             "label": "Rate categories", "kind": "int",
             "cli": [],
             "enabled_when": {"widget": "rate_combo",
                              "in": ["Gamma Distribution (+G)",
                                     "FreeRate Model (+R)"]},
             "notes": "Appended to -m as +G{n} when Gamma is selected or "
                      "+R{n} when FreeRate is selected."},
            {"id": "invariant_sites", "widget": "invar_checkbox",
             "label": "Invariant sites", "kind": "bool",
             "cli": ["+I"],
             "notes": "Appended to -m."},
            {"id": "ascertainment_bias", "widget": "ascertain_bias_checkbox",
             "label": "Ascertainment bias", "kind": "bool",
             "cli": ["+ASC"],
             "notes": "Appended to -m."},
            {"id": "state_freq", "widget": "state_freq_combo",
             "label": "State frequencies", "kind": "choice",
             "choices": {
                 "Estimated": "",
                 "Empirical (+F)": "+F",
                 "ML-optimized (+FO)": "+FO",
                 "Equal (+FQ)": "+FQ",
             },
             "cli": ["{value}"],
             "notes": "Appended to -m; Estimated adds nothing."},
            {"id": "ufboot", "widget": "ufboot_spinbox",
             "label": "UFBoot replicates", "kind": "int",
             "step": 100,
             "cli": ["-bb", "{value}"],
             "enabled_by": "ufboot_checkbox",
             "notes": "Only used when the UFBoot checkbox is enabled."},
            {"id": "abayes", "widget": "abayes_checkbox",
             "label": "aBayes", "kind": "bool",
             "cli": ["-abayes"]},
            {"id": "sh_alrt", "widget": "sh_alrt_spinbox",
             "label": "SH-aLRT replicates", "kind": "int",
             "cli": ["-alrt", "{value}"],
             "enabled_by": "sh_alrt_checkbox",
             "notes": "Only used when the SH-aLRT checkbox is enabled."},
            {"id": "partition_mode", "widget": "",
             "label": "Partition mode", "kind": "choice",
             "choices": {
                 "Edge-linked Proportional (-p)": "-p",
                 "Edge-linked Equal (-q)": "-q",
                 "Edge-unlinked (-Q)": "-Q",
                 "Separate Tree (-S)": "-S",
             },
             "cli": [],
             "notes": "Set in Configure Partitions; the flag consumes the "
                      "partition file."},
            {"id": "partition_bootstrap", "widget": "",
             "label": "Partition bootstrap", "kind": "composite",
             "step": 100,
             "cli": ["-bb", "{value}"],
             "notes": "Partition thread adds -bb N and -abayes; N is the "
                      "replicate count."},
            {"id": "partition_sampling", "widget": "",
             "label": "Partition resampling", "kind": "choice",
             "choices": {
                 "SITE (Site-wise resampling)": "SITE",
                 "GENE (Gene-wise resampling)": "GENE",
                 "GENESITE (Gene then Site resampling)": "GENESITE",
             },
             "cli": [],
             "notes": "Reserved by the partition dialog; not passed by the "
                      "current -bb command builder."},
            {"id": "partition_rcluster", "widget": "",
             "label": "RCluster", "kind": "composite",
             "cli": [],
             "notes": "Optional accelerated partition merging; percentage is "
                      "shown in Configure Partitions."},
        ]

    def cli_preview(self):
        input_file = (self.imported_files[0] if self.imported_files
                      else self.file_path_edit.text().strip() or "<input>")
        parameters = self.get_parameters()
        return [self.tool_path or "iqtree3", "-s", input_file] + parameters

    def get_parameters(self):
        """获取命令行参数"""
        params = []
        
        # 序列类型
        seq_type = self.seq_type_combo.currentText()
        
        if seq_type != "AUTO":
            seq_type_code = {"dna": "DNA", "prot": "AA"}[seq_type.lower()]
            params.extend(["-st", seq_type_code])
        
        # 模型
        model_text = self.model_combo.currentText()
        model = model_text
        
        # 如果选择了带别名的模型显示项，则提取主模型名称
        if " (" in model_text and ")" in model_text:
            model = model_text.split(" (")[0]
        # 如果是别名，则转换为主名称
        elif model_text in self.model_map:
            model = self.model_map[model_text]
        
        # 添加模型扩展参数
        rate_option = self.rate_combo.currentText()
        if rate_option == "Gamma Distribution (+G)":
            model += f"+G{self.rate_categories_spinbox.value()}"
        elif rate_option == "FreeRate Model (+R)":
            model += f"+R{self.rate_categories_spinbox.value()}"
        # Equal 选项不添加任何参数
        
        if self.invar_checkbox.isChecked():
            model += "+I"

        if self.ascertain_bias_checkbox.isChecked():
            model += "+ASC"

        stfreq = self.state_freq_combo.currentText()
        model += {"Estimated": "", "Empirical (+F)": "+F", "ML-optimized (+FO)": "+FO", "Equal (+FQ)": "+FQ"}[stfreq]
        
        if model.split('+')[0] != 'auto':
            params.extend(["-m", model])
        else:
            params.extend(["-m", "MFP"]) # use modelfinder to select best-fit model
        # 分支支持参数
        if self.ufboot_checkbox.isChecked():
            params.extend(["-bb", str(self.ufboot_spinbox.value())])
            
        if self.abayes_checkbox.isChecked():
            params.extend(["-abayes"])
            
        if self.sh_alrt_checkbox.isChecked():
            params.extend(["-alrt", str(self.sh_alrt_spinbox.value())])
        
        params.extend(["-redo"])
        
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
            QMessageBox.critical(self, "Error", "IQ-TREE executable file not found!")
            return
        
        # 如果启用了分区模式但没有定义分区
        if self.partition_mode_enabled and not self.partition_definitions:
            QMessageBox.warning(self, "Warning", "Partition mode is enabled but no partitions are defined. Please configure partitions first.")
            return
            
        # 添加控制台消息
        if self.partition_mode_enabled:
            self.add_console_message(f"Starting IQ-TREE phylogenetic inference with partition mode ({self.partition_mode.value})...", "info")
        else:
            self.add_console_message("Starting IQ-TREE phylogenetic inference...", "info")
        
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
        
        # 根据是否启用分区模式选择不同的线程类
        if self.partition_mode_enabled:
            # 使用分区模式线程
            model_params = {
                'seq_type': self.seq_type_combo.currentText(),
                'model': self.build_model_string(),
                'threads': self.threads_spinbox.value()
            }
            
            self.analysis_thread = IQTreePartitionThread(
                self.tool_path,
                input_files[0] if len(input_files) == 1 else input_files,
                self.partition_definitions,
                self.partition_mode,
                model_params,
                self.partition_bootstrap_enabled,
                self.partition_bootstrap_replicates,
                rcluster_percent=self.partition_rcluster_percent,
            )
        else:
            # 使用常规线程
            self.analysis_thread = IQTreeThread(
                self.tool_path, input_files, self.get_parameters(), self.imported_files,
                workdir=self.workdir,
            )
        
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        if hasattr(self.analysis_thread, "result_ready"):
            self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()
    
    def _on_result(self, payload):
        self._payload = payload

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
        if (self.import_from in ["YR_MPEA", "DATASET_MANAGER"]
                or self.result_has_consumer()):
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
        if hasattr(self, "analysis_thread") and self.analysis_thread.isRunning():
            thread = self.analysis_thread
            for signal_name in ("progress", "result_ready", "finished",
                                "error", "console_output"):
                signal = getattr(thread, signal_name, None)
                if signal is None:
                    continue
                try:
                    signal.disconnect()
                except (TypeError, RuntimeError):
                    pass
            thread.stop()
            thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Phylogenetic inference stopped by user.", "warning")
    
    def display_results(self, output_files):
        treefile = next(
            (path for path in output_files if path.endswith(".treefile")), None
        )
        if treefile:
            try:
                text = Path(treefile).read_text(encoding="utf-8", errors="replace")
            except OSError as exc:
                self.output_info.setText("Cannot read tree: %s" % exc)
                self.output_preview.clear()
                return
            self.output_info.setText("Tree: %s" % Path(treefile).name)
            self.output_preview.setPlainText(text)
        else:
            self.output_info.setText("No .treefile found in the IQ-TREE output.")
            self.output_preview.clear()

    def _build_result_spec(self):
        files = list(getattr(self, "current_output_files", []) or [])
        if not files:
            raise ValueError("No IQ-TREE output is available.")
        items = []
        artifacts = []
        for path in files:
            suffix = Path(path).suffix
            if suffix == ".treefile":
                text = Path(path).read_text(encoding="utf-8", errors="replace").strip()
                if not text:
                    continue
                trees = [line.strip() for line in text.splitlines() if line.strip()]
                if len(trees) > 1:
                    items.append(DatasetItemSpec(
                        loci_name=Path(path).stem,
                        item_type="treeset",
                        file_path=str(path),
                        format="newick",
                        data={"data": [{"content": tree} for tree in trees]},
                        metadata={"source": self.plugin_name, "tool": "IQ-TREE 3"},
                    ))
                else:
                    items.append(DatasetItemSpec(
                        loci_name=Path(path).stem,
                        item_type="phylogeny",
                        file_path=str(path),
                        format="newick",
                        data={"newick": text},
                        metadata={"source": self.plugin_name, "tool": "IQ-TREE 3"},
                    ))
            elif suffix == ".iqtree":
                params = parse_iqtree_file(path) or {}
                items.append(DatasetItemSpec(
                    loci_name=Path(path).stem,
                    item_type="model",
                    file_path=str(path),
                    format="iqtree",
                    data={"parameters": params},
                    metadata={"source": self.plugin_name, "tool": "IQ-TREE 3"},
                ))
            else:
                artifacts.append(str(path))
        if not items:
            raise ValueError("IQ-TREE produced no tree or model report.")
        dataset = DatasetSpec(
            name="IQ-TREE result",
            items=items,
            description="%d tree/model item(s)" % len(items),
            source="iqtree",
            summary={"items": len(items)},
            metadata={
                "input": self.imported_files[0] if len(self.imported_files) == 1 else "",
                "parameters": self.get_parameters()
                if isinstance(self.get_parameters(), list) else {},
            },
        )
        return ResultSpec(
            source="iqtree",
            datasets=[dataset],
            artifacts=artifacts,
        )

    def import_to_platform(self):
        if not getattr(self, "current_output_files", None):
            QMessageBox.warning(self, "No result", "Run IQ-TREE first.")
            return
        if self.result_has_consumer():
            try:
                result = self._build_result_spec()
            except Exception as exc:
                QMessageBox.critical(self, "Export failed", str(exc))
                return
            if self.export_result(result):
                self.import_to_platform_btn.setEnabled(False)
                self.import_to_platform_btn.setText("Sent to workspace")
                QMessageBox.information(
                    self, "Success", "IQ-TREE result sent to the current workspace."
                )
            return
        # Legacy fallback: emit treefiles as phylogenies.
        phylogenies = []
        for path in self.current_output_files:
            if path.endswith(".treefile"):
                phylogenies.append({
                    "filename": os.path.basename(path),
                    "content": Path(path).read_text(encoding="utf-8"),
                    "file_path": path,
                })
        if phylogenies:
            self.import_phylogenies_to_platform(phylogenies)

    def handle_import_data(self, import_data):
        """处理从YR-MPEA或Dataset Manager导入的数据"""
        # 情况1: 从 YR_MPEA 导入的单一序列列表
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

        # 情况2: 从 Dataset Manager 导入的数据
        elif isinstance(import_data, dict):
            dataset_items = import_data.get('dataset_items', [])
            dataset_config = import_data.get('dataset_config', {})

            # 筛选出 alignment 类型的 items；兼容尚未改写类型但已标记为
            # aligned 的数据项。
            from ..platforms.methods.dataset_models import ITEM_TYPE_ALIGNMENT
            selected_items = [
                item for item in dataset_items
                if item.item_type == ITEM_TYPE_ALIGNMENT
                or getattr(item, "is_aligned", False)
            ]

            if len(selected_items) == 0:
                # 没有 alignment：明确提示，而不是静默留空。
                self.import_file = None
                self.imported_files = []
                self.add_console_message(
                    "Warning: the selected dataset contains unaligned "
                    "sequence items only. Align them first (for example "
                    "with MAFFT) and run IQ-TREE on the aligned dataset.",
                    "warning",
                )
                return

            if len(selected_items) == 1:
                # 单一 alignment：优先使用文件；若只有内存序列，写临时文件。
                single = selected_items[0]
                path = single.file_path
                if not path or not os.path.exists(path):
                    if single.sequences:
                        path = self.create_temp_file(suffix=".fas")
                        SeqIO.write(single.sequences, path, "fasta")
                        self.temp_files.append(path)
                        self.add_console_message(
                            "Alignment had no file on disk; wrote its %d "
                            "sequences to a temporary FASTA." % len(single.sequences),
                            "info",
                        )
                    else:
                        path = ""
                self.import_file = path or None
                self.imported_files = [path] if path else []
            else:
                # 多个alignment，合并成超级矩阵（参考ModelFinder的做法）
                try:
                    # 检查所有alignment是否都已对齐
                    unaligned_items = [item for item in selected_items if not item.is_aligned]
                    if unaligned_items:
                        unaligned_names = [item.loci_name for item in unaligned_items]
                        warning_msg = f"The following partitions are not aligned:\n"
                        warning_msg += "\n".join(f"  - {name}" for name in unaligned_names)
                        warning_msg += "\n\nPlease align these partitions first or select only aligned partitions."
                        self.add_console_message(f"Warning: {warning_msg}", "warning")

                    # Get dataset settings
                    topo_linked = dataset_config.get('topo_linked', False)
                    edge_linked = dataset_config.get('edge_linked', False)

                    # Map to partition mode
                    if not topo_linked:
                        self.partition_mode = PartitionMode.TUL  # Topo-unlinked
                    elif edge_linked:
                        self.partition_mode = PartitionMode.TL  # Edge-linked Equal
                    else:
                        self.partition_mode = PartitionMode.EUL  # Edge-unlinked

                    # 合并所有选中的 partition 为 supermatrix
                    supermatrix_sequences = {}
                    partition_definitions = []
                    current_pos = 1

                    # Use the union of taxa across all partitions. Taking the
                    # first partition's list drops taxa that are absent there
                    # but present in other loci (missing-data datasets).
                    all_seq_names, records_by_item = collect_partition_sequences(
                        selected_items
                    )
                    partition_lengths = []
                    for item, records in zip(selected_items, records_by_item):
                        length = item.length
                        if not length and records:
                            length = max(len(str(record.seq)) for record in records)
                        partition_lengths.append(length)

                    # 合并序列
                    for seq_name in all_seq_names:
                        supermatrix_seq = ""
                        for records, length in zip(records_by_item, partition_lengths):
                            seq = next((s for s in records if s.id == seq_name), None)
                            if seq:
                                supermatrix_seq += str(seq.seq)
                            else:
                                # 如果某个 partition 缺少该序列，用 ? 填充
                                supermatrix_seq += "?" * length
                        supermatrix_sequences[seq_name] = supermatrix_seq

                    # 创建 supermatrix 临时文件
                    temp_file = self.create_temp_file(suffix='.fas')

                    with open(temp_file, 'w') as f:
                        for seq_name, seq_content in supermatrix_sequences.items():
                            f.write(f">{seq_name}\n{seq_content}\n")

                    self.import_file = temp_file
                    self.imported_files = [temp_file]

                    # 计算 partition 坐标并创建 partition definitions
                    for item, length in zip(selected_items, partition_lengths):
                        end_pos = current_pos + length - 1
                        partition_def = PartitionDefinition(
                            name=item.loci_name,
                            file_path="",  # 空字符串表示使用 supermatrix 的坐标
                            seq_type="DNA",  # 默认 DNA，后续会检测
                            model_range=f"{current_pos}-{end_pos}",
                            selected_model=None
                        )
                        partition_definitions.append(partition_def)
                        current_pos = end_pos + 1

                    self.partition_definitions = partition_definitions

                    # 检测序列类型冲突
                    self._detect_sequence_type_conflicts(selected_items, partition_definitions)

                    # 启用分区模式
                    self.partition_mode_enabled = True

                    # 更新 UI
                    if hasattr(self, 'file_path_edit') and self.file_path_edit:
                        self.file_path_edit.setText(temp_file)

                    if hasattr(self, 'partition_mode_checkbox'):
                        self.partition_mode_checkbox.setChecked(True)

                    # 更新分区状态显示
                    self.update_partition_status()

                    # 添加控制台消息
                    self.add_console_message(f"Dataset imported: {len(selected_items)} partitions", "info")
                    self.add_console_message(f"Partition mode: {self.partition_mode.value}", "info")

                except Exception as e:
                    self.add_console_message(f"Failed to import dataset: {str(e)}", "error")
                    self.import_file = None
                    self.imported_files = []
                    self.partition_definitions = []

        # 其他情况
        else:
            self.import_file = None
            self.imported_files = []

    def _detect_sequence_type_conflicts(self, alignment_items, partition_definitions):
        """检测序列类型冲突（参考ModelFinder）"""
        from ..platforms.methods.partition_utils import detect_sequence_type_conflicts

        has_conflict, conflict_msg = detect_sequence_type_conflicts(alignment_items)

        if has_conflict:
            self.add_console_message(f"Warning: {conflict_msg}", "warning")
            # 更新分区定义中的序列类型
            for i, item in enumerate(alignment_items):
                from ..platforms.methods.partition_utils import detect_sequence_type
                if item.sequences:
                    seq_type = detect_sequence_type(item.sequences[0])
                    partition_definitions[i].seq_type = seq_type
        else:
            # 没有冲突，统一设置序列类型
            if alignment_items and alignment_items[0].sequences:
                from ..platforms.methods.partition_utils import detect_sequence_type
                seq_type = detect_sequence_type(alignment_items[0].sequences[0])
                for partition_def in partition_definitions:
                    partition_def.seq_type = seq_type

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
        # 更新分区模式
        from .partition_mode import PartitionMode
        mode_map = {
            'EL': PartitionMode.EL,
            'TL': PartitionMode.TL,
            'EUL': PartitionMode.EUL,
            'TUL': PartitionMode.TUL
        }
        self.partition_mode = mode_map.get(partition_model_config.get('partition_mode', 'EL'), PartitionMode.EL)

        # 用模型信息覆盖分区定义
        model_partitions = partition_model_config.get('partitions', [])
        if model_partitions and self.partition_definitions:
            # 确保分区数量匹配
            if len(model_partitions) == len(self.partition_definitions):
                for i, model_part in enumerate(model_partitions):
                    # 用模型信息更新分区定义
                    # 支持两种字段名：'model' 和 'best_model'
                    model_code = model_part.get('model') or model_part.get('best_model')
                    if model_code:
                        self.partition_definitions[i].selected_model = model_code
                        self.add_console_message(f"  Partition {i+1}: {model_part.get('name', '')} -> {model_code}", "info")
                    if model_part.get('seq_type'):
                        self.partition_definitions[i].seq_type = model_part['seq_type']
            else:
                # 分区数量不匹配，警告用户
                self.add_console_message(f"Warning: Model has {len(model_partitions)} partitions but dataset has {len(self.partition_definitions)} partitions. Partition count mismatch!", "warning")

        # 更新分区状态显示
        self.update_partition_status()

        # 添加控制台消息
        self.add_console_message(f"Partition model config applied", "info")
        self.add_console_message(f"  Mode: {self.partition_mode.value}", "info")
        self.add_console_message(f"  Models: {len(model_partitions)} partitions with specified models", "info")

    def on_partition_mode_toggled(self, state):
        """处理分区模式复选框切换"""
        self.partition_mode_enabled = (state == Qt.Checked)
        
        if self.partition_mode_enabled:
            # 启用分区模式
            self.partition_config_btn.setVisible(True)
            self.partition_status_label.setVisible(True)
            self.add_console_message("Partition mode enabled. Please configure partitions.", "info")

            # 禁用Substitution Model Options（分区模式下使用MFP+MERGE）
            self.model_combo.setEnabled(False)
            self.rate_combo.setEnabled(False)
            self.rate_categories_spinbox.setEnabled(False)
            self.invar_checkbox.setEnabled(False)
            self.ascertain_bias_checkbox.setEnabled(False)
            self.state_freq_combo.setEnabled(False)

            # UFBOOT保持启用（分区模式下仍然需要Bootstrap）
            # 不禁用UFBOOT
        else:
            # 禁用分区模式
            self.partition_config_btn.setVisible(False)
            self.partition_status_label.setVisible(False)
            self.partition_definitions = []
            self.add_console_message("Partition mode disabled.", "info")

            # 启用Substitution Model Options
            self.model_combo.setEnabled(True)
            self.rate_combo.setEnabled(True)
            self.on_rate_combo_changed(self.rate_combo.currentText())  # 根据 rate_combo 的选择启用/禁用 spinbox
            self.invar_checkbox.setEnabled(True)
            self.ascertain_bias_checkbox.setEnabled(True)
            self.state_freq_combo.setEnabled(True)

            # UFBOOT保持启用
            self.ufboot_checkbox.setEnabled(True)
            self.ufboot_spinbox.setEnabled(self.ufboot_checkbox.isChecked())
    
    def open_partition_dialog(self):
        """打开分区配置对话框"""
        try:
            # 创建分区对话框
            if self.partition_dialog is None:
                self.partition_dialog = IQTreePartitionDialog(parent=self)
                # 连接信号
                self.partition_dialog.partitions_selected.connect(self.handle_partitions_selected)
            
            # 设置当前分区定义
            if self.partition_definitions:
                self.partition_dialog.set_partitions(self.partition_definitions)
                self.partition_dialog.set_partition_mode(self.partition_mode)
            
            # 显示对话框
            if self.partition_dialog.exec_() == QDialog.Accepted:
                # 用户点击了OK
                partitions = self.partition_dialog.get_partitions()
                mode = self.partition_dialog.get_partition_mode()
                bootstrap_enabled = self.partition_dialog.enable_bootstrap.isChecked()
                bootstrap_replicates = self.partition_dialog.bootstrap_replicates.value()
                rcluster_percent = self.partition_dialog.get_rcluster_percent()

                # 更新分区定义
                self.partition_definitions = partitions
                self.partition_mode = mode
                self.partition_bootstrap_enabled = bootstrap_enabled
                self.partition_bootstrap_replicates = bootstrap_replicates
                self.partition_rcluster_percent = rcluster_percent
                
                # 更新状态显示
                self.update_partition_status()
                
                self.add_console_message(f"Partition configuration updated: {len(partitions)} partitions, mode: {mode.value}", "info")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open partition dialog: {str(e)}")
    
    def handle_partitions_selected(self, partitions: List[PartitionDefinition]):
        """处理分区定义信号"""
        self.partition_definitions = partitions
        self.update_partition_status()
        self.add_console_message(f"Received {len(partitions)} partition definitions", "info")
    
    def update_partition_status(self):
        """更新分区状态显示"""
        if not self.partition_definitions:
            self.partition_status_label.setText("No partitions defined")
        else:
            mode_text = {
                PartitionMode.EL: "Edge-linked",
                PartitionMode.TL: "Edge-linked",
                PartitionMode.EUL: "Edge-unlinked",
                PartitionMode.TUL: "Topo-unlinked"
            }.get(self.partition_mode, "Unknown")
            
            partition_names = [p.name for p in self.partition_definitions]
            if len(partition_names) <= 3:
                names_text = ", ".join(partition_names)
            else:
                names_text = f"{partition_names[0]}, {partition_names[1]}, ... ({len(partition_names)} total)"
            
            bootstrap_text = f", Bootstrap: {self.partition_bootstrap_replicates}" if self.partition_bootstrap_enabled else ""
            self.partition_status_label.setText(f"Mode: {mode_text} | Partitions: {names_text}{bootstrap_text}")
    
    def build_model_string(self) -> str:
        """构建模型字符串（用于分区模式）"""
        # 如果分区中已指定模型，使用MFP+MERGE进行模型选择和合并
        has_models = any(p.selected_model for p in self.partition_definitions)

        if has_models:
            return "MFP+MERGE"
        else:
            # 使用用户指定的模型
            model_text = self.model_combo.currentText()
            model = model_text

            if " (" in model_text and ")" in model_text:
                model = model_text.split(" (")[0]

            # 添加速率异质性参数
            rate_option = self.rate_combo.currentText()
            if rate_option == "Gamma Distribution (+G)":
                model += f"+G{self.rate_categories_spinbox.value()}"
            elif rate_option == "FreeRate Model (+R)":
                model += f"+R{self.rate_categories_spinbox.value()}"
            # Equal 选项不添加任何参数

            if self.invar_checkbox.isChecked():
                model += "+I"

            if self.ascertain_bias_checkbox.isChecked():
                model += "+ASC"

            stfreq = self.state_freq_combo.currentText()
            model += {"Estimated": "", "Empirical (+F)": "+F", "ML-optimized (+FO)": "+FO", "Equal (+FQ)": "+FQ"}[stfreq]

            return model

    def set_partition_config(self, mode: str, partitions: list, alignment_file: str):
        """
        设置分区配置（从平台导入）

        Args:
            mode: 分区模式 (EL/TL/EUL/TUL)
            partitions: 分区定义列表
            alignment_file: 超级矩阵文件路径
        """
        from .partition_mode import PartitionMode, PartitionDefinition

        # 映射模式字符串到枚举
        mode_map = {
            'EL': PartitionMode.EL,
            'TL': PartitionMode.TL,
            'EUL': PartitionMode.EUL,
            'TUL': PartitionMode.TUL
        }

        self.partition_mode = mode_map.get(mode, PartitionMode.EL)

        # 转换分区定义
        self.partition_definitions = []
        for part in partitions:
            part_def = PartitionDefinition(
                name=part.get('name', ''),
                file_path='',  # 使用超级矩阵的坐标，不需要文件路径
                seq_type=part.get('seq_type', 'DNA'),
                model_range=part.get('range', ''),
                selected_model=part.get('model', None)
            )
            self.partition_definitions.append(part_def)

        # 设置超级矩阵文件（设置为导入文件，这样 prepare_input_files 才能找到）
        self.partition_alignment_file = alignment_file
        self.import_file = alignment_file  # 关键：设置为导入文件

        # 启用分区模式
        self.partition_mode_enabled = True
        self.partition_mode_checkbox.setChecked(True)

        # 更新 UI 显示
        self.update_partition_status()

        # 添加控制台消息
        self.add_console_message(f"Partition config imported from platform", "info")
        self.add_console_message(f"  Mode: {self.partition_mode.value}", "info")
        self.add_console_message(f"  Partitions: {len(self.partition_definitions)}", "info")
        self.add_console_message(f"  Alignment: {os.path.basename(alignment_file)}", "info")

        # 显示分区模型信息
        for i, part_def in enumerate(self.partition_definitions, 1):
            model_info = part_def.selected_model if part_def.selected_model else "Auto (MFP+MERGE)"
            self.add_console_message(f"    Partition {i}: {part_def.name} ({part_def.model_range}) - {model_info}", "info")


# 插件入口点
class IQTreePluginEntry:
    """IQ-TREE插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None, **kwargs):
        return IQTreePlugin(import_from=import_from, import_data=import_data, **kwargs)