#
# astral_plugin.py
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

from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit,
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, QDoubleSpinBox, QRadioButton, QButtonGroup)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon, QFont
import tempfile
import os
from typing import List, Optional

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class AstralThread(BaseProcessThread):
    """ASTRAL分析线程类"""

    def __init__(self, tool_path, input_files, parameters, imported_files=None, output_file=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
        self.output_file = output_file

    def get_tool_name(self):
        """返回工具名称"""
        return "ASTRAL"

    def execute_commands(self):
        """执行ASTRAL分析命令"""
        try:
            output_files = []

            # 处理输入文件
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break

                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", "info")

                # 如果指定了输出文件，则使用指定的输出文件，否则创建临时输出文件
                if self.output_file:
                    output_file = self.output_file
                else:
                    output_file = self.create_temp_file(suffix='.tre')

                # 构建命令（input_file 已经在 parameters 中）
                cmd = [
                    self.tool_path,
                    *self.parameters
                ]

                # 执行命令
                result = self.execute_command(cmd)

                if result.returncode != 0:
                    self.error.emit(f"ASTRAL execution failed: {result.stderr}")
                    return

                # 将stdout写入输出文件
                with open(output_file, 'w', encoding='utf-8') as f:
                    f.write(result.stdout)

                output_files.append(output_file)

            self.progress.emit("ASTRAL analysis completed")
            self.finished.emit(output_files, [])

        except Exception as e:
            self.error.emit(f"Analysis exception: {str(e)}")


class AdvancedOptionsDialog(QDialog):
    """高级选项对话框"""
    def __init__(self, astral_version="astral4", parent=None):
        super().__init__(parent)
        self.astral_version = astral_version
        self.setWindowTitle("Advanced Options")
        self.resize(500, 600)

        layout = QVBoxLayout()

        # 创建滚动区域以容纳所有高级参数
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # Constraint tree
        constraint_group = QGroupBox("Constraint Tree")
        constraint_layout = QHBoxLayout()
        self.constraint_edit = QLineEdit()
        constraint_browse_btn = QPushButton("Browse")
        constraint_browse_btn.clicked.connect(lambda: self.browse_file(self.constraint_edit, "Newick files (*.nwk *.newick *.tre *.tree);;All files (*)"))
        constraint_label = QLabel("Newick file containing a binary species tree:")
        constraint_label.setAlignment(Qt.AlignLeft)
        constraint_layout.addWidget(constraint_label)
        constraint_layout.addWidget(self.constraint_edit)
        constraint_layout.addWidget(constraint_browse_btn)
        constraint_group.setLayout(constraint_layout)
        scroll_layout.addWidget(constraint_group)

        # Guide tree
        guide_group = QGroupBox("Guide Tree")
        guide_layout = QHBoxLayout()
        self.guide_edit = QLineEdit()
        guide_browse_btn = QPushButton("Browse")
        guide_browse_btn.clicked.connect(lambda: self.browse_file(self.guide_edit, "Newick files (*.nwk *.newick *.tre *.tree);;All files (*)"))
        guide_label = QLabel("Binary trees as guide trees:")
        guide_label.setAlignment(Qt.AlignLeft)
        guide_layout.addWidget(guide_label)
        guide_layout.addWidget(self.guide_edit)
        guide_layout.addWidget(guide_browse_btn)
        guide_group.setLayout(guide_layout)
        scroll_layout.addWidget(guide_group)

        # Rounds
        round_group = QGroupBox("Initial Rounds")
        round_layout = QHBoxLayout()
        self.round_spin = QSpinBox()
        self.round_spin.setRange(1, 100)
        self.round_spin.setValue(4)
        round_layout.addWidget(self.round_spin)
        round_desc = QLabel("Number of initial rounds of placements")
        round_desc.setAlignment(Qt.AlignLeft)
        round_layout.addWidget(round_desc)
        round_group.setLayout(round_layout)
        scroll_layout.addWidget(round_group)

        # Subsample rounds
        subsample_group = QGroupBox("Subsample Rounds")
        subsample_layout = QHBoxLayout()
        self.subsample_spin = QSpinBox()
        self.subsample_spin.setRange(1, 100)
        self.subsample_spin.setValue(4)
        subsample_layout.addWidget(self.subsample_spin)
        subsample_desc = QLabel("Number of rounds of subsampling per exploration step")
        subsample_desc.setAlignment(Qt.AlignLeft)
        subsample_layout.addWidget(subsample_desc)
        subsample_group.setLayout(subsample_layout)
        scroll_layout.addWidget(subsample_group)

        # Mapping file
        mapping_group = QGroupBox("Mapping File")
        mapping_layout = QHBoxLayout()
        self.mapping_edit = QLineEdit()
        mapping_browse_btn = QPushButton("Browse")
        mapping_browse_btn.clicked.connect(lambda: self.browse_file(self.mapping_edit, "Text files (*.txt);;All files (*)"))
        mapping_label = QLabel("Gene name to taxon name maps:")
        mapping_label.setAlignment(Qt.AlignLeft)
        mapping_layout.addWidget(mapping_label)
        mapping_layout.addWidget(self.mapping_edit)
        mapping_layout.addWidget(mapping_browse_btn)
        mapping_group.setLayout(mapping_layout)
        scroll_layout.addWidget(mapping_group)

        # Root species
        root_group = QGroupBox("Root Species")
        root_layout = QHBoxLayout()
        self.root_edit = QLineEdit()
        self.root_edit.setPlaceholderText("Species name to root at")
        root_label = QLabel("Root at the given species:")
        root_label.setAlignment(Qt.AlignLeft)
        root_layout.addWidget(root_label)
        root_layout.addWidget(self.root_edit)
        root_group.setLayout(root_layout)
        scroll_layout.addWidget(root_group)

        # Proportion
        proportion_group = QGroupBox("Proportion of Taxa")
        proportion_layout = QHBoxLayout()
        self.proportion_spin = QDoubleSpinBox()
        self.proportion_spin.setRange(0.01, 1.0)
        self.proportion_spin.setValue(0.25)
        self.proportion_spin.setSingleStep(0.05)
        proportion_layout.addWidget(self.proportion_spin)
        proportion_desc = QLabel("Proportion of taxa in the subsample in naive algorithm")
        proportion_desc.setAlignment(Qt.AlignLeft)
        proportion_layout.addWidget(proportion_desc)
        proportion_group.setLayout(proportion_layout)
        scroll_layout.addWidget(proportion_group)

        # Seed
        seed_group = QGroupBox("Random Seed")
        seed_layout = QHBoxLayout()
        self.seed_spin = QSpinBox()
        self.seed_spin.setRange(0, 999999)
        self.seed_spin.setValue(233)
        seed_layout.addWidget(self.seed_spin)
        seed_desc = QLabel("Seed for pseudorandomness")
        seed_desc.setAlignment(Qt.AlignLeft)
        seed_layout.addWidget(seed_desc)
        seed_group.setLayout(seed_layout)
        scroll_layout.addWidget(seed_group)

        # Length unit
        length_group = QGroupBox("Length Unit")
        length_layout = QHBoxLayout()
        self.length_combo = QComboBox()
        self.length_combo.addItems(["SULength", "CULength"])
        self.length_combo.setEditable(False)
        length_layout.addWidget(self.length_combo)
        length_desc = QLabel("SULength: substitution-per-site unit; CULength: coalescent unit")
        length_desc.setAlignment(Qt.AlignLeft)
        length_layout.addWidget(length_desc)
        length_group.setLayout(length_layout)
        scroll_layout.addWidget(length_group)

        # Lambda
        lambda_group = QGroupBox("Lambda")
        lambda_layout = QHBoxLayout()
        self.lambda_spin = QDoubleSpinBox()
        self.lambda_spin.setRange(0.01, 10.0)
        self.lambda_spin.setValue(0.5)
        self.lambda_spin.setSingleStep(0.1)
        lambda_layout.addWidget(self.lambda_spin)
        lambda_desc = QLabel("Rate lambda of Yule process under which the species tree is modeled")
        lambda_desc.setAlignment(Qt.AlignLeft)
        lambda_layout.addWidget(lambda_desc)
        lambda_group.setLayout(lambda_layout)
        scroll_layout.addWidget(lambda_group)

        # Gene length
        gene_length_group = QGroupBox("Gene Length")
        gene_length_layout = QHBoxLayout()
        self.gene_length_spin = QDoubleSpinBox()
        self.gene_length_spin.setRange(1.0, 100000.0)
        self.gene_length_spin.setValue(1000.0)
        self.gene_length_spin.setSingleStep(100.0)
        gene_length_layout.addWidget(self.gene_length_spin)
        gene_length_desc = QLabel("Average gene sequence length")
        gene_length_desc.setAlignment(Qt.AlignLeft)
        gene_length_layout.addWidget(gene_length_desc)
        gene_length_group.setLayout(gene_length_layout)
        scroll_layout.addWidget(gene_length_group)

        # Downweight repeat
        downweight_group = QGroupBox("Downweight Repeat")
        downweight_layout = QHBoxLayout()
        self.downweight_spin = QDoubleSpinBox()
        self.downweight_spin.setRange(0.1, 100.0)
        self.downweight_spin.setValue(1.0)
        self.downweight_spin.setSingleStep(0.1)
        downweight_layout.addWidget(self.downweight_spin)
        downweight_desc = QLabel("The number of trees sampled for each locus")
        downweight_desc.setAlignment(Qt.AlignLeft)
        downweight_layout.addWidget(downweight_desc)
        downweight_group.setLayout(downweight_layout)
        scroll_layout.addWidget(downweight_group)

        # Verbose level
        verbose_group = QGroupBox("Verbose Level")
        verbose_layout = QHBoxLayout()
        self.verbose_combo = QComboBox()
        self.verbose_combo.addItems(["1 - Minimum", "2 - Normal"])
        self.verbose_combo.setEditable(False)
        verbose_layout.addWidget(self.verbose_combo)
        verbose_desc = QLabel("Level of logging")
        verbose_desc.setAlignment(Qt.AlignLeft)
        verbose_layout.addWidget(verbose_desc)
        verbose_group.setLayout(verbose_layout)
        scroll_layout.addWidget(verbose_group)

        # ASTRAL-Hybrid specific options
        if self.astral_version == "astral-hybrid":
            hybrid_group = QGroupBox("ASTRAL-Hybrid Specific Options")
            hybrid_layout = QVBoxLayout()

            # Mode
            mode_layout = QHBoxLayout()
            mode_label = QLabel("Weighting Mode:")
            mode_label.setAlignment(Qt.AlignLeft)
            self.mode_combo = QComboBox()
            self.mode_combo.addItems([
                "1 - Hybrid weighting",
                "2 - Support only",
                "3 - Length only",
                "4 - Unweighted"
            ])
            self.mode_combo.setEditable(False)
            mode_layout.addWidget(mode_label)
            mode_layout.addWidget(self.mode_combo)
            hybrid_layout.addLayout(mode_layout)

            # Tree weights file
            weights_layout = QHBoxLayout()
            weights_label = QLabel("Tree Weights File:")
            weights_label.setAlignment(Qt.AlignLeft)
            self.weights_edit = QLineEdit()
            weights_browse_btn = QPushButton("Browse")
            weights_browse_btn.clicked.connect(lambda: self.browse_file(self.weights_edit, "Text files (*.txt);;All files (*)"))
            weights_layout.addWidget(weights_label)
            weights_layout.addWidget(self.weights_edit)
            weights_layout.addWidget(weights_browse_btn)
            hybrid_layout.addLayout(weights_layout)

            # Min support value
            min_layout = QHBoxLayout()
            min_label = QLabel("Min Support Value:")
            min_label.setAlignment(Qt.AlignLeft)
            self.min_spin = QDoubleSpinBox()
            self.min_spin.setRange(-1.0, 1.0)
            self.min_spin.setValue(-1.0)
            self.min_spin.setSingleStep(0.1)
            self.min_spin.setSpecialValueText("Auto-detect")
            min_layout.addWidget(min_label)
            min_layout.addWidget(self.min_spin)
            hybrid_layout.addLayout(min_layout)

            # Max support value
            max_layout = QHBoxLayout()
            max_label = QLabel("Max Support Value:")
            max_label.setAlignment(Qt.AlignLeft)
            self.max_spin = QDoubleSpinBox()
            self.max_spin.setRange(-1.0, 1.0)
            self.max_spin.setValue(-1.0)
            self.max_spin.setSingleStep(0.1)
            self.max_spin.setSpecialValueText("Auto-detect")
            max_layout.addWidget(max_label)
            max_layout.addWidget(self.max_spin)
            hybrid_layout.addLayout(max_layout)

            # Default support value
            default_layout = QHBoxLayout()
            default_label = QLabel("Default Support Value:")
            default_label.setAlignment(Qt.AlignLeft)
            self.default_spin = QDoubleSpinBox()
            self.default_spin.setRange(0.0, 1.0)
            self.default_spin.setValue(0.0)
            self.default_spin.setSingleStep(0.1)
            default_layout.addWidget(default_label)
            default_layout.addWidget(self.default_spin)
            hybrid_layout.addLayout(default_layout)

            hybrid_group.setLayout(hybrid_layout)
            scroll_layout.addWidget(hybrid_group)

        # ASTRAL-Pro3 specific options
        if self.astral_version == "astral-pro3":
            pro3_group = QGroupBox("ASTRAL-Pro3 Specific Options")
            pro3_layout = QVBoxLayout()

            # Exit mode
            exit_layout = QHBoxLayout()
            exit_label = QLabel("Exit Mode:")
            exit_label.setAlignment(Qt.AlignLeft)
            self.exit_combo = QComboBox()
            self.exit_combo.addItems([
                "0 - Print warning when polytomies",
                "1 - Resolve polytomies",
                "2 - Print rooted/tagged and exit"
            ])
            self.exit_combo.setEditable(False)
            exit_layout.addWidget(exit_label)
            exit_layout.addWidget(self.exit_combo)
            pro3_layout.addLayout(exit_layout)

            pro3_group.setLayout(pro3_layout)
            scroll_layout.addWidget(pro3_group)

        scroll_layout.addStretch()
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)

        layout.addWidget(scroll_area)

        # 对话框按钮
        button_layout = QHBoxLayout()
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addStretch()
        button_layout.addWidget(ok_btn)
        button_layout.addWidget(cancel_btn)

        layout.addLayout(button_layout)
        self.setLayout(layout)

    def browse_file(self, target_edit, filter_text):
        """通用文件浏览函数"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select File",
            "",
            filter_text
        )
        if file_path:
            target_edit.setText(file_path)


class AstralPlugin(BasePlugin):
    """ASTRAL插件：用于物种树推断的基于基因树四叶树采样的工具"""

    # 类常量：支持的 ASTRAL 版本
    ASTRAL_VERSIONS = ["astral4", "astral-pro3", "wastral", "astral-hybrid"]

    def __init__(self, import_from=None, import_data=None):
        super().__init__(import_from, import_data)
        
        # 如果提供了import_data且它是包含多棵基因树的newick文件路径，则设置输入文件
        if import_data and isinstance(import_data, str) and os.path.exists(import_data):
            # 延迟设置，因为file_path_edit在setup_input_tab中创建
            self._preimport_gene_trees_file = import_data
        # 如果提供了import_data且它是来自dataset_selection_manager的TreeSet数据
        elif import_from == "DATASET_MANAGER" and import_data and isinstance(import_data, dict):
            dataset_items = import_data.get('dataset_items', [])
            # 查找 ITEM_TYPE_TREESET 类型的数据项
            from ..platforms.methods.dataset_models import ITEM_TYPE_TREESET
            treeset_items = [item for item in dataset_items if item.item_type == ITEM_TYPE_TREESET]
            if treeset_items:
                # 创建临时文件，包含所有树的 Newick 字符串
                import tempfile
                trees = []
                for tree_data in treeset_items[0].data.get('data', []):
                    if 'content' in tree_data:
                        trees.append(tree_data['content'])
                    elif 'file_path' in tree_data:
                        try:
                            with open(tree_data['file_path'], 'r') as f:
                                trees.append(f.read().strip())
                        except:
                            pass
                
                if trees:
                    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.trees', delete=False)
                    temp_file.write('\n'.join(trees))
                    temp_file.close()
                    self._preimport_gene_trees_file = temp_file.name

    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "ASTRAL"
        self.tool_name = "astral"
        self.citation = [
            "Chao Zhang, Rasmus Nielsen, Siavash Mirarab, ASTER: A Package for Large-Scale Phylogenomic Reconstructions, Molecular Biology and Evolution, Volume 42, Issue 8, August 2025, msaf172, https://doi.org/10.1093/molbev/msaf172"
        ]
        self.input_types = {
            "Newick": ["nwk", "newick", "tre", "tree"]
        }
        self.output_types = {
            "Newick": ".tre"
        }
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)

        # ASTRAL 版本选择
        version_group = QGroupBox("ASTRAL Version")
        version_layout = QHBoxLayout()
        version_group.setLayout(version_layout)

        version_label = QLabel("Select ASTRAL Version:")
        version_label.setAlignment(Qt.AlignLeft)
        version_layout.addWidget(version_label)

        self.version_combo = QComboBox()
        self.version_combo.addItems([
            "ASTRAL-IV (Standard species tree inference)",
            "ASTRAL-Pro3 (Paralogs and orthologs)",
            "Weighted ASTRAL (wastral)",
            "ASTRAL-Hybrid (Hybrid weighting)"
        ])
        self.version_combo.setEditable(False)
        self.version_combo.currentIndexChanged.connect(self.on_version_changed)
        version_layout.addWidget(self.version_combo)

        layout.addWidget(version_group)

        # 文件选择组
        file_group = QGroupBox("Input File")
        file_layout = QFormLayout()
        file_group.setLayout(file_layout)

        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select input file...")
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_input_file)
        import_dataset_btn = QPushButton("Import from Selected TreeSet")
        import_dataset_btn.clicked.connect(self.import_from_selected_treeset)
        file_hbox = QHBoxLayout()
        file_hbox.addWidget(self.file_path_edit)
        file_hbox.addWidget(browse_btn)
        file_hbox.addWidget(import_dataset_btn)
        file_layout.addRow("Input File:", file_hbox)

        layout.addWidget(file_group)

        # 参数设置组
        params_group = QGroupBox("Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)

        # 线程数
        self.thread_spin = QSpinBox()
        self.thread_spin.setRange(1, 64)
        self.thread_spin.setValue(1)
        params_layout.addRow("Threads:", self.thread_spin)

        # 支持选项
        self.support_combo = QComboBox()
        self.support_combo.addItems([
            "0 - No branch or support",
            "1 - Length and support only",
            "2 - Detailed",
            "3 - freqQuad.csv"
        ])
        self.support_combo.setEditable(False)
        self.support_combo.setCurrentIndex(1)  # 默认选择"1 - Length and support only"
        params_layout.addRow("Support Option:", self.support_combo)

        # 输出文件
        self.output_path_edit = QLineEdit()
        self.output_path_edit.setPlaceholderText("Specify output file (optional)...")
        output_browse_btn = QPushButton("Browse")
        output_browse_btn.clicked.connect(self.browse_output_file)
        output_hbox = QHBoxLayout()
        output_hbox.addWidget(self.output_path_edit)
        output_hbox.addWidget(output_browse_btn)
        params_layout.addRow("Output File:", output_hbox)

        layout.addWidget(params_group)

        # ASTRAL-Hybrid 预设选项组
        self.hybrid_preset_group = QGroupBox("ASTRAL-Hybrid Preset Options")
        hybrid_preset_layout = QHBoxLayout()
        self.hybrid_preset_group.setLayout(hybrid_preset_layout)

        self.bayes_radio = QRadioButton("Bayesian (-B)")
        self.bayes_radio.setToolTip("Probability (abayes) support value mode")
        self.bootstrap_radio = QRadioButton("Bootstrap (-S)")
        self.bootstrap_radio.setToolTip("Bootstrap support value mode (default)")
        self.bootstrap_radio.setChecked(True)
        self.lrt_radio = QRadioButton("LRT (-L)")
        self.lrt_radio.setToolTip("Likelihood (alrt) support value mode")

        hybrid_preset_layout.addWidget(self.bayes_radio)
        hybrid_preset_layout.addWidget(self.bootstrap_radio)
        hybrid_preset_layout.addWidget(self.lrt_radio)
        hybrid_preset_layout.addStretch()

        self.hybrid_preset_group.setVisible(False)  # 默认隐藏，只在选择 Hybrid 版本时显示
        layout.addWidget(self.hybrid_preset_group)

        # 预设选项组
        preset_group = QGroupBox("Preset Options")
        preset_layout = QHBoxLayout()
        preset_group.setLayout(preset_layout)

        self.scoring_checkbox = QCheckBox("Scoring Mode (-C/--scoring)")
        self.moreround_checkbox = QCheckBox("More Rounds (-R/--moreround)")
        preset_layout.addWidget(self.scoring_checkbox)
        preset_layout.addWidget(self.moreround_checkbox)
        preset_layout.addStretch()

        layout.addWidget(preset_group)

        # 高级选项按钮
        advanced_btn_layout = QHBoxLayout()
        advanced_btn = QPushButton("Advanced Options...")
        advanced_btn.clicked.connect(self.show_advanced_options)
        advanced_btn_layout.addStretch()
        advanced_btn_layout.addWidget(advanced_btn)
        layout.addLayout(advanced_btn_layout)

        layout.addStretch()

        # 初始化高级参数对话框
        self.advanced_dialog = AdvancedOptionsDialog("astral4", self)
        
        # 如果有预导入的基因树文件，设置到输入文件框中
        if hasattr(self, '_preimport_gene_trees_file') and os.path.exists(self._preimport_gene_trees_file):
            self.file_path_edit.setText(self._preimport_gene_trees_file)
            self.import_file = self._preimport_gene_trees_file
            self.imported_files = [self._preimport_gene_trees_file]
            self.add_console_message(f"Loaded gene trees from: {os.path.basename(self._preimport_gene_trees_file)}", "info")

    def setup_output_tab(self):
        """设置输出标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)

        # 添加输出信息标签
        self.output_info = QLabel("No tree available yet.")
        self.output_info.setAlignment(Qt.AlignCenter)
        self.output_info.setStyleSheet("color: #6c757d; padding: 20px;")
        layout.addWidget(self.output_info)

        # 添加输出预览文本框
        self.output_preview = QTextEdit()
        self.output_preview.setReadOnly(True)
        self.output_preview.setFont(QFont("Courier", 10))
        layout.addWidget(QLabel("Phylogenetic Tree Preview (Newick Format):"))
        layout.addWidget(self.output_preview)

    def on_version_changed(self, index):
        """当ASTRAL版本改变时更新UI"""
        # 更新高级选项对话框
        selected_version = self.ASTRAL_VERSIONS[index]
        self.advanced_dialog = AdvancedOptionsDialog(selected_version, self)

        # 显示/隐藏 Hybrid 预设选项
        self.hybrid_preset_group.setVisible(selected_version == "astral-hybrid")

    def show_advanced_options(self):
        """显示高级选项对话框"""
        result = self.advanced_dialog.exec_()
        if result == QDialog.Accepted:
            self.add_console_message("Advanced options configured", "info")

    def browse_input_file(self):
        """浏览选择输入文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Input File",
            "",
            "Newick Files (*.nwk *.newick *.tre *.tree);;All Files (*)"
        )
        if file_path:
            self.file_path_edit.setText(file_path)
            if not self.import_file:
                self.import_file = file_path
                self.imported_files = [file_path]

    def browse_output_file(self):
        """浏览选择输出文件"""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Output File",
            "",
            "Newick Files (*.nwk *.newick *.tre *.tree);;All Files (*)"
        )
        if file_path:
            self.output_path_edit.setText(file_path)
    
    def import_from_selected_treeset(self):
        """从选中的 TreeSet 导入基因树到 ASTRAL"""
        # 获取 dataset_selection_manager
        dataset_selection_manager = getattr(self, '_dataset_selection_manager', None)
        if not dataset_selection_manager:
            QMessageBox.warning(self, "Warning", "No dataset selection manager available")
            return
        
        # 获取选中的 TreeSet
        from ..platforms.methods.dataset_models import ITEM_TYPE_TREESET, SELECTION_STATE_GREEN
        selected_items = dataset_selection_manager.get_selected_items()
        treeset_items = [item for item in selected_items if item.item_type == ITEM_TYPE_TREESET]
        
        if not treeset_items:
            QMessageBox.warning(self, "Warning", "Please select a TreeSet first (click on it to select)")
            return
        
        # 提取所有树的 Newick 字符串
        trees = []
        for tree_data in treeset_items[0].data.get('data', []):
            if 'content' in tree_data:
                trees.append(tree_data['content'])
            elif 'file_path' in tree_data:
                try:
                    with open(tree_data['file_path'], 'r') as f:
                        trees.append(f.read().strip())
                except:
                    pass
        
        if not trees:
            QMessageBox.warning(self, "Warning", "No trees found in the selected TreeSet")
            return
        
        # 创建临时文件，所有树用换行符分隔
        import tempfile
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.trees', delete=False)
        temp_file.write('\n'.join(trees))
        temp_file.close()
        
        # 设置输入文件
        self.file_path_edit.setText(temp_file.name)
        if not self.import_file:
            self.import_file = temp_file.name
            self.imported_files = [temp_file.name]
        
        self.add_console_message(f"Imported {len(trees)} trees from selected TreeSet", "info")

    def run_analysis(self):
        """运行ASTRAL分析"""
        # 调试信息
        print(f"[AstralPlugin] run_analysis called")
        print(f"[AstralPlugin] tool_path: {self.tool_path}")
        print(f"[AstralPlugin] tool_path exists: {os.path.exists(self.tool_path) if self.tool_path else 'None'}")
        print(f"[AstralPlugin] plugin_path: {self.plugin_path}")

        # 检查工具路径
        if not self.tool_path or not os.path.exists(self.tool_path):
            error_msg = f"Tool not found at: {self.tool_path}\n\n"
            if self.tool_path:
                error_msg += f"Please check if the file exists at this location.\n"
            else:
                error_msg += f"Tool path is not configured.\n"
            error_msg += "Please configure the tool path in config.json"
            QMessageBox.warning(self, "Warning", error_msg)
            return

        input_file = self.file_path_edit.text()
        if not input_file or not os.path.exists(input_file):
            QMessageBox.warning(self, "Warning", "Please select a valid input file!")
            return

        # 获取选择的 ASTRAL 版本
        selected_version = self.ASTRAL_VERSIONS[self.version_combo.currentIndex()]

        # 构建参数列表，第一个参数是 ASTRAL 版本
        params = [selected_version]

        # 线程数
        params.extend(["-t", str(self.thread_spin.value())])

        output_file = self.output_path_edit.text()
        if output_file:
            params.extend(["-o", output_file])

        # 提取支持选项的数值部分
        support_val = self.support_combo.currentText().split(" ")[0]
        params.extend(["-u", support_val])

        # 从高级选项对话框获取参数（只在非默认值时添加）
        if self.advanced_dialog.constraint_edit.text():
            params.extend(["-c", self.advanced_dialog.constraint_edit.text()])

        if self.advanced_dialog.guide_edit.text():
            params.extend(["-g", self.advanced_dialog.guide_edit.text()])

        if self.advanced_dialog.round_spin.value() != 4:  # 默认值为4
            params.extend(["-r", str(self.advanced_dialog.round_spin.value())])

        if self.advanced_dialog.subsample_spin.value() != 4:  # 默认值为4
            params.extend(["-s", str(self.advanced_dialog.subsample_spin.value())])

        if self.advanced_dialog.mapping_edit.text():
            params.extend(["-a", self.advanced_dialog.mapping_edit.text()])

        if self.advanced_dialog.root_edit.text():
            params.extend(["--root", self.advanced_dialog.root_edit.text()])

        # 只在非默认值时添加这些参数
        if self.advanced_dialog.proportion_spin.value() != 0.25:  # 默认值为0.25
            params.extend(["--proportion", f"{self.advanced_dialog.proportion_spin.value():.2f}"])

        if self.advanced_dialog.seed_spin.value() != 233:  # 默认值为233
            params.extend(["--seed", str(self.advanced_dialog.seed_spin.value())])

        # --length 参数只适用于 astral4 和 astral-pro3
        if selected_version in ["astral4", "astral-pro3"]:
            if self.advanced_dialog.length_combo.currentIndex() != 0:  # 默认值为SULength
                params.extend(["--length", self.advanced_dialog.length_combo.currentText()])

        if self.advanced_dialog.lambda_spin.value() != 0.5:  # 默认值为0.5
            params.extend(["-l", f"{self.advanced_dialog.lambda_spin.value():.1f}"])

        # --genelength 参数只适用于 astral4 和 astral-pro3
        if selected_version in ["astral4", "astral-pro3"]:
            if self.advanced_dialog.gene_length_spin.value() != 1000.0:  # 默认值为1000
                params.extend(["--genelength", f"{self.advanced_dialog.gene_length_spin.value():.0f}"])

        if self.advanced_dialog.downweight_spin.value() != 1.0:  # 默认值为1.0
            params.extend(["-w", f"{self.advanced_dialog.downweight_spin.value():.1f}"])

        if self.advanced_dialog.verbose_combo.currentIndex() != 1:  # 默认值为2 (Normal)
            verbose_val = self.advanced_dialog.verbose_combo.currentText().split(" ")[0]
            params.extend(["-v", verbose_val])

        # ASTRAL-Hybrid 特定参数
        selected_version = self.ASTRAL_VERSIONS[self.version_combo.currentIndex()]

        if selected_version == "astral-hybrid":
            # 预设选项
            if self.bayes_radio.isChecked():
                params.append("-B")
            elif self.bootstrap_radio.isChecked():
                params.append("-S")
            elif self.lrt_radio.isChecked():
                params.append("-L")

            # 高级选项
            if self.advanced_dialog.weights_edit.text():
                params.extend(["--treeweights", self.advanced_dialog.weights_edit.text()])

            if self.advanced_dialog.mode_combo.currentIndex() != 0:  # 默认值为1 (Hybrid)
                mode_val = self.advanced_dialog.mode_combo.currentText().split(" ")[0]
                params.extend(["--mode", mode_val])

            if self.advanced_dialog.min_spin.value() != -1.0:  # 默认值为-1 (Auto-detect)
                params.extend(["-n", f"{self.advanced_dialog.min_spin.value():.3f}"])

            if self.advanced_dialog.max_spin.value() != -1.0:  # 默认值为-1 (Auto-detect)
                params.extend(["-x", f"{self.advanced_dialog.max_spin.value():.3f}"])

            if self.advanced_dialog.default_spin.value() != 0.0:  # 默认值为0
                params.extend(["-d", f"{self.advanced_dialog.default_spin.value():.3f}"])

        # ASTRAL-Pro3 特定参数
        if selected_version == "astral-pro3":
            if self.advanced_dialog.exit_combo.currentIndex() != 0:  # 默认值为0
                exit_val = self.advanced_dialog.exit_combo.currentText().split(" ")[0]
                params.extend(["-e", exit_val])

        # 预设选项
        if self.scoring_checkbox.isChecked():
            params.append("-C")

        if self.moreround_checkbox.isChecked():
            params.append("-R")

        # 添加输入文件到最后
        params.append(input_file)

        # 创建并启动线程
        self.analysis_thread = AstralThread(
            self.tool_path,
            [input_file],
            params,
            self.imported_files,
            output_file if output_file else None
        )

        self.analysis_thread.progress.connect(self.update_progress)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.error.connect(self.handle_error)
        self.analysis_thread.finished.connect(self.handle_finished)

        self.analysis_thread.start()
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)

    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'analysis_thread') and self.analysis_thread:
            self.analysis_thread.stop()
            self.is_running = False
            self.run_button.setEnabled(True)
            self.stop_button.setEnabled(False)
            self.progress_bar.setVisible(False)
            self.add_console_message("Analysis stopped by user", "info")

    def handle_error(self, error_msg):
        """处理错误"""
        self.add_console_message(error_msg, "error")
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

    def display_results(self, output_files):
        """显示结果到输出标签页"""
        if not output_files:
            QMessageBox.warning(self, "Error", "No output files generated")
            return

        treefile = None
        for file in output_files:
            # 支持多种树文件扩展名
            if file.endswith('.nwk') or file.endswith('.tre') or file.endswith('.tree') or file.endswith('.newick'):
                treefile = file
                break

        if treefile:
            try:
                with open(treefile, 'r', encoding='utf-8') as f:
                    tree_content = f.read().strip()

                if not tree_content:
                    QMessageBox.warning(self, "Error", "Tree file is empty")
                    return

                # 导入并创建IcyTree插件
                from .icytree import IcyTreePlugin
                plugin_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '')
                icytree_plugin = IcyTreePlugin(plugin_path=plugin_path)

                # 连接IcyTree状态变化信号
                icytree_plugin.status_changed.connect(lambda status: self._on_icytree_status_changed(status, icytree_plugin))

                # 设置newick字符串
                icytree_plugin.set_newick_string(tree_content)

                # 清理output_tab布局中的现有组件（保留output_info）
                output_layout = self.output_tab.layout()
                if output_layout:
                    for i in reversed(range(output_layout.count())):
                        widget = output_layout.itemAt(i).widget()
                        if widget and widget != self.output_info:
                            widget.setParent(None)

                # 添加IcyTree插件到输出标签页
                output_layout.addWidget(icytree_plugin)

            except ImportError:
                QMessageBox.warning(self, "Error", "IcyTree plugin not available")
                # 如果IcyTree不可用，显示原始文本
                if hasattr(self, 'output_preview'):
                    with open(treefile, 'r', encoding='utf-8') as f:
                        self.output_preview.setText(f.read())
                    if hasattr(self, 'output_info'):
                        self.output_info.setParent(None)
                        self.output_info = None
            except Exception as e:
                error_msg = f"Error processing tree file: {str(e)}"
                QMessageBox.warning(self, "Error", error_msg)
                self.add_console_message(error_msg, "error")
        else:
            QMessageBox.warning(self, "Error", f"No tree file found. Generated {len(output_files)} file(s).")

    def _on_icytree_status_changed(self, status, icytree_plugin):
        """处理IcyTree状态变化"""
        if status == "Tree loaded to IcyTree":
            # 当树加载到IcyTree时，移除"No tree available yet."标签
            if hasattr(self, 'output_info') and self.output_info:
                self.output_info.setParent(None)
                self.output_info = None

    def handle_finished(self, output_files, reports):
        """处理分析完成事件"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

        # 更新输出预览
        self.reports = output_files
        self.current_output_files = output_files

        # 显示结果并切换到输出标签页
        self.display_results(output_files)
        if hasattr(self, 'tab_widget'):
            self.tab_widget.setCurrentIndex(1)

        self.add_console_message(f"ASTRAL analysis completed! Output saved to: {', '.join(output_files)}", "info")
        QMessageBox.information(self, "Completed", f"ASTRAL analysis completed! Output saved to: {', '.join(output_files)}")

    def update_progress(self, progress_text):
        """更新进度"""
        self.progress_bar.setFormat(progress_text)
        self.progress_bar.setValue(0)  # 无限进度条


class AstralPluginEntry:
    def run(self, import_from=None, import_data=None):
        return AstralPlugin(import_from=import_from, import_data=import_data)
