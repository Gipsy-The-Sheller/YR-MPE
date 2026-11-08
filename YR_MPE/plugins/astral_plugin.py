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

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.nexus_tree import read_tree_file
import subprocess


def normalize_weighted_support_labels(input_path, output_path):
    """Collapse composite gene-tree support labels to one numeric value.

    IQ-TREE can write ``SH-aLRT/aBayes/UFBoot`` labels such as ``/1/100``.
    ASTER's wASTRAL and ASTRAL-Hybrid parse node labels with ``std::stod`` and
    abort on the leading slash (``std::invalid_argument: stod``). Keep the
    first numeric field (aBayes when present, otherwise aLRT/UFBoot) and drop
    the rest. Returns ``(trees, changed)``.
    """
    from Bio import Phylo

    trees = list(Phylo.parse(input_path, "newick"))
    if not trees:
        raise ValueError("no Newick trees found")

    changed = 0
    for tree in trees:
        for clade in tree.get_nonterminals():
            label = (clade.name or "").strip()
            if not label:
                continue
            value = None
            for token in label.split("/"):
                token = token.strip()
                if not token:
                    continue
                try:
                    float(token)
                    value = token
                    break
                except ValueError:
                    continue
            if value != clade.name:
                clade.name = value
                changed += 1

    with open(output_path, "w", encoding="utf-8") as handle:
        Phylo.write(trees, handle, "newick")
    return len(trees), changed


class AstralThread(BaseProcessThread):
    """ASTRAL分析线程类"""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters, imported_files=None,
                 output_file=None, workdir=None):
        super().__init__(tool_path, input_files, parameters, imported_files,
                         workdir=workdir)
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

                # ASTER writes the tree to -o; if the dialog did not specify
                # one, insert a temporary -o before the input file.
                output_file = self.output_file
                cmd = [self.tool_path, *self.parameters]
                if not output_file:
                    output_file = self.create_temp_file(suffix='.tre')
                    cmd = cmd[:-1] + ["-o", output_file, cmd[-1]]

                result = self.execute_command(
                    cmd, cwd=self.workdir if self.workdir else None
                )

                if result.returncode != 0:
                    detail = (result.stderr or "").strip()
                    raise RuntimeError(
                        "ASTRAL failed%s" % (": " + detail if detail else "")
                    )

                # Prefer the file ASTER wrote. stdout may contain only logs
                # when -o is used, so never overwrite the tree with it.
                if not (os.path.exists(output_file) and os.path.getsize(output_file) > 0):
                    with open(output_file, 'w', encoding='utf-8') as f:
                        f.write(result.stdout or "")
                if not (os.path.exists(output_file) and os.path.getsize(output_file) > 0):
                    raise RuntimeError("ASTRAL produced no tree output.")

                self._ensure_terminal_lengths(output_file)
                output_files.append(output_file)

            self.progress.emit("ASTRAL analysis completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])

        except Exception as e:
            self.error.emit(f"Analysis exception: {str(e)}")

    def _ensure_terminal_lengths(self, path):
        """Add dummy terminal lengths for viewers that reject missing tips.

        ASTRAL-IV/CASTER SU output already has terminal lengths. wASTRAL
        (coalescent units) intentionally leaves them empty; ape and some
        viewers fail on that. The dummy is tiny and display-only.
        """
        try:
            from Bio import Phylo
            tree = Phylo.read(path, "newick")
            changed = False
            for tip in tree.get_terminals():
                if tip.branch_length is None:
                    tip.branch_length = 1e-5
                    changed = True
            if changed:
                Phylo.write(tree, path, "newick")
                self.console_output.emit(
                    "Added display-only terminal branch lengths (1e-5) "
                    "because this ASTER mode leaves them empty.", "info")
        except Exception as exc:
            self.console_output.emit(
                "Could not add terminal branch lengths: %s" % exc, "warning")


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
    ASTRAL_VERSIONS = ["astral4", "wastral", "astral-hybrid"]

    def __init__(self, import_from=None, import_data=None, **kwargs):
        super().__init__(import_from, import_data, **kwargs)
        
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

        # setup_input_tab() runs during super().__init__, before the preimport
        # attribute above exists. Apply it again once the UI is ready.
        self._apply_preimport_gene_trees()

    def _apply_preimport_gene_trees(self):
        path = getattr(self, "_preimport_gene_trees_file", "")
        if not path or not os.path.exists(path):
            return
        if not hasattr(self, "file_path_edit") or self.file_path_edit is None:
            return
        self.file_path_edit.setText(path)
        self.import_file = path
        self.imported_files = [path]
        self.add_console_message(
            "Loaded gene trees from: %s" % os.path.basename(path), "info"
        )

    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "ASTRAL"
        self.tool_name = "astral"
        self.citation = []  # citation text lives in YR_MPE/citations/*.bib
        self.citation_keys = ["zhang2025aster", "zhang2022wastral"]
        self.input_types = {
            "Newick": ["nwk", "newick", "tre", "tree"]
        }
        self.output_types = {
            "Newick": ".tre"
        }
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

    def control_bindings(self):
        """GUI controls mapped to the ASTER command line."""
        bindings = [
            {"id": "input", "widget": "file_path_edit",
             "label": "Input gene trees", "kind": "path", "cli": [],
             "notes": "Positional input, appended after all options."},
            {"id": "version", "widget": "version_combo",
             "label": "ASTRAL version", "kind": "choice",
             "choices": {
                 "ASTRAL-IV (Standard species tree inference)": "astral4",
                 "Weighted ASTRAL (wastral)": "wastral",
                 "ASTRAL-Hybrid (Hybrid weighting)": "astral-hybrid",
             },
             "cli": ["{value}"],
             "notes": "First argv token; the wrapper selects the matching "
                      "ASTER executable."},
            {"id": "threads", "widget": "thread_spin",
             "label": "Threads", "kind": "int",
             "cli": ["-t", "{value}"]},
            {"id": "output", "widget": "output_path_edit",
             "label": "Output tree", "kind": "path",
             "cli": ["-o", "{value}"],
             "notes": "Optional; if empty the plugin writes a temporary .tre."},
            {"id": "support", "widget": "support_combo",
             "label": "Support option", "kind": "choice",
             "choices": {
                 "0 - No branch or support": "0",
                 "1 - Length and support only": "1",
                 "2 - Detailed": "2",
                 "3 - freqQuad.csv": "3",
             },
             "cli": ["-u", "{value}"]},
            {"id": "constraint", "widget": "advanced_dialog.constraint_edit",
             "label": "Constraint tree", "kind": "path",
             "cli": ["-c", "{value}"]},
            {"id": "guide", "widget": "advanced_dialog.guide_edit",
             "label": "Guide trees", "kind": "path",
             "cli": ["-g", "{value}"]},
            {"id": "rounds", "widget": "advanced_dialog.round_spin",
             "label": "Rounds", "kind": "int",
             "cli": ["-r", "{value}"], "omit_values": [4]},
            {"id": "subsample", "widget": "advanced_dialog.subsample_spin",
             "label": "Subsample rounds", "kind": "int",
             "cli": ["-s", "{value}"], "omit_values": [4]},
            {"id": "mapping", "widget": "advanced_dialog.mapping_edit",
             "label": "Gene-to-taxon map", "kind": "path",
             "cli": ["-a", "{value}"]},
            {"id": "root", "widget": "advanced_dialog.root_edit",
             "label": "Root taxon", "kind": "text",
             "cli": ["--root", "{value}"]},
            {"id": "proportion", "widget": "advanced_dialog.proportion_spin",
             "label": "Subsample proportion", "kind": "float",
             "cli": ["--proportion", "{value}"], "omit_values": [0.25]},
            {"id": "seed", "widget": "advanced_dialog.seed_spin",
             "label": "Seed", "kind": "int",
             "cli": ["--seed", "{value}"], "omit_values": [233]},
            {"id": "lambda", "widget": "advanced_dialog.lambda_spin",
             "label": "Yule lambda", "kind": "float",
             "cli": ["-l", "{value}"], "omit_values": [0.5]},
            {"id": "downweight", "widget": "advanced_dialog.downweight_spin",
             "label": "Downweight repeats", "kind": "float",
             "cli": ["-w", "{value}"], "omit_values": [1.0]},
            {"id": "verbose", "widget": "advanced_dialog.verbose_combo",
             "label": "Verbose level", "kind": "choice",
             "choices": {"1 - Minimum": "1", "2 - Normal": "2"},
             "cli": ["-v", "{value}"], "omit_values": ["2"]},
            {"id": "scoring", "widget": "scoring_checkbox",
             "label": "Scoring only", "kind": "bool", "cli": ["-C"]},
            {"id": "moreround", "widget": "moreround_checkbox",
             "label": "More rounds", "kind": "bool", "cli": ["-R"]},
            {"id": "length", "widget": "advanced_dialog.length_combo",
             "label": "Branch lengths", "kind": "choice",
             "choices": {"SULength": "SULength", "CULength": "CULength"},
             "cli": ["--length", "{value}"], "omit_values": ["SULength"],
             "enabled_when": {"widget": "version_combo",
                              "equals": "ASTRAL-IV (Standard species tree inference)"},
             "notes": "ASTRAL-IV only."},
            {"id": "gene_length", "widget": "advanced_dialog.gene_length_spin",
             "label": "Gene length", "kind": "float",
             "cli": ["--genelength", "{value}"], "omit_values": [1000.0],
             "enabled_when": {"widget": "version_combo",
                              "equals": "ASTRAL-IV (Standard species tree inference)"},
             "notes": "ASTRAL-IV only."},
        ]

        # Hybrid-only widgets exist only when the advanced dialog was built
        # for astral-hybrid.
        if hasattr(self.advanced_dialog, "mode_combo"):
            bindings.extend([
                {"id": "hybrid_bayes", "widget": "bayes_radio",
                 "label": "Bayesian support preset", "kind": "bool",
                 "cli": ["-B"]},
                {"id": "hybrid_bootstrap", "widget": "bootstrap_radio",
                 "label": "Bootstrap support preset", "kind": "bool",
                 "cli": ["-S"]},
                {"id": "hybrid_lrt", "widget": "lrt_radio",
                 "label": "LRT support preset", "kind": "bool",
                 "cli": ["-L"]},
                {"id": "treeweights", "widget": "advanced_dialog.weights_edit",
                 "label": "Tree weights", "kind": "path",
                 "cli": ["--treeweights", "{value}"]},
                {"id": "hybrid_mode", "widget": "advanced_dialog.mode_combo",
                 "label": "Weighting mode", "kind": "choice",
                 "choices": {
                     "1 - Hybrid weighting": "1",
                     "2 - Support only": "2",
                     "3 - Length only": "3",
                     "4 - Unweighted": "4",
                 },
                 "cli": ["--mode", "{value}"], "omit_values": ["1"]},
                {"id": "min_support", "widget": "advanced_dialog.min_spin",
                 "label": "Minimum support", "kind": "float",
                 "cli": ["-n", "{value}"], "omit_values": [-1.0]},
                {"id": "max_support", "widget": "advanced_dialog.max_spin",
                 "label": "Maximum support", "kind": "float",
                 "cli": ["-x", "{value}"], "omit_values": [-1.0]},
                {"id": "default_support", "widget": "advanced_dialog.default_spin",
                 "label": "Default support", "kind": "float",
                 "cli": ["-d", "{value}"], "omit_values": [0.0]},
            ])
        return bindings

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
        self.import_to_platform_btn = QPushButton("Send to workspace")
        self.import_to_platform_btn.setVisible(False)
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        layout.addWidget(self.import_to_platform_btn)

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
        """Import selected phylogenies and/or TreeSets into the tree input."""
        dataset_selection_manager = getattr(self, '_dataset_selection_manager', None)
        if not dataset_selection_manager:
            QMessageBox.warning(self, "Warning", "No dataset selection manager available")
            return

        from ..platforms.methods.dataset_models import (
            ITEM_TYPE_PHYLOGENY, ITEM_TYPE_TREESET,
        )
        selected_items = dataset_selection_manager.get_selected_items()
        trees = []
        for item in selected_items:
            if item.item_type == ITEM_TYPE_TREESET:
                for tree_data in (item.data or {}).get('data', []):
                    if isinstance(tree_data, str):
                        trees.append(tree_data)
                    elif 'content' in tree_data:
                        trees.append(tree_data['content'])
                    elif tree_data.get('file_path'):
                        try:
                            with open(tree_data['file_path'], 'r', encoding='utf-8') as f:
                                trees.append(f.read().strip())
                        except OSError:
                            pass
            elif item.item_type == ITEM_TYPE_PHYLOGENY:
                data = item.data or {}
                content = ""
                if isinstance(data, dict):
                    content = data.get('content') or data.get('newick') or data.get('tree') or ""
                elif isinstance(data, str):
                    content = data
                if not content and item.file_path:
                    try:
                        with open(item.file_path, 'r', encoding='utf-8') as f:
                            content = f.read().strip()
                    except OSError:
                        content = ""
                if content:
                    trees.extend(line.strip() for line in content.splitlines()
                                 if line.strip())

        trees = [tree for tree in trees if tree]
        if not trees:
            QMessageBox.warning(
                self, "Warning",
                "Please select a Phylogeny or TreeSet first."
            )
            return

        import tempfile
        handle = tempfile.NamedTemporaryFile(
            mode='w', suffix='.trees', delete=False, encoding='utf-8'
        )
        handle.write('\n'.join(trees))
        handle.close()

        self.file_path_edit.setText(handle.name)
        self.import_file = handle.name
        self.imported_files = [handle.name]
        self.add_console_message(
            "Imported %d tree(s) from the current selection." % len(trees),
            "info",
        )

    def run_analysis(self):
        """运行ASTRAL分析"""
        # 调试信息

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

        # --length 参数只适用于 astral4
        if selected_version == "astral4":
            if self.advanced_dialog.length_combo.currentIndex() != 0:  # 默认值为SULength
                params.extend(["--length", self.advanced_dialog.length_combo.currentText()])

        if self.advanced_dialog.lambda_spin.value() != 0.5:  # 默认值为0.5
            params.extend(["-l", f"{self.advanced_dialog.lambda_spin.value():.1f}"])

        # --genelength 参数只适用于 astral4
        if selected_version == "astral4":
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

        # 预设选项
        if self.scoring_checkbox.isChecked():
            params.append("-C")

        if self.moreround_checkbox.isChecked():
            params.append("-R")

        # Weighted modes use ASTER's std::stod parser, which aborts on
        # composite labels like "/1/100" from IQ-TREE. Normalize first.
        if selected_version in ("wastral", "astral-hybrid"):
            try:
                sanitized = self.create_temp_file(
                    suffix="_%s.trees" % selected_version
                )
                tree_count, changed = normalize_weighted_support_labels(
                    input_file, sanitized
                )
                if changed:
                    input_file = sanitized
                    self.add_console_message(
                        "Normalized %d composite support label(s) across %d "
                        "gene tree(s); kept one numeric value per node."
                        % (changed, tree_count),
                        "info",
                    )
            except Exception as exc:
                self.add_console_message(
                    "Could not normalize support labels (%s); passing the "
                    "input unchanged." % exc,
                    "warning",
                )

        # 添加输入文件到最后
        params.append(input_file)

        # 创建并启动线程
        self.analysis_thread = AstralThread(
            self.tool_path,
            [input_file],
            params,
            self.imported_files,
            output_file if output_file else None,
            workdir=self.workdir,
        )

        self.analysis_thread.progress.connect(self.update_progress)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.error.connect(self.handle_error)
        self.analysis_thread.finished.connect(self.handle_finished)

        self.analysis_thread.start()
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)

    def _on_result(self, payload):
        self._result_payload = payload

    def stop_analysis(self):
        """停止分析"""
        thread = getattr(self, 'analysis_thread', None)
        if thread is not None and thread.isRunning():
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
            self.add_console_message("Analysis stopped by user", "info")
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

    def handle_error(self, error_msg):
        """处理错误"""
        self.add_console_message(error_msg, "error")
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

    def display_results(self, output_files):
        """Show the inferred species tree without embedding a web view."""
        treefile = None
        for path in output_files:
            if path.endswith((".nwk", ".tre", ".tree", ".newick")):
                treefile = path
                break
        if treefile is None:
            self.output_info.setText(
                "No tree file found. Generated %d file(s)." % len(output_files)
            )
            self.output_preview.clear()
            return
        try:
            trees = read_tree_file(treefile)
        except OSError as exc:
            self.output_info.setText("Cannot read %s: %s" % (treefile, exc))
            self.output_preview.clear()
            return
        if not trees:
            self.output_info.setText(
                "No tree could be parsed from %s." % os.path.basename(treefile)
            )
            self.output_preview.clear()
            return
        self.output_info.setText(
            "Species tree: %s · %d output file(s)"
            % (os.path.basename(treefile), len(output_files))
        )
        self.output_preview.setPlainText("\n".join(trees))

    def _build_result_spec(self):
        files = list(getattr(self, "current_output_files", []) or [])
        if not files:
            raise ValueError("No ASTRAL output is available.")
        items = []
        artifacts = []
        for path in files:
            name = os.path.basename(path)
            if name.endswith((".nwk", ".tre", ".tree", ".newick")):
                trees = read_tree_file(path)
                if not trees:
                    raise ValueError("No tree could be parsed from %s." % name)
                if len(trees) > 1:
                    items.append(DatasetItemSpec(
                        loci_name=os.path.splitext(name)[0],
                        item_type="treeset",
                        file_path=str(path),
                        format="newick",
                        data={"data": [{"content": tree} for tree in trees]},
                        metadata={"source": self.plugin_name, "tool": "ASTRAL"},
                    ))
                else:
                    items.append(DatasetItemSpec(
                        loci_name=os.path.splitext(name)[0],
                        item_type="phylogeny",
                        file_path=str(path),
                        format="newick",
                        data={"newick": trees[0]},
                        metadata={"source": self.plugin_name, "tool": "ASTRAL"},
                    ))
            else:
                artifacts.append(str(path))
        if not items:
            raise ValueError("ASTRAL produced no tree item.")
        dataset = DatasetSpec(
            name="ASTRAL species tree",
            items=items,
            description="%d tree item(s)" % len(items),
            source="astral",
            summary={"items": len(items)},
            metadata={
                "input": self.file_path_edit.text(),
                "version": self.ASTRAL_VERSIONS[self.version_combo.currentIndex()],
            },
        )
        return ResultSpec(source="astral", datasets=[dataset], artifacts=artifacts)

    def import_to_platform(self):
        if not getattr(self, "current_output_files", None):
            QMessageBox.warning(self, "No result", "Run ASTRAL first.")
            return
        if not self.result_has_consumer():
            QMessageBox.warning(
                self, "No workspace",
                "This window is not connected to a YR-MPE workspace."
            )
            return
        try:
            result = self._build_result_spec()
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        if self.export_result(result):
            self.import_to_platform_btn.setEnabled(False)
            self.import_to_platform_btn.setText("Sent to workspace")
            QMessageBox.information(
                self, "Success", "ASTRAL species tree sent to the current workspace."
            )

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

        if (self.import_from in ["YR_MPEA", "DATASET_MANAGER"]
                or self.result_has_consumer()):
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)

        self.add_console_message(f"ASTRAL analysis completed! Output saved to: {', '.join(output_files)}", "info")
        QMessageBox.information(self, "Completed", f"ASTRAL analysis completed! Output saved to: {', '.join(output_files)}")

    def update_progress(self, progress_text):
        """更新进度"""
        self.progress_bar.setFormat(progress_text)
        self.progress_bar.setValue(0)  # 无限进度条


class AstralPluginEntry:
    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return AstralPlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )
