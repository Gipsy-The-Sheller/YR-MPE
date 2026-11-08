from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit,
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, QDoubleSpinBox)
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


class CasterSiteThread(BaseProcessThread):
    """CASTER-site分析线程类"""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters, imported_files=None,
                 output_file=None, workdir=None):
        super().__init__(tool_path, input_files, parameters, imported_files,
                         workdir=workdir)
        self.output_file = output_file
    
    def get_tool_name(self):
        """返回工具名称"""
        return "CASTER-site"
        
    def execute_commands(self):
        """执行CASTER-site分析命令"""
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
                result = self.execute_command(
                    cmd, cwd=self.workdir if self.workdir else None
                )

                if result.returncode != 0:
                    detail = (result.stderr or "").strip()
                    raise RuntimeError(
                        "CASTER-site failed%s" % (": " + detail if detail else "")
                    )

                # ASTER writes the tree to stdout.
                with open(output_file, 'w', encoding='utf-8') as f:
                    f.write(result.stdout or "")
                if not (os.path.exists(output_file) and os.path.getsize(output_file) > 0):
                    raise RuntimeError("CASTER-site produced no tree output.")

                output_files.append(output_file)

            self.progress.emit("CASTER-site analysis completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Analysis exception: {str(e)}")


class AdvancedOptionsDialog(QDialog):
    """高级选项对话框"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Advanced Options")
        self.resize(500, 600)
        
        layout = QVBoxLayout()
        
        # 创建滚动区域以容纳所有高级参数
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
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
        
        # Chunk size
        chunk_group = QGroupBox("Chunk Size")
        chunk_layout = QHBoxLayout()
        self.chunk_spin = QSpinBox()
        self.chunk_spin.setRange(1000, 100000)
        self.chunk_spin.setValue(10000)
        chunk_layout.addWidget(self.chunk_spin)
        chunk_desc = QLabel("The chunk size of each local region for parameter estimation")
        chunk_desc.setAlignment(Qt.AlignLeft)
        chunk_layout.addWidget(chunk_desc)
        chunk_group.setLayout(chunk_layout)
        scroll_layout.addWidget(chunk_group)
        
        # Ambiguity handling
        ambiguity_group = QGroupBox("Ambiguity Handling")
        ambiguity_layout = QHBoxLayout()
        self.ambiguity_combo = QComboBox()
        self.ambiguity_combo.addItems(["0", "1"])  # 0: treat as N, 1: treat as diploid
        self.ambiguity_combo.setEditable(False)
        ambiguity_layout.addWidget(self.ambiguity_combo)
        ambiguity_desc = QLabel("0: ambiguity codes as N, 1: ambiguity codes as diploid unphased sites")
        ambiguity_desc.setAlignment(Qt.AlignLeft)
        ambiguity_layout.addWidget(ambiguity_desc)
        ambiguity_group.setLayout(ambiguity_layout)
        scroll_layout.addWidget(ambiguity_group)
        
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


class CasterSitePlugin(BasePlugin):
    """CASTER-site插件：用于物种树估计的coalescence感知的基于比对的工具"""
    
    def __init__(self, import_from=None, import_data=None, **kwargs):
        super().__init__(import_from, import_data, **kwargs)

    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "CASTER-site"
        self.tool_name = "caster-site"
        self.citation = []  # citation text lives in YR_MPE/citations/*.bib
        self.citation_keys = ["zhang2025aster", "zhang2025caster"]
        self.input_types = {
            "FASTA": ["fas", "fna", "fa", "fasta"],
            "List": ["txt"],
            "PHYLIP": ["phy", "phylip"]
        }
        self.output_types = {
            "Newick": ".tre"
        }
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')


    def control_bindings(self):
        """GUI controls mapped to the CASTER-site command line."""
        return [
            {"id": "input", "widget": "file_path_edit",
             "label": "Input alignment", "kind": "path", "cli": [],
             "notes": "Positional input, appended last."},
            {"id": "format", "widget": "format_combo",
             "label": "Input format", "kind": "choice",
             "choices": {"auto": "", "fasta": "fasta", "list": "list",
                         "phylip": "phylip"},
             "cli": ["--format", "{value}"],
             "notes": "auto lets CASTER detect the format."},
            {"id": "threads", "widget": "thread_spin",
             "label": "Threads", "kind": "int",
             "cli": ["--thread", "{value}"]},
            {"id": "output", "widget": "output_path_edit",
             "label": "Output tree", "kind": "path",
             "cli": ["--output", "{value}"],
             "notes": "Optional; a temporary tree is used if empty."},
            {"id": "support", "widget": "support_combo",
             "label": "Support option", "kind": "choice",
             "choices": {
                 "0 - No branch or support": "0",
                 "1 - Length and support only": "1",
                 "2 - Detailed": "2",
             },
             "cli": ["--support", "{value}"]},
            {"id": "subsample", "widget": "advanced_dialog.subsample_spin",
             "label": "Subsample rounds", "kind": "int",
             "cli": ["--subsample", "{value}"], "omit_values": [4]},
            {"id": "constraint", "widget": "advanced_dialog.constraint_edit",
             "label": "Constraint tree", "kind": "path",
             "cli": ["--constraint", "{value}"]},
            {"id": "guide", "widget": "advanced_dialog.guide_edit",
             "label": "Guide trees", "kind": "path",
             "cli": ["--guide", "{value}"]},
            {"id": "rounds", "widget": "advanced_dialog.round_spin",
             "label": "Rounds", "kind": "int",
             "cli": ["--round", "{value}"], "omit_values": [4]},
            {"id": "mapping", "widget": "advanced_dialog.mapping_edit",
             "label": "Gene-to-taxon map", "kind": "path",
             "cli": ["--mapping", "{value}"]},
            {"id": "root", "widget": "advanced_dialog.root_edit",
             "label": "Root taxon", "kind": "text",
             "cli": ["--root", "{value}"]},
            {"id": "proportion", "widget": "advanced_dialog.proportion_spin",
             "label": "Subsample proportion", "kind": "float",
             "cli": ["--proportion", "{value}"], "omit_values": [0.25]},
            {"id": "seed", "widget": "advanced_dialog.seed_spin",
             "label": "Seed", "kind": "int",
             "cli": ["--seed", "{value}"], "omit_values": [233]},
            {"id": "chunk", "widget": "advanced_dialog.chunk_spin",
             "label": "Chunk size", "kind": "int",
             "cli": ["--chunk", "{value}"], "omit_values": [10000]},
            {"id": "ambiguity", "widget": "advanced_dialog.ambiguity_combo",
             "label": "Ambiguity handling", "kind": "choice",
             "choices": {"0": "0", "1": "1"},
             "cli": ["--ambiguity", "{value}"], "omit_values": ["0"],
             "notes": "0 treats ambiguity codes as N; 1 treats them as diploid."},
            {"id": "length", "widget": "advanced_dialog.length_combo",
             "label": "Branch lengths", "kind": "choice",
             "choices": {"SULength": "SULength", "CULength": "CULength"},
             "cli": ["--length", "{value}"], "omit_values": ["SULength"]},
            {"id": "scoring", "widget": "scoring_checkbox",
             "label": "Scoring only", "kind": "bool", "cli": ["--scoring"]},
            {"id": "moreround", "widget": "moreround_checkbox",
             "label": "More rounds", "kind": "bool", "cli": ["--moreround"]},
        ]

    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # 文件选择组
        file_group = QGroupBox("Input File")
        file_layout = QFormLayout()
        file_group.setLayout(file_layout)
        
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select input file...")
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_input_file)
        file_hbox = QHBoxLayout()
        file_hbox.addWidget(self.file_path_edit)
        file_hbox.addWidget(browse_btn)
        file_layout.addRow("Input File:", file_hbox)
        
        # 输入格式选择
        self.format_combo = QComboBox()
        self.format_combo.addItems(["auto", "fasta", "list", "phylip"])
        file_layout.addRow("Format:", self.format_combo)
        
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
            "2 - Detailed"
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
        self.advanced_dialog = AdvancedOptionsDialog(self)

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
            "Supported Files (*.fas *.fna *.fa *.fasta *.txt *.phy *.phylip);;FASTA Files (*.fas *.fna *.fa *.fasta);;List Files (*.txt);;PHYLIP Files (*.phy *.phylip);;All Files (*)"
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

    def run_analysis(self):
        """运行CASTER-site分析"""
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

        # 构建参数列表
        params = []

        # 基本参数
        format_val = self.format_combo.currentText()
        if format_val != "auto":
            params.extend(["--format", format_val])

        params.extend(["--thread", str(self.thread_spin.value())])

        output_file = self.output_path_edit.text()
        if output_file:
            params.extend(["--output", output_file])

        # 提取支持选项的数值部分
        support_val = self.support_combo.currentText().split(" ")[0]
        params.extend(["--support", support_val])

        # 从高级选项对话框获取参数（只在非默认值时添加）
        if self.advanced_dialog.subsample_spin.value() != 4:  # 默认值为4
            params.extend(["--subsample", str(self.advanced_dialog.subsample_spin.value())])

        if self.advanced_dialog.constraint_edit.text():
            params.extend(["--constraint", self.advanced_dialog.constraint_edit.text()])

        if self.advanced_dialog.guide_edit.text():
            params.extend(["--guide", self.advanced_dialog.guide_edit.text()])

        if self.advanced_dialog.round_spin.value() != 4:  # 默认值为4
            params.extend(["--round", str(self.advanced_dialog.round_spin.value())])

        if self.advanced_dialog.mapping_edit.text():
            params.extend(["--mapping", self.advanced_dialog.mapping_edit.text()])

        if self.advanced_dialog.root_edit.text():
            params.extend(["--root", self.advanced_dialog.root_edit.text()])

        # 只在非默认值时添加这些参数
        if self.advanced_dialog.proportion_spin.value() != 0.25:  # 默认值为0.25
            params.extend(["--proportion", f"{self.advanced_dialog.proportion_spin.value():.2f}"])

        if self.advanced_dialog.seed_spin.value() != 233:  # 默认值为233
            params.extend(["--seed", str(self.advanced_dialog.seed_spin.value())])

        if self.advanced_dialog.chunk_spin.value() != 10000:  # 默认值为10000
            params.extend(["--chunk", str(self.advanced_dialog.chunk_spin.value())])

        if self.advanced_dialog.ambiguity_combo.currentIndex() != 0:  # 默认值为0
            params.extend(["--ambiguity", str(self.advanced_dialog.ambiguity_combo.currentText())])

        if self.advanced_dialog.length_combo.currentIndex() != 0:  # 默认值为SULength
            params.extend(["--length", self.advanced_dialog.length_combo.currentText()])

        # 预设选项
        if self.scoring_checkbox.isChecked():
            params.append("--scoring")

        if self.moreround_checkbox.isChecked():
            params.append("--moreround")

        # 添加输入文件到最后
        params.append(input_file)

        # 创建并启动线程
        self.analysis_thread = CasterSiteThread(
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
            raise ValueError("No CASTER-site output is available.")
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
                        metadata={"source": self.plugin_name, "tool": "CASTER-site"},
                    ))
                else:
                    items.append(DatasetItemSpec(
                        loci_name=os.path.splitext(name)[0],
                        item_type="phylogeny",
                        file_path=str(path),
                        format="newick",
                        data={"newick": trees[0]},
                        metadata={"source": self.plugin_name, "tool": "CASTER-site"},
                    ))
            else:
                artifacts.append(str(path))
        if not items:
            raise ValueError("CASTER-site produced no tree item.")
        dataset = DatasetSpec(
            name="CASTER-site species tree",
            items=items,
            description="%d tree item(s)" % len(items),
            source="caster_site",
            summary={"items": len(items)},
            metadata={"input": self.file_path_edit.text()},
        )
        return ResultSpec(source="caster_site", datasets=[dataset], artifacts=artifacts)

    def import_to_platform(self):
        if not getattr(self, "current_output_files", None):
            QMessageBox.warning(self, "No result", "Run CASTER-site first.")
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
                self, "Success", "CASTER-site species tree sent to the current workspace."
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

        self.add_console_message(f"CASTER-site analysis completed! Output saved to: {', '.join(output_files)}", "info")
        QMessageBox.information(self, "Completed", f"CASTER-site analysis completed! Output saved to: {', '.join(output_files)}")

    def update_progress(self, progress_text):
        """更新进度"""
        self.progress_bar.setFormat(progress_text)
        self.progress_bar.setValue(0)  # 无限进度条

class CasterSitePluginEntry:
    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return CasterSitePlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )