from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit,
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, QDoubleSpinBox)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon, QFont
import tempfile
import os
from typing import List, Optional

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class CasterSiteThread(BaseProcessThread):
    """CASTER-site分析线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None, output_file=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
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
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"CASTER-site execution failed: {result.stderr}")
                    return
                
                # 将stdout写入输出文件
                with open(output_file, 'w', encoding='utf-8') as f:
                    f.write(result.stdout)
                
                output_files.append(output_file)
            
            self.progress.emit("CASTER-site analysis completed")
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
    
    def __init__(self, import_from=None, import_data=None):
        super().__init__(import_from, import_data)
    
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "CASTER-site"
        self.tool_name = "caster-site"
        self.citation = [
            "Chao Zhang, Rasmus Nielsen, Siavash Mirarab, ASTER: A Package for Large-Scale Phylogenomic Reconstructions, Molecular Biology and Evolution, Volume 42, Issue 8, August 2025, msaf172, https://doi.org/10.1093/molbev/msaf172",
            "Chao Zhang et al., CASTER: Direct species tree inference from whole-genome alignments. Science 387, eadk9688 (2025). DOI: 10.1126/science.adk9688"
        ]
        self.input_types = {
            "FASTA": ["fas", "fna", "fa", "fasta"],
            "List": ["txt"],
            "PHYLIP": ["phy", "phylip"]
        }
        self.output_types = {
            "Newick": ".tre"
        }
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')


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
        print(f"[CasterSitePlugin] run_analysis called")
        print(f"[CasterSitePlugin] tool_path: {self.tool_path}")
        print(f"[CasterSitePlugin] tool_path exists: {os.path.exists(self.tool_path) if self.tool_path else 'None'}")
        print(f"[CasterSitePlugin] plugin_path: {self.plugin_path}")

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

        self.add_console_message(f"CASTER-site analysis completed! Output saved to: {', '.join(output_files)}", "info")
        QMessageBox.information(self, "Completed", f"CASTER-site analysis completed! Output saved to: {', '.join(output_files)}")

    def update_progress(self, progress_text):
        """更新进度"""
        self.progress_bar.setFormat(progress_text)
        self.progress_bar.setValue(0)  # 无限进度条

class CasterSitePluginEntry:
    def run(self, import_from=None, import_data=None):
        return CasterSitePlugin(import_from=import_from, import_data=import_data)