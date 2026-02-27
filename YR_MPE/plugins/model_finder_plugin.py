# model_finder_plugin.py
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
                             QSpinBox, QCheckBox, QLabel, QComboBox, QTableWidget, QTableWidgetItem,
                             QTextEdit, QToolButton, QApplication, QFrame, QWidget, QDialog)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon, QFont
import tempfile
import os
import re
from typing import List, Optional, Dict
from Bio import SeqIO

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess

# 导入分区模式相关模块
from .partition_mode import PartitionDefinition, PartitionMode, create_partition_from_dict
from .model_finder_partition_ui import ModelFinderPartitionDialog
from .model_finder_partition_thread import ModelFinderPartitionThread


class ModelFinderThread(BaseProcessThread):
    """ModelFinder线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "IQ-TREE ModelFinder"
        
    def execute_commands(self):
        """执行ModelFinder命令"""
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
                    "-m", "TESTONLY",
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"IQ-TREE ModelFinder execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的.iqtree文件
                iqtree_file = input_file + ".iqtree"
                if os.path.exists(iqtree_file):
                    output_files.append(iqtree_file)
                else:
                    # 尝试其他可能的命名方式
                    alternate_iqtree = input_file + ".phy.iqtree"
                    if os.path.exists(alternate_iqtree):
                        output_files.append(alternate_iqtree)
                    else:
                        self.console_output.emit(f"Warning: Could not find .iqtree file for {input_file}", "warning")
            
            self.progress.emit("Model finding completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Model finding exception: {str(e)}")


class ModelFinderPlugin(BasePlugin):
    """IQ-TREE ModelFinder插件类"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    export_model_result_signal = pyqtSignal(dict)  # 导出模型结果信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化ModelFinder插件"""
        super().__init__(import_from, import_data)

        # 初始化分区模式相关变量
        self.partition_mode_enabled = False
        self.partition_definitions: List[PartitionDefinition] = []
        self.partition_mode = PartitionMode.EL
        self.partition_rcluster = False
        self.partition_rcluster_percent = None
        self.partition_dialog = None

        # 处理不同来源的导入数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
        elif import_from == "DATASET_MANAGER" and import_data is not None:
            self.handle_dataset_import(import_data)
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "IQ-TREE ModelFinder"
        self.tool_name = "IQ-Tree 3"
        self.citation = [
            """Bui Quang Minh, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., in press. https://doi.org/10.1038/s41586-020-2176-9""",
            """Subha Kalyaanamoorthy, Bui Quang Minh, Thomas KF Wong, Arndt von Haeseler, and Lars S Jermiin (2017) ModelFinder: Fast model selection for accurate phylogenetic estimates. Nature Methods, 14:587–589. https://doi.org/10.1038/nmeth.4285"""
        ]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "PHYLIP": ["phy"]}
        self.output_types = {"IQ-TREE": ".iqtree"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

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

    def handle_dataset_import(self, import_data):
        """处理从DatasetManager导入的数据"""
        
        try:
            # import_data 应该是一个包含 dataset items 和配置信息的字典
            if not isinstance(import_data, dict):
                QMessageBox.warning(self, "Import Error", "Invalid dataset import data format")
                return

            # 获取 dataset items 和配置
            dataset_items = import_data.get('dataset_items', [])
            dataset_config = import_data.get('dataset_config', {})

            # dataset_items 应该已经是选中的 items（由 _prepare_import_data 过滤）
            selected_items = dataset_items
            if not selected_items:
                QMessageBox.warning(self, "No Selection", "No partitions selected from dataset. Please select at least one partition.")
                return
            
            # 检查所有选中的 partition 是否已比对
            unaligned_items = [item for item in selected_items if not item.is_aligned]
            if unaligned_items:
                unaligned_names = [item.loci_name for item in unaligned_items]
                warning_msg = "The following partitions are not aligned:\n" + "\n".join(f"  - {name}" for name in unaligned_names)
                warning_msg += "\n\nPlease align these partitions first or select only aligned partitions."
                QMessageBox.warning(self, "Alignment Required", warning_msg)
                # 强制打断自动导入，剩余一个空界面
                return
            
            # Get dataset settings (topo-linked/topo-unlinked, edge-linked/edge-unlinked)
            topo_linked = dataset_config.get('topo_linked', False)
            edge_linked = dataset_config.get('edge_linked', False)
            
            # Map to partition mode
            # topo-unlinked -> Topo-unlinked (TUL)
            # edge-linked -> Edge-linked Equal (TL)
            # edge-unlinked -> Edge-unlinked (EUL)
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
            
            # 获取所有序列名称（假设所有 partition 有相同的序列）
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
            temp_file = self.create_temp_file(suffix='.fas')
            
            with open(temp_file, 'w') as f:
                for seq_name, seq_content in supermatrix_sequences.items():
                    f.write(f">{seq_name}\n{seq_content}\n")
            
            self.import_file = temp_file
            self.imported_files = [temp_file]
            
            # 计算 partition 坐标并创建 partition definitions
            for item in selected_items:
                end_pos = current_pos + item.length - 1
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
            
            # 默认保持 rcluster 打开
            self.partition_rcluster = True
            self.partition_rcluster_percent = 0.95
            
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
            self.add_console_message(f"Rcluster enabled: {self.partition_rcluster}", "info")
            
        except Exception as e:
            QMessageBox.critical(self, "Import Error", f"Failed to import dataset: {str(e)}")
            # 强制打断自动导入，剩余一个空界面
            self.import_file = None
            self.imported_files = []
            self.partition_definitions = []
            self.partition_mode_enabled = False

    def _detect_sequence_type_conflicts(self, dataset_items, partition_definitions):
        """检测序列类型冲突"""
        # 检查每个 partition 的序列类型
        conflicts = []
        seq_types = []
        
        for item, partition_def in zip(dataset_items, partition_definitions):
            # 检测序列类型
            seq_type = self._detect_sequence_type(item.sequences)
            seq_types.append(seq_type)
            partition_def.seq_type = seq_type
        
        # 检查是否有冲突
        if len(set(seq_types)) > 1:
            conflict_types = set(seq_types)
            conflict_msg = f"Sequence type conflict detected: {', '.join(conflict_types)}\n"
            conflict_msg += "Partitions have different sequence types:\n"
            for item, seq_type in zip(dataset_items, seq_types):
                conflict_msg += f"  - {item.loci_name}: {seq_type}\n"
            conflict_msg += "\nThis may cause issues during model selection. Proceed with caution."
            
            QMessageBox.warning(self, "Sequence Type Conflict", conflict_msg)
            self.add_console_message("WARNING: Sequence type conflict detected", "warning")
        
        return seq_types
    
    def _detect_sequence_type(self, sequences):
        """检测序列类型"""
        if not sequences:
            return "DNA"
        
        # 获取第一个序列的字符串
        first_seq = str(sequences[0].seq).upper()
        
        # 检查是否包含蛋白质特殊字符
        protein_chars = set("ACDEFGHIKLMNPQRSTVWY*")
        dna_chars = set("ACGTNRYSWKMBDHV")
        
        # 统计字符
        total_chars = len(first_seq)
        protein_count = sum(1 for c in first_seq if c in protein_chars)
        dna_count = sum(1 for c in first_seq if c in dna_chars)
        
        # 如果包含大量的蛋白质特殊字符，则为蛋白质
        if protein_count > total_chars * 0.5:
            return "AA"
        # 否则为 DNA
        else:
            return "DNA"

    def setup_input_tab(self):
        """设置输入标签页"""
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
        self.file_path_edit.setPlaceholderText("Select sequence files...")
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
        self.sequence_text.setPlaceholderText("Or paste sequence in FASTA format...")
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
        
        # 参数组
        params_group = QGroupBox("ModelFinder Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # 序列类型
        self.seq_type_combo = QComboBox()
        self.seq_type_combo.addItems(["AUTO", "DNA", "PROT"])
        self.seq_type_combo.setCurrentText("AUTO")
        params_layout.addRow("Sequence Type:", self.seq_type_combo)
        
        # 准则
        self.criterion_combo = QComboBox()
        self.criterion_combo.addItems(["BIC", "AICc", "AIC"])
        self.criterion_combo.setCurrentText("BIC")
        params_layout.addRow("Criterion:", self.criterion_combo)
        
        # 线程数
        self.threads_spinbox = QSpinBox()
        self.threads_spinbox.setRange(1, 5)
        self.threads_spinbox.setValue(1)
        self.threads_spinbox.setSpecialValueText("AUTO")
        params_layout.addRow("Threads:", self.threads_spinbox)

        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets
    
    def browse_files(self):
        """浏览选择文件"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select sequence files", "", "Sequence files (*.fas *.fna *.fa *.fasta *.phy);;All files (*)"
        )
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
                background-color: #dc3545;
                color: white;
                border: none;
                border-radius: 10px;
                font-weight: bold;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_file_tag(file_path, tag_widget))
        tag_layout.addWidget(close_btn)
        
        # Add to layout
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append((file_path, tag_widget))
        
        # Show container if it was hidden
        self.file_tags_container.setVisible(True)
    
    def remove_file_tag(self, file_path, tag_widget):
        """Remove a file tag"""
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        
        # Remove from tags list
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        
        # Remove widget
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Update file path display
        if not self.imported_files:
            self.file_path_edit.setText("")
            self.file_tags_container.setVisible(False)
            # Show text input when no files
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
        elif len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        else:
            self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def get_display_name(self, file_path):
        """Get display name for file, handling duplicates"""
        filename = os.path.basename(file_path)
        
        # Check for duplicates
        duplicate_count = 0
        for existing_path in self.imported_files:
            if os.path.basename(existing_path) == filename and existing_path != file_path:
                duplicate_count += 1
        
        if duplicate_count > 0:
            # Use last directory name for duplicates
            parent_dir = os.path.basename(os.path.dirname(file_path))
            return f"{filename} ({parent_dir})"
        else:
            return filename
    
    def on_text_changed(self):
        """Handle text input changes"""
        text = self.sequence_text.toPlainText().strip()
        if text:
            # Hide file input when text is present
            self.file_path_edit.setVisible(False)
            self.file_browse_btn.setVisible(False)
            self.file_tags_container.setVisible(False)
            # Clear imported files
            self.clear_all_file_tags()
        else:
            # Show file input when text is empty
            self.file_path_edit.setVisible(True)
            self.file_browse_btn.setVisible(True)
            if self.imported_files:
                self.file_tags_container.setVisible(True)
    
    def clear_all_file_tags(self):
        """Clear all file tags"""
        for file_path, tag_widget in self.file_tags:
            self.file_tags_layout.removeWidget(tag_widget)
            tag_widget.deleteLater()
        
        self.imported_files.clear()
        self.file_tags.clear()
        self.file_tags_container.setVisible(False)
    
    def on_partition_mode_toggled(self, state):
        """处理分区模式复选框切换"""
        self.partition_mode_enabled = (state == Qt.Checked)
        
        if self.partition_mode_enabled:
            # 启用分区模式
            self.partition_config_btn.setVisible(True)
            self.partition_status_label.setVisible(True)
            self.add_console_message("Partition mode enabled. Please configure partitions.", "info")
        else:
            # 禁用分区模式
            self.partition_config_btn.setVisible(False)
            self.partition_status_label.setVisible(False)
            self.partition_definitions = []
            self.add_console_message("Partition mode disabled.", "info")
    
    def open_partition_dialog(self):
        """打开分区配置对话框"""
        try:
            # 创建分区对话框
            if self.partition_dialog is None:
                self.partition_dialog = ModelFinderPartitionDialog(parent=self)
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
                rcluster = self.partition_dialog.get_rcluster_enabled()
                rcluster_percent = self.partition_dialog.get_rcluster_percent()
                
                # 更新分区定义
                self.partition_definitions = partitions
                self.partition_mode = mode
                self.partition_rcluster = rcluster
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
        """Update partition status display"""
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
            
            self.partition_status_label.setText(f"Mode: {mode_text} | Partitions: {names_text}")

    def setup_output_tab(self):
        """设置输出预览标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 说明文字
        info_label = QLabel("Model selection results will be displayed here after analysis completes.")
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        
        # 结果表格
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(8)
        self.results_table.setHorizontalHeaderLabels([
            "Model", "LogL", "AIC", "w-AIC", "AICc", "w-AICc", "BIC", "w-BIC"
        ])
        self.results_table.setSortingEnabled(True)
        layout.addWidget(self.results_table)
        
        # 输出文件信息
        self.output_info = QLabel("Output file information will be displayed here")
        layout.addWidget(self.output_info)
        
        # # 添加导出到YR-MPEA按钮
        # self.export_to_mpea_btn = QPushButton("Export Best Model to YR-MPEA")
        # self.export_to_mpea_btn.clicked.connect(self.export_to_mpea)
        # self.export_to_mpea_btn.setVisible(False)  # 初始隐藏
        # layout.addWidget(self.export_to_mpea_btn)
        
        # # 添加导出完整模型表到YR-MPEA按钮
        # self.export_table_to_mpea_btn = QPushButton("Export Full Model Table to YR-MPEA")
        # self.export_table_to_mpea_btn.clicked.connect(self.export_table_to_mpea)
        # self.export_table_to_mpea_btn.setVisible(False)  # 初始隐藏
        # layout.addWidget(self.export_table_to_mpea_btn)
        
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.export_table_to_mpea)
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
                QMessageBox.warning(self, "Warning", "Please input sequence text!")
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
        
        # 如果启用了分区模式
        if self.partition_mode_enabled and self.partition_definitions:
            # 使用分区模式线程类处理
            # 这里返回特殊标记，让run_analysis方法使用分区线程
            return {"partition_mode": True}
        
        # 常规模型查找参数
        # 序列类型
        seq_type = self.seq_type_combo.currentText()
        
        if seq_type != "AUTO":
            seq_type_code = {"dna": "DNA", "prot": "AA"}[seq_type.lower()]
            params.extend(["-st", seq_type_code])
        
        # 准则
        criterion = self.criterion_combo.currentText()
        if criterion == "BIC":
            pass
        elif criterion == "AICc":
            params.extend(["-AICc"])
        elif criterion == "AIC":
            params.extend(["-AIC"])
        else:
            pass
        params.extend(["-mset", "ALL"])  # 测试所有模型
        
        # 线程数
        threads = self.threads_spinbox.value()
        if threads > 1:
            params.extend(["-nt", str(threads)])
        
        return params
    
    def parse_iqtree_file(self, iqtree_file_path):
        """解析.iqtree文件并提取模型信息"""
        models_data = []
        
        try:
            with open(iqtree_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 查找模型列表部分
            model_section_pattern = r"List of models sorted by.*?\n\n(.*?)\n\n"
            model_section_match = re.search(model_section_pattern, content, re.DOTALL)
            
            if model_section_match:
                model_lines = model_section_match.group(1).strip().split('\n')
                if len(model_lines) > 1:
                    # 跳过标题行
                    for line in model_lines[1:]:
                        # 使用正则表达式提取模型数据
                        # 格式: Model LogL AIC w-AIC AICc w-AICc BIC w-BIC
                        parts = line.split()
                        if len(parts) >= 10:
                            model_data = {
                                'Model': parts[0],
                                'LogL': parts[1],
                                'AIC': parts[2],
                                'w-AIC': parts[4],
                                'AICc': parts[5],
                                'w-AICc': parts[7],
                                'BIC': parts[8],
                                'w-BIC': parts[10]
                            }
                            models_data.append(model_data)
            
            return models_data
        except Exception as e:
            self.add_console_message(f"Error parsing .iqtree file: {str(e)}", "error")
            return []
    
    def parse_partition_results(self, iqtree_file_path):
        """解析分区模式的.iqtree文件并提取分区模型信息"""
        from .model_finder_partition_thread import PartitionResultParser
        import os
        
        self.add_console_message(f"Starting to parse partition results from: {os.path.basename(iqtree_file_path)}", "info")
        
        try:
            # 先从 .iqtree 文件获取统计信息和分区模型
            self.add_console_message(f"Parsing statistics from .iqtree file...", "info")
            iqt_results = PartitionResultParser.parse_partition_results(iqtree_file_path)
            
            # 检查 .iqtree 文件解析结果
            if not iqt_results:
                self.add_console_message("Failed to parse .iqtree file", "error")
                return None
            
            self.add_console_message(f"IQ-TREE parsed: {len(iqt_results.get('partition_models', {}))} partitions, LogL={iqt_results.get('log_likelihood', 'N/A')}", "info")
            
            # 然后尝试从 .best_scheme.nex 文件获取分区定义（如果存在）
            best_scheme_file = iqtree_file_path.replace('.iqtree', '.best_scheme.nex')
            if os.path.exists(best_scheme_file):
                self.add_console_message(f"Parsing partition definitions from .best_scheme.nex file...", "info")
                try:
                    bs_results = PartitionResultParser.parse_best_scheme_nex(best_scheme_file)
                    
                    if bs_results:
                        # 合并结果：使用 .best_scheme.nex 的分区定义，.iqtree 的统计信息
                        bs_results['log_likelihood'] = iqt_results.get('log_likelihood', 0.0)
                        bs_results['aic'] = iqt_results.get('aic', 0.0)
                        bs_results['aicc'] = iqt_results.get('aicc', 0.0)
                        bs_results['bic'] = iqt_results.get('bic', 0.0)
                        bs_results['aic_score'] = iqt_results.get('aic_score', 0.0)
                        bs_results['aicc_score'] = iqt_results.get('aicc_score', 0.0)
                        bs_results['bic_score'] = iqt_results.get('bic_score', 0.0)
                        self.add_console_message(f"Combined results: {len(bs_results.get('partition_models', {}))} partitions", "info")
                        return bs_results
                    else:
                        self.add_console_message("parse_best_scheme_nex returned None, using .iqtree results only", "warning")
                except Exception as e:
                    import traceback
                    self.add_console_message(f"Failed to parse .best_scheme.nex: {str(e)}", "warning")
                    self.add_console_message(f"Traceback: {traceback.format_exc()}", "warning")
            
            # 使用 .iqtree 文件的结果
            self.add_console_message(f"Using .iqtree results: {len(iqt_results.get('partition_models', {}))} partitions", "info")
            return iqt_results
            
        except Exception as e:
            import traceback
            self.add_console_message(f"Error parsing partition results: {str(e)}", "error")
            self.add_console_message(f"Traceback: {traceback.format_exc()}", "error")
            return None
    
    def _match_partition_by_range(self, partition_range: str, charset_ranges: Dict[str, List[str]]) -> str:
        """
        通过位点范围匹配分区到 charset 名称
        
        Args:
            partition_range: 分区的位点范围（如 "1-744" 或 "gene1.fas:1-1000"）
            charset_ranges: charset 名称到位点范围列表的映射
            
        Returns:
            匹配的 charset 名称，如果没有匹配则返回 None
        """
        try:
            # 提取纯数字范围（去除文件前缀）
            if ':' in partition_range:
                # 格式: "file:range"，提取range部分
                pure_range = partition_range.split(':', 1)[1]
            else:
                # 格式: "range"
                pure_range = partition_range
            
            # 解析分区范围
            start, end = map(int, pure_range.split('-'))
            
            # 检查每个 charset 的范围
            for charset_name, ranges in charset_ranges.items():
                for range_str in ranges:
                    # 解析 charset 范围
                    cs_start, cs_end = map(int, range_str.split('-'))
                    # 检查是否在范围内
                    if start >= cs_start and end <= cs_end:
                        return charset_name
            
            return None
        except:
            return None
    
    def _get_partition_model_by_range(self, partition_name: str, partition_range: str,
                                          partition_results: Dict) -> str:
            """
            通过位点范围获取分区的模型
            
            Args:
                partition_name: 原始分区名称
                partition_range: 分区的位点范围
                partition_results: 分区结果字典
                
            Returns:
                模型字符串
            """
            partition_models = partition_results.get('partition_models', {})
            charset_ranges = partition_results.get('charset_ranges', {})
            
            self.add_console_message(f"_get_partition_model_by_range: partition_name={partition_name}, partition_range={partition_range}", "info")
            self.add_console_message(f"_get_partition_model_by_range: partition_models keys={list(partition_models.keys())}", "info")
            self.add_console_message(f"_get_partition_model_by_range: charset_ranges={charset_ranges}", "info")
            
            # 直接匹配分区名称
            if partition_name in partition_models:
                model = partition_models[partition_name]['model']
                self.add_console_message(f"_get_partition_model_by_range: Direct match found! model={model}", "info")
                return model
            
            # 通过位点范围匹配
            if charset_ranges:
                matched_charset = self._match_partition_by_range(partition_range, charset_ranges)
                self.add_console_message(f"_get_partition_model_by_range: Matched charset={matched_charset}", "info")
                if matched_charset and matched_charset in partition_models:
                    model = partition_models[matched_charset]['model']
                    self.add_console_message(f"_get_partition_model_by_range: Range match found! model={model}", "info")
                    return model
            
            # 如果没有匹配，返回空字符串
            self.add_console_message(f"_get_partition_model_by_range: No match found, returning empty string", "info")
            return ''    
    def display_results(self, iqtree_files):
        """在输出标签页中显示结果"""
        if not iqtree_files:
            self.results_table.setRowCount(0)
            return
        
        # 根据是否启用分区模式选择不同的显示方式
        if self.partition_mode_enabled:
            # 分区模式显示
            self.display_partition_results(iqtree_files[0])
        else:
            # 常规模型显示
            self.display_single_model_results(iqtree_files[0])
    
    def display_single_model_results(self, iqtree_file):
        """显示单一模型结果"""
        # 解析结果
        models_data = self.parse_iqtree_file(iqtree_file)
        
        # 更新表格
        self.results_table.setRowCount(len(models_data))
        
        for row, model_data in enumerate(models_data):
            for col, (key, value) in enumerate(model_data.items()):
                item = QTableWidgetItem(value)
                if key not in ['Model']:
                    # 数值列右对齐
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.results_table.setItem(row, col, item)
        
        # 调整列宽
        self.results_table.resizeColumnsToContents()
        
        # 更新输出信息
        self.output_info.setText(f"Results from: {os.path.basename(iqtree_file)}")
    
    def display_partition_results(self, iqtree_file):
        """显示分区模型结果"""
        # 解析分区结果
        partition_results = self.parse_partition_results(iqtree_file)
        
        if not partition_results:
            self.results_table.setRowCount(1)
            self.results_table.setColumnCount(1)
            item = QTableWidgetItem("Failed to parse partition results")
            self.results_table.setItem(0, 0, item)
            self.output_info.setText(f"Results from: {os.path.basename(iqtree_file)} (Parse Error)")
            return
        
        # 调试：输出 partition_results 的所有键
        self.add_console_message(f"Partition results keys: {list(partition_results.keys())}", "info")
        self.add_console_message(f"Log-likelihood value: {partition_results.get('log_likelihood', 'N/A')}", "info")
        self.add_console_message(f"AICc values: aicc_score={partition_results.get('aicc_score', 'N/A')}, aicc={partition_results.get('aicc', 'N/A')}", "info")
        self.add_console_message(f"BIC values: bic_score={partition_results.get('bic_score', 'N/A')}, bic={partition_results.get('bic', 'N/A')}", "info")
        
        # 调试：输出 partition_models 数据
        partition_models = partition_results.get('partition_models', {})
        self.add_console_message(f"Partition models: {partition_models}", "info")
        
        # 调试：输出 partition_definitions 数据
        self.add_console_message(f"Partition definitions count: {len(self.partition_definitions)}", "info")
        for p in self.partition_definitions:
            self.add_console_message(f"  Partition: name={p.name}, range={p.model_range}, display_range={p.get_display_range()}", "info")
        
        # 修改表格为分区显示格式
        self.results_table.setColumnCount(4)
        self.results_table.setHorizontalHeaderLabels([
            "Partition", "Range", "Best Model", "LogL"
        ])
        
        # 先设置表格行数（在填充项目之前）
        total_rows = len(self.partition_definitions)
        self.results_table.setRowCount(total_rows)
        
        # 填充分区数据
        row = 0
        
        # 遍历分区定义
        for partition in self.partition_definitions:
            partition_name = partition.name
            partition_range = partition.get_display_range()
            
            self.add_console_message(f"Processing partition: name={partition_name}, model_range={partition.model_range}, display_range={partition_range}", "info")
            
            # 通过位点范围匹配获取模型
            best_model = self._get_partition_model_by_range(
                partition_name, partition.model_range, partition_results
            )
            
            self.add_console_message(f"Best model for partition '{partition_name}': {best_model}", "info")
            
            if not best_model:
                best_model = 'Unknown'
            
            # 设置表格项
            self.results_table.setItem(row, 0, QTableWidgetItem(partition_name))
            self.results_table.setItem(row, 1, QTableWidgetItem(partition_range))
            self.results_table.setItem(row, 2, QTableWidgetItem(best_model))
            self.results_table.setItem(row, 3, QTableWidgetItem('N/A'))
            
            row += 1
        
        # 设置行高
        self.results_table.resizeRowsToContents()
        
        # 调整列宽
        self.results_table.resizeColumnsToContents()
        self.results_table.setAlternatingRowColors(True)
        
        # 获取 info_label（第一个子组件）
        info_label = self.output_tab.layout().itemAt(0).widget()
        
        # 构建分区模式信息
        mode_text = {
            "EL": "Edge-linked partition model",
            "TL": "Edge-linked partition model",
            "EUL": "Edge-unlinked partition model",
            "TUL": "Topo-unlinked partition model"
        }.get(self.partition_mode.value, "Unknown partition model")
        
        # 构建统计信息 - 优先使用 score 后缀的键，否则使用无后缀的键
        logl = partition_results.get('log_likelihood', 
               partition_results.get('tree_log_likelihood', 'N/A'))
        aicc = partition_results.get('aicc_score', 
              partition_results.get('aicc', 'N/A'))
        bic = partition_results.get('bic_score', 
             partition_results.get('bic', 'N/A'))
        
        # 如果仍然是 'N/A' 或空值，检查是否有其他键名
        if aicc == 'N/A' or aicc == '':
            # 尝试从 partition_models 中获取第一个分区的 AICc
            partition_models = partition_results.get('partition_models', {})
            if partition_models:
                first_partition = list(partition_models.values())[0]
                self.add_console_message(f"First partition data: {first_partition}", "info")
                if 'aic' in first_partition:
                    aicc = first_partition['aic']
        
        if bic == 'N/A' or bic == '':
            # 尝试从 partition_models 中获取第一个分区的 BIC
            partition_models = partition_results.get('partition_models', {})
            if partition_models:
                first_partition = list(partition_models.values())[0]
                if 'bic' in first_partition:
                    bic = first_partition['bic']
        
        # 格式化数值显示
        if logl != 'N/A' and isinstance(logl, float):
            logl = f"{logl:.4f}"
        if aicc != 'N/A' and isinstance(aicc, float):
            aicc = f"{aicc:.4f}"
        if bic != 'N/A' and isinstance(bic, float):
            bic = f"{bic:.4f}"
        
        # 调试：输出最终值
        self.add_console_message(f"Final display values - LogL: {logl}, AICc: {aicc}, BIC: {bic}", "info")
        
        # 更新 info_label 显示分区模式
        info_label.setText(f"{mode_text}")
        
        # 更新 output_info 显示统计信息
        self.output_info.setText(f"Log-likelihood: {logl} | AICc: {aicc} | BIC: {bic}")
    
    def run_analysis(self):
        """运行ModelFinder分析"""
        if self.is_running:
            return
            
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
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
            self.add_console_message(f"Starting IQ-TREE ModelFinder with partition mode ({self.partition_mode.value})...", "info")
        else:
            self.add_console_message("Starting IQ-TREE ModelFinder...", "info")
        
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
        params = self.get_parameters()
        if isinstance(params, dict) and params.get("partition_mode"):
            # 使用分区模式线程
            self.analysis_thread = ModelFinderPartitionThread(
                self.tool_path,
                input_files[0] if len(input_files) == 1 else input_files,
                self.partition_definitions,
                self.partition_mode,
                self.partition_rcluster,
                self.partition_rcluster_percent,
                {
                    'seq_type': self.seq_type_combo.currentText(),
                    'criterion': self.criterion_combo.currentText(),
                    'threads': self.threads_spinbox.value()
                }
            )
        else:
            # 使用常规线程
            self.analysis_thread = ModelFinderThread(
                self.tool_path, input_files, params, self.imported_files
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
        
        # 显示结果
        self.display_results(output_files)
        
        # 切换到输出标签页
        self.tab_widget.setCurrentIndex(1)
        
        # 添加控制台消息
        self.add_console_message(f"ModelFinder completed successfully! Found {len(output_files)} result file(s)", "info")
        
        # 根据是否启用分区模式发送不同的数据格式
        if self.partition_mode_enabled and output_files:
            # 解析分区结果并发送分区模型数据
            partition_results = self.parse_partition_results(output_files[0])
            if partition_results:
                # 统一使用一种键名格式
                aicc_value = partition_results.get('aicc_score', partition_results.get('aicc', 0.0))
                bic_value = partition_results.get('bic_score', partition_results.get('bic', 0.0))
                
                # 构建分区模型数据
                model_data = {
                    "type": "partitioned",
                    "partition_mode": self.partition_mode.value,
                    "best_scheme": partition_results.get('overall_model', ''),
                    "partitions": [],
                    "statistics": {
                        "logL": partition_results.get('log_likelihood', 0.0),
                        "aicc": aicc_value,
                        "bic": bic_value
                    }
                }
                
                # 添加分区信息
                for partition in self.partition_definitions:
                    partition_name = partition.name
                    partition_range = partition.model_range  # 使用原始位点范围
                    
                    # 通过位点范围匹配获取模型
                    best_model = self._get_partition_model_by_range(
                        partition_name, partition_range, partition_results
                    )
                    
                    model_data["partitions"].append({
                        "name": partition_name,
                        "range": partition.get_display_range(),
                        "best_model": best_model,
                        "logL": 0.0  # IQ-TREE 输出中可能没有每个分区的 logL
                    })
                
                # 发送信号
                self.export_model_result_signal.emit(model_data)
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from == "YR_MPEA":
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
            
        # 保存输出文件路径供导入使用
        self.alignment_output_files = output_files
        
        QMessageBox.information(self, "Completed", "Model selection completed!")
    
    def analysis_error(self, error_message):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 添加控制台消息
        self.add_console_message(f"ModelFinder failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"ModelFinder failed: {error_message}")
    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "ModelFinder has been aborted.")

    def import_alignment_to_platform(self, sequences):
        """将比对结果导入到平台的工作区"""
        # 发送信号将数据导入到平台
        self.import_alignment_signal.emit(sequences)

    def export_to_mpea(self):
        """导出最佳模型到YR-MPEA"""
        if self.results_table.rowCount() == 0:
            QMessageBox.warning(self, "Warning", "No model results to export.")
            return
            
        # 获取第一行（最佳模型）的数据
        best_model_data = {}
        for col in range(self.results_table.columnCount()):
            header = self.results_table.horizontalHeaderItem(col).text()
            value = self.results_table.item(0, col).text()
            best_model_data[header] = value
            
        # 发送信号将模型结果导出到YR-MPEA
        self.export_model_result_signal.emit(best_model_data)
        
        QMessageBox.information(self, "Success", f"Best model '{best_model_data['Model']}' exported to YR-MPEA successfully!")

    def export_table_to_mpea(self):
        """导出完整模型表到YR-MPEA"""
        if self.results_table.rowCount() == 0:
            QMessageBox.warning(self, "Warning", "No model results to export.")
            return
            
        # 获取完整的模型表数据
        model_table_data = []
        headers = []
        
        # 获取表头
        for col in range(self.results_table.columnCount()):
            headers.append(self.results_table.horizontalHeaderItem(col).text())
        
        # 获取所有行数据
        for row in range(self.results_table.rowCount()):
            row_data = {}
            for col in range(self.results_table.columnCount()):
                header = headers[col]
                value = self.results_table.item(row, col).text()
                row_data[header] = value
            model_table_data.append(row_data)
            
        # 发送信号将完整模型表导出到YR-MPEA
        table_export_data = {
            "type": "model_table",
            "headers": headers,
            "data": model_table_data
        }
        self.export_model_result_signal.emit(table_export_data)
        
        QMessageBox.information(self, "Success", "Full model table exported to YR-MPEA successfully!")

    def copy_citation(self):
        """复制引用"""
        clipboard = QApplication.clipboard()
        clipboard.setText("\n".join(self.citation))
        QMessageBox.information(self, "Copied", "Citation copied to clipboard!")


class ModelFinderPluginEntry:
    """ModelFinder插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return ModelFinderPlugin(import_from, import_data)