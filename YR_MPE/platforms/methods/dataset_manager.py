"""
Dataset Manager Module
实现多序列数据集的管理和批量处理功能
"""
import os
import csv
from typing import List, Dict, Optional
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QTableWidget, 
                            QTableWidgetItem, QPushButton, QCheckBox, QLabel, 
                            QFileDialog, QMessageBox, QMenuBar, QMenu, QAction,
                            QRadioButton, QButtonGroup, QGroupBox, QWidget)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor
from Bio import SeqIO


class DatasetItem:
    """Dataset数据项模型"""
    
    def __init__(self):
        self.selected = False          # 是否选中
        self.loci_name = ""           # 位点名称
        self.length = 0               # 序列长度
        self.sequence_count = 0       # 序列数量
        self.is_aligned = False       # 是否已比对（可手动修改）
        self.file_path = ""           # 原始文件路径
        self.sequences = []           # 序列数据
        
    def __str__(self):
        return f"DatasetItem(loci_name={self.loci_name}, length={self.length}, count={self.sequence_count})"


class DatasetManager(QDialog):
    """Dataset管理对话框"""
    
    # 信号定义
    dataset_processed = pyqtSignal(list)  # 处理完成的dataset列表
    
    def __init__(self, dataset_name: str = "Default Dataset", plugin_factory=None):
        super().__init__()
        self.dataset_name = dataset_name
        self.dataset_items: List[DatasetItem] = []
        self.plugin_factory = plugin_factory  # 添加plugin_factory引用
        self.plugin_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self.init_ui()
        
    def init_ui(self):
        """初始化UI"""
        self.setWindowTitle(f"Dataset Manager - {self.dataset_name}")
        self.setMinimumSize(800, 600)
        
        # 主布局
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # 菜单栏
        menubar = QMenuBar()
        file_menu = menubar.addMenu("&File")
        
        # 导出选项
        export_fasta_action = QAction("Export (to multiple FASTA files)", self)
        export_fasta_action.triggered.connect(self.export_to_multiple_fasta)
        file_menu.addAction(export_fasta_action)
        
        export_nexus_action = QAction("Export (to partitioned NEXUS)", self)
        export_nexus_action.triggered.connect(self.export_to_partitioned_nexus)
        file_menu.addAction(export_nexus_action)
        
        export_summary_action = QAction("Export (Summary)", self)
        export_summary_action.triggered.connect(self.export_summary)
        file_menu.addAction(export_summary_action)
        
        # Batch Processing 菜单（独立顶级菜单）
        batch_menu = menubar.addMenu("&Batch Processing")
        
        align_menu = QMenu("Align by...", self)
        batch_menu.addMenu(align_menu)
        
        # 添加比对工具选项
        clustal_omega_action = QAction("Clustal Omega", self)
        clustal_omega_action.triggered.connect(self.batch_align_clustal_omega)
        align_menu.addAction(clustal_omega_action)
        
        mafft_action = QAction("MAFFT", self)
        mafft_action.triggered.connect(self.batch_align_mafft)
        align_menu.addAction(mafft_action)
        
        muscle5_action = QAction("Muscle 5", self)
        muscle5_action.triggered.connect(self.batch_align_muscle5)
        align_menu.addAction(muscle5_action)
        
        macse2_action = QAction("MACSE 2", self)
        macse2_action.triggered.connect(self.batch_align_macse2)
        align_menu.addAction(macse2_action)
        
        trim_menu = QMenu("Trim by...", self)
        batch_menu.addMenu(trim_menu)
        
        # 添加修剪工具选项
        trimal_action = QAction("TrimAl", self)
        trimal_action.triggered.connect(self.batch_trim_trimal)
        trim_menu.addAction(trimal_action)
        
        gblocks_action = QAction("GBlocks", self)
        gblocks_action.triggered.connect(self.batch_trim_gblocks)
        trim_menu.addAction(gblocks_action)
        
        main_layout.setMenuBar(menubar)
        
        # 基础设置区域
        settings_group = QGroupBox("Settings")
        settings_layout = QHBoxLayout()
        settings_group.setLayout(settings_layout)
        
        # Topology设置
        topo_group = QButtonGroup(self)
        self.topo_linked_radio = QRadioButton("Topo linked")
        self.topo_unlinked_radio = QRadioButton("Topo unlinked")
        self.topo_unlinked_radio.setChecked(True)
        topo_group.addButton(self.topo_linked_radio)
        topo_group.addButton(self.topo_unlinked_radio)
        
        settings_layout.addWidget(QLabel("Topology:"))
        settings_layout.addWidget(self.topo_linked_radio)
        settings_layout.addWidget(self.topo_unlinked_radio)
        
        # Edge设置（仅在Topo linked时启用）
        edge_group = QButtonGroup(self)
        self.edge_linked_radio = QRadioButton("Edge linked")
        self.edge_unlinked_radio = QRadioButton("Edge unlinked")
        self.edge_unlinked_radio.setChecked(True)
        edge_group.addButton(self.edge_linked_radio)
        edge_group.addButton(self.edge_unlinked_radio)
        
        self.edge_settings_widget = QWidget()
        edge_layout = QHBoxLayout()
        self.edge_settings_widget.setLayout(edge_layout)
        edge_layout.addWidget(QLabel("Edge:"))
        edge_layout.addWidget(self.edge_linked_radio)
        edge_layout.addWidget(self.edge_unlinked_radio)
        edge_layout.addStretch()
        
        # 默认禁用Edge设置
        self.edge_settings_widget.setEnabled(False)
        
        # 连接Topology选择信号
        self.topo_linked_radio.toggled.connect(self.on_topo_linked_toggled)
        
        settings_layout.addWidget(self.edge_settings_widget)
        settings_layout.addStretch()
        
        main_layout.addWidget(settings_group)
        
        # 表格区域
        table_layout = QVBoxLayout()
        
        # "+"按钮
        add_button = QPushButton("+")
        add_button.clicked.connect(self.add_datasets)
        table_layout.addWidget(add_button)
        
        # 数据表格
        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels([
            "Selected", "Loci Name", "Length", "Sequence Count", 
            "Aligned", "View"
        ])
        self.table.horizontalHeader().setStretchLastSection(True)
        # 启用右键菜单
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.show_table_context_menu)
        table_layout.addWidget(self.table)
        
        main_layout.addLayout(table_layout)
        
    def on_topo_linked_toggled(self, checked: bool):
        """Topology链接状态切换时的处理"""
        self.edge_settings_widget.setEnabled(checked)
        
    def add_datasets(self):
        """添加新的序列数据集"""
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFiles)
        file_dialog.setNameFilter(
            "Sequence files (*.fas *.fasta *.fa *.fna *.phy *.nex *.nexus)"
        )
        
        if file_dialog.exec_():
            selected_files = file_dialog.selectedFiles()
            for file_path in selected_files:
                try:
                    dataset_item = self.parse_sequence_file(file_path)
                    if dataset_item:
                        self.dataset_items.append(dataset_item)
                        self.add_dataset_to_table(dataset_item)
                except Exception as e:
                    QMessageBox.warning(
                        self, "Error", 
                        f"Failed to parse file {os.path.basename(file_path)}: {str(e)}"
                    )
                    
    def parse_sequence_file(self, file_path: str) -> Optional[DatasetItem]:
        """解析序列文件并创建DatasetItem"""
        # 确定文件格式
        file_ext = os.path.splitext(file_path)[1].lower()
        format_map = {
            '.fas': 'fasta', '.fasta': 'fasta', '.fa': 'fasta', '.fna': 'fasta',
            '.phy': 'phylip', '.nex': 'nexus', '.nexus': 'nexus'
        }
        
        if file_ext not in format_map:
            raise ValueError(f"Unsupported file format: {file_ext}")
            
        file_format = format_map[file_ext]
        
        # 读取序列
        sequences = list(SeqIO.parse(file_path, file_format))
        if not sequences:
            raise ValueError("No sequences found in file")
            
        # 创建DatasetItem
        item = DatasetItem()
        item.file_path = file_path
        item.loci_name = os.path.splitext(os.path.basename(file_path))[0]
        item.sequences = sequences
        item.sequence_count = len(sequences)
        item.length = len(str(sequences[0].seq))
        
        # 检查是否可能已比对（所有序列长度相同）
        lengths = [len(str(seq.seq)) for seq in sequences]
        all_same_length = len(set(lengths)) == 1
        
        # 如果所有序列长度相同，询问用户是否已比对
        if all_same_length:
            reply = QMessageBox.question(
                self, "Alignment Status", 
                f"All sequences in '{item.loci_name}' have the same length.\n"
                "Is this dataset already aligned?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
            )
            if reply == QMessageBox.Yes:
                item.is_aligned = True
            elif reply == QMessageBox.No:
                item.is_aligned = False
            else:  # Cancel
                return None
        else:
            # 长度不同，肯定未比对
            item.is_aligned = False
        
        return item
        
    def add_dataset_to_table(self, item: DatasetItem):
        """将DatasetItem添加到表格中"""
        row = self.table.rowCount()
        self.table.insertRow(row)
        
        # Selected checkbox
        selected_checkbox = QCheckBox()
        selected_checkbox.setChecked(item.selected)
        selected_checkbox.stateChanged.connect(
            lambda state, r=row: self.on_selected_changed(r, state)
        )
        self.table.setCellWidget(row, 0, selected_checkbox)
        
        # Loci Name
        loci_label = QLabel(item.loci_name)
        self.table.setCellWidget(row, 1, loci_label)
        
        # Length
        length_label = QLabel(str(item.length))
        self.table.setCellWidget(row, 2, length_label)
        
        # Sequence Count (带颜色)
        count_label = QLabel(str(item.sequence_count))
        self.colorize_sequence_count(count_label, item.sequence_count)
        self.table.setCellWidget(row, 3, count_label)
        
        # Aligned status
        aligned_label = QLabel("✓" if item.is_aligned else "✗")
        self.table.setCellWidget(row, 4, aligned_label)
        
        # View button
        view_button = QPushButton("View")
        view_button.clicked.connect(lambda _, r=row: self.view_dataset(r))
        self.table.setCellWidget(row, 5, view_button)
        
    def colorize_sequence_count(self, label: QLabel, count: int):
        """根据序列数量着色"""
        if not self.dataset_items:
            return
            
        # 计算平均值和标准差
        counts = [item.sequence_count for item in self.dataset_items]
        mean_count = sum(counts) / len(counts)
        std_count = (sum((x - mean_count) ** 2 for x in counts) / len(counts)) ** 0.5
        
        # 设置颜色
        if std_count == 0:
            # 所有值相同
            label.setStyleSheet("color: black;")
        elif count < mean_count - std_count:
            # 异常低值
            label.setStyleSheet("background-color: #ffcccc; color: black;")
        elif count > mean_count + std_count:
            # 异常高值
            label.setStyleSheet("background-color: #cce6ff; color: black;")
        else:
            # 正常范围
            label.setStyleSheet("color: black;")
            
    def on_selected_changed(self, row: int, state: int):
        """选中状态改变时的处理"""
        if 0 <= row < len(self.dataset_items):
            self.dataset_items[row].selected = (state == Qt.Checked)
            
    def view_dataset(self, row: int):
        """查看数据集"""
        try:
            if 0 <= row < len(self.dataset_items):
                dataset = self.dataset_items[row]
                # 使用PluginFactory获取序列查看器
                seq_viewer = self.plugin_factory.get_sequence_viewer()
                # 序列查看器不需要参数，直接显示
                viewer_widget = seq_viewer.run()
                
                # 创建对话框显示查看器
                from PyQt5.QtWidgets import QDialog
                dialog = QDialog()
                dialog.setWindowTitle(f"View Dataset: {dataset.name}")
                dialog.setMinimumSize(800, 600)
                dialog.setLayout(QVBoxLayout())
                dialog.layout().addWidget(viewer_widget)
                dialog.exec_()
        except Exception as e:
            QMessageBox.critical(None, "Error", f"Failed to open sequence viewer:\n{str(e)}")
            
    def export_to_multiple_fasta(self):
        """导出为多个FASTA文件"""
        # 获取选中的datasets
        selected_items = [item for item in self.dataset_items if item.selected]
        if not selected_items:
            QMessageBox.warning(self, "Warning", "No datasets selected for export.")
            return
            
        # 选择导出目录
        export_dir = QFileDialog.getExistingDirectory(
            self, "Select Export Directory"
        )
        if not export_dir:
            return
            
        try:
            for item in selected_items:
                export_path = os.path.join(export_dir, f"{item.loci_name}.fasta")
                with open(export_path, 'w') as f:
                    SeqIO.write(item.sequences, f, 'fasta')
                    
            QMessageBox.information(
                self, "Success", 
                f"Successfully exported {len(selected_items)} datasets to {export_dir}"
            )
        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to export datasets: {str(e)}"
            )
            
    def export_to_partitioned_nexus(self):
        """导出为分区NEXUS格式"""
        # 获取选中的datasets
        selected_items = [item for item in self.dataset_items if item.selected]
        if not selected_items:
            QMessageBox.warning(self, "Warning", "No datasets selected for export.")
            return
            
        # 选择导出文件
        export_path, _ = QFileDialog.getSaveFileName(
            self, "Save Partitioned NEXUS File", "", "NEXUS Files (*.nex)"
        )
        if not export_path:
            return
            
        try:
            # 合并所有序列
            all_sequences = []
            charsets = []
            current_pos = 1
            
            for item in selected_items:
                # 为每个dataset的序列添加前缀以避免名称冲突
                for seq in item.sequences:
                    new_seq = seq.__class__()
                    new_seq.id = f"{item.loci_name}_{seq.id}"
                    new_seq.name = f"{item.loci_name}_{seq.name}"
                    new_seq.description = seq.description
                    new_seq.seq = seq.seq
                    all_sequences.append(new_seq)
                    
                # 记录charset
                end_pos = current_pos + item.length - 1
                charsets.append(f"    CHARSET {item.loci_name} = {current_pos}-{end_pos};")
                current_pos = end_pos + 1
                
            # 写入NEXUS文件
            with open(export_path, 'w') as f:
                f.write("#NEXUS\n\n")
                f.write("BEGIN DATA;\n")
                f.write(f"    DIMENSIONS NTAX={len(all_sequences)} NCHAR={current_pos - 1};\n")
                f.write("    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
                f.write("    MATRIX\n")
                
                for seq in all_sequences:
                    f.write(f"    {seq.id}    {str(seq.seq)}\n")
                    
                f.write("    ;\n")
                f.write("END;\n\n")
                
                f.write("BEGIN SETS;\n")
                for charset in charsets:
                    f.write(charset + "\n")
                f.write("END;\n")
                
            QMessageBox.information(
                self, "Success", 
                f"Successfully exported partitioned NEXUS file to {export_path}"
            )
        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to export partitioned NEXUS: {str(e)}"
            )
            
    def export_summary(self):
        """导出摘要信息"""
        if not self.dataset_items:
            QMessageBox.warning(self, "Warning", "No datasets to export summary.")
            return
            
        export_path, _ = QFileDialog.getSaveFileName(
            self, "Save Summary CSV", "", "CSV Files (*.csv)"
        )
        if not export_path:
            return
            
        try:
            with open(export_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow([
                    "Loci Name", "Length", "Sequence Count", 
                    "Aligned", "File Path"
                ])
                
                for item in self.dataset_items:
                    writer.writerow([
                        item.loci_name,
                        item.length,
                        item.sequence_count,
                        "Yes" if item.is_aligned else "No",
                        item.file_path
                    ])
                    
            QMessageBox.information(
                self, "Success", 
                f"Successfully exported summary to {export_path}"
            )
        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"Failed to export summary: {str(e)}"
            )
            
    def batch_align_clustal_omega(self):
        """批量执行Clustal Omega比对"""
        self._batch_align_with_tool("clustal_omega")
        
    def batch_align_mafft(self):
        """批量执行MAFFT比对"""
        self._batch_align_with_tool("mafft")
        
    def batch_align_muscle5(self):
        """批量执行Muscle 5比对"""
        self._batch_align_with_tool("muscle5")
        
    def batch_align_macse2(self):
        """批量执行MACSE 2比对"""
        self._batch_align_with_tool("macse2")
        
    def batch_trim_trimal(self):
        """批量执行TrimAl修剪"""
        self._batch_trim_with_tool("trimal")
        
    def batch_trim_gblocks(self):
        """批量执行GBlocks修剪"""
        self._batch_trim_with_tool("gblocks")
        
    def _batch_align_with_tool(self, tool_name: str):
        """使用指定工具批量执行比对"""
        # 获取选中的且未比对的数据集
        selected_items = [item for item in self.dataset_items if item.selected and not item.is_aligned]
        if not selected_items:
            QMessageBox.warning(self, "Warning", f"No unaligned datasets selected for {tool_name} alignment.")
            return
            
        try:
            # 获取插件管理器
            from ..methods import PluginManager, WorkspaceManager
            workspace_manager = WorkspaceManager()
            plugin_manager = PluginManager(workspace_manager)
            plugin_manager.register_all_plugins()
            
            # 创建插件实例
            plugin_instance = plugin_manager.create_plugin_instance(tool_name)
            if not plugin_instance:
                QMessageBox.critical(self, "Error", f"Plugin {tool_name} not available.")
                return
                
            # TODO: 实现批量比对逻辑
            # 这里需要调用插件的执行方法，传入所有选中的序列数据
            QMessageBox.information(
                self, "Batch Alignment", 
                f"Starting batch alignment with {tool_name} for {len(selected_items)} datasets.\n"
                "Note: Full implementation requires integration with PluginExecutor."
            )
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start batch alignment: {str(e)}")
            
    def _batch_trim_with_tool(self, tool_name: str):
        """使用指定工具批量执行修剪"""
        # 获取选中的且已比对的数据集
        selected_items = [item for item in self.dataset_items if item.selected and item.is_aligned]
        if not selected_items:
            QMessageBox.warning(self, "Warning", f"No aligned datasets selected for {tool_name} trimming.")
            return
            
        try:
            # 获取插件管理器
            from ..methods import PluginManager, WorkspaceManager
            workspace_manager = WorkspaceManager()
            plugin_manager = PluginManager(workspace_manager)
            plugin_manager.register_all_plugins()
            
            # 创建插件实例
            plugin_instance = plugin_manager.create_plugin_instance(tool_name)
            if not plugin_instance:
                QMessageBox.critical(self, "Error", f"Plugin {tool_name} not available.")
                return
                
            # TODO: 实现批量修剪逻辑
            QMessageBox.information(
                self, "Batch Trimming", 
                f"Starting batch trimming with {tool_name} for {len(selected_items)} datasets.\n"
                "Note: Full implementation requires integration with PluginExecutor."
            )
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start batch trimming: {str(e)}")

    def show_table_context_menu(self, position):
        """显示表格右键菜单"""
        row = self.table.rowAt(position.y())
        if row < 0 or row >= len(self.dataset_items):
            return
            
        menu = QMenu()
        item = self.dataset_items[row]
        
        if item.is_aligned:
            mark_unaligned_action = QAction("Mark as Unaligned", self)
            mark_unaligned_action.triggered.connect(lambda: self.mark_dataset_unaligned(row))
            menu.addAction(mark_unaligned_action)
        else:
            mark_aligned_action = QAction("Mark as Aligned", self)
            mark_aligned_action.triggered.connect(lambda: self.mark_dataset_aligned(row))
            menu.addAction(mark_aligned_action)
            
        remove_action = QAction("Remove Dataset", self)
        remove_action.triggered.connect(lambda: self.remove_dataset(row))
        menu.addAction(remove_action)
        
        menu.exec_(self.table.mapToGlobal(position))
        
    def mark_dataset_aligned(self, row: int):
        """标记数据集为已比对"""
        if 0 <= row < len(self.dataset_items):
            self.dataset_items[row].is_aligned = True
            aligned_label = QLabel("✓")
            self.table.setCellWidget(row, 4, aligned_label)
            
    def mark_dataset_unaligned(self, row: int):
        """标记数据集为未比对"""
        if 0 <= row < len(self.dataset_items):
            self.dataset_items[row].is_aligned = False
            aligned_label = QLabel("✗")
            self.table.setCellWidget(row, 4, aligned_label)
            
    def remove_dataset(self, row: int):
        """移除数据集"""
        if 0 <= row < len(self.dataset_items):
            del self.dataset_items[row]
            self.table.removeRow(row)
