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

# 添加新的导入
from .export_partitioned_nexus import export_partitioned_nexus
from .import_partitioned_nexus import import_partitioned_nexus
from .dataset_models import SELECTION_STATE_GREEN, SELECTION_STATE_NONE


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
        self.name = ""                # 添加name属性以兼容Sequence Viewer
        
    def __str__(self):
        return f"DatasetItem(loci_name={self.loci_name}, length={self.length}, count={self.sequence_count})"
        
    def set_name(self, name):
        """设置name属性"""
        self.name = name
        self.loci_name = name  # 保持loci_name同步
        
    def get_name(self):
        """获取name属性"""
        return self.name


class DatasetManager(QDialog):
    """Dataset管理对话框"""
    
    # 信号定义
    dataset_processed = pyqtSignal(list)  # 处理完成的dataset列表
    
    def __init__(self, dataset_name: str = "Default Dataset", plugin_factory=None, workspace=None):
        super().__init__()
        self.dataset_name = dataset_name
        self.dataset_items: List[DatasetItem] = []
        self.plugin_factory = plugin_factory  # 添加plugin_factory引用
        self.workspace = workspace  # 添加workspace引用
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
        
        # 导入选项
        import_nexus_action = QAction("Import (from partitioned NEXUS)", self)
        import_nexus_action.triggered.connect(self.import_from_partitioned_nexus)
        file_menu.addAction(import_nexus_action)
        
        file_menu.addSeparator()
        
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
            "Aligned", "Options"
        ])
        self.table.horizontalHeader().setStretchLastSection(True)
        # 启用右键菜单
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.show_table_context_menu)
        table_layout.addWidget(self.table)
        
        main_layout.addLayout(table_layout)
        
        # 恢复保存的设置
        self._restore_settings()
        
    def _restore_settings(self):
        """从 dataset.settings 恢复保存的设置"""
        print("[DEBUG] _restore_settings called")
        print(f"[DEBUG] Before restore: dataset_items count = {len(self.dataset_items)}")
        if not self.workspace or not hasattr(self.workspace, 'current_dataset_id'):
            print("[DEBUG] No workspace or current_dataset_id")
            return
        
        dataset_id = self.workspace.current_dataset_id
        print(f"[DEBUG] dataset_id: {dataset_id}")
        if not dataset_id or not self.workspace.dataset_selection_manager:
            print("[DEBUG] No dataset_id or dataset_selection_manager")
            return
        
        dataset = self.workspace.dataset_selection_manager.get_dataset(dataset_id)
        print(f"[DEBUG] dataset: {dataset}")
        if not dataset:
            print("[DEBUG] Dataset not found")
            return
        
        # 清空现有的dataset_items和表格
        print(f"[DEBUG] Clearing existing {len(self.dataset_items)} items and table")
        self.dataset_items.clear()
        self.table.setRowCount(0)
        
        # 恢复 topo_linked 和 edge_linked 设置
        if 'topo_linked' in dataset.settings:
            self.topo_linked_radio.setChecked(dataset.settings['topo_linked'])
            self.topo_unlinked_radio.setChecked(not dataset.settings['topo_linked'])
        
        if 'edge_linked' in dataset.settings:
            self.edge_linked_radio.setChecked(dataset.settings['edge_linked'])
            self.edge_unlinked_radio.setChecked(not dataset.settings['edge_linked'])
        
        # 恢复 dataset_items
        if 'dataset_items' in dataset.settings:
            saved_items = dataset.settings['dataset_items']
            print(f"[DEBUG] Found {len(saved_items)} saved items")
            print(f"[DEBUG] Dataset items before restore: {len(dataset.items)}")
            
            # 清空dataset.items
            dataset.items.clear()
            
            # 禁用自动保存，避免在恢复过程中保存不完整的数据
            original_disable_auto_save = self.workspace.dataset_selection_manager.disable_auto_save
            self.workspace.dataset_selection_manager.disable_auto_save = True
            
            from .dataset_models import DatasetItem as NewDatasetItem, ITEM_TYPE_ALIGNMENT
            
            for item_data in saved_items:
                print(f"[DEBUG] Restoring item: {item_data.get('loci_name')}, is_aligned: {item_data.get('is_aligned')}, selected: {item_data.get('selected')}")
                # 创建 DatasetItem（本地类）
                item = DatasetItem()
                item.selected = item_data.get('selected', False)
                item.loci_name = item_data.get('loci_name', '')
                item.length = item_data.get('length', 0)
                item.sequence_count = item_data.get('sequence_count', 0)
                item.is_aligned = item_data.get('is_aligned', False)
                item.file_path = item_data.get('file_path', '')
                item.name = item.loci_name
                
                # 处理序列数据
                sequences_data = item_data.get('sequences')
                if sequences_data is None:
                    # 从文件加载
                    if item.file_path and os.path.exists(item.file_path):
                        try:
                            from Bio import SeqIO
                            sequences = list(SeqIO.parse(item.file_path, 'fasta'))
                            if sequences:
                                item.sequences = sequences
                        except Exception as e:
                            print(f"Failed to load sequences from {item.file_path}: {e}")
                elif sequences_data:
                    # 从保存的序列数据恢复
                    item.sequences = []
                    from Bio.SeqRecord import SeqRecord
                    for seq_dict in sequences_data:
                        seq = SeqRecord(
                            seq_dict['sequence'],
                            id=seq_dict['id'],
                            description=seq_dict.get('description', '')
                        )
                        item.sequences.append(seq)
                
                self.dataset_items.append(item)
                self.add_dataset_to_table(item)
                
                # 创建新格式的DatasetItem并添加到Dataset Selection Manager
                new_item = NewDatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
                new_item.dataset_id = dataset_id
                new_item.loci_name = item.loci_name
                new_item.file_path = item.file_path
                new_item.length = item.length
                new_item.sequence_count = item.sequence_count
                new_item.is_aligned = item.is_aligned
                
                # 转换序列数据
                if hasattr(item, 'sequences') and item.sequences:
                    new_item.sequences = item.sequences
                
                # 设置选择状态
                if item.selected:
                    new_item.selection_state = SELECTION_STATE_GREEN
                else:
                    new_item.selection_state = SELECTION_STATE_NONE
                
                # 添加到 Dataset Selection Manager（使用add_item来正确更新dataset.items）
                success = self.workspace.dataset_selection_manager.add_item(new_item, dataset_id)
                
                # 如果 item 被选中，添加到 selected_items 集合
                if item.selected and success:
                    self.workspace.dataset_selection_manager.selected_items.add(new_item.id)
                    print(f"[DEBUG] Added to selected_items: {new_item.loci_name} ({new_item.id})")
            
            # 恢复自动保存设置
            self.workspace.dataset_selection_manager.disable_auto_save = original_disable_auto_save
            
            # 手动保存一次完整的状态
            self.workspace.dataset_selection_manager._save_state()
            
            print(f"[DEBUG] Dataset items after restore: {len(dataset.items)}")
        
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
                        
                        # 同时创建DatasetItem对象并添加到Dataset Selection Manager
                        if self.workspace and self.workspace.dataset_selection_manager and hasattr(self.workspace, 'current_dataset_id'):
                            dataset_id = self.workspace.current_dataset_id
                            if dataset_id:
                                from .dataset_models import DatasetItem as NewDatasetItem, ITEM_TYPE_ALIGNMENT
                                new_item = NewDatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
                                new_item.dataset_id = dataset_id
                                new_item.loci_name = dataset_item.loci_name
                                new_item.file_path = dataset_item.file_path
                                new_item.length = dataset_item.length
                                new_item.sequence_count = dataset_item.sequence_count
                                new_item.is_aligned = dataset_item.is_aligned
                                new_item.sequences = dataset_item.sequences
                                new_item.selection_state = SELECTION_STATE_NONE  # 初始状态为NONE
                                
                                # 添加到manager
                                success = self.workspace.dataset_selection_manager.add_item(new_item, dataset_id)
                                print(f"[DEBUG] add_datasets: Created DatasetItem {dataset_item.loci_name} ({new_item.id}), success={success}")
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
        print(f"[DEBUG] add_dataset_to_table called for: {item.loci_name}, current dataset_items count = {len(self.dataset_items)}")
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
        
        # Options layout with View and Remove buttons
        options_widget = QWidget()
        options_layout = QHBoxLayout()
        options_layout.setContentsMargins(2, 2, 2, 2)
        options_layout.setSpacing(4)
        
        view_button = QPushButton("View")
        view_button.clicked.connect(lambda _, r=row: self.view_dataset(r))
        options_layout.addWidget(view_button)
        
        remove_button = QPushButton("Remove")
        remove_button.clicked.connect(lambda _, r=row: self.remove_dataset(r))
        options_layout.addWidget(remove_button)
        
        options_widget.setLayout(options_layout)
        self.table.setCellWidget(row, 5, options_widget)
        
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
            item = self.dataset_items[row]
            is_selected = (state == Qt.Checked)
            item.selected = is_selected
            
            print(f"[DEBUG] on_selected_changed: row={row}, item={item.loci_name}, selected={is_selected}")
            
            # 同步到Dataset Selection Manager
            if self.workspace and self.workspace.dataset_selection_manager and hasattr(self.workspace, 'current_dataset_id'):
                dataset_id = self.workspace.current_dataset_id
                print(f"[DEBUG] dataset_id={dataset_id}")
                if dataset_id:
                    dataset = self.workspace.dataset_selection_manager.get_dataset(dataset_id)
                    print(f"[DEBUG] dataset={dataset}")
                    if dataset:
                        # 检查是否有任何item被选中
                        has_selected = any(i.selected for i in self.dataset_items)
                        if has_selected:
                            dataset.selection_state = SELECTION_STATE_GREEN
                        else:
                            dataset.selection_state = SELECTION_STATE_NONE
                        print(f"[DEBUG] Dataset selection_state updated to: {dataset.selection_state}")
                    
                    # 查找对应的DatasetItem
                    print(f"[DEBUG] Searching for DatasetItem with loci_name={item.loci_name}, dataset_id={dataset_id}")
                    print(f"[DEBUG] Total items in manager: {len(self.workspace.dataset_selection_manager.items)}")
                    
                    found = False
                    for ds_item in self.workspace.dataset_selection_manager.items.values():
                        print(f"[DEBUG] Checking item: {ds_item.loci_name} (dataset_id={ds_item.dataset_id})")
                        if ds_item.loci_name == item.loci_name and ds_item.dataset_id == dataset_id:
                            if is_selected:
                                ds_item.selection_state = SELECTION_STATE_GREEN
                                self.workspace.dataset_selection_manager.selected_items.add(ds_item.id)
                                print(f"[DEBUG] Added to selected_items: {ds_item.loci_name} ({ds_item.id})")
                            else:
                                ds_item.selection_state = SELECTION_STATE_NONE
                                if ds_item.id in self.workspace.dataset_selection_manager.selected_items:
                                    self.workspace.dataset_selection_manager.selected_items.remove(ds_item.id)
                                print(f"[DEBUG] Removed from selected_items: {ds_item.loci_name} ({ds_item.id})")
                            found = True
                            break
                    
                    if not found:
                        print(f"[DEBUG] WARNING: No matching DatasetItem found for {item.loci_name}!")
                        print(f"[DEBUG] This explains why data can't be imported until the dialog is closed and reopened.")
            
    def view_dataset(self, row: int):
        """查看数据集"""
        try:
            if 0 <= row < len(self.dataset_items):
                dataset = self.dataset_items[row]
                # 确保dataset有name属性
                if not dataset.name:
                    dataset.name = dataset.loci_name
                
                # 转换SeqRecord对象为SequenceAlignmentViewer期望的字典格式
                sequences_for_viewer = []
                for seq_record in dataset.sequences:
                    seq_dict = {
                        'header': getattr(seq_record, 'id', getattr(seq_record, 'name', 'Unknown')),
                        'sequence': str(seq_record.seq) if hasattr(seq_record, 'seq') else str(seq_record)
                    }
                    sequences_for_viewer.append(seq_dict)
                
                # 使用PluginFactory获取序列查看器
                from YR_MPE.sequence_editor import SequenceAlignmentViewer
                viewer = SequenceAlignmentViewer(sequences_for_viewer)
                viewer.show()
                
                # 保存viewer引用以防被垃圾回收
                if not hasattr(self, 'viewers'):
                    self.viewers = []
                self.viewers.append(viewer)
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

    def import_from_partitioned_nexus(self):
        """从分区 NEXUS 文件导入数据"""
        # 选择导入文件
        import_path, _ = QFileDialog.getOpenFileName(
            self, "Select Partitioned NEXUS File", "", "NEXUS Files (*.nex *.nexus)"
        )
        if not import_path:
            return

        try:
            # 调用导入函数
            dataset_items, partition_scheme, summary = import_partitioned_nexus(import_path)

            # 询问用户是替换还是追加
            if self.dataset_items:
                reply = QMessageBox.question(
                    self, "Import Mode",
                    "Do you want to replace existing datasets or append to them?",
                    QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
                    QMessageBox.No  # 默认追加
                )

                if reply == QMessageBox.Yes:
                    # 替换现有数据
                    self.dataset_items = dataset_items
                    self.table.setRowCount(0)  # 清空表格
                    for item in dataset_items:
                        self.add_dataset_to_table(item)
                        
                        # 同时创建DatasetItem对象并添加到Dataset Selection Manager
                        if self.workspace and self.workspace.dataset_selection_manager and hasattr(self.workspace, 'current_dataset_id'):
                            dataset_id = self.workspace.current_dataset_id
                            if dataset_id:
                                from .dataset_models import DatasetItem as NewDatasetItem, ITEM_TYPE_ALIGNMENT
                                new_item = NewDatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
                                new_item.dataset_id = dataset_id
                                new_item.loci_name = item.loci_name
                                new_item.file_path = item.file_path
                                new_item.length = item.length
                                new_item.sequence_count = item.sequence_count
                                new_item.is_aligned = item.is_aligned
                                new_item.sequences = item.sequences
                                new_item.selection_state = SELECTION_STATE_NONE
                                
                                # 添加到manager
                                success = self.workspace.dataset_selection_manager.add_item(new_item, dataset_id)
                                print(f"[DEBUG] import_from_partitioned_nexus (replace): Created DatasetItem {item.loci_name} ({new_item.id}), success={success}")
                elif reply == QMessageBox.No:
                    # 追加到现有数据
                    for item in dataset_items:
                        self.dataset_items.append(item)
                        self.add_dataset_to_table(item)
                        
                        # 同时创建DatasetItem对象并添加到Dataset Selection Manager
                        if self.workspace and self.workspace.dataset_selection_manager and hasattr(self.workspace, 'current_dataset_id'):
                            dataset_id = self.workspace.current_dataset_id
                            if dataset_id:
                                from .dataset_models import DatasetItem as NewDatasetItem, ITEM_TYPE_ALIGNMENT
                                new_item = NewDatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
                                new_item.dataset_id = dataset_id
                                new_item.loci_name = item.loci_name
                                new_item.file_path = item.file_path
                                new_item.length = item.length
                                new_item.sequence_count = item.sequence_count
                                new_item.is_aligned = item.is_aligned
                                new_item.sequences = item.sequences
                                new_item.selection_state = SELECTION_STATE_NONE
                                
                                # 添加到manager
                                success = self.workspace.dataset_selection_manager.add_item(new_item, dataset_id)
                                print(f"[DEBUG] import_from_partitioned_nexus (append): Created DatasetItem {item.loci_name} ({new_item.id}), success={success}")
                else:
                    # 取消导入
                    return
            else:
                # 没有现有数据，直接添加
                self.dataset_items = dataset_items
                for item in dataset_items:
                    self.add_dataset_to_table(item)
                    
                    # 同时创建DatasetItem对象并添加到Dataset Selection Manager
                    if self.workspace and self.workspace.dataset_selection_manager and hasattr(self.workspace, 'current_dataset_id'):
                        dataset_id = self.workspace.current_dataset_id
                        if dataset_id:
                            from .dataset_models import DatasetItem as NewDatasetItem, ITEM_TYPE_ALIGNMENT
                            new_item = NewDatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
                            new_item.dataset_id = dataset_id
                            new_item.loci_name = item.loci_name
                            new_item.file_path = item.file_path
                            new_item.length = item.length
                            new_item.sequence_count = item.sequence_count
                            new_item.is_aligned = item.is_aligned
                            new_item.sequences = item.sequences
                            new_item.selection_state = SELECTION_STATE_NONE
                            
                            # 添加到manager
                            success = self.workspace.dataset_selection_manager.add_item(new_item, dataset_id)
                            print(f"[DEBUG] import_from_partitioned_nexus: Created DatasetItem {item.loci_name} ({new_item.id}), success={success}")

            # 显示成功消息
            message = (
                f"Successfully imported partitioned NEXUS file from {import_path}\n"
                f"Taxa: {summary['total_taxa']}, "
                f"Partitions: {summary['partition_count']}, "
                f"Total length: {summary['nchar']}"
            )
            QMessageBox.information(self, "Import Success", message)

            # 发送信号通知数据已更新
            self.dataset_processed.emit(self.dataset_items)

        except FileNotFoundError as e:
            QMessageBox.critical(self, "File Not Found", str(e))
        except ValueError as e:
            QMessageBox.critical(self, "Parse Error", str(e))
        except Exception as e:
            QMessageBox.critical(self, "Import Error", f"Failed to import NEXUS file: {str(e)}")

    def export_to_partitioned_nexus(self):
        """导出为分区NEXUS格式 - 使用改进的实现"""
        # 获取选中的datasets
        selected_items = [item for item in self.dataset_items if item.selected]
        if not selected_items:
            QMessageBox.warning(self, "Warning", "No loci selected for export.")
            return
            
        # 选择导出文件
        export_path, _ = QFileDialog.getSaveFileName(
            self, "Save Partitioned NEXUS File", "", "NEXUS Files (*.nex)"
        )
        if not export_path:
            return
            
        try:
            # 准备数据格式：转换为新函数需要的格式
            loci = []
            loci_names = []
            
            for item in selected_items:
                locus = []
                for seq_record in item.sequences:
                    # 获取序列名称和序列内容
                    seq_name = getattr(seq_record, 'id', getattr(seq_record, 'name', 'Unknown'))
                    seq_content = str(seq_record.seq) if hasattr(seq_record, 'seq') else str(seq_record)
                    locus.append([seq_name, seq_content])
                loci.append(locus)
                loci_names.append(item.loci_name)
            
            # 调用新的export_partitioned_nexus函数
            nexus_content, partition_scheme, missing_info = export_partitioned_nexus(loci, loci_names)
            
            # 写入NEXUS文件
            with open(export_path, 'w') as f:
                f.write(nexus_content)
                
            # 显示成功消息
            total_taxa = len(missing_info) if missing_info else 0
            total_loci = len(selected_items)
            QMessageBox.information(
                self, "Success", 
                f"Successfully exported partitioned NEXUS file to {export_path}\n"
                f"Taxa: {total_taxa}, Loci: {total_loci}"
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
            
    def batch_align_mafft(self):
        """批量执行MAFFT比对"""
        self._batch_align_with_tool("mafft")
    
    def batch_align_clustal_omega(self):
        """批量执行Clustal Omega比对"""
        self._batch_align_with_tool("clustal_omega")
    
    def batch_align_muscle5(self):
        """批量执行MUSCLE5比对"""
        self._batch_align_with_tool("muscle5")
    
    def batch_align_macse2(self):
        """批量执行MACSE2比对"""
        self._batch_align_with_tool("macse2")
    
    def batch_trim_gblocks(self):
        """批量执行GBlocks修剪"""
        self._batch_align_with_tool("gblocks")
    
    def batch_trim_trimal(self):
        """批量执行TrimAl修剪"""
        self._batch_align_with_tool("trimal")
        
    def _batch_align_with_tool(self, tool_name: str):
        """使用指定工具批量执行比对"""
        # 获取选中的数据集（不限制是否已比对，由用户决定）
        selected_items = [item for item in self.dataset_items if item.selected]
        if not selected_items:
            QMessageBox.warning(self, "Warning", f"No datasets selected for {tool_name} alignment.")
            return
            
        try:
            # 准备批量输入数据：[[SeqRecord, ...], [SeqRecord, ...]]
            batch_input_data = []
            for item in selected_items:
                batch_input_data.append(item.sequences)
                
            # 使用PluginFactory获取插件实例
            if tool_name == "mafft":
                plugin_instance = self.plugin_factory.get_mafft_plugin()
            elif tool_name == "clustal_omega":
                plugin_instance = self.plugin_factory.get_clustal_omega_plugin()
            elif tool_name == "muscle5":
                plugin_instance = self.plugin_factory.get_muscle5_plugin()
            elif tool_name == "trimal":
                plugin_instance = self.plugin_factory.get_trimal_plugin()
            elif tool_name == "gblocks":
                plugin_instance = self.plugin_factory.get_gblocks_plugin()
            elif tool_name == "macse2":
                plugin_instance = self.plugin_factory.get_macse_plugin()
            else:
                QMessageBox.critical(self, "Error", f"Plugin {tool_name} not available.")
                return
                
            # 创建插件对话框并传递批量数据
            dialog = QDialog()
            dialog.setWindowTitle(f"Batch {tool_name.upper()} Alignment")
            dialog.setMinimumSize(800, 600)
            dialog.setLayout(QVBoxLayout())
            
            # 运行插件（关键：传递批量数据，来源标识为DATASET_MANAGER）
            plugin_widget = plugin_instance.run(
                import_from="DATASET_MANAGER", 
                import_data=batch_input_data
            )
            
            # 连接批量信号（关键：连接到Dataset Manager自身的方法）
            if hasattr(plugin_widget, 'batch_import_alignment_signal'):
                plugin_widget.batch_import_alignment_signal.connect(
                    self.handle_batch_alignment_results
                )
                
            # 连接单文件信号（向后兼容，但批量场景不会触发）
            plugin_widget.import_alignment_signal.connect(
                self.handle_single_alignment_result
            )
            
            dialog.layout().addWidget(plugin_widget)
            dialog.exec_()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start batch alignment: {str(e)}")
    
    def handle_batch_alignment_results(self, batch_results):
        """处理批量比对结果（在Dataset Manager内部完成）"""
        if not batch_results or not isinstance(batch_results, list):
            return
            
        try:
            # 更新Dataset Items的状态
            selected_items = [item for item in self.dataset_items if item.selected]
            if len(batch_results) != len(selected_items):
                QMessageBox.warning(self, "Warning", "Result count mismatch!")
                return
                
            for i, (item, result_sequences) in enumerate(zip(selected_items, batch_results)):
                # 更新序列数据
                item.sequences = result_sequences
                item.is_aligned = True
                item.length = len(result_sequences[0].seq) if result_sequences else 0
                item.sequence_count = len(result_sequences)
                
            # 刷新表格显示
            self.refresh_table_display()
            
            # 注意：这里不再发送信号到主平台！
            # 用户需要手动选择"Export to Platform"来导出数据
            
            QMessageBox.information(
                self, "Success", 
                f"Successfully completed batch alignment for {len(batch_results)} datasets!"
            )
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to process batch results: {str(e)}")

    def handle_single_alignment_result(self, single_result):
        """处理单文件比对结果（向后兼容，但批量场景不会调用）"""
        pass
    
    def refresh_table_display(self):
        """刷新表格显示"""
        # 清空现有行
        self.table.setRowCount(0)
        
        # 重新填充数据
        for item in self.dataset_items:
            self.add_dataset_to_table(item)
            
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
            
    def closeEvent(self, event):
        """关闭事件处理 - 保存设置"""
        print("[DEBUG] DatasetManager.closeEvent called")
        if self.workspace and hasattr(self.workspace, 'current_dataset_id'):
            dataset_id = self.workspace.current_dataset_id
            print(f"[DEBUG] Saving settings for dataset_id: {dataset_id}")

            if dataset_id and self.workspace.dataset_selection_manager:
                dataset = self.workspace.dataset_selection_manager.get_dataset(dataset_id)

                if dataset:
                    print(f"[DEBUG] Dataset found: {dataset.name}")
                    print(f"[DEBUG] Current dataset_items count: {len(self.dataset_items)}")
                    
                    # 保存 topo_linked 和 edge_linked 设置
                    dataset.settings['topo_linked'] = self.topo_linked_radio.isChecked()
                    dataset.settings['edge_linked'] = self.edge_linked_radio.isChecked()

                    # 只有在有数据时才保存
                    if self.dataset_items:
                        # 保存到 settings（保留向后兼容）
                        dataset.settings['dataset_items'] = []
                        for item in self.dataset_items:
                            print(f"[DEBUG] Saving item: {item.loci_name}, is_aligned: {item.is_aligned}")
                            # 如果有文件路径且文件存在，只保存文件路径
                            # 否则保存序列数据
                            if item.file_path and os.path.exists(item.file_path):
                                item_data = {
                                    'selected': item.selected,
                                    'loci_name': item.loci_name,
                                    'length': item.length,
                                    'sequence_count': item.sequence_count,
                                    'is_aligned': item.is_aligned,
                                    'file_path': item.file_path,
                                    'sequences': None  # 从文件加载
                                }
                            else:
                                # 没有文件路径或文件不存在，保存序列数据
                                sequences_dict = []
                                if hasattr(item, 'sequences') and item.sequences:
                                    for seq in item.sequences:
                                        sequences_dict.append({
                                            'id': getattr(seq, 'id', getattr(seq, 'name', 'Unknown')),
                                            'description': getattr(seq, 'description', ''),
                                            'sequence': str(seq.seq) if hasattr(seq, 'seq') else str(seq)
                                        })

                                item_data = {
                                    'selected': item.selected,
                                    'loci_name': item.loci_name,
                                    'length': item.length,
                                    'sequence_count': item.sequence_count,
                                    'is_aligned': item.is_aligned,
                                    'file_path': item.file_path,
                                    'sequences': sequences_dict  # 保存序列数据
                                }

                            dataset.settings['dataset_items'].append(item_data)
                        
                        print(f"[DEBUG] Saved {len(dataset.settings['dataset_items'])} items to settings")
                    else:
                        print("[DEBUG] No dataset_items to save, preserving existing settings")

                    # 保存状态
                    if self.workspace.dataset_selection_manager.workspace_path:
                        self.workspace.dataset_selection_manager._save_state()
                        print("[DEBUG] State saved to file")
                else:
                    print("[DEBUG] Dataset not found")
            else:
                print("[DEBUG] No dataset_id or dataset_selection_manager")
        else:
            print("[DEBUG] No workspace or current_dataset_id")

        # 继续正常的关闭流程
        super().closeEvent(event)
