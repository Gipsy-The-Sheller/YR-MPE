#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Advanced Sequence Editor for YR-MPE
A modern, feature-rich sequence editor for molecular phylogenetics analysis
"""

import os
import re
import json
from typing import Dict, List, Optional, Tuple, Union
from collections import OrderedDict
from dataclasses import dataclass
from enum import Enum

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

try:
    from .sequence_models import SequenceModel, SequenceItem
    from .sequence_commands import *
    from .sequence_highlighter import SequenceHighlighter
    from .sequence_utils import SequenceUtils, SequenceValidator
except ImportError:
    # Fallback for direct execution
    from sequence_models import SequenceModel, SequenceItem
    from sequence_commands import *
    from sequence_highlighter import SequenceHighlighter
    from sequence_utils import SequenceUtils, SequenceValidator


class EditorMode(Enum):
    """编辑器模式"""
    NORMAL = "normal"
    ALIGNMENT = "alignment"
    CONSENSUS = "consensus"


@dataclass
class EditorSettings:
    """编辑器设置"""
    font_family: str = "Consolas"
    font_size: int = 12
    show_line_numbers: bool = True
    show_amino_acid_numbers: bool = True
    wrap_text: bool = False
    syntax_highlighting: bool = True
    auto_completion: bool = True
    show_consensus: bool = False
    color_scheme: str = "default"


class SequenceEditor(QMainWindow):
    """高级序列编辑器主窗口"""
    
    # 信号定义
    sequence_changed = pyqtSignal(str)  # 序列改变信号
    file_opened = pyqtSignal(str)       # 文件打开信号
    file_saved = pyqtSignal(str)        # 文件保存信号
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.current_file = None
        self.sequence_model = SequenceModel()
        self.settings = EditorSettings()
        self.undo_stack = QUndoStack()
        self.search_results = []
        self.current_search_index = -1
        
        self.init_ui()
        self.setup_connections()
        self.load_settings()
        
    def init_ui(self):
        """初始化用户界面"""
        self.setWindowTitle("YR-MPE Sequence Editor")
        self.setMinimumSize(1000, 700)
        
        # 创建中央部件
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # 主布局
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        
        # 创建工具栏
        self.create_toolbar()
        main_layout.addWidget(self.toolbar)
        
        # 创建状态栏
        self.create_status_bar()
        
        # 创建主编辑区域
        self.create_editor_area(main_layout)
        
        # 创建侧边栏
        self.create_sidebar(main_layout)
        
        # 创建菜单栏
        self.create_menu_bar()
        
    def create_toolbar(self):
        """创建工具栏"""
        self.toolbar = QToolBar("Main Toolbar")
        self.toolbar.setMovable(False)
        
        # 文件操作
        self.toolbar.addAction(self.create_action("Open", "open", "Open File", self.open_file))
        self.toolbar.addAction(self.create_action("Save", "save", "Save File", self.save_file))
        self.toolbar.addAction(self.create_action("Save As", "save_as", "Save As", self.save_as_file))
        self.toolbar.addSeparator()
        
        # 编辑操作
        self.toolbar.addAction(self.create_action("Undo", "undo", "Undo", self.undo_stack.undo))
        self.toolbar.addAction(self.create_action("Redo", "redo", "Redo", self.undo_stack.redo))
        self.toolbar.addSeparator()
        
        # 序列操作
        self.toolbar.addAction(self.create_action("Reverse", "reverse", "Reverse Sequence", self.reverse_sequence))
        self.toolbar.addAction(self.create_action("Complement", "complement", "Complement Sequence", self.complement_sequence))
        self.toolbar.addAction(self.create_action("Reverse Complement", "rev_comp", "Reverse Complement", self.reverse_complement))
        self.toolbar.addSeparator()
        
        # 搜索功能
        self.toolbar.addAction(self.create_action("Find", "find", "Find", self.show_find_dialog))
        self.toolbar.addAction(self.create_action("Replace", "replace", "Replace", self.show_replace_dialog))
        self.toolbar.addSeparator()
        
        # 视图选项
        self.toolbar.addAction(self.create_action("Settings", "settings", "Settings", self.show_settings))
        
    def create_status_bar(self):
        """创建状态栏"""
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        
        # 位置信息
        self.position_label = QLabel("Line: 1, Column: 1")
        self.status_bar.addWidget(self.position_label)
        
        # 序列信息
        self.sequence_info_label = QLabel("Sequences: 0, Length: 0")
        self.status_bar.addPermanentWidget(self.sequence_info_label)
        
        # 模式指示器
        self.mode_label = QLabel("Normal Mode")
        self.status_bar.addPermanentWidget(self.mode_label)
        
    def create_editor_area(self, parent_layout):
        """创建编辑区域"""
        # 创建分割器
        splitter = QSplitter(Qt.Horizontal)
        parent_layout.addWidget(splitter)
        
        # 左侧编辑区域
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        
        # 创建标签页控件
        self.tab_widget = QTabWidget()
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self.close_tab)
        self.tab_widget.currentChanged.connect(self.on_tab_changed)
        
        left_layout.addWidget(self.tab_widget)
        splitter.addWidget(left_widget)
        
        # 右侧信息面板
        self.info_panel = self.create_info_panel()
        splitter.addWidget(self.info_panel)
        
        # 设置分割器比例
        splitter.setSizes([800, 200])
        
    def create_info_panel(self):
        """创建信息面板"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        
        # 序列统计
        stats_group = QGroupBox("Sequence Statistics")
        stats_layout = QVBoxLayout(stats_group)
        
        self.stats_tree = QTreeWidget()
        self.stats_tree.setHeaderLabels(["Property", "Value"])
        stats_layout.addWidget(self.stats_tree)
        
        layout.addWidget(stats_group)
        
        # 序列质量
        quality_group = QGroupBox("Sequence Quality")
        quality_layout = QVBoxLayout(quality_group)
        
        self.quality_text = QTextEdit()
        self.quality_text.setMaximumHeight(150)
        self.quality_text.setReadOnly(True)
        quality_layout.addWidget(self.quality_text)
        
        layout.addWidget(quality_group)
        
        # 搜索历史
        search_group = QGroupBox("Search History")
        search_layout = QVBoxLayout(search_group)
        
        self.search_list = QListWidget()
        self.search_list.setMaximumHeight(100)
        search_layout.addWidget(self.search_list)
        
        layout.addWidget(search_group)
        
        return panel
        
    def create_sidebar(self, parent_layout):
        """创建侧边栏"""
        # 这里可以添加侧边栏功能，如文件浏览器、序列列表等
        pass
        
    def create_menu_bar(self):
        """创建菜单栏"""
        menubar = self.menuBar()
        
        # 文件菜单
        file_menu = menubar.addMenu("File")
        file_menu.addAction(self.create_action("New", "new", "New File", self.new_file))
        file_menu.addAction(self.create_action("Open", "open", "Open File", self.open_file))
        file_menu.addSeparator()
        file_menu.addAction(self.create_action("Save", "save", "Save File", self.save_file))
        file_menu.addAction(self.create_action("Save As", "save_as", "Save As", self.save_as_file))
        file_menu.addSeparator()
        file_menu.addAction(self.create_action("Export", "export", "Export", self.export_sequences))
        file_menu.addSeparator()
        file_menu.addAction(self.create_action("Exit", "exit", "Exit", self.close))
        
        # 编辑菜单
        edit_menu = menubar.addMenu("Edit")
        edit_menu.addAction(self.create_action("Undo", "undo", "Undo", self.undo_stack.undo))
        edit_menu.addAction(self.create_action("Redo", "redo", "Redo", self.undo_stack.redo))
        edit_menu.addSeparator()
        edit_menu.addAction(self.create_action("Cut", "cut", "Cut", self.cut_sequence))
        edit_menu.addAction(self.create_action("Copy", "copy", "Copy", self.copy_sequence))
        edit_menu.addAction(self.create_action("Paste", "paste", "Paste", self.paste_sequence))
        edit_menu.addSeparator()
        edit_menu.addAction(self.create_action("Find", "find", "Find", self.show_find_dialog))
        edit_menu.addAction(self.create_action("Replace", "replace", "Replace", self.show_replace_dialog))
        
        # 序列菜单
        sequence_menu = menubar.addMenu("Sequence")
        sequence_menu.addAction(self.create_action("Reverse", "reverse", "Reverse", self.reverse_sequence))
        sequence_menu.addAction(self.create_action("Complement", "complement", "Complement", self.complement_sequence))
        sequence_menu.addAction(self.create_action("Reverse Complement", "rev_comp", "Reverse Complement", self.reverse_complement))
        sequence_menu.addSeparator()
        sequence_menu.addAction(self.create_action("Translate", "translate", "Translate to Protein", self.translate_sequence))
        sequence_menu.addAction(self.create_action("Find ORFs", "orf", "Find Open Reading Frames", self.find_orfs))
        
        # 视图菜单
        view_menu = menubar.addMenu("View")
        view_menu.addAction(self.create_action("Settings", "settings", "Settings", self.show_settings))
        view_menu.addAction(self.create_action("Zoom In", "zoom_in", "Zoom In", self.zoom_in))
        view_menu.addAction(self.create_action("Zoom Out", "zoom_out", "Zoom Out", self.zoom_out))
        
        # 工具菜单
        tools_menu = menubar.addMenu("Tools")
        tools_menu.addAction(self.create_action("Align", "align", "Align Sequences", self.align_sequences))
        tools_menu.addAction(self.create_action("Consensus", "consensus", "Generate Consensus", self.generate_consensus))
        tools_menu.addAction(self.create_action("Statistics", "stats", "Sequence Statistics", self.show_statistics))
        
    def create_action(self, text, icon_name, tooltip, callback):
        """创建动作"""
        action = QAction(text, self)
        action.setToolTip(tooltip)
        action.triggered.connect(callback)
        return action
        
    def setup_connections(self):
        """设置信号连接"""
        self.undo_stack.canUndoChanged.connect(self.update_undo_action)
        self.undo_stack.canRedoChanged.connect(self.update_redo_action)
        
    def load_settings(self):
        """加载设置"""
        # 这里可以从配置文件加载设置
        pass
        
    def save_settings(self):
        """保存设置"""
        # 这里可以保存设置到配置文件
        pass
        
    # 文件操作方法
    def new_file(self):
        """新建文件"""
        self.create_new_tab("Untitled")
        
    def open_file(self):
        """打开文件"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Sequence File", "",
            "All Supported (*.fasta *.fas *.fa *.phy *.phylip *.nex *.nexus *.gb *.genbank);;"
            "FASTA (*.fasta *.fas *.fa);;"
            "Phylip (*.phy *.phylip);;"
            "Nexus (*.nex *.nexus);;"
            "GenBank (*.gb *.genbank)"
        )
        
        if file_path:
            self.load_file(file_path)
            
    def load_file(self, file_path):
        """加载文件"""
        try:
            sequences = SequenceUtils.load_sequences(file_path)
            if sequences:
                self.sequence_model.load_sequences(sequences)
                self.create_new_tab(os.path.basename(file_path), file_path)
                self.file_opened.emit(file_path)
                self.update_sequence_info()
            else:
                QMessageBox.warning(self, "Warning", "No sequences found in file.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file: {str(e)}")
            
    def save_file(self):
        """保存文件"""
        if self.current_file:
            self.save_file_to_path(self.current_file)
        else:
            self.save_as_file()
            
    def save_as_file(self):
        """另存为文件"""
        file_path, file_type = QFileDialog.getSaveFileName(
            self, "Save Sequence File", "",
            "FASTA (*.fasta);;"
            "Phylip (*.phy);;"
            "Nexus (*.nex);;"
            "All Files (*.*)"
        )
        
        if file_path:
            self.save_file_to_path(file_path)
            
    def save_file_to_path(self, file_path):
        """保存文件到指定路径"""
        try:
            sequences = self.sequence_model.get_sequences()
            SequenceUtils.save_sequences(sequences, file_path)
            self.current_file = file_path
            self.file_saved.emit(file_path)
            self.status_bar.showMessage(f"File saved: {file_path}", 3000)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save file: {str(e)}")
            
    def export_sequences(self):
        """导出序列"""
        # 实现序列导出功能
        pass
        
    # 编辑操作方法
    def cut_sequence(self):
        """剪切序列"""
        # 实现剪切功能
        pass
        
    def copy_sequence(self):
        """复制序列"""
        # 实现复制功能
        pass
        
    def paste_sequence(self):
        """粘贴序列"""
        # 实现粘贴功能
        pass
        
    # 序列操作方法
    def reverse_sequence(self):
        """反向序列"""
        # 实现反向功能
        pass
        
    def complement_sequence(self):
        """互补序列"""
        # 实现互补功能
        pass
        
    def reverse_complement(self):
        """反向互补序列"""
        # 实现反向互补功能
        pass
        
    def translate_sequence(self):
        """翻译序列"""
        # 实现翻译功能
        pass
        
    def find_orfs(self):
        """查找开放阅读框"""
        # 实现ORF查找功能
        pass
        
    # 搜索和替换方法
    def show_find_dialog(self):
        """显示查找对话框"""
        # 实现查找对话框
        pass
        
    def show_replace_dialog(self):
        """显示替换对话框"""
        # 实现替换对话框
        pass
        
    # 视图方法
    def show_settings(self):
        """显示设置对话框"""
        # 实现设置对话框
        pass
        
    def zoom_in(self):
        """放大"""
        # 实现放大功能
        pass
        
    def zoom_out(self):
        """缩小"""
        # 实现缩小功能
        pass
        
    # 工具方法
    def align_sequences(self):
        """对齐序列"""
        # 实现序列对齐功能
        pass
        
    def generate_consensus(self):
        """生成一致性序列"""
        # 实现一致性序列生成
        pass
        
    def show_statistics(self):
        """显示统计信息"""
        # 实现统计信息显示
        pass
        
    # 辅助方法
    def create_new_tab(self, title, file_path=None):
        """创建新标签页"""
        # 创建新的标签页内容
        tab_widget = QWidget()
        layout = QVBoxLayout(tab_widget)
        
        # 创建序列表格视图
        try:
            from .sequence_models import SequenceTableModel
        except ImportError:
            from sequence_models import SequenceTableModel
        table_model = SequenceTableModel(self.sequence_model.get_sequences())
        table_view = QTableView()
        
        table_view.setModel(table_model)
        table_view.setAlternatingRowColors(True)
        table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        
        # 设置表格样式
        table_view.setStyleSheet("""
            QTableView {
                gridline-color: #d0d0d0;
                background-color: white;
                alternate-background-color: #f8f8f8;
            }
            QTableView::item:selected {
                background-color: #3daee9;
                color: white;
            }
        """)

        layout.addWidget(table_view)
        
        # 添加到标签页
        tab_index = self.tab_widget.addTab(tab_widget, title)
        self.tab_widget.setCurrentIndex(tab_index)
        
        # 保存文件路径
        if file_path:
            self.tab_widget.setTabToolTip(tab_index, file_path)
            self.current_file = file_path
        
    def close_tab(self, index):
        """关闭标签页"""
        if index >= 0 and index < self.tab_widget.count():
            self.tab_widget.removeTab(index)
            if self.tab_widget.count() == 0:
                self.current_file = None
        
    def on_tab_changed(self, index):
        """标签页改变事件"""
        if index >= 0 and index < self.tab_widget.count():
            file_path = self.tab_widget.tabToolTip(index)
            if file_path:
                self.current_file = file_path
                self.setWindowTitle(f"YR-MPE Sequence Editor - {os.path.basename(file_path)}")
            self.update_sequence_info()
        
    def update_sequence_info(self):
        """更新序列信息"""
        sequences = self.sequence_model.get_sequences()
        if sequences:
            total_sequences = len(sequences)
            total_length = sum(seq.get_length() for seq in sequences)
            self.sequence_info_label.setText(f"Sequences: {total_sequences}, Total Length: {total_length}")
            
            # 更新统计信息
            self.stats_tree.clear()
            self.stats_tree.setHeaderLabels(["Property", "Value"])
            
            # 添加统计信息
            stats_items = [
                ("Total Sequences", str(total_sequences)),
                ("Total Length", str(total_length)),
                ("Average Length", f"{total_length // total_sequences if total_sequences > 0 else 0}"),
                ("Sequence Types", ", ".join(set(seq.sequence_type.value for seq in sequences)))
            ]
            
            for prop, value in stats_items:
                item = QTreeWidgetItem([prop, value])
                self.stats_tree.addTopLevelItem(item)
        else:
            self.sequence_info_label.setText("Sequences: 0, Total Length: 0")
            self.stats_tree.clear()
        
    def update_undo_action(self, can_undo):
        """更新撤销动作状态"""
        # 更新工具栏中的撤销按钮状态
        for action in self.toolbar.actions():
            if action.text() == "Undo":
                action.setEnabled(can_undo)
                break
        
    def update_redo_action(self, can_redo):
        """更新重做动作状态"""
        # 更新工具栏中的重做按钮状态
        for action in self.toolbar.actions():
            if action.text() == "Redo":
                action.setEnabled(can_redo)
                break
        
    def closeEvent(self, event):
        """关闭事件"""
        self.save_settings()
        event.accept()


class SequenceEditorEntry:
    """序列编辑器入口类"""
    def run(self):
        return SequenceEditor()


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    editor = SequenceEditor()
    editor.show()
    sys.exit(app.exec_())
