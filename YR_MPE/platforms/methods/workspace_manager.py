import os
import json
from datetime import datetime
from typing import Dict, List, Any, Optional
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QLabel, QToolButton, 
                            QGridLayout, QMessageBox)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon

class WorkspaceManager:
    """
    工作区管理器 - 负责管理当前工作区状态、历史记录和文件组织
    """
    
    def __init__(self):
        self.workspace_path: Optional[str] = None
        self.current_state: Dict[str, Any] = {
            "created_at": None,
            "last_modified": None,
            "version": "1.0",
            "current_items": {
                "alignments": [],
                "models": [],
                "distances": [],
                "phylogenies": [],
                "variants": [],
                "coalescent": []
            }
        }
        self.plugin_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self._workspace_widget = None
        self._workspace_hint = None
        self._main_layout = None
        
    def _init_workspace_widget(self):
        """初始化工作区UI组件（延迟初始化）"""
        if self._workspace_widget is None:
            self._workspace_widget = QWidget()
            self._main_layout = QVBoxLayout()
            self._workspace_widget.setLayout(self._main_layout)
            
            # 添加提示标签
            self._workspace_hint = QLabel("Single Gene Workspace\nAdd an alignment or drag and drop a file here to start")
            self._workspace_hint.setAlignment(Qt.AlignCenter)
            self._workspace_hint.setStyleSheet("color: #555555;")
            self._main_layout.addWidget(self._workspace_hint)
            
            # 添加Dataset和Single Sequence按钮
            dataset_button = QToolButton()
            dataset_button.setText("Dataset")
            dataset_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/dataset.svg")))
            dataset_button.setCheckable(True)  # 设置为可选中状态
            dataset_button.setChecked(False)  # 初始未选中
            dataset_button.clicked.connect(self._on_dataset_button_clicked)
            self._main_layout.addWidget(dataset_button)
            
            single_sequence_button = QToolButton()
            single_sequence_button.setText("Single Sequence")
            single_sequence_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/single_sequence.svg")))
            single_sequence_button.setCheckable(True)  # 设置为可选中状态
            single_sequence_button.setChecked(True)  # 初始选中
            single_sequence_button.clicked.connect(self._on_single_sequence_button_clicked)
            self._main_layout.addWidget(single_sequence_button)
            
            # 存储按钮引用以便后续管理
            self.dataset_button = dataset_button
            self.single_sequence_button = single_sequence_button
            
            # 初始化数据项列表
            self.datasets = []
            self.single_sequences = []
            
    def _on_dataset_button_clicked(self):
        """Dataset按钮点击事件处理"""
        # 如果是单击，则切换选中状态
        if not self.dataset_button.isChecked():
            # 单击：选中Dataset功能模式
            self.dataset_button.setChecked(True)
            self.single_sequence_button.setChecked(False)
            
            # 打开Dataset管理器
            from .dataset_manager import DatasetManagerDialog
            dialog = DatasetManagerDialog()
            dialog.exec_()
            
            # 保存创建的dataset
            if hasattr(dialog, 'created_datasets'):
                self.datasets.extend(dialog.created_datasets)
                
    def _on_single_sequence_button_clicked(self):
        """Single Sequence按钮点击事件处理"""
        # 如果是单击，则切换选中状态
        if not self.single_sequence_button.isChecked():
            # 单击：选中Single Sequence功能模式
            self.single_sequence_button.setChecked(True)
            self.dataset_button.setChecked(False)
            
            # 可以在这里添加单序列相关的功能
            pass
            
    def get_workspace_widget(self) -> QWidget:
        """获取工作区UI组件"""
        self._init_workspace_widget()
        return self._workspace_widget
        
    def set_workspace_path(self, path: str):
        """设置工作区路径并初始化结构"""
        self.workspace_path = path
        self._initialize_workspace_structure()
        self.current_state["created_at"] = datetime.now().isoformat()
        self.current_state["last_modified"] = datetime.now().isoformat()
        self._save_workspace_state()
        
    def _initialize_workspace_structure(self):
        """初始化工作区目录结构"""
        if not self.workspace_path:
            return
            
        # 创建必要的子目录
        subdirs = ["history", "temp", "results"]
        for subdir in subdirs:
            os.makedirs(os.path.join(self.workspace_path, subdir), exist_ok=True)
            
    def _save_workspace_state(self):
        """保存工作区状态到JSON文件"""
        if not self.workspace_path:
            return
            
        state_file = os.path.join(self.workspace_path, "workspace_state.json")
        with open(state_file, 'w', encoding='utf-8') as f:
            json.dump(self.current_state, f, indent=2, ensure_ascii=False)
            
    def add_files_to_history(self, file_paths: List[str], operation_type: str = "unknown"):
        """将文件添加到历史记录"""
        if not self.workspace_path or not file_paths:
            return
            
        # 创建时间戳目录
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        history_dir = os.path.join(self.workspace_path, "history", f"{timestamp}_{operation_type}")
        os.makedirs(history_dir, exist_ok=True)
        
        # 复制文件到历史目录
        for file_path in file_paths:
            if os.path.exists(file_path):
                filename = os.path.basename(file_path)
                dest_path = os.path.join(history_dir, filename)
                # TODO: 实现文件复制逻辑
                
    def refresh_workspace_display(self):
        """刷新工作区显示"""
        # TODO: 实现工作区UI刷新逻辑
        pass
        
    def get_current_alignments(self) -> List[Any]:
        """获取当前比对结果"""
        return self.current_state["current_items"]["alignments"]
        
    def get_current_models(self) -> List[Any]:
        """获取当前模型结果"""
        return self.current_state["current_items"]["models"]
        
    def get_current_phylogenies(self) -> List[Any]:
        """获取当前系统发育树结果"""
        return self.current_state["current_items"]["phylogenies"]
        
    def add_alignment(self, alignment_data: Any):
        """添加比对结果到工作区"""
        self.current_state["current_items"]["alignments"].append(alignment_data)
        self.current_state["last_modified"] = datetime.now().isoformat()
        self._save_workspace_state()
        
    def add_model(self, model_data: Any):
        """添加模型结果到工作区"""
        self.current_state["current_items"]["models"].append(model_data)
        self.current_state["last_modified"] = datetime.now().isoformat()
        self._save_workspace_state()
        
    def add_phylogeny(self, phylogeny_data: Any):
        """添加系统发育树结果到工作区"""
        self.current_state["current_items"]["phylogenies"].append(phylogeny_data)
        self.current_state["last_modified"] = datetime.now().isoformat()
        self._save_workspace_state()