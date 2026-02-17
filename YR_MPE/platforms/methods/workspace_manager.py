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
                "coalescent": [],
                "chains": []
            }
        }
        self.plugin_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        self._workspace_widget = None
        self._workspace_hint = None
        self._main_layout = None
        
        # 工作区历史相关
        self.max_history = 20  # 支持至少20条历史记录
        self.workspace_history: List[str] = []
        self.settings_file = os.path.join(os.path.dirname(self.plugin_path), "settings.json")
        
        # 加载工作区历史
        self._load_workspace_history()
        
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
        
    def _initialize_workspace_structure(self):
        """初始化工作区目录结构

        Note:
            使用 exist_ok=True 避免重复创建导致的警告
        """
        if not self.workspace_path:
            return

        # 创建必要的子目录
        subdirs = ["history", "temp", "results"]
        for subdir in subdirs:
            dir_path = os.path.join(self.workspace_path, subdir)
            os.makedirs(dir_path, exist_ok=True)

            # 检查是否是新建的目录
            if not os.listdir(dir_path):
                print(f"Created new directory: {dir_path}")
            
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
        import shutil
        copied_files = []
        for file_path in file_paths:
            if os.path.exists(file_path):
                filename = os.path.basename(file_path)
                dest_path = os.path.join(history_dir, filename)
                try:
                    shutil.copy2(file_path, dest_path)
                    copied_files.append(dest_path)
                    print(f"Copied to history: {filename}")
                except Exception as e:
                    print(f"Failed to copy {filename} to history: {e}")
        
        if copied_files:
            print(f"Added {len(copied_files)} files to history: {operation_type}")
        else:
            print(f"No files added to history for operation: {operation_type}")
                
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
    
    def add_chain(self, chain_data: Any):
        """添加MCMC链文件到工作区"""
        self.current_state["current_items"]["chains"].append(chain_data)
        self.current_state["last_modified"] = datetime.now().isoformat()
        self._save_workspace_state()
    
    def get_current_chains(self) -> List[Any]:
        """获取当前MCMC链文件"""
        return self.current_state["current_items"]["chains"]
    
    def _load_workspace_history(self):
        """从配置文件加载工作区历史"""
        try:
            if os.path.exists(self.settings_file):
                with open(self.settings_file, 'r', encoding='utf-8') as f:
                    settings = json.load(f)
                    workspace_config = settings.get('workspace', {})
                    
                    # 加载历史记录
                    loaded_history = workspace_config.get('history', [])
                    # 过滤掉不存在的路径
                    self.workspace_history = [
                        path for path in loaded_history 
                        if path and os.path.exists(path)
                    ]
                    
                    # 注释掉自动加载当前工作区，确保程序启动时进入临时模式
                    # 用户需要手动选择工作区才会加载
                    # current_workspace = workspace_config.get('current', None)
                    # if current_workspace and os.path.exists(current_workspace):
                    #     self.workspace_path = current_workspace
                    #     self._load_workspace_state()
                    
                    # 更新max_history（如果配置中有，且大于等于1）
                    config_max = workspace_config.get('max_history', 20)
                    if config_max >= 1:
                        self.max_history = config_max
                    
                    # 如果配置中的max_history小于默认值20，更新为20
                    if self.max_history < 20:
                        print(f"Updating max_history from {self.max_history} to 20")
                        self.max_history = 20
                        self._save_workspace_history()
        except json.JSONDecodeError as e:
            print(f"Failed to parse settings file: {e}")
            self.workspace_history = []
        except Exception as e:
            print(f"Failed to load workspace history: {e}")
            self.workspace_history = []
    
    def _save_workspace_history(self):
        """保存工作区历史到配置文件"""
        try:
            # 读取现有配置
            settings = {}
            if os.path.exists(self.settings_file):
                try:
                    with open(self.settings_file, 'r', encoding='utf-8') as f:
                        settings = json.load(f)
                except json.JSONDecodeError:
                    # 如果文件损坏，从空开始
                    settings = {}
            
            # 更新工作区配置
            if 'workspace' not in settings:
                settings['workspace'] = {}
            
            settings['workspace']['current'] = self.workspace_path
            settings['workspace']['history'] = self.workspace_history
            settings['workspace']['max_history'] = self.max_history
            settings['workspace']['last_updated'] = datetime.now().isoformat()
            
            # 保存配置
            with open(self.settings_file, 'w', encoding='utf-8') as f:
                json.dump(settings, f, indent=2, ensure_ascii=False)
                
            print(f"Workspace history saved: {len(self.workspace_history)} entries")
        except Exception as e:
            print(f"Failed to save workspace history: {e}")
    
    def add_to_history(self, workspace_path: str):
        """添加工作区到历史记录"""
        if not workspace_path or not os.path.exists(workspace_path):
            print(f"Invalid workspace path: {workspace_path}")
            return False
        
        # 如果已经在历史中，先移除（这样会移动到开头）
        if workspace_path in self.workspace_history:
            self.workspace_history.remove(workspace_path)
        
        # 添加到开头
        self.workspace_history.insert(0, workspace_path)
        
        # 限制历史记录数量
        if len(self.workspace_history) > self.max_history:
            removed_paths = self.workspace_history[self.max_history:]
            self.workspace_history = self.workspace_history[:self.max_history]
            print(f"Removed from history (limit reached): {removed_paths}")
        
        # 保存历史
        self._save_workspace_history()
        return True
    
    def remove_from_history(self, workspace_path: str):
        """从历史记录中移除工作区"""
        if workspace_path in self.workspace_history:
            self.workspace_history.remove(workspace_path)
            
            # 如果移除的是当前工作区，也要更新current
            if self.workspace_path == workspace_path:
                self.workspace_path = None
                
            self._save_workspace_history()
            return True
        return False
    
    def clear_history(self):
        """清空所有历史记录（保留当前工作区）"""
        self.workspace_history = []
        self._save_workspace_history()
        print("Workspace history cleared")
    
    def get_history_list(self) -> List[str]:
        """获取工作区历史列表（过滤不存在的路径）"""
        # 返回前过滤掉不存在的路径
        valid_history = [
            path for path in self.workspace_history 
            if path and os.path.exists(path)
        ]
        
        # 如果有无效路径，更新并保存
        if len(valid_history) != len(self.workspace_history):
            self.workspace_history = valid_history
            self._save_workspace_history()
        
        return valid_history.copy()
    
    def set_max_history(self, max_count: int):
        """设置最大历史记录数量"""
        if max_count > 0:
            self.max_history = max_count
            # 如果当前历史超过新限制，截断
            if len(self.workspace_history) > self.max_history:
                self.workspace_history = self.workspace_history[:self.max_history]
            self._save_workspace_history()
    
    def get_history_info(self) -> Dict[str, Any]:
        """获取历史记录信息"""
        valid_history = self.get_history_list()
        return {
            "total_count": len(self.workspace_history),
            "valid_count": len(valid_history),
            "max_count": self.max_history,
            "current_workspace": self.workspace_path,
            "history_list": valid_history,
            "last_updated": self._get_last_updated()
        }
    
    def _get_last_updated(self) -> Optional[str]:
        """获取历史记录最后更新时间"""
        try:
            if os.path.exists(self.settings_file):
                with open(self.settings_file, 'r', encoding='utf-8') as f:
                    settings = json.load(f)
                    workspace_config = settings.get('workspace', {})
                    return workspace_config.get('last_updated')
        except:
            pass
        return None
    
    def validate_history_paths(self) -> int:
        """验证历史记录中的路径，返回无效路径的数量"""
        invalid_count = 0
        valid_history = []
        
        for path in self.workspace_history:
            if path and os.path.exists(path):
                valid_history.append(path)
            else:
                invalid_count += 1
                print(f"Invalid workspace path in history: {path}")
        
        if invalid_count > 0:
            self.workspace_history = valid_history
            self._save_workspace_history()
            print(f"Removed {invalid_count} invalid paths from history")
        
        return invalid_count
    
    def get_workspace_info(self, workspace_path: str = None) -> Optional[Dict[str, Any]]:
        """获取指定工作区的信息"""
        target_path = workspace_path or self.workspace_path
        
        if not target_path or not os.path.exists(target_path):
            return None
        
        state_file = os.path.join(target_path, "workspace_state.json")
        
        if not os.path.exists(state_file):
            return {
                "path": target_path,
                "has_state": False,
                "created_at": None,
                "last_modified": None,
                "items_count": {}
            }
        
        try:
            with open(state_file, 'r', encoding='utf-8') as f:
                state = json.load(f)
                
            items = state.get("current_items", {})
            return {
                "path": target_path,
                "has_state": True,
                "created_at": state.get("created_at"),
                "last_modified": state.get("last_modified"),
                "version": state.get("version"),
                "items_count": {
                    "alignments": len(items.get("alignments", [])),
                    "models": len(items.get("models", [])),
                    "distances": len(items.get("distances", [])),
                    "phylogenies": len(items.get("phylogenies", [])),
                    "datasets": len(items.get("datasets", [])),
                }
            }
        except Exception as e:
            print(f"Failed to read workspace state: {e}")
            return None
    
    def truncate_path(self, path: str, max_length: int = 50) -> str:
        """截断路径显示，保留目录名和路径（截断中间部分）"""
        if not path:
            return ""
        
        if len(path) <= max_length:
            return path
        
        # 获取目录名
        dirname = os.path.basename(path)
        dirname_len = len(dirname)
        
        # 如果目录名本身就超过了最大长度，直接返回目录名
        if dirname_len > max_length - 5:
            return dirname[:max_length-3] + "..."
        
        # 保留目录名和部分路径
        available_for_path = max_length - dirname_len - 7  # 减去" [... "和目录名
        if available_for_path < 10:
            # 可用空间太少，只返回目录名
            return dirname
        
        # 截断路径开头部分
        path_prefix = path[:available_for_path]
        return f"{dirname} [{path_prefix}...]"
    
    def set_workspace_path(self, path: str):
        """设置工作区路径并初始化结构

        Args:
            path: 工作区路径

        Note:
            此方法只在必要时创建目录结构，避免重复创建
        """
        if not path:
            raise ValueError("Workspace path cannot be empty")

        # 验证路径
        if not os.path.exists(path):
            try:
                os.makedirs(path, exist_ok=True)
            except Exception as e:
                raise ValueError(f"Failed to create workspace directory: {e}")

        self.workspace_path = path

        # 只在目录不存在时才初始化结构
        self._initialize_workspace_structure()

        # 如果是新的工作区（没有状态文件），设置创建时间
        state_file = os.path.join(path, "workspace_state.json")
        if not os.path.exists(state_file):
            self.current_state["created_at"] = datetime.now().isoformat()

        self.current_state["last_modified"] = datetime.now().isoformat()
        self._save_workspace_state()

        # 添加到历史记录（自动持久化）
        self.add_to_history(path)

        print(f"Workspace set: {path}")
    
    def switch_workspace(self, new_path: str, overwrite: bool = False):
        """
        切换工作区
        
        Args:
            new_path: 新工作区路径
            overwrite: 是否用当前数据覆盖新工作区
        """
        if not new_path:
            raise ValueError("New workspace path cannot be empty")
        
        # 如果切换到同一个工作区，不做任何操作
        if new_path == self.workspace_path:
            print(f"Already in workspace: {new_path}")
            return
        
        # 验证新路径
        if not os.path.exists(new_path):
            try:
                os.makedirs(new_path, exist_ok=True)
            except Exception as e:
                raise ValueError(f"Failed to create new workspace directory: {e}")
        
        old_path = self.workspace_path
        
        if overwrite and old_path:
            # 覆盖模式：将当前数据保存到新工作区
            self._save_current_data_to_workspace(new_path)
        else:
            # 非覆盖模式：清空当前数据，加载新工作区数据
            self.clear_current_state()
            self.workspace_path = new_path
            self._load_workspace_state()
        
        # 添加到历史记录
        self.add_to_history(new_path)
        
        print(f"Switched workspace from {old_path} to {new_path}")
    
    def _save_current_data_to_workspace(self, target_path: str):
        """保存当前数据到指定工作区"""
        if not target_path:
            return
        
        # 确保目标工作区目录存在
        os.makedirs(target_path, exist_ok=True)
        
        # 创建子目录
        subdirs = ["history", "temp", "results"]
        for subdir in subdirs:
            os.makedirs(os.path.join(target_path, subdir), exist_ok=True)
        
        # 临时保存当前工作区状态到新工作区
        temp_path = self.workspace_path
        self.workspace_path = target_path
        
        # 更新修改时间
        self.current_state["last_modified"] = datetime.now().isoformat()
        self._save_workspace_state()
        
        # 恢复原路径（因为set_workspace_path会再次调用）
        self.workspace_path = temp_path
    
    def _load_workspace_state(self):
        """从工作区加载状态"""
        if not self.workspace_path:
            return
            
        state_file = os.path.join(self.workspace_path, "workspace_state.json")
        if os.path.exists(state_file):
            try:
                with open(state_file, 'r', encoding='utf-8') as f:
                    self.current_state = json.load(f)
            except Exception as e:
                print(f"Failed to load workspace state: {e}")
    
    def has_data(self) -> bool:
        """检查工作区是否有数据"""
        if not self.current_state:
            return False
        
        items = self.current_state.get("current_items", {})
        for item_type, item_list in items.items():
            if item_list:  # 如果有任何非空列表
                return True
        
        return False
    
    def clear_current_state(self):
        """清空当前工作区状态"""
        self.current_state = {
            "created_at": None,
            "last_modified": None,
            "version": "1.0",
            "current_items": {
                "alignments": [],
                "models": [],
                "distances": [],
                "phylogenies": [],
                "variants": [],
                "coalescent": [],
                "chains": []
            }
        }
        if self.workspace_path:
            self._save_workspace_state()