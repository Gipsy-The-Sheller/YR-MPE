# YR-MPEA
#
# Copyright (c) 2026 Zhi-Jie Xu
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

import os
import uuid
from typing import Optional, Dict, List, Any
from Bio import SeqIO
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout,
QMenuBar, QToolBar, QToolButton, QGroupBox, QLabel,
QAction, QMenu, QSizePolicy, QGridLayout, QFileDialog, QMessageBox, QDialog, QPushButton)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon

# 导入selection状态常量
from .methods.dataset_models import SELECTION_STATE_GREEN, SELECTION_STATE_NONE

# 导入新的模块架构
from .methods import (
    WorkspaceManager, PluginManager, PluginExecutor,
    FileOperations, UIHelpers
)
from .methods.import_partitioned_nexus import import_partitioned_nexus

# 导入数据集选择管理器和按钮组件
from .methods.dataset_selection_manager import DatasetSelectionManager
from .methods.dataset_item_button import DatasetItemButton, DatasetButton
from .methods.dataset_models import (
    DatasetItem, DatasetInfo,
    ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT, ITEM_TYPE_MODEL,
    ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY, ITEM_TYPE_TREESET,
    ITEM_TYPE_CHAIN, ITEM_TYPE_VARIANT, ITEM_TYPE_COALESCENT,
    ITEM_TYPE_CLOCK, SELECTION_STATE_NONE, SELECTION_STATE_GREEN,
    SELECTION_STATE_BLUE, SELECTION_STATE_RED
)

# 导入工厂模式
from ..factories import ResourceFactory, PluginFactory

class YR_MPEA_Widget(QWidget):
    def __init__(self):
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.mode = "single_gene"  # 添加缺失的mode属性
        
        # 初始化核心管理器
        self.workspace_manager = WorkspaceManager()
        self.plugin_manager = PluginManager(self.workspace_manager)
        self.plugin_executor = PluginExecutor(self.plugin_manager, self.workspace_manager)
        self.file_ops = FileOperations()
        self.ui_helpers = UIHelpers()
        
        # 初始化工厂
        self.resource_factory = ResourceFactory()
        self.plugin_factory = PluginFactory()
        
        # 初始化数据集选择管理器
        self.dataset_selection_manager = DatasetSelectionManager()
        
        # 注册所有插件
        self.plugin_manager.register_all_plugins()
        
        self.init_ui()
        
    def select_workspace_folder(self):
        """选择工作区文件夹"""
        folder_path = QFileDialog.getExistingDirectory(
            self, 
            "Select Workspace Folder", 
            os.path.expanduser("~")  # 默认打开用户主目录
        )
        if folder_path:
            self._switch_workspace(folder_path, is_new=True)
    
    def _switch_workspace(self, new_workspace_path: str, is_new: bool = False):
        """切换工作区"""
        try:
            # 检查当前工作区是否有数据
            current_has_data = self.workspace_manager.has_data() if self.workspace_manager.workspace_path else False
            
            if current_has_data and self.workspace_manager.workspace_path != new_workspace_path:
                # 显示切换选项对话框
                reply = QMessageBox.question(
                    self,
                    "Workspace Switch",
                    f"Current workspace contains data.\n\n"
                    f"Current: {self.workspace_manager.workspace_path}\n"
                    f"Target: {new_workspace_path}\n\n"
                    f"Choose an option:\n"
                    f"• Yes: Overwrite target workspace with current data\n"
                    f"• No: Clear current data and load target workspace\n"
                    f"• Cancel: Keep current workspace",
                    QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
                    QMessageBox.Cancel
                )
                
                if reply == QMessageBox.Cancel:
                    return  # 取消切换
                elif reply == QMessageBox.Yes:
                    # 选项1: 用当前数据覆盖目标工作区
                    self._save_current_data_to_workspace(new_workspace_path)
                elif reply == QMessageBox.No:
                    # 选项2: 清空当前数据并加载目标工作区
                    self.workspace_manager.clear_current_state()
                    self.workspace_manager.workspace_path = new_workspace_path
                    self.workspace_manager._load_workspace_state()
                    # 加载数据集选择管理器的状态
                    self.dataset_selection_manager.workspace_path = new_workspace_path
                    self.dataset_selection_manager._load_state()
            
            # 设置新工作区
            self.workspace_manager.set_workspace_path(new_workspace_path)
            
            # 设置数据集选择管理器的工作区路径
            self.dataset_selection_manager.set_workspace_path(new_workspace_path)
            
            # 加载数据集选择管理器的状态（如果存在的话）
            self.dataset_selection_manager._load_state()
            
            # 刷新UI
            self._update_workspace_display()
            
            # 显示成功消息
            display_name = self.workspace_manager.truncate_path(new_workspace_path)
            QMessageBox.information(
                self,
                "Workspace Selected",
                f"Workspace set to:\n{new_workspace_path}\n\n"
                f"History, temp, and results directories have been created."
            )
            
        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to switch workspace:\n{str(e)}"
            )
    
    def _save_current_data_to_workspace(self, workspace_path: str):
        """保存当前数据到指定工作区"""
        # 确保工作区目录存在
        os.makedirs(workspace_path, exist_ok=True)
        
        # 创建子目录
        subdirs = ["history", "temp", "results"]
        for subdir in subdirs:
            os.makedirs(os.path.join(workspace_path, subdir), exist_ok=True)
        
        # 保存当前工作区状态到新工作区
        temp_path = self.workspace_manager.workspace_path
        self.workspace_manager.workspace_path = workspace_path
        self.workspace_manager._save_workspace_state()
        
        # 保存数据集选择管理器的状态到新工作区
        temp_dataset_path = self.dataset_selection_manager.workspace_path
        self.dataset_selection_manager.workspace_path = workspace_path
        self.dataset_selection_manager._save_state()
        self.dataset_selection_manager.workspace_path = temp_dataset_path
        
        # 恢复原路径（因为set_workspace_path会再次调用）
        self.workspace_manager.workspace_path = temp_path
    
    def _update_workspace_display(self):
        """更新工作区显示"""
        if self.workspace_manager.workspace_path:
            display_name = self.workspace_manager.truncate_path(
                self.workspace_manager.workspace_path
            )
            # 更新workspace按钮的工具提示
            if hasattr(self, 'workspace_button'):
                self.workspace_button.setToolTip(f"Current Workspace: {self.workspace_manager.workspace_path}")
        else:
            # 临时模式
            if hasattr(self, 'workspace_button'):
                self.workspace_button.setToolTip("Temporary Mode - No workspace selected")
    
    def switch_to_temporary(self):
        """切换到临时模式（无工作区）"""
        try:
            # 检查当前是否有工作区且有数据
            current_has_data = self.workspace_manager.has_data() if self.workspace_manager.workspace_path else False

            if current_has_data:
                # 提示用户数据会丢失
                reply = QMessageBox.warning(
                    self,
                    "Switch to Temporary Mode",
                    "Switching to temporary mode will clear all current workspace data.\n\n"
                    "This action cannot be undone.\n\n"
                    "Do you want to continue?",
                    QMessageBox.Yes | QMessageBox.No,
                    QMessageBox.No
                )

                if reply == QMessageBox.No:
                    return  # 用户取消

            # 清空工作区路径和数据
            self.workspace_manager.workspace_path = None
            self.workspace_manager.clear_current_state()

            # 保存配置（清空current）
            self.workspace_manager._save_workspace_history()

            # 刷新UI
            self._update_workspace_display()
            
            # 刷新工作区菜单（更新Temporary选项状态）
            self._create_workspace_menu(self.workspace_button)
            
            # 显示成功消息
            QMessageBox.information(
                self,
                "Temporary Mode",
                "Switched to temporary mode.\n\n"
                "All workspace data has been cleared.\n"
                "Select a workspace folder to save your work."
            )
            
        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to switch to temporary mode:\n{str(e)}"
            )

    def get_workdir(self):
        """获取工作区 temp 目录路径

        Returns:
            str: 工作区 temp 目录路径，如果没有工作区则返回 None
        """
        if self.workspace_manager.workspace_path:
            return os.path.join(self.workspace_manager.workspace_path, "temp")
        return None

    def _prepare_import_data(self, workspace_item_types=None):
        """准备插件导入数据（抽象层，统一处理 dataset 和 SingleGeneWorkspace）
        
        Args:
            workspace_item_types: 工作区中的数据项类型列表，例如 ["alignments"]
                                   如果为 None，则不检查工作区数据
        
        Returns:
            tuple: (import_from, import_data)
                   import_from: "DATASET_MANAGER" 或 "YR_MPEA" 或 None
                   import_data: 导入数据的字典或单个数据项
        """
        import_from = None
        import_data = None
        
        # 1. 优先检查是否有选中的 dataset items
        if self.dataset_selection_manager:
            # 获取状态为 green 的 dataset
            from .methods.dataset_models import SELECTION_STATE_GREEN
            green_datasets = [ds for ds in self.dataset_selection_manager.get_all_datasets()
                             if ds.selection_state == SELECTION_STATE_GREEN]

            
            if green_datasets:
                # 获取第一个 green dataset
                green_dataset = green_datasets[0]

                # 获取该 dataset 中所有 selected 的 items（使用 selected_items 集合）
                selected_items = []
                for item_id in green_dataset.items:
                    # 检查 item_id 是否在 selected_items 集合中
                    if item_id in self.dataset_selection_manager.selected_items:
                        item = self.dataset_selection_manager.get_item(item_id)
                        if item:
                            selected_items.append(item)

                
                if selected_items:
                    import_from = "DATASET_MANAGER"
                    import_data = {
                        'dataset_items': selected_items,
                        'dataset_config': {
                            'topo_linked': green_dataset.settings.get('topo_linked', False),
                            'edge_linked': green_dataset.settings.get('edge_linked', False)
                        }
                    }
                    return import_from, import_data
            else:
                pass
        
        # 2. 如果没有选中的 dataset items，则检查 SingleGeneWorkspace
        if workspace_item_types:
            workspace_type = type(self.workspace).__name__
            
            if workspace_type == "SingleGeneWorkspace":
                for item_type in workspace_item_types:
                    if len(self.workspace.items.get(item_type, [])) >= 1:
                        import_from = "YR_MPEA"
                        import_data = self.workspace.items[item_type][0]
                        break
        return import_from, import_data
    
    def _get_dataset_config(self):
        """获取当前 dataset 的配置信息
        
        Returns:
            dict: 包含 topo_linked 和 edge_linked 的配置字典
        """
        if not self.dataset_selection_manager:
            return {}
        
        # 获取当前选中的 dataset
        selected_items = self.dataset_selection_manager.get_selected_items()
        if not selected_items:
            return {}
        
        # 获取 dataset ID
        dataset_id = selected_items[0].dataset_id
        if not dataset_id:
            return {}
        
        dataset = self.dataset_selection_manager.get_dataset(dataset_id)
        if not dataset:
            return {}
        
        # 返回配置信息
        return {
            'topo_linked': dataset.settings.get('topo_linked', False),
            'edge_linked': dataset.settings.get('edge_linked', False)
        }

    def _create_workspace_menu(self, workspace_button):
        """创建工作区菜单"""
        # 添加Select Workspace子菜单
        select_workspace_action = QAction("Select Workspace", workspace_button)
        workspace_menu = QMenu(workspace_button)
        
        # "Select Folder..." 选项
        select_folder_action = QAction("Select Folder...", workspace_menu)
        select_folder_action.triggered.connect(self.select_workspace_folder)
        workspace_menu.addAction(select_folder_action)
        
        # 添加分隔符
        separator = QAction(workspace_menu)
        separator.setSeparator(True)
        workspace_menu.addAction(separator)
        
        # 添加历史记录
        history_list = self.workspace_manager.get_history_list()
        if history_list:
            for i, ws_path in enumerate(history_list):
                display_name = self.workspace_manager.truncate_path(ws_path)
                
                # 标记当前工作区
                if ws_path == self.workspace_manager.workspace_path:
                    display_name += " (Current)"
                
                history_action = QAction(display_name, workspace_menu)
                history_action.setData(ws_path)
                history_action.triggered.connect(
                    lambda checked, path=ws_path: self._switch_workspace(path)
                )
                workspace_menu.addAction(history_action)
        else:
            # 如果没有历史记录，添加提示
            no_history_action = QAction("No recent workspaces", workspace_menu)
            no_history_action.setEnabled(False)
            workspace_menu.addAction(no_history_action)
        
        # 添加分隔符
        separator2 = QAction(workspace_menu)
        separator2.setSeparator(True)
        workspace_menu.addAction(separator2)
        
        # 添加Temporary选项
        temporary_action = QAction("Temporary", workspace_menu)
        # 如果已经在临时模式，标记为当前
        if not self.workspace_manager.workspace_path:
            temporary_action.setText("Temporary (Current)")
            temporary_action.setEnabled(False)
        else:
            temporary_action.triggered.connect(self.switch_to_temporary)
        workspace_menu.addAction(temporary_action)
        
        select_workspace_action.setMenu(workspace_menu)
        workspace_button.addAction(select_workspace_action)
            
    def init_ui(self):
        self.setWindowTitle("YR_MPEA")
        
        main_layout = QVBoxLayout()
        # no margins
        main_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(main_layout)
        # A window inspired by MEGA (Molecular Evolutionary Genetics Analysis)

        # 1. Main menu bar (On the top). Items: Align / Models / Distance / Phylogeny
        main_toolbar = QToolBar()
        main_layout.addWidget(main_toolbar, alignment = Qt.AlignTop)
        main_toolbar.setIconSize(QSize(45, 45))
        # vertically compact
        main_toolbar.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        # main_toolbar.layout.setSpacing(20)

        # 添加Workspace按钮（在所有其他按钮之前）
        workspace_button = QToolButton()
        workspace_button.setText("WORKSPACE")
        workspace_button.setIcon(self.resource_factory.get_icon("workspace.svg"))
        workspace_button.setPopupMode(QToolButton.InstantPopup)
        workspace_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        main_toolbar.addWidget(workspace_button)
        
        # 保存workspace_button引用
        self.workspace_button = workspace_button
        
        # 创建工作区菜单
        self._create_workspace_menu(workspace_button)
        
        # Dataset子菜单 - 使用正确的Qt API
        dataset_action = QAction("Dataset", workspace_button)
        dataset_menu = QMenu(workspace_button)
        
        new_dataset_action = QAction("New", dataset_menu)
        new_dataset_action.setIcon(self.resource_factory.get_icon("new.svg"))
        new_dataset_action.triggered.connect(self.create_new_dataset)
        dataset_menu.addAction(new_dataset_action)
        
        import_nexus_action = QAction("Import from partitioned NEXUS file", dataset_menu)
        import_nexus_action.setIcon(self.resource_factory.get_icon("partition.svg"))
        import_nexus_action.triggered.connect(self.import_dataset_from_nexus)
        dataset_menu.addAction(import_nexus_action)
        
        create_seqmatrix_action = QAction("Create by SeqMatrix", dataset_menu)
        create_seqmatrix_action.setIcon(self.resource_factory.get_icon("seqmatrix/seqmatrix.svg"))
        create_seqmatrix_action.triggered.connect(self.create_dataset_by_seqmatrix)
        dataset_menu.addAction(create_seqmatrix_action)
        
        create_seqdbg_action = QAction("Create by SeqDBG", dataset_menu)
        create_seqdbg_action.setIcon(self.resource_factory.get_icon("seqdbg/seqdbg.svg"))
        create_seqdbg_action.triggered.connect(self.create_dataset_by_seqdbg)
        dataset_menu.addAction(create_seqdbg_action)
        
        dataset_action.setMenu(dataset_menu)
        workspace_button.addAction(dataset_action)
        
        align_button = QToolButton()
        align_button.setText("ALIGN")
        # align_button.setStyleSheet("color: #555555;")
        align_button.setIcon(self.resource_factory.get_icon("align.svg"))
        align_button.setPopupMode(QToolButton.InstantPopup)
        align_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        # align_button.setFixedSize(60, 60)
        main_toolbar.addWidget(align_button)


        open_action = QAction("Open Sequence Files", align_button)
        open_action.setIcon(self.resource_factory.get_icon("open.svg"))
        open_action.triggered.connect(self.open_sequence_files)
        # build_action = QAction("New Alignment", align_button)
        # build_action.setIcon(self.resource_factory.get_icon("new.svg"))
        align_action = QAction("Align by...", align_button)
        align_action.setIcon(self.resource_factory.get_icon("alignby.svg"))
        trim_action = QAction("Trim Alignment by...", align_button)
        trim_action.setIcon(self.resource_factory.get_icon("trim.svg"))
        
        align_button.addAction(open_action)
        # align_button.addAction(build_action)
        align_button.addAction(align_action)
        align_button.addAction(trim_action)

        # aligners submenu for align_action
        aligners_menu = QMenu()
        # aligners_menu.addAction(QAction("Clustal Omega", align_button))
        # aligners_menu.addAction(QAction("MAFFT", align_button))
        clustal_omega_action = QAction("Clustal Omega", align_button)
        clustal_omega_action.setIcon(self.resource_factory.get_icon("software/clustalo.svg"))
        clustal_omega_action.triggered.connect(self.open_clustal_omega_wrapper)
        aligners_menu.addAction(clustal_omega_action)
        mafft_action = QAction("MAFFT", align_button)
        mafft_action.setIcon(self.resource_factory.get_icon("software/mafft.svg"))
        mafft_action.triggered.connect(self.open_mafft_wrapper)
        aligners_menu.addAction(mafft_action)
        muscle5_action = QAction("Muscle 5", align_button)
        muscle5_action.setIcon(self.resource_factory.get_icon("software/muscle.svg"))
        muscle5_action.triggered.connect(self.open_muscle5_wrapper)
        aligners_menu.addAction(muscle5_action)
        macse_action = QAction("MACSE", align_button)
        macse_action.setIcon(self.resource_factory.get_icon("software/macse.svg"))
        macse_action.triggered.connect(self.open_macse_wrapper)
        aligners_menu.addAction(macse_action)
        align_action.setMenu(aligners_menu)

        trim_menu = QMenu()
        trimal_action = QAction("Trimal", trim_action)
        trimal_action.setIcon(self.resource_factory.get_icon("software/trimal.svg"))
        trimal_action.triggered.connect(self.open_trimal_wrapper)
        trim_menu.addAction(trimal_action)
        gblocks_action = QAction("GBlocks", trim_action)
        gblocks_action.setIcon(self.resource_factory.get_icon("software/gblocks.svg"))
        gblocks_action.triggered.connect(self.open_gblocks_wrapper)
        trim_menu.addAction(gblocks_action)
        trim_action.setMenu(trim_menu)

        model_button = QToolButton()
        model_button.setText("MODEL")
        model_button.setIcon(self.resource_factory.get_icon("model.svg"))
        model_button.setPopupMode(QToolButton.InstantPopup)
        model_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        # model_button.setFixedSize(60, 60)
        main_toolbar.addWidget(model_button)

        models_menu = QMenu()
        find_best_model_action = QAction("Find Best Substitution Model (ML)", model_button)
        find_best_model_action.setIcon(self.resource_factory.get_icon("find_model.svg"))
        find_best_model_action.triggered.connect(self.open_modelfinder_wrapper)
        models_menu.addAction(find_best_model_action)
        
        estimate_model_params_action = QAction("Estimate Model Parameters (DNA, ML)", model_button)
        estimate_model_params_action.setIcon(self.resource_factory.get_icon("find_model.svg"))
        estimate_model_params_action.triggered.connect(self.open_model_parameter_wrapper)
        models_menu.addAction(estimate_model_params_action)
        
        model_button.setMenu(models_menu)

        distance_button = QToolButton()
        distance_button.setText("DISTANCE")
        distance_button.setIcon(self.resource_factory.get_icon("distance.svg"))
        distance_button.setPopupMode(QToolButton.InstantPopup)
        distance_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        # distance_button.setFixedSize(60, 60)
        main_toolbar.addWidget(distance_button)

        # 添加"Compute Pairwise Distances"动作
        self.comp_dist_action = QAction("Compute Pairwise Distances (ML)", self)
        self.comp_dist_action.setIcon(self.resource_factory.get_icon("dist.svg"))
        self.comp_dist_action.triggered.connect(self.open_ml_distance_wrapper)
        distance_button.addAction(self.comp_dist_action)
        
        # 添加"Compute Uncorrected Distances"动作
        self.comp_uncorr_dist_action = QAction("Compute Uncorrected Distances (p-distance)", self)
        self.comp_uncorr_dist_action.setIcon(self.resource_factory.get_icon("dist.svg"))
        self.comp_uncorr_dist_action.triggered.connect(self.open_uncorrected_distance_wrapper)
        distance_button.addAction(self.comp_uncorr_dist_action)
        
        separator = QAction(workspace_button)
        separator.setSeparator(True)
        distance_button.addAction(separator)
        
        # 添加"Calculate & Plot Substitution Saturation"动作
        self.saturation_action = QAction("Calculate && Plot Substitution Saturation", self)
        self.saturation_action.setIcon(self.resource_factory.get_icon("saturation.svg"))
        self.saturation_action.triggered.connect(self.open_substitution_saturation_wrapper)
        distance_button.addAction(self.saturation_action)
        
        # 添加"Compute Overall Mean Distances"动作
        # self.comp_mdist_action = QAction("Compute Overall Mean Distances (ML)", self)
        # self.comp_mdist_action.setIcon(self.resource_factory.get_icon("mdist.svg"))
        # TODO: 实现平均距离计算功能
        # self.comp_mdist_action.triggered.connect(self.compute_mean_distances)
        # distance_button.addAction(self.comp_mdist_action)

        phylogeny_button = QToolButton()
        phylogeny_button.setText("PHYLOGENY")
        phylogeny_button.setIcon(self.resource_factory.get_icon("phylogeny.svg"))
        phylogeny_button.setPopupMode(QToolButton.InstantPopup)
        phylogeny_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        phylogeny_button_menu = QMenu()
        phylogeny_button.setMenu(phylogeny_button_menu)
        # phylogeny_button.setFixedSize(60, 60)
        main_toolbar.addWidget(phylogeny_button)

        cons_ml_action = QAction("Maximum Likelihood (ML)", phylogeny_button)
        cons_ml_action.setIcon(self.resource_factory.get_icon("ml.svg"))

        ml_menu = QMenu()
        cons_ml_action.setMenu(ml_menu)
        
        # Add IQ-TREE action
        iqtree_action = QAction("IQ-Tree 3", phylogeny_button)
        iqtree_action.setIcon(self.resource_factory.get_icon("software/iqtree.svg"))
        iqtree_action.triggered.connect(self.open_iqtree_wrapper)

        ml_menu.addAction(iqtree_action)
        
        # TODO: ML Programs: IQ-TREE 3 / FastTree

        cons_bi_action = QAction("Bayesian Inference (BI)", phylogeny_button)
        cons_bi_action.setIcon(self.resource_factory.get_icon("bi.svg"))

        # 创建Bayesian Inference子菜单
        bi_menu = QMenu()
        cons_bi_action.setMenu(bi_menu)
        
        # 添加MrBayes选项
        mrbayes_action = QAction("MrBayes3-MPI-BEAGLE", phylogeny_button)
        mrbayes_action.setIcon(self.resource_factory.get_icon("software/mrbayes.svg"))
        mrbayes_action.triggered.connect(self.open_mrbayes_wrapper)
        bi_menu.addAction(mrbayes_action)
        
        # 添加PhyloBayes选项
        phylobayes_action = QAction("PhyloBayes-MPI", phylogeny_button)
        phylobayes_action.setIcon(self.resource_factory.get_icon("software/phylobayes.svg"))
        phylobayes_action.triggered.connect(self.open_phylobayes_wrapper)
        bi_menu.addAction(phylobayes_action)
        
        # 添加miniTracer选项
        minitracer_action = QAction("MCMC Diagnostics (MiniTracer)", phylogeny_button)
        minitracer_action.setIcon(self.resource_factory.get_icon("software/minitracer.svg"))
        minitracer_action.triggered.connect(self.open_minitracer_wrapper)
        bi_menu.addAction(minitracer_action)

        # TODO: BI Programs: MrBayes

        # cons_mp_action = QAction("Cladistics - Maximum Parsimony Phylogenies (MP)", phylogeny_button)
        # cons_mp_action.setIcon(self.resource_factory.get_icon("mp.svg"))

        # TODO: MP Programs: TNT

        # 添加Maximum Parsimony (MPBoot)选项
        cons_mp_action = QAction("Maximum Parsimony (MPBoot)", phylogeny_button)
        cons_mp_action.setIcon(self.resource_factory.get_icon("mp.svg"))
        cons_mp_action.triggered.connect(self.open_mpboot_wrapper)

        cons_distance_action = QAction("Distance Methods (DecentTree)", phylogeny_button)
        cons_distance_action.setIcon(self.resource_factory.get_icon("bionj.svg"))
        cons_distance_action.triggered.connect(self.open_decenttree_wrapper)

        # cons_rand_action = QAction("Simulation - Random Trees (Evolver)", phylogeny_button)
        # cons_rand_action.setIcon(self.resource_factory.get_icon("randtree.svg"))

        # cons_seq_action = QAction("Simulation - Simulated Sequences (Evolver)", phylogeny_button)
        # cons_seq_action.setIcon(self.resource_factory.get_icon("randseq.svg"))

        tree_viewer_action = QAction("Tree Viewer (IcyTree)", phylogeny_button)
        tree_viewer_action.setIcon(self.resource_factory.get_icon("software/icytree.svg"))
        tree_viewer_action.triggered.connect(self.open_icytree_wrapper)

        phylogeny_button_menu.addAction(cons_ml_action)
        phylogeny_button_menu.addAction(cons_bi_action)
        phylogeny_button_menu.addSeparator()
        # phylogeny_button_menu.addAction(cons_mp_action)
        phylogeny_button_menu.addAction(cons_mp_action)
        phylogeny_button_menu.addSeparator()
        phylogeny_button_menu.addAction(cons_distance_action)
        phylogeny_button_menu.addSeparator()
        # phylogeny_button_menu.addAction(cons_rand_action)
        # phylogeny_button_menu.addAction(cons_seq_action)
        phylogeny_button_menu.addSeparator()
        phylogeny_button_menu.addAction(tree_viewer_action)
        # phylogeny_button_menu.addSeparator()

        # variants_button = QToolButton()
        # variants_button.setText("VARIANTS")
        # variants_button.setIcon(self.resource_factory.get_icon("variants.svg"))
        # variants_button.setPopupMode(QToolButton.InstantPopup)
        # variants_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        # variants_button_menu = QMenu()
        # variants_button.setMenu(variants_button_menu)
        # # phylogeny_button.setFixedSize(60, 60)
        # main_toolbar.addWidget(variants_button)

        coalescent_button = QToolButton()
        coalescent_button.setText("COALESCENCE")
        coalescent_button.setIcon(self.resource_factory.get_icon("coalescent.svg"))
        coalescent_button.setPopupMode(QToolButton.InstantPopup)
        coalescent_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        coalescent_button_menu = QMenu()
        coalescent_button.setMenu(coalescent_button_menu)
        main_toolbar.addWidget(coalescent_button)

        # 添加CASTER-site功能到COALESCENT菜单
        caster_site_action = QAction("CASTER-site", coalescent_button)
        caster_site_action.setIcon(self.resource_factory.get_icon("software/caster.svg"))
        caster_site_action.triggered.connect(self.open_caster_site_wrapper)
        coalescent_button_menu.addAction(caster_site_action)

        # 添加ASTRAL功能到COALESCENT菜单
        astral_action = QAction("ASTRAL", coalescent_button)
        astral_action.setIcon(self.resource_factory.get_icon("software/astral.svg"))
        astral_action.triggered.connect(self.open_astral_wrapper)
        coalescent_button_menu.addAction(astral_action)

        clock_button = QToolButton()
        clock_button.setText("CLOCK")
        clock_button.setIcon(self.resource_factory.get_icon("clock.svg"))
        clock_button.setPopupMode(QToolButton.InstantPopup)
        clock_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        clock_button_menu = QMenu()
        clock_button.setMenu(clock_button_menu)
        main_toolbar.addWidget(clock_button)

        # TODO: Molecular Clock: LSD2 (IQ-TREE 3) / Bayesian Inference (MrBayes / PAML-mcmctree) / Penalized Likelihood (r8s)

        # Add PD-Guide action
        pdguide_action = QAction("PD-Guide", clock_button)
        pdguide_action.setIcon(self.resource_factory.get_icon("software/pdguide.svg"))
        pdguide_action.triggered.connect(self.open_pdguide_wrapper)
        clock_button_menu.addAction(pdguide_action)

        # mainworkspace_group = QGroupBox("Workspace")
        # # mainworkspace_layout = QGridLayout(10,4)
        # # # forbid stretch of rows and columns
        # # for i in range(10):
        # #     for j in range(4):
        # #         mainworkspace_layout.setColumnStretch(j, 0)
        # #         mainworkspace_layout.setRowStretch(i, 0)
        # # mainworkspace_group.setLayout(mainworkspace_layout)
        # mainworkspace_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # main_layout.addWidget(mainworkspace_group)

        # 创建SingleGeneWorkspace实例
        self.workspace = SingleGeneWorkspace(
            resource_factory=self.resource_factory,
            dataset_selection_manager=self.dataset_selection_manager,
            parent=self
        )
        self.workspace.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(self.workspace)

    def handle_dataset_click(self, button):
        """处理Dataset按钮点击事件"""
        if not button.doubleClicked:
            # 单击：选中Dataset功能模式
            button.setChecked(True)
    
    def open_dataset_manager(self):
        """打开Dataset管理对话框"""
        try:
            # 使用PluginFactory获取DatasetManager插件
            dataset_manager = self.plugin_factory.get_dataset_manager()
            
            # 保存引用防止被垃圾回收
            if not hasattr(self, 'dataset_managers'):
                self.dataset_managers = []
            self.dataset_managers.append(dataset_manager)
            
            # 设置必要的参数
            dataset_manager.dataset_name = "Dataset Manager"
            dataset_manager.plugin_factory = self.plugin_factory
            dataset_manager.workspace = self
            
            # 连接信号以便在Dataset管理器中查看序列
            if hasattr(dataset_manager, 'view_sequence_signal'):
                dataset_manager.view_sequence_signal.connect(self.view_sequence_in_viewer)
            
            dataset_manager.show()  # 改为非模态显示
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Dataset Manager: {str(e)}")
    
    def view_sequence_in_viewer(self, sequences):
        """使用SeqViewer插件查看序列"""
        try:
            from YR_MPE.sequence_editor import SequenceAlignmentViewer
            viewer = SequenceAlignmentViewer(sequences)
            viewer.show()
            # 保存viewer引用以防被垃圾回收
            if not hasattr(self, 'viewers'):
                self.viewers = []
            self.viewers.append(viewer)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open sequence viewer: {str(e)}")
    
    def open_muscle5_wrapper(self):
        # create a QDialog to open the muscle5_wrapper
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("MUSCLE5")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/muscle.svg"))
        dialog.setLayout(QVBoxLayout())
        
        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # use QDialog to open the muscle5_wrapper
        muscle5_entry = self.plugin_factory.get_muscle5_plugin()
        muscle5_wrapper = muscle5_entry.run(import_from=import_from, import_data=import_data)
        muscle5_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(muscle5_wrapper)
        dialog.exec_()

    
    def add_alignment_to_workspace(self, sequences):
        """将比对结果添加到工作区"""
        self.workspace.add_sequence(sequences)

    def add_model_to_workspace(self, model_data):
        """将模型选择结果添加到工作区"""
        self.workspace.add_model(model_data)

    def add_distance_matrix_to_workspace(self, distance_matrix):
        """将距离矩阵添加到工作区"""
        # Extract the data part from the signal
        if isinstance(distance_matrix, dict) and 'data' in distance_matrix:
            # Store only the distance matrix data, not the wrapper
            for dist_data in distance_matrix['data']:
                self.workspace.add_distance(dist_data)
        else:
            # Fallback: store as-is
            self.workspace.add_distance(distance_matrix)

    def add_phylogeny_to_workspace(self, phylogeny):
        """将树形图添加到工作区"""
        self.workspace.add_phylogeny(phylogeny)

        # 如果有文件路径，保存到 history
        if isinstance(phylogeny, dict) and 'data' in phylogeny:
            for tree_data in phylogeny['data']:
                if isinstance(tree_data, dict) and 'file_path' in tree_data:
                    file_path = tree_data['file_path']
                    if os.path.exists(file_path):
                        # 查找所有相关的输出文件（.treefile, .iqtree, .log 等）
                        output_files = []
                        base_path = os.path.splitext(file_path)[0]
                        for ext in ['.treefile', '.iqtree', '.log']:
                            if os.path.exists(base_path + ext):
                                output_files.append(base_path + ext)

                        # 添加到 history
                        if output_files:
                            self.workspace_manager.add_files_to_history(
                                output_files,
                                operation_type="iqtree"
                            )

    def add_chain_to_workspace(self, chain_item):
        """将MCMC链文件添加到工作区"""
        self.workspace.add_chain(chain_item)

    def open_caster_site_wrapper(self):
        """打开CASTER-site插件"""
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("CASTER-site - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/caster.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]

        # use PluginFactory to get the plugin
        caster_site_entry = self.plugin_factory.get_caster_site_plugin()
        caster_site_wrapper = caster_site_entry.run(import_from=import_from, import_data=import_data)
        # 连接信号，如果插件发出新序列或结果，添加到工作区
        # 注意：根据插件实际情况决定是否需要连接信号
        dialog.layout().addWidget(caster_site_wrapper)
        dialog.exec_()

    def open_astral_wrapper(self, import_data=None):
        """打开ASTRAL插件
        
        Args:
            import_data: 可选，包含多棵基因树的newick文件路径
                        如果提供，将直接使用这个文件作为ASTRAL的输入
        """
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("ASTRAL - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/astral.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # use PluginFactory to get the plugin
        astral_entry = self.plugin_factory.get_astral_plugin()
        
        # 如果提供了import_data（包含多棵基因树的newick文件），则直接使用
        if import_data and os.path.exists(import_data):
            import_from = "YR_MPEA"
            astral_wrapper = astral_entry.run(import_from=import_from, import_data=import_data)
        else:
            # 否则打开空的ASTRAL插件
            astral_wrapper = astral_entry.run()
        
        # 传递 dataset_selection_manager 给 ASTRAL 插件
        if hasattr(astral_wrapper, '_dataset_selection_manager') is False:
            astral_wrapper._dataset_selection_manager = self.dataset_selection_manager
        
        # 连接信号，如果插件发出新序列或结果，添加到工作区
        # 注意：根据插件实际情况决定是否需要连接信号
        dialog.layout().addWidget(astral_wrapper)
        dialog.exec_()

    def open_clustal_omega_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("Clustal Omega - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/clustalo.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # use PluginFactory to get the plugin
        clustal_omega_entry = self.plugin_factory.get_clustal_omega_plugin()
        clustalo_wrapper = clustal_omega_entry.run(import_from=import_from, import_data=import_data)
        clustalo_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(clustalo_wrapper)
        dialog.exec_()
    
    def open_mafft_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("MAFFT - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/mafft.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # use PluginFactory to get the plugin
        mafft_entry = self.plugin_factory.get_mafft_plugin()
        mafft_wrapper = mafft_entry.run(import_from=import_from, import_data=import_data)
        mafft_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(mafft_wrapper)
        dialog.exec_()
    
    def open_macse_wrapper(self):
        """打开MACSE插件"""
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("MACSE - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/macse.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # use PluginFactory to get the plugin
        macse_entry = self.plugin_factory.get_macse_plugin()
        macse_wrapper = macse_entry.run(import_from=import_from, import_data=import_data)
        macse_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(macse_wrapper)
        dialog.exec_()
    
    def open_modelfinder_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("ModelFinder - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("find_model.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])

        # use PluginFactory to get the plugin
        model_finder_entry = self.plugin_factory.get_model_finder_plugin()
        modelfinder_wrapper = model_finder_entry.run(import_from=import_from, import_data=import_data)
        modelfinder_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        modelfinder_wrapper.export_model_result_signal.connect(self.add_model_to_workspace)
        dialog.layout().addWidget(modelfinder_wrapper)
        dialog.exec_()
    
    def open_model_parameter_wrapper(self):
        """打开Model Parameter Estimation进行模型参数估计"""
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("Model Parameter Estimation - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("find_model.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # use PluginFactory to get the plugin
        model_parameter_entry = self.plugin_factory.get_model_parameter_plugin()
        model_parameter_wrapper = model_parameter_entry.run(import_from=import_from, import_data=import_data)
        model_parameter_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        model_parameter_wrapper.export_model_result_signal.connect(self.add_model_to_workspace)
        dialog.layout().addWidget(model_parameter_wrapper)
        dialog.exec_()
    
    def open_minitracer_wrapper(self, chain_item=None):
        """打开MiniTracer进行MCMC诊断分析
        
        Args:
            chain_item: 可选的ChainItem对象，如果提供则自动加载该chain文件
        """
        from PyQt5.QtWidgets import QDialog
        try:
            dialog = QDialog()
            dialog.setWindowTitle("MiniTracer - MCMC Diagnostics - YR-MPEA")
            dialog.setWindowIcon(self.resource_factory.get_icon("software/minitracer.svg"))
            # dialog.setMinimumSize(600, 500)

            dlayout = QVBoxLayout()
            dlayout.setContentsMargins(0, 0, 0, 0)
            dialog.setLayout(dlayout)
            
            # 使用PluginFactory获取MiniTracer插件
            plugin_entry = self.plugin_factory.get_minitracer_plugin()
            minitracer_wrapper = plugin_entry.run()
            
            # 如果提供了chain_item，自动加载该chain文件
            if chain_item:
                # 检查是旧的file_path还是新的file_paths
                if hasattr(chain_item, 'file_paths') and chain_item.file_paths:
                    # 新的ChainItem数据结构，支持多个链文件
                    for file_path in chain_item.file_paths:
                        if file_path not in minitracer_wrapper.mcmc_files:
                            minitracer_wrapper.mcmc_files.append(file_path)
                            try:
                                data = minitracer_wrapper.parse_trace_file(file_path)
                                minitracer_wrapper.mcmc_data[file_path] = data
                            except Exception as e:
                                QMessageBox.warning(dialog, "Error", f"Failed to load chain file {file_path}: {str(e)}")
                                if file_path in minitracer_wrapper.mcmc_files:
                                    minitracer_wrapper.mcmc_files.remove(file_path)
                    
                    # 更新界面
                    minitracer_wrapper.update_file_table()
                    # 选择第一个文件
                    if minitracer_wrapper.mcmc_files:
                        minitracer_wrapper.file_table.selectRow(0)
                        minitracer_wrapper.update_stats_table()
                        minitracer_wrapper.update_analysis_tabs()
                        
                elif hasattr(chain_item, 'file_path') and chain_item.file_path:
                    # 旧的ChainItem数据结构，单个链文件
                    file_path = chain_item.file_path
                    if file_path not in minitracer_wrapper.mcmc_files:
                        minitracer_wrapper.mcmc_files.append(file_path)
                        try:
                            data = minitracer_wrapper.parse_trace_file(file_path)
                            minitracer_wrapper.mcmc_data[file_path] = data
                            minitracer_wrapper.update_file_table()
                            # 自动选择刚添加的文件
                            file_index = minitracer_wrapper.mcmc_files.index(file_path)
                            minitracer_wrapper.file_table.selectRow(file_index)
                            minitracer_wrapper.update_stats_table()
                            minitracer_wrapper.update_analysis_tabs()
                        except Exception as e:
                            QMessageBox.warning(dialog, "Error", f"Failed to load chain file: {str(e)}")
                            if file_path in minitracer_wrapper.mcmc_files:
                                minitracer_wrapper.mcmc_files.remove(file_path)
            
            dialog.layout().addWidget(minitracer_wrapper)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open MiniTracer: {str(e)}")

    def open_icytree_wrapper(self, phylogeny_data=None):
        """打开IcyTree查看系统发育树
        
        Args:
            phylogeny_data: 可选的系统发育树数据字典，包含file_path等信息
        """
        from PyQt5.QtWidgets import QDialog
        try:
            dialog = QDialog()
            dialog.setWindowTitle("IcyTree - YR-MPEA")
            dialog.setWindowIcon(self.resource_factory.get_icon("software/icytree.svg"))
            dialog.setMinimumSize(800, 600)
            dialog.setLayout(QVBoxLayout())
            
            # 使用PluginFactory获取IcyTree插件
            plugin_entry = self.plugin_factory.get_icytree_plugin()
            icytree_wrapper = plugin_entry.run()
            
            # 如果提供了phylogeny_data，自动加载该树文件
            if phylogeny_data and isinstance(phylogeny_data, dict):
                tree_content = None
                
                # 方式1: 从 file_path 读取
                if 'file_path' in phylogeny_data:
                    try:
                        with open(phylogeny_data['file_path'], 'r') as f:
                            tree_content = f.read().strip()
                    except Exception as e:
                        QMessageBox.warning(dialog, "Error", f"Failed to load tree file: {str(e)}")
                
                # 方式2: 从 data[0]['content'] 读取（旧版platform_multigene兼容）
                elif 'data' in phylogeny_data and len(phylogeny_data['data']) > 0:
                    tree_data = phylogeny_data['data'][0]
                    if 'content' in tree_data:
                        tree_content = tree_data['content'].strip()
                
                # 方式3: 从直接 content 读取（TreeSet 多棵树）
                elif 'content' in phylogeny_data:
                    tree_content = phylogeny_data['content'].strip()
                
                # 设置树的content
                if tree_content:
                    icytree_wrapper.set_newick_string(tree_content)
            
            dialog.layout().addWidget(icytree_wrapper)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open IcyTree: {str(e)}")

    def open_trimal_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("TrimAl - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/trimal.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])

        # use PluginFactory to get the plugin
        trimal_entry = self.plugin_factory.get_trimal_plugin()
        trimal_wrapper = trimal_entry.run(import_from=import_from, import_data=import_data)
        trimal_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(trimal_wrapper)
        dialog.exec_()
        
    def open_gblocks_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("GBlocks - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/gblocks.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])

        # use PluginFactory to get the plugin
        gblocks_entry = self.plugin_factory.get_gblocks_plugin()
        gblocks_wrapper = gblocks_entry.run(import_from=import_from, import_data=import_data)
        gblocks_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(gblocks_wrapper)
        dialog.exec_()
    
    def open_sequence_files(self):
        file_dialog = QFileDialog()

        # supported formats: FASTA / Phylip / NEXUS
        # single select
        file_dialog.setFileMode(QFileDialog.ExistingFile)        
        file_dialog.setNameFilter("FASTA files (*.fas *.fasta *.fa *.fna);;Phylip files (*.phy);;NEXUS files (*.nex *.nexus)")
        file_types = {
            "fas": "fasta",
            "fasta": "fasta",
            "fa": "fasta",
            "fna": "fasta",
            "phy": "phylip",
            "nex": "nexus",
            "nexus": "nexus",
        }
        file_dialog.exec_()
        files = file_dialog.selectedFiles()
        if len(files) == 0:
            return
        file = files[0]
        try:
            # Use SeqIO.parse() to handle multiple sequences
            file_format = file_types[file.split(".")[-1].lower()]
            sequences = list(SeqIO.parse(file, file_format))
            
            if not sequences:
                QMessageBox.warning(self, "Error", "No sequences found in the file")
                return
            
            # Add all sequences to workspace
            self.workspace.add_sequence(sequences)
                
            # Show success message
            QMessageBox.information(self, "Success", f"Successfully loaded {len(sequences)} sequence(s) from {file}")
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error opening file: {e}")
            return
    
    def open_ml_distance_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        from .methods.dataset_models import ITEM_TYPE_MODEL, SELECTION_STATE_GREEN
        dialog = QDialog()
        dialog.setWindowTitle(f"Distance Calculator [implemented from IQ-TREE] - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("dist.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())
    
        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # 获取最佳模型（只检查green状态的model）
        best_model = ""
        # 从dataset_selection_manager获取green状态的model
        if self.dataset_selection_manager:
            green_models = [item for item in self.dataset_selection_manager.get_selected_items() 
                          if item.item_type == ITEM_TYPE_MODEL]
            if green_models:
                model_data = green_models[0].data
                # 处理模型表或单个模型
                if isinstance(model_data, dict):
                    if "type" in model_data and model_data["type"] == "model_table" and model_data["data"]:
                        # 取模型表中的第一个（最佳）模型
                        best_model = model_data["data"][0]['Model']
                    elif "Model" in model_data:
                        # 直接使用单个模型数据
                        best_model = model_data
    
        # Use PluginFactory to get the plugin
        ml_distance_entry = self.plugin_factory.get_ml_distance_plugin()
        plugin_wrapper = ml_distance_entry.run(import_from=import_from, import_data=import_data)
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_distance_result_signal.connect(self.add_distance_matrix_to_workspace)
    
        # parse model
        model_entries = best_model.split("+")

        model_noalias = ['JC69 (JC)', 'F81', 'K2P (K80)', 'HKY85 (HKY)', 'TNe', 'TN93 (TN)', 'K3P (K81)', 
                        'K81u', 'TPM2', 'TPM2u', 'TPM3', 'TPM3u', 'TIM', 'TIMe', 'TIM2', 'TIM2e', 'TIM3', 
                        'TIM3e', 'TVM', 'TVMe', 'SYM', 'GTR', 'Blosum62', 'cpREV', 
                        'Dayhoff', 'DCMut', 'EAL', 'ELM', 'FLAVI', 'FLU', 'GTR20', 'HIVb', 'HIVw', 'JTT', 
                        'JTTDCMut', 'LG', 'mtART', 'mtMAM', 'mtREV', 'mtZOA', 'mtMet', 'mtVer', 'mtInv', 
                        'NQ.bird', 'NQ.insect', 'NQ.mammal', 'NQ.pfam', 'NQ.plant', 'NQ.yeast', 'Poisson', 
                        'PMB', 'Q.bird', 'Q.insect', 'Q.mammal', 'Q.pfam', 'Q.plant', 'Q.yeast', 'rtREV', 'VT', 'WAG']
        model_alias = {'JC': 'JC69 (JC)', 'JC69': 'JC69 (JC)', 'K80': 'K2P (K80)', 'K2P': 'K2P (K80)', 
                       'TN': 'TN93 (TN)', 'TN93': 'TN93 (TN)', 'K81': 'K3P (K81)', 'K3P': 'K3P (K81)',
                       'K81u': 'K81u (K3Pu)', 'K3Pu': 'K81u (K3Pu)', 'HKY': 'HKY85 (HKY)', 'HKY85': 'HKY85 (HKY)'}
        if model_entries[0] in model_noalias:
            plugin_wrapper.model_combo.setCurrentText(model_entries[0])
        elif model_entries[0] in model_alias.keys():
            plugin_wrapper.model_combo.setCurrentText(model_alias[model_entries[0]])
        else:
            plugin_wrapper.model_combo.setCurrentText("auto")

        # Invariable sites?
        if "I" in model_entries:
            plugin_wrapper.invar_checkbox.setChecked(True)
        
        # empirical?
        if "F" in model_entries:
            # plugin_wrapper.empirical_checkbox.setChecked(True)
            plugin_wrapper.state_freq_combo.setCurrentText("Empirical (+F)")
        
        elif "FO" in model_entries:
            plugin_wrapper.state_freq_combo.setCurrentText("ML-optimized (+FO)")

        elif "FQ" in model_entries:
            plugin_wrapper.state_freq_combo.setCurrentText("Equal (+FQ)")
        
        # FreeRate?
        if "R" in model_entries:
            plugin_wrapper.freerate_checkbox.setChecked(True)
        
        if "ASC" in model_entries:
            plugin_wrapper.ascertain_checkbox.setChecked(True)
        
        # Gamma Caterories [identify Gx]
        for item in model_entries[1:]:
            if item.startswith("G"):
                plugin_wrapper.gamma_checkbox.setChecked(True)
                plugin_wrapper.gamma_spinbox.setValue(int(item[1:]))
    
        dialog.layout().addWidget(plugin_wrapper)
        dialog.exec_()        
    def open_uncorrected_distance_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle(f"Uncorrected Distance Calculator - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("dist.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])

        # Use PluginFactory to get the plugin
        p_distance_entry = self.plugin_factory.get_p_distance_plugin()
        plugin_wrapper = p_distance_entry.run(import_from=import_from, import_data=import_data)
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_distance_result_signal.connect(self.add_distance_matrix_to_workspace)

        dialog.layout().addWidget(plugin_wrapper)
        dialog.exec_()
    
    def open_substitution_saturation_wrapper(self):
        """打开替代饱和度分析插件"""
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("Substitution Saturation Analysis - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("saturation.svg"))
        dialog.setMinimumSize(900, 700)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])

        # Use PluginFactory to get the plugin
        saturation_entry = self.plugin_factory.get_substitution_saturation_plugin()
        plugin_wrapper = saturation_entry.run(import_from=import_from, import_data=import_data)
        
        dialog.layout().addWidget(plugin_wrapper)
        dialog.exec_()

    def _extract_partition_model_info(self, model_item):
        """
        提取分区模型信息

        Args:
            model_item: 分区模型 DatasetItem

        Returns:
            dict: {
                'partition_mode': str,  # EL/TL/EUL/TUL
                'partitions': list,    # 分区定义列表（包含模型信息）
                'model_count': int     # 分区数量
            }
        """
        # 获取分区模式
        partition_mode = model_item.partition_mode or "EL"

        # 获取分区定义
        model_data = model_item.data
        if isinstance(model_data, dict):
            partitions = model_data.get('partitions', [])
        else:
            partitions = []

        return {
            'partition_mode': partition_mode,
            'partitions': partitions,
            'model_count': len(partitions)
        }

    def open_iqtree_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        from .methods.dataset_models import ITEM_TYPE_MODEL, SELECTION_STATE_GREEN
        dialog = QDialog()
        dialog.setWindowTitle(f"IQ-TREE 3 - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/iqtree.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])

        # 检查是否有选中的分区模型
        partition_model_config = None
        if self.dataset_selection_manager:
            green_models = [item for item in self.dataset_selection_manager.get_selected_items()
                          if item.item_type == ITEM_TYPE_MODEL]
            if green_models:
                model_item = green_models[0]
                # 检查是否为分区模型
                if model_item.model_sub_type == "partitioned":
                    partition_model_config = self._extract_partition_model_info(model_item)

        # 获取最佳单一模型（如果没有分区模型）
        best_model = ""
        if not partition_model_config:
            # 从dataset_selection_manager获取green状态的model
            if self.dataset_selection_manager:
                green_models = [item for item in self.dataset_selection_manager.get_selected_items()
                              if item.item_type == ITEM_TYPE_MODEL]
                if green_models:
                    model_data = green_models[0].data
                    # 处理模型表或单个模型
                    if isinstance(model_data, dict):
                        if "type" in model_data and model_data["type"] == "model_table" and model_data["data"]:
                            # 取模型表中的第一个（最佳）模型
                            best_model = model_data["data"][0]['Model']
                        elif "Model" in model_data:
                            # 直接使用单个模型数据
                            best_model = model_data

        # Use PluginFactory to get the plugin entry
        iqtree_entry = self.plugin_factory.get_iqtree_plugin()

        # 获取工作目录
        workdir = self.get_workdir()

        plugin_wrapper = iqtree_entry.run(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir  # 添加工作目录参数
        )
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_model_result_signal.connect(self.add_model_to_workspace)
        plugin_wrapper.export_phylogeny_result_signal.connect(self.add_phylogeny_to_workspace)

        # 如果有分区模型配置，应用它
        if partition_model_config:
            plugin_wrapper.apply_partition_model_config(partition_model_config)
        else:
            # 没有分区模型，解析单一模型
            # parse model
            model_entries = best_model.split("+")

            model_noalias = ['JC69 (JC)', 'F81', 'K2P (K80)', 'HKY85 (HKY)', 'TNe', 'TN93 (TN)', 'K3P (K81)',
                            'K81u', 'TPM2', 'TPM2u', 'TPM3', 'TPM3u', 'TIM', 'TIMe', 'TIM2', 'TIM2e', 'TIM3',
                            'TIM3e', 'TVM', 'TVMe', 'SYM', 'GTR', 'Blosum62', 'cpREV',
                            'Dayhoff', 'DCMut', 'EAL', 'ELM', 'FLAVI', 'FLU', 'GTR20', 'HIVb', 'HIVw', 'JTT',
                            'JTTDCMut', 'LG', 'mtART', 'mtMAM', 'mtREV', 'mtZOA', 'mtMet', 'mtVer', 'mtInv',
                            'NQ.bird', 'NQ.insect', 'NQ.mammal', 'NQ.pfam', 'NQ.plant', 'NQ.yeast', 'Poisson',
                            'PMB', 'Q.bird', 'Q.insect', 'Q.mammal', 'Q.pfam', 'Q.plant', 'Q.yeast', 'rtREV', 'VT', 'WAG']
            model_alias = {'JC': 'JC69 (JC)', 'JC69': 'JC69 (JC)', 'K80': 'K2P (K80)', 'K2P': 'K2P (K80)',
                           'TN': 'TN93 (TN)', 'TN93': 'TN93 (TN)', 'K81': 'K3P (K81)', 'K3P': 'K3P (K81)',
                           'K81u': 'K81u (K3Pu)', 'K3Pu': 'K81u (K3Pu)', 'HKY': 'HKY85 (HKY)', 'HKY85': 'HKY85 (HKY)'}
            if model_entries[0] in model_noalias:
                plugin_wrapper.model_combo.setCurrentText(model_entries[0])
            elif model_entries[0] in model_alias.keys():
                plugin_wrapper.model_combo.setCurrentText(model_alias[model_entries[0]])
            else:
                plugin_wrapper.model_combo.setCurrentText("auto")

            # Invariable sites?
            if "I" in model_entries:
                plugin_wrapper.invar_checkbox.setChecked(True)

            # empirical?
            if "F" in model_entries:
                # plugin_wrapper.empirical_checkbox.setChecked(True)
                plugin_wrapper.state_freq_combo.setCurrentText("Empirical (+F)")

            elif "FO" in model_entries:
                plugin_wrapper.state_freq_combo.setCurrentText("ML-optimized (+FO)")

            elif "FQ" in model_entries:
                plugin_wrapper.state_freq_combo.setCurrentText("Equal (+FQ)")

            # Rate Heterogeneity?
            if "R" in model_entries:
                plugin_wrapper.rate_combo.setCurrentText("FreeRate Model (+R)")
            else:
                # 检查是否有 Gamma 参数
                for item in model_entries[1:]:
                    if item.startswith("G"):
                        plugin_wrapper.rate_combo.setCurrentText("Gamma Distribution (+G)")
                        plugin_wrapper.rate_categories_spinbox.setValue(int(item[1:]))
                        break

        dialog.layout().addWidget(plugin_wrapper)
        dialog.exec_()

    def open_decenttree_wrapper(self):
        """打开 DecentTree 插件"""
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("DecentTree Distance Methods - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("bionj.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        
        if workspace_type == "SingleGeneWorkspace":
            # 优先导入距离矩阵数据
            if len(self.workspace.items["distances"]) >= 1:
                import_from = "YR_MPEA"
                import_data = {
                    'type': 'distance_matrix',
                    'data': self.workspace.items["distances"]
                }
            # 如果没有距离矩阵，可以提示用户先计算距离
            else:
                QMessageBox.information(
                    self, 
                    "Distance Matrix Required",
                    "Please compute ML distances first using the DISTANCE menu before building a tree."
                )
                return

        # Use PluginFactory to get the plugin
        decenttree_entry = self.plugin_factory.get_decenttree_plugin()
        decenttree_wrapper = decenttree_entry.run(import_from=import_from, import_data=import_data)
        decenttree_wrapper.export_phylogeny_result_signal.connect(self.add_phylogeny_to_workspace)

        dialog.layout().addWidget(decenttree_wrapper)
        dialog.exec_()

    def open_mpboot_wrapper(self):
        """打开 MPBoot 插件"""
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("Maximum Parsimony (MPBoot) - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("mp.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data (使用抽象层，支持两种导入方式)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # Use PluginFactory to get the plugin
        mpboot_entry = self.plugin_factory.get_mpboot_plugin()
        
        # 获取工作目录
        workdir = self.get_workdir()
        
        mpboot_wrapper = mpboot_entry.run(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir  # 添加工作目录参数
        )
        
        mpboot_wrapper.export_phylogeny_result_signal.connect(self.add_phylogeny_to_workspace)
        
        dialog.layout().addWidget(mpboot_wrapper)
        dialog.exec_()

    def open_pdguide_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        from ..plugins.pdguide_plugin import PdGuidePluginEntry
        dialog = QDialog()
        dialog.setWindowTitle(f"PD-Guide - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/pdguide.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data
        import_from = None
        import_data = None
        
        if len(self.workspace.items["alignments"]) > 0:
            import_from = "alignment"
            import_data = self.workspace.items["alignments"][0]
            
        plugin_entry = PdGuidePluginEntry()
        plugin_wrapper = plugin_entry.run(import_from=import_from, import_data=import_data)
        dialog.layout().addWidget(plugin_wrapper)
        dialog.exec_()

    def open_seqdbg_wrapper(self):
        """Open SeqDBG plugin for annotation graph visualization"""
        from PyQt5.QtWidgets import QDialog
        from ..plugins.seqdbg_plugin import SeqDBGPluginEntry
        
        # Create SeqDBG plugin instance using the entry point
        plugin_entry = SeqDBGPluginEntry()
        plugin_window = plugin_entry.run(import_from="YR_MPEA")
        
        # Connect import signal to workspace
        plugin_window.import_dataset_signal.connect(self.add_seqdbg_export_to_workspace)
        
        # Show the window (it's a QMainWindow, not a dialog)
        plugin_window.show()
        
        # Keep reference to prevent garbage collection
        if not hasattr(self, 'seqdbg_windows'):
            self.seqdbg_windows = []
        self.seqdbg_windows.append(plugin_window)
    
    def add_seqdbg_export_to_workspace(self, dataset_items):
        """Add exported items from SeqDBG to workspace as a Dataset"""
        import logging
        from PyQt5.QtCore import QTimer
        logger = logging.getLogger(__name__)
        
        if not dataset_items:
            logger.warning("add_seqdbg_export_to_workspace: dataset_items is empty")
            return
        
        try:
            logger.info(f"add_seqdbg_export_to_workspace: received {len(dataset_items)} items")
            
            # 创建Dataset对象（按照SeqMatrix的方式）
            class DatasetObject:
                def __init__(self, name, items=None, partition_scheme=None):
                    self.dataset_name = name
                    self.items = items if items else []
                    self.partition_scheme = partition_scheme
                    self.summary = {}
            
            # 收集所有有效的dataset_items
            valid_items = []
            for i, item in enumerate(dataset_items):
                try:
                    logger.info(f"Processing item {i}: item_type={getattr(item, 'item_type', 'N/A')}, "
                               f"loci_name={getattr(item, 'loci_name', 'N/A')}, "
                               f"sequence_count={getattr(item, 'sequence_count', 'N/A')}")
                    
                    if item.item_type == "sequence" and item.sequences:
                        logger.info(f"Adding {len(item.sequences)} sequences to dataset")
                        valid_items.append(item)
                    else:
                        logger.warning(f"Item {i} skipped: no sequences or wrong type")
                except Exception as e:
                    logger.error(f"Error processing item {i}: {e}", exc_info=True)
            
            if not valid_items:
                logger.warning("No valid items to export")
                QTimer.singleShot(100, lambda: QMessageBox.warning(
                    self, "Export Warning",
                    "No valid sequences found to export."
                ))
                return
            
            # 创建Dataset
            dataset_name = "SeqDBG_Export"
            dataset = DatasetObject(
                name=dataset_name,
                items=valid_items
            )
            
            # 计算摘要信息
            total_sequences = sum(item.sequence_count for item in valid_items)
            total_loci = len(valid_items)
            dataset.summary = {
                'total_sequences': total_sequences,
                'total_loci': total_loci,
                'loci_names': [item.loci_name for item in valid_items]
            }
            
            # 添加到工作区
            self.workspace.add_dataset(dataset)
            
            logger.info(f"Successfully added dataset with {total_loci} loci and {total_sequences} sequences")
            
            # 延迟显示成功消息
            QTimer.singleShot(100, lambda: self._show_export_success_message(dataset))
            
        except Exception as e:
            logger.error(f"Error in add_seqdbg_export_to_workspace: {e}", exc_info=True)
            QTimer.singleShot(100, lambda: self._show_export_error_message(str(e)))
    
    def _show_export_success_message(self, dataset):
        import logging
        logger = logging.getLogger(__name__)        
        """显示导出成功消息"""
        try:
            summary = getattr(dataset, 'summary', {})
            message = f"Exported dataset from SeqDBG to workspace.\n\n"
            message += f"Dataset Name: {getattr(dataset, 'dataset_name', 'Unnamed')}\n"
            message += f"Loci: {summary.get('total_loci', 0)}\n"
            message += f"Total Sequences: {summary.get('total_sequences', 0)}\n"
            
            loci_names = summary.get('loci_names', [])
            if loci_names:
                message += f"Loci Names: {', '.join(loci_names[:5])}"
                if len(loci_names) > 5:
                    message += f" ... ({len(loci_names)} total)"
            
            message += "\n\nThe dataset is now available for alignment and phylogenetic analysis."
            
            QMessageBox.information(self, "Export Successful", message)
        except Exception as e:
            logger.error(f"Error showing success message: {e}", exc_info=True)
    
    def _show_export_error_message(self, error_msg):
        import logging
        logger = logging.getLogger(__name__)  
        """显示导出错误消息"""
        try:
            QMessageBox.critical(
                self,
                "Export Error",
                f"Failed to add SeqDBG export to workspace:\n{error_msg}"
            )
        except Exception as e:
            logger.error(f"Error showing error message: {e}", exc_info=True)
    
    def _show_export_error_message(self, error_msg):
        import logging
        logger = logging.getLogger(__name__)  
        """显示导出错误消息"""
        try:
            QMessageBox.critical(
                self,
                "Export Error",
                f"Failed to add SeqDBG export to workspace:\n{error_msg}"
            )
        except Exception as e:
            logger.error(f"Error showing error message: {e}", exc_info=True)

    def create_new_dataset(self):
        """创建新的Dataset - 只创建空对象和ToolButton，不打开对话框"""
        try:
            # 创建Dataset对象用于工作区显示
            class DatasetObject:
                def __init__(self, name):
                    self.dataset_name = name
                    self.items = []  # 用于存储实际的dataset items
            
            dataset = DatasetObject("New Dataset")
            
            # 添加到工作区（创建QToolButton）
            self.workspace.add_dataset(dataset)
                    
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create new dataset: {str(e)}")
    
    def import_dataset_from_nexus(self):
        """从分区NEXUS文件导入Dataset"""
        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        file_dialog.setNameFilter("NEXUS files (*.nex *.nexus)")

        if file_dialog.exec_():
            selected_files = file_dialog.selectedFiles()
            if selected_files:
                nexus_file = selected_files[0]
                try:
                    # 调用导入函数解析NEXUS文件
                    dataset_items, partition_scheme, summary = import_partitioned_nexus(nexus_file)

                    # 创建Dataset对象
                    class DatasetObject:
                        def __init__(self, name, items=None, partition_scheme=None):
                            self.dataset_name = name
                            self.items = items if items else []
                            self.partition_scheme = partition_scheme
                            self.summary = summary if summary else {}

                    # 使用文件名（不含扩展名）作为数据集名称
                    dataset_name = os.path.splitext(os.path.basename(nexus_file))[0]
                    dataset = DatasetObject(
                        name=dataset_name,
                        items=dataset_items,
                        partition_scheme=partition_scheme
                    )

                    # 添加到工作区
                    self.workspace.add_dataset(dataset)

                    # 显示成功消息
                    message = (
                        f"Successfully imported dataset from:\n{nexus_file}\n\n"
                        f"Dataset Name: {dataset_name}\n"
                        f"Taxa: {summary.get('ntax', 0)}\n"
                        f"Total Length: {summary.get('nchar', 0)}\n"
                        f"Partitions: {summary.get('partition_count', 0)}\n"
                        f"Partition Names: {', '.join(summary.get('partition_names', []))}"
                    )
                    QMessageBox.information(self, "Import Success", message)

                except FileNotFoundError as e:
                    QMessageBox.critical(self, "File Not Found", str(e))
                except ValueError as e:
                    QMessageBox.critical(self, "Parse Error", str(e))
                except Exception as e:
                    QMessageBox.critical(self, "Import Error", f"Failed to import dataset from NEXUS: {str(e)}")
    
    def create_dataset_by_seqmatrix(self):
        """通过SeqMatrix创建Dataset"""
        try:
            # 获取 SeqMatrix 插件
            seqmatrix_entry = self.plugin_factory.get_seqmatrix_plugin()
            seqmatrix_window = seqmatrix_entry.run(import_from="YR_MPEA")
            
            # 连接信号
            seqmatrix_window.import_dataset_signal.connect(self.handle_seqmatrix_import)
            
            # 保存引用防止被垃圾回收
            if not hasattr(self, 'seqmatrix_windows'):
                self.seqmatrix_windows = []
            self.seqmatrix_windows.append(seqmatrix_window)
            
            # 显示窗口
            seqmatrix_window.show()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open SeqMatrix: {str(e)}")
    
    def handle_seqmatrix_import(self, sequences_list, partition_definitions):
        """处理从 SeqMatrix 导入的数据
        
        Args:
            sequences_list: 序列列表 [{'name': str, 'seq': str}, ...]
            partition_definitions: 分区定义列表 [{'name': str, 'start': int, 'end': int}, ...]
        """
        try:
            from YR_MPE.platforms.methods.import_partitioned_nexus import TempDatasetItem
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord

            # 创建 DatasetItem 列表
            dataset_items = []

            for partition in partition_definitions:
                partition_name = partition['name']
                start_pos = partition['start'] - 1  # 转换为0-based索引
                end_pos = partition['end']  # 包含结束位置

                # 创建 DatasetItem
                dataset_item = TempDatasetItem()
                dataset_item.loci_name = partition_name
                dataset_item.file_path = "SeqMatrix"
                dataset_item.format = 'nexus'

                # 提取该分区的序列
                item_sequences = []
                for seq_data in sequences_list:
                    full_seq = seq_data['seq']
                    partition_seq = full_seq[start_pos:end_pos]

                    # 创建 SeqRecord
                    seq_record = SeqRecord(
                        Seq(partition_seq),
                        id=seq_data['name'],
                        description=f"{partition_name} partition"
                    )
                    item_sequences.append(seq_record)

                dataset_item.sequences = item_sequences
                dataset_item.length = end_pos - start_pos
                dataset_item.sequence_count = len(item_sequences)

                # 检查是否已比对
                all_chars = ''.join(str(record.seq) for record in item_sequences)
                has_missing = '?' in all_chars or '-' in all_chars
                dataset_item.is_aligned = not has_missing or (has_missing and len(set(all_chars) - set('?-\n')) > 1)

                dataset_items.append(dataset_item)

            # 生成分区方案字符串
            partition_scheme = "begin sets;\n"
            for partition in partition_definitions:
                partition_scheme += f"    charset {partition['name']} = {partition['start']}-{partition['end']};\n"
            partition_scheme += "end;"

            # 创建 Dataset 对象
            class DatasetObject:
                def __init__(self, name, items=None, partition_scheme=None):
                    self.dataset_name = name
                    self.items = items if items else []
                    self.partition_scheme = partition_scheme
                    self.summary = {}

            dataset = DatasetObject(
                name="SeqMatrix Dataset",
                items=dataset_items,
                partition_scheme=partition_scheme
            )

            # 添加到工作区
            self.workspace.add_dataset(dataset)

            # 显示成功消息
            partition_names = [p['name'] for p in partition_definitions]
            message = (
                f"Successfully imported dataset from SeqMatrix\n\n"
                f"Taxa: {len(sequences_list)}\n"
                f"Total Length: {len(sequences_list[0]['seq']) if sequences_list else 0}\n"
                f"Partitions: {len(partition_definitions)}\n"
                f"Partition Names: {', '.join(partition_names)}"
            )
            QMessageBox.information(self, "Import Success", message)

        except Exception as e:
            QMessageBox.critical(self, "Import Error", f"Failed to import dataset from SeqMatrix: {str(e)}")
    
    def create_dataset_by_seqdbg(self):
        """通过SeqDBG创建Dataset"""
        try:
            # 调用SeqDBG插件
            self.open_seqdbg_wrapper()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create dataset by SeqDBG: {str(e)}")
    
    def _check_workspace_required(self, plugin_name: str) -> bool:
        """
        检查是否需要工作区，如果没有则提示用户选择
        
        Returns:
            bool: True表示可以继续，False表示用户取消
        """
        if not self.workspace_manager.workspace_path:
            reply = QMessageBox.warning(
                self,
                "Workspace Required",
                f"{plugin_name} requires a workspace folder to store\n"
                f"intermediate files and run state.\n\n"
                f"Bayesian inference can take several hours to days,\n"
                f"and it's important to save all intermediate files\n"
                f"for analysis and potential recovery.\n\n"
                f"Do you want to select a workspace folder now?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes
            )
            
            if reply == QMessageBox.Yes:
                # 让用户选择工作区
                folder_path = QFileDialog.getExistingDirectory(
                    self,
                    "Select Workspace Folder",
                    os.path.expanduser("~")
                )
                if folder_path:
                    self._switch_workspace(folder_path, is_new=True)
                    return True
                else:
                    return False  # 用户取消了文件夹选择
            else:
                return False  # 用户拒绝选择工作区
        
        return True  # 已有工作区，可以继续
    
    def open_mrbayes_wrapper(self):
        """打开MrBayes插件"""
        # 检查工作区
        if not self._check_workspace_required("MrBayes"):
            return
        
        # 如果已有对话框，先删除
        if hasattr(self, 'mrbayes_dialog') and self.mrbayes_dialog:
            self.mrbayes_dialog.close()
            self.mrbayes_dialog.deleteLater()
        
        # 创建MrBayes对话框（非模态窗口，允许后台运行）
        from PyQt5.QtWidgets import QDialog, QMessageBox
        from PyQt5.QtCore import Qt
        from .methods.dataset_models import ITEM_TYPE_MODEL, SELECTION_STATE_GREEN
        from ..plugins.partition_mode import MrBayesModelConverter
        
        dialog = QDialog(self)  # 设置父窗口
        dialog.setWindowTitle("MrBayes3-MPI-BEAGLE - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/mrbayes.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())
        dialog.setWindowModality(Qt.NonModal)  # 设置为非模态窗口
        dialog.setAttribute(Qt.WA_DeleteOnClose, False)  # 禁止自动删除
        
        # 准备导入数据 (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # 获取最佳模型（只检查green状态的model）
        best_model = ""
        seq_type = "DNA"  # 默认DNA
        # 从dataset_selection_manager获取green状态的model
        if self.dataset_selection_manager:
            green_models = [item for item in self.dataset_selection_manager.get_selected_items() 
                          if item.item_type == ITEM_TYPE_MODEL]
            if green_models:
                model_data = green_models[0].data
                # 处理模型表或单个模型
                if isinstance(model_data, dict):
                    if "type" in model_data and model_data["type"] == "model_table" and model_data["data"]:
                        # 取模型表中的第一个（最佳）模型
                        best_model = model_data["data"][0]['Model']
                        # 尝试获取序列类型
                        if "Data" in model_data["data"][0]:
                            seq_type = model_data["data"][0]['Data']
                    elif "Model" in model_data:
                        # 直接使用单个模型数据
                        best_model = model_data
                        # 尝试获取序列类型
                        if "Data" in model_data:
                            seq_type = model_data['Data']
        
        # 处理蛋白质模型的转换和警告
        model_conversion_result = None
        if best_model:
            if seq_type.upper() == "PROTEIN":
                model_conversion_result, warnings = MrBayesModelConverter.convert_model_to_mrbayes(
                    best_model, seq_type
                )
                
                if warnings:
                    # 显示警告对话框
                    warning_text = "\n".join(warnings)
                    warning_text += "\n\n请选择："
                    
                    reply = QMessageBox.warning(
                        self,
                        "模型转换警告",
                        warning_text,
                        QMessageBox.Yes | QMessageBox.No,
                        QMessageBox.No
                    )
                    
                    if reply == QMessageBox.Yes:
                        # 选项1: 使用GTR
                        model_conversion_result = "gtr"
                    else:
                        # 选项2: 使用mixed
                        model_conversion_result = "mixed"
        
        # 使用PluginFactory获取MrBayes插件
        mrbayes_entry = self.plugin_factory.get_mrbayes_plugin()
        workdir = self.get_workdir()
        mrbayes_wrapper = mrbayes_entry.run(
            import_from=import_from, 
            import_data=import_data,
            workdir=workdir,
            imported_model=best_model,
            model_conversion_result=model_conversion_result,
            seq_type=seq_type
        )
        
        # 连接信号
        mrbayes_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        mrbayes_wrapper.export_phylogeny_result_signal.connect(self.add_phylogeny_to_workspace)
        mrbayes_wrapper.export_chain_result_signal.connect(self.add_chain_to_workspace)
        
        dialog.layout().addWidget(mrbayes_wrapper)
        # 保存到实例变量，防止被垃圾回收
        self.mrbayes_dialog = dialog
        dialog.show()  # 使用show()而非exec_()，显示非模态窗口
    
    def open_phylobayes_wrapper(self):
        """打开PhyloBayes插件"""
        # 检查工作区
        if not self._check_workspace_required("PhyloBayes"):
            return
        
        # 如果已有对话框，先删除
        if hasattr(self, 'phylobayes_dialog') and self.phylobayes_dialog:
            self.phylobayes_dialog.close()
            self.phylobayes_dialog.deleteLater()
        
        # 创建PhyloBayes对话框（非模态窗口，允许后台运行）
        from PyQt5.QtWidgets import QDialog
        from PyQt5.QtCore import Qt
        dialog = QDialog(self)  # 设置父窗口
        dialog.setWindowTitle("PhyloBayes-MPI - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/phylobayes.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())
        dialog.setWindowModality(Qt.NonModal)  # 设置为非模态窗口
        dialog.setAttribute(Qt.WA_DeleteOnClose, False)  # 禁止自动删除
        
        # 准备导入数据 (使用抽象层)
        import_from, import_data = self._prepare_import_data(workspace_item_types=["alignments"])
        
        # 使用PluginFactory获取PhyloBayes插件
        phylobayes_entry = self.plugin_factory.get_phylobayes_plugin()
        workdir = self.get_workdir()
        phylobayes_wrapper = phylobayes_entry.run(
            import_from=import_from, 
            import_data=import_data,
            workdir=workdir
        )
        
        # 连接信号
        phylobayes_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        phylobayes_wrapper.export_phylogeny_result_signal.connect(self.add_phylogeny_to_workspace)
        phylobayes_wrapper.export_chain_result_signal.connect(self.add_chain_to_workspace)
        
        dialog.layout().addWidget(phylobayes_wrapper)
        # 保存到实例变量，防止被垃圾回收
        self.phylobayes_dialog = dialog
        dialog.show()  # 使用show()而非exec_()，显示非模态窗口

class SingleGeneWorkspace(QWidget):
    def __init__(self, resource_factory=None, dataset_selection_manager=None, parent=None):
        super().__init__()
        self.resource_factory = resource_factory  # 添加resource_factory引用
        self.parent_window = parent  # 添加parent引用
        self.dataset_selection_manager = dataset_selection_manager  # 添加数据集选择管理器
        
        # 保留旧的 items 结构以保持向后兼容
        self.items = {
            "alignments": [],
            "models": [],
            "distances": [],
            "phylogenies": [],
            "chains": [],  # 添加chains项
            "datasets": [],  # 添加datasets项
            # "variants": [], TODO
            # "coalescent": [], TODO
            # "clock": [], TODO
        }
        
        # 新的按钮存储，按类型分组
        self.item_buttons = {
            "alignments": [],
            "models": [],
            "distances": [],
            "phylogenies": [],
            "treesets": [],
            "chains": [],
            "datasets": []
        }
        
        # 当前选中的数据集ID
        self.current_dataset_id = None
        
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.init_ui()

    def init_ui(self):
        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)
        # background color: white
        # self.setStyleSheet("background-color: white;")
        # add a label: "Single Gene Workspace"
        # self.workspace_hint = QLabel("Single Gene Workspace\nAdd an alignment or drag and drop a file here to start")
        self.workspace_hint = QLabel("Add a sequence file or create a new dataset to start.") 
        self.workspace_hint.setAlignment(Qt.AlignCenter)
        self.workspace_hint.setStyleSheet("color: #555555;")
        self.main_layout.addWidget(self.workspace_hint)

        # TODO: drag and drop event
    
    def add_sequence(self, sequences):
        """添加序列到工作区"""
        # 如果还没有数据集，创建一个默认数据集
        if not self.current_dataset_id:
            if self.dataset_selection_manager:
                self.current_dataset_id = self.dataset_selection_manager.create_dataset(
                    name="Default Dataset",
                    description="Auto-created default dataset",
                    is_multigene=False
                )
        
        # 确保网格布局存在
        if not hasattr(self, 'grid_layout'):
            self.workspace_hint.setVisible(False)
            self.grid_widget = QWidget()
            self.grid_layout = QGridLayout()
            self.grid_widget.setLayout(self.grid_layout)
            self.main_layout.addWidget(self.grid_widget)
            self.grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        
        # 创建 DatasetItem
        from datetime import datetime
        dataset_item = DatasetItem(item_type=ITEM_TYPE_SEQUENCE)
        dataset_item.dataset_id = self.current_dataset_id
        dataset_item.loci_name = f"Sequence_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        dataset_item.name = dataset_item.loci_name
        
        # 设置所有权 UUID（初始数据生成新的 UUID）
        dataset_item.ownership_uuid = str(uuid.uuid4())
        
        # 存储序列数据
        dataset_item.sequences = sequences
        dataset_item.sequence_count = len(sequences)
        if sequences and len(sequences) > 0:
            if hasattr(sequences[0], 'seq'):
                try:
                    dataset_item.length = len(str(sequences[0].seq))
                except:
                    dataset_item.length = 0
            else:
                dataset_item.length = 0
        dataset_item.is_aligned = False  # 未比对
        
        # 添加到管理器
        if self.dataset_selection_manager:
            success = self.dataset_selection_manager.add_item(dataset_item, self.current_dataset_id)
            if not success:
                QMessageBox.warning(self, "Warning", "Failed to add item to dataset")
                return
            # 使用 dataset_item 的 ID 来获取数据项（因为 add_item 会设置 ID）
            dataset_item = self.dataset_selection_manager.get_item(dataset_item.id)
            if not dataset_item:
                QMessageBox.warning(self, "Warning", "Failed to retrieve added item")
                return
        
        # 保留旧的数据结构（向后兼容）
        self.items["alignments"].append(sequences)
        
        # 创建 DatasetItemButton
        item_button = DatasetItemButton(dataset_item, parent=self)
        item_button.clicked_single.connect(self._on_item_single_click)
        item_button.clicked_double.connect(self._on_item_double_click)
        
        # 设置图标
        if self.resource_factory:
            icon = self.resource_factory.get_icon("file/sequence.svg")
            if icon:
                item_button.setIcon(icon)
                item_button.setIconSize(QSize(45, 45))
        
        # 存储按钮引用
        self.item_buttons["alignments"].append(item_button)
        
        # 添加到网格布局
        # alignments 总是从第 0 列开始，依次向右排列
        self.grid_layout.addWidget(item_button, 0, len(self.item_buttons["alignments"])-1)
        
        # 如果有 dataset 按钮，需要将它们向右移动以避免覆盖
        if "datasets" in self.item_buttons and self.item_buttons["datasets"]:
            self._rearrange_row_0_buttons()
        
        # 自动选中该数据项
        if self.dataset_selection_manager:
            self.dataset_selection_manager.select_item(dataset_item.id)
            # 更新所有按钮的样式（select_item 会改变多个项目的状态）
            self._update_all_button_styles()
    
    def _on_item_single_click(self, item_id: str):
        """处理数据项单击事件"""
        if not self.dataset_selection_manager:
            return
        
        # 清除所有数据集的高亮
        self._clear_dataset_highlights()
        
        # 切换选择状态
        item = self.dataset_selection_manager.get_item(item_id)
        if item and item.selection_state == SELECTION_STATE_GREEN:
            # 如果已经选中，只取消这个项目的选择，不影响其他项目
            item.selection_state = SELECTION_STATE_NONE
            item.selection_reason = ""
        else:
            # 否则选中该项（select_item会处理同一类型和不同UUID的deactivate逻辑）
            self.dataset_selection_manager.select_item(item_id)
        
        # 更新所有按钮的样式
        self._update_all_button_styles()
    
    def _clear_dataset_highlights(self):
        """清除所有数据集的高亮"""
        if not self.dataset_selection_manager:
            return
        
        for dataset in self.dataset_selection_manager.get_all_datasets():
            dataset.selection_state = SELECTION_STATE_NONE
        
        # 保存状态
        self.dataset_selection_manager._save_state()
    
    def _on_item_double_click(self, item_id: str):
        """处理数据项双击事件"""
        if not self.dataset_selection_manager:
            return
        
        item = self.dataset_selection_manager.get_item(item_id)
        if not item:
            return
        
        # 根据数据类型执行相应的查看操作
        if item.item_type == ITEM_TYPE_SEQUENCE or item.item_type == ITEM_TYPE_ALIGNMENT:
            self.view_alignment(item.sequences)
        elif item.item_type == ITEM_TYPE_MODEL:
            self.view_model_result(item.data)
        elif item.item_type == ITEM_TYPE_DISTANCE:
            self.view_distance_matrix(item.data)
        elif item.item_type == ITEM_TYPE_PHYLOGENY:
            self.view_phylogeny(item.data)
        elif item.item_type == ITEM_TYPE_TREESET:
            # TreeSet 双击打开树集合查看器
            self._on_treeset_click(item.data)
        elif item.item_type == ITEM_TYPE_CHAIN:
            # Chain 数据的特殊处理
            if self.parent_window:
                self.parent_window.open_minitracer_wrapper(chain_item=item.data)
    
    def _update_all_button_styles(self):
        """更新所有按钮的样式"""
        for button_list in self.item_buttons.values():
            for button in button_list:
                if isinstance(button, DatasetItemButton):
                    button.update_style()
                elif isinstance(button, DatasetButton):
                    button.update_style()
    
    def _rearrange_row_0_buttons(self):
        """重新排列第 0 行的所有按钮，确保 alignments 在左，datasets 在右"""
        # 首先移除所有第 0 行的按钮
        alignments_count = len(self.item_buttons.get("alignments", []))
        datasets_count = len(self.item_buttons.get("datasets", []))
        
        # 从右向左移除，避免索引变化
        for i in range(max(alignments_count, datasets_count) - 1, -1, -1):
            # 移除 alignments
            if i < alignments_count:
                button = self.item_buttons["alignments"][i]
                self.grid_layout.removeWidget(button)
            # 移除 datasets
            if i < datasets_count:
                button = self.item_buttons["datasets"][i]
                self.grid_layout.removeWidget(button)
        
        # 重新添加所有按钮，alignments 在左，datasets 在右
        # alignments 从第 0 列开始
        for i, button in enumerate(self.item_buttons.get("alignments", [])):
            self.grid_layout.addWidget(button, 0, i)
        
        # datasets 从 alignments 的数量之后开始
        start_col = len(self.item_buttons.get("alignments", []))
        for i, button in enumerate(self.item_buttons.get("datasets", [])):
            self.grid_layout.addWidget(button, 0, start_col + i)
        
        # 重新排列后更新所有按钮的样式
        self._update_all_button_styles()
    
    def _create_item_button(self, dataset_item: DatasetItem, item_type_key: str, row: int, col: int, icon_name: str):
        """创建数据项按钮的通用方法
        
        Args:
            dataset_item: 数据项对象
            item_type_key: 数据类型键值（如 "models", "distances"）
            row: 网格布局行号
            col: 网格布局列号
            icon_name: 图标文件名
        """
        # 确保网格布局存在
        if not hasattr(self, 'grid_layout'):
            self.workspace_hint.setVisible(False)
            self.grid_widget = QWidget()
            self.grid_layout = QGridLayout()
            self.grid_widget.setLayout(self.grid_layout)
            self.main_layout.addWidget(self.grid_widget)
            self.grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        
        # 创建 DatasetItemButton
        item_button = DatasetItemButton(dataset_item, parent=self)
        item_button.clicked_single.connect(self._on_item_single_click)
        item_button.clicked_double.connect(self._on_item_double_click)
        
        # 设置图标
        if self.resource_factory:
            icon = self.resource_factory.get_icon(icon_name)
            if icon:
                item_button.setIcon(icon)
                item_button.setIconSize(QSize(45, 45))
        
        # 存储按钮引用
        if item_type_key not in self.item_buttons:
            self.item_buttons[item_type_key] = []
        self.item_buttons[item_type_key].append(item_button)
        
        # 添加到网格布局
        self.grid_layout.addWidget(item_button, row, col)
        
        # 检查新项目是否应该显示为蓝色高亮
        # 如果新项目继承了当前激活的 ownership UUID，将其设置为蓝色
        if self.dataset_selection_manager and dataset_item.ownership_uuid:
            active_uuid = self._get_active_ownership_uuid()
            if active_uuid and dataset_item.ownership_uuid == active_uuid:
                # 将新项目设置为蓝色高亮（上下文高亮）
                dataset_item.selection_state = SELECTION_STATE_BLUE
                dataset_item.selection_reason = "Same ownership as active items"
                item_button.update_style()
        
        return item_button
    
    def _get_active_ownership_uuid(self) -> Optional[str]:
        """获取当前激活的序列/alignment 的 ownership_uuid"""
        if not self.dataset_selection_manager:
            return None
        
        # 获取当前激活（绿色）的数据项
        green_items = self.dataset_selection_manager.get_items_by_state(SELECTION_STATE_GREEN)
        
        # 优先返回 sequence 或 alignment 类型的 UUID
        for item in green_items:
            if item.item_type in [ITEM_TYPE_SEQUENCE, ITEM_TYPE_ALIGNMENT]:
                return item.ownership_uuid
        
        # 如果没有 sequence/alignment，返回第一个绿色项的 UUID
        if green_items:
            return green_items[0].ownership_uuid
        
        return None
    
    def add_model(self, model):
        """添加模型到工作区"""
        # 如果还没有数据集，创建一个默认数据集
        if not self.current_dataset_id:
            if self.dataset_selection_manager:
                self.current_dataset_id = self.dataset_selection_manager.create_dataset(
                    name="Default Dataset",
                    description="Auto-created default dataset",
                    is_multigene=False
                )
        
        # 创建 DatasetItem
        from datetime import datetime
        dataset_item = DatasetItem(item_type=ITEM_TYPE_MODEL)
        dataset_item.dataset_id = self.current_dataset_id
        dataset_item.loci_name = f"Model_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        dataset_item.name = dataset_item.loci_name
        
        # 继承所有权 UUID（如果存在）
        ownership_uuid = self._get_active_ownership_uuid()
        if ownership_uuid:
            dataset_item.ownership_uuid = ownership_uuid
        else:
            # 如果没有可继承的 UUID，生成新的
            dataset_item.ownership_uuid = str(uuid.uuid4())
        
        # 存储模型数据
        if isinstance(model, dict):
            dataset_item.data = model
            
            # 检测模型类型（单一模型 vs 分区模型）
            model_type = model.get("type", "single")
            if model_type == "partitioned":
                # 分区模型
                dataset_item.model_sub_type = "partitioned"
                dataset_item.partition_mode = model.get("partition_mode", "")
                model_name = f"PartitionedModel_{model.get('partition_mode', '')}"
            else:
                # 单一模型
                dataset_item.model_sub_type = "single"
                model_name = model.get("name", model.get("model", "Unknown Model"))
        else:
            # 兼容旧格式
            dataset_item.data = {"model": getattr(model, 'model_name', 'Unknown Model')}
            dataset_item.model_sub_type = "single"
            model_name = getattr(model, 'model_name', 'Unknown Model')
        
        # 保留旧的数据结构（向后兼容）
        self.items["models"].append(model)
        
        # 添加到管理器
        if self.dataset_selection_manager:
            success = self.dataset_selection_manager.add_item(dataset_item, self.current_dataset_id)
            if not success:
                QMessageBox.warning(self, "Warning", "Failed to add model to dataset")
                return
            dataset_item = self.dataset_selection_manager.get_item(dataset_item.id)
            if not dataset_item:
                QMessageBox.warning(self, "Warning", "Failed to retrieve added model")
                return
        
        # 根据模型类型选择图标
        if dataset_item.model_sub_type == "partitioned":
            # 分区模型：根据分区模式选择图标
            icon_map = {
                "EL": "file/modelset-el.svg",
                "TL": "file/modelset-el.svg",
                "EUL": "file/modelset-eul.svg",
                "TUL": "file/modelset-tul.svg"
            }
            icon_name = icon_map.get(dataset_item.partition_mode, "file/model.svg")
        else:
            # 单一模型
            icon_name = "file/model.svg"
        
        # 创建按钮
        self._create_item_button(
            dataset_item=dataset_item,
            item_type_key="models",
            row=1,
            col=len(self.item_buttons["models"]),
            icon_name=icon_name
        )
    
    def add_distance(self, distance):
        """添加距离矩阵到工作区"""
        # 如果还没有数据集，创建一个默认数据集
        if not self.current_dataset_id:
            if self.dataset_selection_manager:
                self.current_dataset_id = self.dataset_selection_manager.create_dataset(
                    name="Default Dataset",
                    description="Auto-created default dataset",
                    is_multigene=False
                )
        
        # 创建 DatasetItem
        from datetime import datetime
        dataset_item = DatasetItem(item_type=ITEM_TYPE_DISTANCE)
        dataset_item.dataset_id = self.current_dataset_id
        dataset_item.loci_name = f"Distance_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        dataset_item.name = dataset_item.loci_name
        
        # 继承所有权 UUID（如果存在）
        ownership_uuid = self._get_active_ownership_uuid()
        if ownership_uuid:
            dataset_item.ownership_uuid = ownership_uuid
        else:
            # 如果没有可继承的 UUID，生成新的
            dataset_item.ownership_uuid = str(uuid.uuid4())
        
        # 存储距离数据
        if isinstance(distance, dict):
            dataset_item.data = distance
        else:
            dataset_item.data = {"content": distance}
        
        # 保留旧的数据结构（向后兼容）
        self.items["distances"].append(distance)
        
        # 添加到管理器
        if self.dataset_selection_manager:
            success = self.dataset_selection_manager.add_item(dataset_item, self.current_dataset_id)
            if not success:
                QMessageBox.warning(self, "Warning", "Failed to add distance to dataset")
                return
            dataset_item = self.dataset_selection_manager.get_item(dataset_item.id)
            if not dataset_item:
                QMessageBox.warning(self, "Warning", "Failed to retrieve added distance")
                return
        
        # 创建按钮
        self._create_item_button(
            dataset_item=dataset_item,
            item_type_key="distances",
            row=2,
            col=len(self.item_buttons["distances"]),
            icon_name="file/distance.svg"
        )
    
    def add_phylogeny(self, phylogeny):
        """添加系统发育树到工作区"""
        # 检测是否有多个树
        tree_count = 1
        if isinstance(phylogeny, dict) and 'data' in phylogeny:
            tree_count = len(phylogeny['data'])
        
        # 如果有多个树，创建 treeset 对象
        if tree_count > 1:
            self.add_treeset(phylogeny)
        else:
            # 单个树，创建 phylogeny 对象
            self._add_single_phylogeny(phylogeny)
    
    def _add_single_phylogeny(self, phylogeny):
        """添加单个系统发育树到工作区"""
        # 如果还没有数据集，创建一个默认数据集
        if not self.current_dataset_id:
            if self.dataset_selection_manager:
                self.current_dataset_id = self.dataset_selection_manager.create_dataset(
                    name="Default Dataset",
                    description="Auto-created default dataset",
                    is_multigene=False
                )
        
        # 创建 DatasetItem
        from datetime import datetime
        dataset_item = DatasetItem(item_type=ITEM_TYPE_PHYLOGENY)
        dataset_item.dataset_id = self.current_dataset_id
        dataset_item.loci_name = f"Phylogeny_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        dataset_item.name = dataset_item.loci_name
        
        # 继承所有权 UUID（如果存在）
        ownership_uuid = self._get_active_ownership_uuid()
        if ownership_uuid:
            dataset_item.ownership_uuid = ownership_uuid
        else:
            # 如果没有可继承的 UUID，生成新的
            dataset_item.ownership_uuid = str(uuid.uuid4())
        
        # 存储树数据
        if isinstance(phylogeny, dict):
            dataset_item.data = phylogeny
        else:
            dataset_item.data = {"content": phylogeny}
        
        # 保留旧的数据结构（向后兼容）
        self.items["phylogenies"].append(phylogeny)
        
        # 添加到管理器
        if self.dataset_selection_manager:
            success = self.dataset_selection_manager.add_item(dataset_item, self.current_dataset_id)
            if not success:
                QMessageBox.warning(self, "Warning", "Failed to add phylogeny to dataset")
                return
            dataset_item = self.dataset_selection_manager.get_item(dataset_item.id)
            if not dataset_item:
                QMessageBox.warning(self, "Warning", "Failed to retrieve added phylogeny")
                return
        
        # 创建按钮
        self._create_item_button(
            dataset_item=dataset_item,
            item_type_key="phylogenies",
            row=3,
            col=len(self.item_buttons.get("phylogenies", [])),
            icon_name="file/phylogeny.svg"
        )
    
    def add_treeset(self, treeset_data):
        """添加基因树集合（treeset）到工作区"""
        # 如果还没有数据集，创建一个默认数据集
        if not self.current_dataset_id:
            if self.dataset_selection_manager:
                self.current_dataset_id = self.dataset_selection_manager.create_dataset(
                    name="Default Dataset",
                    description="Auto-created default dataset",
                    is_multigene=False
                )
        
        # 创建 DatasetItem
        from datetime import datetime
        dataset_item = DatasetItem(item_type=ITEM_TYPE_TREESET)
        dataset_item.dataset_id = self.current_dataset_id
        dataset_item.loci_name = f"TreeSet_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        dataset_item.name = dataset_item.loci_name
        
        # 继承所有权 UUID（如果存在）
        ownership_uuid = self._get_active_ownership_uuid()
        if ownership_uuid:
            dataset_item.ownership_uuid = ownership_uuid
        else:
            # 如果没有可继承的 UUID，生成新的
            dataset_item.ownership_uuid = str(uuid.uuid4())
        
        # 存储树集合数据
        if isinstance(treeset_data, dict):
            dataset_item.data = treeset_data
        else:
            dataset_item.data = {"content": treeset_data}
        
        # 保留旧的数据结构（向后兼容）
        if 'treesets' not in self.items:
            self.items['treesets'] = []
        self.items["treesets"].append(treeset_data)
        
        # 添加到管理器
        if self.dataset_selection_manager:
            success = self.dataset_selection_manager.add_item(dataset_item, self.current_dataset_id)
            if not success:
                QMessageBox.warning(self, "Warning", "Failed to add treeset to dataset")
                return
            dataset_item = self.dataset_selection_manager.get_item(dataset_item.id)
            if not dataset_item:
                QMessageBox.warning(self, "Warning", "Failed to retrieve added treeset")
                return
        
        # 创建按钮
        self._create_item_button(
            dataset_item=dataset_item,
            item_type_key="treesets",
            row=5,
            col=len(self.item_buttons.get("treesets", [])),
            icon_name="file/treeset.svg"
        )
    
    def _on_treeset_click(self, treeset_data):
        """处理 treeset 按钮点击事件 - 直接用 IcyTree 呈现所有树"""
        # 从 treeset_data 中提取所有树的 Newick 字符串
        trees = []
        if isinstance(treeset_data, dict) and 'data' in treeset_data:
            for tree_data in treeset_data['data']:
                if 'content' in tree_data:
                    trees.append(tree_data['content'])
                elif 'file_path' in tree_data:
                    # 从文件读取
                    try:
                        with open(tree_data['file_path'], 'r') as f:
                            trees.append(f.read().strip())
                    except:
                        pass
        elif isinstance(treeset_data, str):
            # 如果直接是字符串
            trees.append(treeset_data)
        
        if not trees:
            QMessageBox.warning(self, "Warning", "No trees found in the tree set")
            return
        
        # 用 IcyTree 呈现所有树
        if self.parent_window:
            # 将所有树的 Newick 字符串用换行符连接
            combined_newick = '\n'.join(trees)
            self.parent_window.open_icytree_wrapper(phylogeny_data={'content': combined_newick, 'treeset': True})
    
    def _import_treeset_to_astral(self, treeset_data):
        """导入树集合到 ASTRAL"""
        # 将树集合格式化为 ASTRAL 可接受的格式
        trees = []
        if isinstance(treeset_data, dict) and 'data' in treeset_data:
            for tree_data in treeset_data['data']:
                if 'content' in tree_data:
                    trees.append(tree_data['content'])
                elif 'file_path' in tree_data:
                    # 从文件读取
                    try:
                        with open(tree_data['file_path'], 'r') as f:
                            trees.append(f.read().strip())
                    except:
                        pass
        
        if trees:
            # 创建临时文件，所有树用换行符分隔
            import tempfile
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.trees', delete=False).name
            with open(temp_file, 'w') as f:
                f.write('\n'.join(trees))
            
            # 打开 ASTRAL 插件
            # 注意：这里需要访问 parent_window 的 astral 相关方法
            # 由于访问限制，我们需要通过 parent_window 来调用
            # 先检查 parent_window 是否有 astral 相关方法
            if hasattr(self.parent_window, 'open_astral_wrapper'):
                self.parent_window.open_astral_wrapper(import_data=temp_file)
            else:
                QMessageBox.information(self, "Info", "ASTRAL plugin integration not available")
        else:
            QMessageBox.warning(self, "Warning", "No trees found in the tree set")
    
    def add_chain(self, chain_item):
        """添加MCMC链文件到工作区"""
        # 确保grid_layout存在
        if not hasattr(self, 'grid_layout'):
            # disable hint label
            self.workspace_hint.setVisible(False)
            # add a grid layout
            self.grid_widget = QWidget()
            self.grid_layout = QGridLayout()
            self.grid_widget.setLayout(self.grid_layout)
            self.main_layout.addWidget(self.grid_widget)
            self.grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        
        # add a chain to items["chains"]
        self.items["chains"].append(chain_item)
        # add a chain icon to workspace
        chain_icon = self.resource_factory.get_icon("file/chain.svg")
        chain_button = QToolButton()
        chain_button.setIcon(chain_icon)
        chain_button.setIconSize(QSize(45, 45))
        
        # 设置样式，背景透明（与其他未选中的按钮保持一致）
        chain_button.setStyleSheet("""
            QToolButton {
                background-color: transparent;
            }
        """)
        
        # 创建tooltip，显示链信息
        tooltip_text = f"MCMC Chains (Run {chain_item.run_number}, {chain_item.chain_count} chain(s), Tool: {chain_item.tool})"
        chain_button.setToolTip(tooltip_text)
        # 点击打开miniTracer（通过parent_window）
        chain_button.clicked.connect(lambda: self.parent_window.open_minitracer_wrapper(chain_item=chain_item))
        self.grid_layout.addWidget(chain_button, 4, len(self.items["chains"])-1)
    
    def add_dataset(self, dataset):
        """添加数据集到工作区"""
        # add a dataset to items["datasets"]
        self.items["datasets"].append(dataset)
        
        # 确保grid_layout存在
        if not hasattr(self, 'grid_layout'):
            # disable hint label
            self.workspace_hint.setVisible(False)
            # add a grid layout
            self.grid_widget = QWidget()
            self.grid_layout = QGridLayout()
            self.grid_widget.setLayout(self.grid_layout)
            self.main_layout.addWidget(self.grid_widget)
            self.grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        
        # 创建或获取 DatasetInfo
        dataset_name = getattr(dataset, 'dataset_name', 'Unnamed Dataset')
        
        # 检查是否已经有这个数据集在管理器中
        dataset_info = None
        dataset_id = None
        if self.dataset_selection_manager:
            # 查找是否已存在同名数据集
            for ds in self.dataset_selection_manager.get_all_datasets():
                if ds.name == dataset_name:
                    dataset_info = ds
                    dataset_id = ds.id
                    break
        
        # 如果不存在，创建新的数据集
        if not dataset_info and self.dataset_selection_manager:
            is_multigene = hasattr(dataset, 'items') and len(dataset.items) > 1
            dataset_id = self.dataset_selection_manager.create_dataset(
                name=dataset_name,
                description=f"Dataset created from {dataset_name}",
                is_multigene=is_multigene
            )
            dataset_info = self.dataset_selection_manager.get_dataset(dataset_id)
            
            # 如果有 items，添加到数据集
            if hasattr(dataset, 'items') and dataset.items:
                for item in dataset.items:
                    # 这里需要将旧格式的 item 转换为 DatasetItem
                    # 暂时跳过，因为需要更复杂的转换逻辑
                    pass
        
        # 保存 dataset_id 到 dataset 对象中，以便后续使用
        if dataset_id:
            if isinstance(dataset, dict):
                dataset['dataset_id'] = dataset_id
            else:
                setattr(dataset, 'dataset_id', dataset_id)
        
        # 创建 DatasetButton
        if dataset_info:
            dataset_button = DatasetButton(dataset_info, parent=self)
            dataset_button.clicked_single.connect(self._on_dataset_click)
            # 双击事件需要适配参数：DatasetButton 传递 dataset_id，但 _on_dataset_double_click 期望 button, event
            dataset_button.clicked_double.connect(lambda dataset_id: self._on_dataset_double_click_by_id(dataset_id))
        else:
            # 如果没有数据集信息，使用普通按钮（向后兼容）
            dataset_button = QToolButton()
        
        # 设置图标
        if self.resource_factory:
            icon = self.resource_factory.get_icon("file/dataset.svg")
            if icon:
                dataset_button.setIcon(icon)
                dataset_button.setIconSize(QSize(45, 45))
        
        # 设置 tooltip
        if hasattr(dataset, 'dataset_name'):
            dataset_button.setToolTip(f"Dataset: {dataset.dataset_name}")
        
        # Dataset数据项与 alignments 放在同一行（第0行），因为它们在所有权机制上是并列的
        dataset_row = 0
        
        # 存储按钮引用
        if "datasets" not in self.item_buttons:
            self.item_buttons["datasets"] = []
        self.item_buttons["datasets"].append(dataset_button)
        
        # 添加到网格布局
        # Dataset 的列号从 alignments 数量开始，避免与 alignments 重叠
        dataset_col = len(self.item_buttons.get("alignments", [])) + len(self.item_buttons["datasets"]) - 1
        self.grid_layout.addWidget(dataset_button, dataset_row, dataset_col)
        
        # 存储dataset引用到按钮上，便于后续访问
        dataset_button.dataset_ref = dataset
    
    def _on_dataset_click(self, dataset_id: str):
        """处理数据集单击事件"""
        if not self.dataset_selection_manager:
            return
        
        # 获取数据集信息
        dataset_info = self.dataset_selection_manager.get_dataset(dataset_id)
        if not dataset_info:
            return
        
        # 切换数据集的选择状态
        if dataset_info.selection_state == SELECTION_STATE_GREEN:
            # 如果当前是绿色，取消选择（设为无色）
            dataset_info.selection_state = SELECTION_STATE_NONE
            # 清除所有数据项的选择
            self.dataset_selection_manager.selection_engine.clear_all_selections()
        else:
            # 如果当前是无色或蓝色，设置为绿色
            dataset_info.selection_state = SELECTION_STATE_GREEN
            
            # 从 dataset.settings['dataset_items'] 恢复选中状态，而不是将所有设为蓝色
            if 'dataset_items' in dataset_info.settings:
                for item_data in dataset_info.settings['dataset_items']:
                    if item_data.get('selected', False):
                        loci_name = item_data.get('loci_name', '')
                        for ds_item in self.dataset_selection_manager.items.values():
                            if ds_item.loci_name == loci_name and ds_item.dataset_id == dataset_id:
                                ds_item.selection_state = SELECTION_STATE_GREEN
                                self.dataset_selection_manager.selected_items.add(ds_item.id)
                                break
            else:
                # 如果没有保存的状态，才将所有数据项设置为蓝色
                self.dataset_selection_manager.selection_engine.select_dataset(dataset_id)
        
        # 保存状态
        self.dataset_selection_manager._save_state()
        
        # 更新所有按钮的样式
        self._update_all_button_styles()
    
    def _on_dataset_double_click_by_id(self, dataset_id: str):
        """处理数据集双击事件（通过 dataset_id）"""
        if not self.dataset_selection_manager:
            return
        
        # 获取数据集信息
        dataset_info = self.dataset_selection_manager.get_dataset(dataset_id)
        if dataset_info:
            # 构造 dataset 对象（兼容旧的接口）
            class DatasetObject:
                def __init__(self, dataset_info, dataset_id):
                    self.dataset_name = dataset_info.name
                    self.dataset_id = dataset_id  # 添加 dataset_id 属性
                    self.items = []
                    self.partition_scheme = None
                    self.summary = {}
            
            dataset = DatasetObject(dataset_info, dataset_id)
            self.open_dataset_manager_for_dataset(dataset)
    
    def refresh_workspace_layout(self):
        """刷新工作区布局，重新创建所有按钮"""
        # 清除现有的grid_layout内容
        if hasattr(self, 'grid_widget'):
            self.main_layout.removeWidget(self.grid_widget)
            self.grid_widget.deleteLater()
            del self.grid_widget
            del self.grid_layout
        
        # 重新创建grid_layout
        self.grid_widget = QWidget()
        self.grid_layout = QGridLayout()
        self.grid_widget.setLayout(self.grid_layout)
        self.main_layout.addWidget(self.grid_widget)
        self.grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        
        # 重新添加所有项目
        # 添加alignments
        for i, alignment in enumerate(self.items["alignments"]):
            alignment_icon = self.resource_factory.get_icon("file/sequence.svg")
            alignment_button = QToolButton()
            alignment_button.setIcon(alignment_icon)
            alignment_button.setIconSize(QSize(45, 45))
            alignment_button.clicked.connect(lambda checked, seq=alignment: self.view_alignment(seq))
            self.grid_layout.addWidget(alignment_button, 0, i)
        
        # 添加models
        for i, model in enumerate(self.items["models"]):
            model_icon = self.resource_factory.get_icon("file/model.svg")
            model_button = QToolButton()
            model_button.setIcon(model_icon)
            model_button.setIconSize(QSize(45, 45))
            model_button.setToolTip(f"Substitution Model: {getattr(model, 'model_name', 'Unknown Model')}")
            model_button.clicked.connect(lambda checked, m=model: self.view_model_result(m))
            self.grid_layout.addWidget(model_button, 1, i)
        
        # 添加distances
        for i, distance in enumerate(self.items["distances"]):
            distance_icon = self.resource_factory.get_icon("file/distance.svg")
            distance_button = QToolButton()
            distance_button.setIcon(distance_icon)
            distance_button.setIconSize(QSize(45, 45))
            distance_button.setToolTip(f"Pairwise Distance Matrix")
            distance_button.clicked.connect(lambda checked, d=distance: self.view_distance_matrix(d))
            self.grid_layout.addWidget(distance_button, 2, i)
        
        # 添加phylogenies
        for i, phylogeny in enumerate(self.items["phylogenies"]):
            phylogeny_icon = self.resource_factory.get_icon("file/phylogeny.svg")
            phylogeny_button = QToolButton()
            phylogeny_button.setIcon(phylogeny_icon)
            phylogeny_button.setIconSize(QSize(45, 45))
            phylogeny_button.setToolTip(f"Phylogenetic Tree")
            phylogeny_button.clicked.connect(self.open_icytree_wrapper)
            self.grid_layout.addWidget(phylogeny_button, 3, i)
        
        # 添加chains
        for i, chain in enumerate(self.items["chains"]):
            chain_icon = self.resource_factory.get_icon("file/chain.svg")
            chain_button = QToolButton()
            chain_button.setIcon(chain_icon)
            chain_button.setIconSize(QSize(45, 45))
            tooltip_text = f"MCMC Chain (Run {chain.run_number}, Chain {chain.chain_number}, Tool: {chain.tool})"
            chain_button.setToolTip(tooltip_text)
            chain_button.clicked.connect(lambda checked, c=chain: self.open_minitracer_wrapper(chain_item=c))
            self.grid_layout.addWidget(chain_button, 4, i)
        
        # 添加datasets
        for i, dataset in enumerate(self.items["datasets"]):
            dataset_icon = self.resource_factory.get_icon("file/dataset.svg")
            dataset_button = QToolButton()
            dataset_button.setIcon(dataset_icon)
            dataset_button.setIconSize(QSize(45, 45))
            dataset_button.setToolTip(f"Dataset: {getattr(dataset, 'dataset_name', 'Unnamed Dataset')}")
            dataset_row = 0  # Dataset与alignments放在同一行（第0行）
            
            # 实现单击选中，双击查看的交互逻辑
            dataset_button.doubleClicked = False
            dataset_button.mousePressEvent = lambda event, btn=dataset_button: self._on_dataset_mouse_press(btn, event)
            dataset_button.mouseReleaseEvent = lambda event, btn=dataset_button, ds=dataset: self._on_dataset_mouse_release(btn, ds, event)
            dataset_button.mouseDoubleClickEvent = lambda event, btn=dataset_button: self._on_dataset_double_click(btn, event)
            
            # Dataset 的列号从 alignments 数量开始，避免与 alignments 重叠
            dataset_col = len(self.items["alignments"]) + i
            self.grid_layout.addWidget(dataset_button, dataset_row, dataset_col)
            # 存储dataset引用到按钮上，便于后续访问
            dataset_button.dataset_ref = dataset
    
    def _on_dataset_mouse_press(self, button, event):
        """Dataset按钮鼠标按下事件"""
        button.doubleClicked = False
        from PyQt5.QtWidgets import QToolButton
        QToolButton.mousePressEvent(button, event)
    
    def _on_dataset_mouse_release(self, button, dataset, event):
        """Dataset按钮鼠标释放事件"""
        if not button.doubleClicked:
            # 单击：选中Dataset功能模式（这里可以添加选中逻辑）
            pass
        from PyQt5.QtWidgets import QToolButton
        QToolButton.mouseReleaseEvent(button, event)
    
    def _on_dataset_double_click(self, button, event):
        """Dataset按钮双击事件"""
        button.doubleClicked = True
        # 直接使用存储的dataset引用
        if hasattr(button, 'dataset_ref'):
            dataset = button.dataset_ref
            self.open_dataset_manager_for_dataset(dataset)
        from PyQt5.QtWidgets import QToolButton
        QToolButton.mouseDoubleClickEvent(button, event)
    
    def open_dataset_manager_for_dataset(self, dataset):
        """打开Dataset Manager查看特定数据集"""
        try:
            # 直接导入并实例化 DatasetManager 类
            from .methods.dataset_manager import DatasetManager

            # 获取 dataset_id
            dataset_id = None
            if isinstance(dataset, dict):
                dataset_id = dataset.get('dataset_id')
            else:
                dataset_id = getattr(dataset, 'dataset_id', None)

            # 先设置当前 dataset_id（这样DatasetManager初始化时就能找到它）
            if dataset_id:
                self.current_dataset_id = dataset_id

            # 禁用自动保存和selection_engine，防止它们在对话框操作期间修改状态
            original_disable_auto_save = self.dataset_selection_manager.disable_auto_save
            self.dataset_selection_manager.disable_auto_save = True

            # 创建 DatasetManager 实例，传递必要的参数
            dataset_manager = DatasetManager(
                dataset_name=getattr(dataset, 'dataset_name', 'Dataset'),
                plugin_factory=self.parent_window.plugin_factory,
                workspace=self
            )

            # 设置dataset的selection_state为GREEN，表示当前正在编辑
            if dataset_id and self.dataset_selection_manager:
                ds = self.dataset_selection_manager.get_dataset(dataset_id)
                if ds:
                    ds.selection_state = SELECTION_STATE_GREEN

            # 保存引用防止被垃圾回收
            if not hasattr(self.parent_window, 'dataset_managers'):
                self.parent_window.dataset_managers = []
            self.parent_window.dataset_managers.append(dataset_manager)

            dialog = dataset_manager

            # 执行对话框
            dialog.show()  # 改为非模态显示

            # 对话框关闭后恢复自动保存并手动保存一次
            def on_dialog_closed():
                self.dataset_selection_manager.disable_auto_save = original_disable_auto_save
                self.dataset_selection_manager._save_state()
            
            dialog.finished.connect(on_dialog_closed)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Dataset Manager: {str(e)}")
    
    def view_alignment(self, sequences):
        """查看序列比对结果 - 修复序列格式转换问题"""
        try:
            from YR_MPE.sequence_editor import SequenceAlignmentViewer
            
            # 确保序列数据格式正确
            # 如果是Bio.SeqRecord对象列表，需要转换为字典格式
            if sequences and hasattr(sequences[0], 'seq'):
                # 转换Bio.SeqRecord对象为字典格式
                converted_sequences = []
                for seq_record in sequences:
                    seq_dict = {
                        'header': getattr(seq_record, 'id', getattr(seq_record, 'name', 'Unknown')),
                        'sequence': str(seq_record.seq)
                    }
                    converted_sequences.append(seq_dict)
                sequences = converted_sequences
            
            viewer = SequenceAlignmentViewer(sequences)
            viewer.show()
            # 保存viewer引用以防被垃圾回收
            if not hasattr(self, 'viewers'):
                self.viewers = []
            self.viewers.append(viewer)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open alignment viewer: {str(e)}")
    
    def view_model_result(self, model_result):
        """查看模型选择结果 - 增强版支持替换矩阵可视化和分区模型显示"""
        try:
            from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextEdit, QPushButton, QTabWidget, QWidget, QLabel, QTableWidget, QTableWidgetItem, QHeaderView
            from PyQt5.QtCore import Qt
            from PyQt5.QtGui import QColor, QFont
            import re
            
            dialog = QDialog(self)
            dialog.setWindowTitle("Model Selection Result")
            dialog.resize(800, 600)
            
            layout = QVBoxLayout()
            
            # 创建标签页
            tab_widget = QTabWidget()
            layout.addWidget(tab_widget)
            
            # 检测是否为分区模型
            if isinstance(model_result, dict) and model_result.get("type") == "partitioned":
                # 显示分区模型结果
                partition_mode = model_result.get("partition_mode", "")
                partitions = model_result.get("partitions", [])
                statistics = model_result.get("statistics", {})
                
                # 分区模型标签页
                partition_tab = QWidget()
                partition_layout = QVBoxLayout()
                
                # 添加分区模式信息
                mode_text = {
                    "EL": "Edge-linked partition model",
                    "TL": "Edge-linked partition model",
                    "EUL": "Edge-unlinked partition model",
                    "TUL": "Topo-unlinked partition model"
                }.get(partition_mode, "Unknown partition model")
                
                mode_label = QLabel(mode_text)
                partition_layout.addWidget(mode_label)
                
                # 创建分区表格
                partition_table = QTableWidget()
                partition_table.setColumnCount(4)
                partition_table.setHorizontalHeaderLabels(["Partition", "Range", "Best Model", "LogL"])
                partition_table.setRowCount(len(partitions))
                
                # 填充分区数据
                for row, partition in enumerate(partitions):
                    partition_table.setItem(row, 0, QTableWidgetItem(partition.get("name", "")))
                    partition_table.setItem(row, 1, QTableWidgetItem(partition.get("range", "")))
                    partition_table.setItem(row, 2, QTableWidgetItem(partition.get("best_model", "")))
                    partition_table.setItem(row, 3, QTableWidgetItem(str(partition.get("logL", "N/A"))))
                
                # 调整列宽
                partition_table.resizeColumnsToContents()
                partition_table.setAlternatingRowColors(True)
                partition_layout.addWidget(partition_table)
                
                # 添加统计信息
                stats_text = f"Log-likelihood: {statistics.get('logL', 'N/A')} | AICc: {statistics.get('aicc', 'N/A')} | BIC: {statistics.get('bic', 'N/A')}"
                
                stats_label = QLabel(stats_text)
                partition_layout.addWidget(stats_label)
                
                partition_tab.setLayout(partition_layout)
                tab_widget.addTab(partition_tab, "Partition Models")
                
            elif isinstance(model_result, dict) and "type" in model_result and model_result["type"] == "model_table":
                # 显示模型表 - 使用QTableWidget而不是QTextEdit
                headers = model_result.get("headers", ["Model", "AIC", "BIC"])
                model_data = model_result.get("data", [])
                
                # 模型表标签页
                model_table_tab = QWidget()
                model_table_layout = QVBoxLayout()
                
                # 创建表格
                table = QTableWidget()
                table.setColumnCount(len(headers))
                table.setHorizontalHeaderLabels(headers)
                table.setRowCount(len(model_data))
                
                # 填充数据
                for row, model_info in enumerate(model_data):
                    for col, header in enumerate(headers):
                        value = model_info.get(header, "N/A")
                        item = QTableWidgetItem(str(value))
                        # 数值列右对齐
                        if header != "Model":
                            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                        table.setItem(row, col, item)
                
                # 调整列宽
                table.resizeColumnsToContents()
                table.setSortingEnabled(True)
                
                model_table_layout.addWidget(table)
                model_table_tab.setLayout(model_table_layout)
                tab_widget.addTab(model_table_tab, "Model Table")
                
                # 尝试解析替换矩阵（如果可用）
                self.add_substitution_matrix_tab(tab_widget, model_result)
                
            else:
                # 显示单个模型
                content = f"Selected Model: {getattr(model_result, 'model_name', 'Unknown')}\n"
                content += f"Details: {str(model_result)}"
                
                # 模型详情标签页
                model_detail_tab = QWidget()
                model_detail_layout = QVBoxLayout()
                text_edit = QTextEdit()
                text_edit.setReadOnly(True)
                text_edit.setText(content)
                model_detail_layout.addWidget(text_edit)
                model_detail_tab.setLayout(model_detail_layout)
                tab_widget.addTab(model_detail_tab, "Model Details")
                
                # 尝试解析替换矩阵（如果可用）
                self.add_substitution_matrix_tab(tab_widget, model_result)
            
            close_button = QPushButton("Close")
            close_button.clicked.connect(dialog.accept)
            layout.addWidget(close_button)
            
            dialog.setLayout(layout)
            dialog.exec_()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to display model result: {str(e)}")
    
    def add_substitution_matrix_tab(self, tab_widget, model_result):
        """添加替换矩阵标签页"""
        try:
            from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem, QHeaderView
            from PyQt5.QtCore import Qt
            from PyQt5.QtGui import QColor
            
            # 尝试从模型结果中提取替换矩阵数据
            substitution_matrix = self.extract_substitution_matrix(model_result)
            
            if substitution_matrix is not None:
                # 创建替换矩阵标签页
                sub_matrix_tab = QWidget()
                sub_matrix_layout = QVBoxLayout()
                
                # 添加说明标签
                info_label = QLabel("Substitution Matrix (Relative rates)")
                info_label.setAlignment(Qt.AlignCenter)
                sub_matrix_layout.addWidget(info_label)
                
                # 创建表格
                matrix_size = len(substitution_matrix['labels'])
                table = QTableWidget(matrix_size, matrix_size)
                table.setEditTriggers(QTableWidget.NoEditTriggers)
                
                # 设置表头
                for i in range(matrix_size):
                    table.setHorizontalHeaderItem(i, QTableWidgetItem(substitution_matrix['labels'][i]))
                    table.setVerticalHeaderItem(i, QTableWidgetItem(substitution_matrix['labels'][i]))
                
                # 找到最大值用于着色
                max_value = 0.0
                for row in substitution_matrix['matrix']:
                    for val in row:
                        if val > max_value:
                            max_value = val
                
                # 填充数据
                for i in range(matrix_size):
                    for j in range(matrix_size):
                        if i == j:
                            # 对角线通常为0或1
                            item = QTableWidgetItem("0.000")
                        else:
                            value = substitution_matrix['matrix'][i][j]
                            # 格式化为小数点后三位
                            formatted_value = f"{value:.3f}"
                            item = QTableWidgetItem(formatted_value)
                            
                            # 应用单元格着色 - 数值越大，颜色越深
                            if max_value > 0:
                                intensity = min(1.0, value / max_value)
                                # 使用蓝色系着色（替换矩阵常用蓝色）
                                r = int(200 * (1.0 - intensity))  # 红色分量随强度减小
                                g = int(200 * (1.0 - intensity))  # 绿色分量随强度减小  
                                b = int(255 * intensity)          # 蓝色分量随强度增加
                                
                                bg_color = QColor(r, g, b)
                                item.setBackground(bg_color)
                        
                        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                        table.setItem(i, j, item)
                
                # 设置表格属性
                table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
                table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
                table.setSizeAdjustPolicy(QTableWidget.AdjustToContents)
                
                sub_matrix_layout.addWidget(table)
                sub_matrix_tab.setLayout(sub_matrix_layout)
                tab_widget.addTab(sub_matrix_tab, "Substitution Matrix")
                
        except Exception as e:
            # 如果无法解析替换矩阵，不显示该标签页
            pass
    
    def extract_substitution_matrix(self, model_result):
        """从模型结果中提取替换矩阵数据"""
        try:
            # 尝试从字典格式的模型结果中提取
            if isinstance(model_result, dict):
                # 检查是否包含替换矩阵信息
                if 'substitution_matrix' in model_result:
                    return model_result['substitution_matrix']
                elif 'data' in model_result and isinstance(model_result['data'], list):
                    # 检查第一个模型是否包含替换矩阵
                    first_model = model_result['data'][0] if model_result['data'] else {}
                    if 'substitution_matrix' in first_model:
                        return first_model['substitution_matrix']
            
            # 从模型名称推断替换矩阵类型
            model_name = ""
            if isinstance(model_result, dict) and 'Model' in model_result:
                model_name = model_result['Model']
            elif hasattr(model_result, 'model_name'):
                model_name = getattr(model_result, 'model_name', '')
            elif isinstance(model_result, str):
                model_name = model_result
            
            # 解析模型名称
            if model_name:
                # 移除参数部分（+I, +G等）
                base_model = model_name.split('+')[0]
                
                # 核苷酸替换矩阵
                if base_model in ['JC69', 'JC']:
                    # Jukes-Cantor: 所有替换率相等
                    labels = ['A', 'C', 'G', 'T']
                    matrix = [
                        [0.0, 1.0, 1.0, 1.0],
                        [1.0, 0.0, 1.0, 1.0],
                        [1.0, 1.0, 0.0, 1.0],
                        [1.0, 1.0, 1.0, 0.0]
                    ]
                    return {'labels': labels, 'matrix': matrix}
                
                elif base_model in ['K2P', 'K80']:
                    # Kimura 2-parameter: 转换和颠换
                    labels = ['A', 'C', 'G', 'T']
                    # 假设转换率=2.0, 颠换率=1.0
                    matrix = [
                        [0.0, 1.0, 2.0, 1.0],
                        [1.0, 0.0, 1.0, 2.0],
                        [2.0, 1.0, 0.0, 1.0],
                        [1.0, 2.0, 1.0, 0.0]
                    ]
                    return {'labels': labels, 'matrix': matrix}
                
                elif base_model == 'HKY85':
                    # Hasegawa-Kishino-Yano: 类似K2P但有碱基频率
                    labels = ['A', 'C', 'G', 'T']
                    matrix = [
                        [0.0, 1.0, 2.0, 1.0],
                        [1.0, 0.0, 1.0, 2.0],
                        [2.0, 1.0, 0.0, 1.0],
                        [1.0, 2.0, 1.0, 0.0]
                    ]
                    return {'labels': labels, 'matrix': matrix}
                
                elif base_model == 'GTR':
                    # General Time Reversible: 6个不同参数
                    labels = ['A', 'C', 'G', 'T']
                    # 使用示例值
                    matrix = [
                        [0.0, 1.2, 2.5, 1.8],
                        [1.2, 0.0, 3.1, 2.4],
                        [2.5, 3.1, 0.0, 1.5],
                        [1.8, 2.4, 1.5, 0.0]
                    ]
                    return {'labels': labels, 'matrix': matrix}
                
                # 氨基酸替换矩阵（简化显示）
                elif base_model in ['JTT', 'WAG', 'LG', 'Blosum62', 'Dayhoff']:
                    # 显示20种氨基酸的简化矩阵
                    labels = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
                    # 创建对称矩阵（示例数据）
                    size = len(labels)
                    matrix = [[0.0 for _ in range(size)] for _ in range(size)]
                    for i in range(size):
                        for j in range(i+1, size):
                            # 使用简单的距离值
                            value = abs(i - j) * 0.5 + 1.0
                            matrix[i][j] = value
                            matrix[j][i] = value
                    return {'labels': labels, 'matrix': matrix}
            
            # 如果无法确定模型，返回None
            return None
            
        except Exception:
            return None
    
    def view_distance_matrix(self, distance_matrix):
        """查看距离矩阵 - 增强版支持单元格着色和精确到小数点后三位"""
        try:
            from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, QHeaderView, QMenuBar, QMenu, QAction, QFileDialog, QMessageBox
            from PyQt5.QtCore import Qt
            from PyQt5.QtGui import QColor, QFont
            import tempfile
            import os
            
            dialog = QDialog(self)
            dialog.setWindowTitle("Distance Matrix Visualization")
            dialog.resize(800, 600)
            layout = QVBoxLayout(dialog)
            
            # 创建菜单栏
            menubar = QMenuBar(dialog)
            layout.setMenuBar(menubar)
            
            # 创建Export菜单
            export_menu = QMenu("&Export", dialog)
            menubar.addMenu(export_menu)
            
            # 获取距离矩阵数据
            if isinstance(distance_matrix, dict):
                if 'content' in distance_matrix:
                    # 新格式：dist_data 直接包含 content（由 add_distance_matrix_to_workspace 提取）
                    content = distance_matrix['content']
                elif 'data' in distance_matrix:
                    # 旧格式：data[0]['content']
                    content = distance_matrix['data'][0]['content']
                else:
                    content = str(distance_matrix)
            else:
                content = str(distance_matrix)
            
            # 解析距离矩阵内容
            lines = content.strip().split('\n')
            if not lines:
                QMessageBox.warning(self, "Warning", "Invalid distance matrix format.")
                return
                
            # 第一行可能是序列数量或直接是数据
            sequence_names = []
            matrix_values = []
            
            # 尝试解析格式
            if len(lines) > 1:
                # 检查第一行是否为数字（序列数量）
                try:
                    num_sequences = int(lines[0].strip())
                    data_lines = lines[1:num_sequences+1]
                except ValueError:
                    # 第一行不是数字，直接作为数据行处理
                    data_lines = lines
                    num_sequences = len(data_lines)
                
                # 解析数据行
                for i, line in enumerate(data_lines):
                    parts = line.split()
                    if len(parts) == 0:
                        continue
                        
                    # 第一个部分是序列名称
                    seq_name = parts[0]
                    sequence_names.append(seq_name)
                    
                    # 剩余部分是距离值
                    row_values = []
                    for j, val_str in enumerate(parts[1:]):
                        try:
                            val = float(val_str)
                            row_values.append(val)
                        except ValueError:
                            row_values.append(0.0)
                    matrix_values.append(row_values)
            
            if not sequence_names or not matrix_values:
                QMessageBox.warning(self, "Warning", "Could not parse distance matrix data.")
                return
            
            # 确保矩阵是方阵
            num_sequences = len(sequence_names)
            for i in range(len(matrix_values)):
                while len(matrix_values[i]) < num_sequences:
                    matrix_values[i].append(0.0)
                matrix_values[i] = matrix_values[i][:num_sequences]
            
            # 找到最大距离值用于着色
            max_distance = 0.0
            for row in matrix_values:
                for val in row:
                    if val > max_distance:
                        max_distance = val
            
            # 创建表格控件
            table = QTableWidget()
            table.setEditTriggers(QTableWidget.NoEditTriggers)
            table.setRowCount(num_sequences)
            table.setColumnCount(num_sequences)
            
            # 设置表头
            for i in range(num_sequences):
                table.setHorizontalHeaderItem(i, QTableWidgetItem(sequence_names[i]))
                table.setVerticalHeaderItem(i, QTableWidgetItem(sequence_names[i]))
            
            # 填充数据并应用着色
            color_base = "#ba3e45"  # 默认颜色
            for i in range(num_sequences):
                for j in range(num_sequences):
                    if i == j:
                        # 对角线为0
                        item = QTableWidgetItem("0.000")
                    else:
                        if i < len(matrix_values) and j < len(matrix_values[i]):
                            value = matrix_values[i][j]
                            # 格式化为小数点后三位
                            formatted_value = f"{value:.3f}"
                            item = QTableWidgetItem(formatted_value)
                            
                            # 应用单元格着色 - 数值越大，颜色越深
                            if max_distance > 0:
                                intensity = min(1.0, value / max_distance)
                                # 将十六进制颜色转换为RGB
                                r = int(color_base[1:3], 16) - 255
                                g = int(color_base[3:5], 16) - 255
                                b = int(color_base[5:7], 16) - 255
                                
                                # 计算着色后的颜色（保持色调，调整亮度）
                                # 更亮的颜色表示更小的值，更暗的颜色表示更大的值
                                brightness_factor = intensity  # 0.3到1.0的范围
                                new_r = min(255, 255 + int(r * brightness_factor))
                                new_g = min(255, 255 + int(g * brightness_factor))
                                new_b = min(255, 255 + int(b * brightness_factor))
                                
                                bg_color = QColor(new_r, new_g, new_b)
                                item.setBackground(bg_color)
                        else:
                            item = QTableWidgetItem("0.000")
                    
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    table.setItem(i, j, item)
            
            # 设置表格属性
            table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
            table.setSizeAdjustPolicy(QTableWidget.AdjustToContents)
            
            layout.addWidget(table)
            
            # 添加导出选项
            export_csv_action = QAction("&To CSV", dialog)
            export_csv_action.triggered.connect(lambda: self.export_distance_matrix_to_csv(table, sequence_names))
            export_menu.addAction(export_csv_action)
            
            export_xlsx_action = QAction("To &XLSX", dialog)
            export_xlsx_action.triggered.connect(lambda: self.export_distance_matrix_to_xlsx(table, sequence_names))
            export_menu.addAction(export_xlsx_action)
            
            dialog.exec_()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to display distance matrix: {str(e)}")
    
    def export_distance_matrix_to_csv(self, table, sequence_names):
        """导出距离矩阵到CSV文件"""
        file_path, _ = QFileDialog.getSaveFileName(None, "Save Distance Matrix to CSV", "", "CSV Files (*.csv)")
        if file_path:
            try:
                with open(file_path, 'w', encoding='utf-8') as f:
                    # 写入表头
                    f.write("," + ",".join(sequence_names) + "\n")
                    
                    # 写入数据行
                    for i in range(table.rowCount()):
                        row_data = [sequence_names[i]]
                        for j in range(table.columnCount()):
                            item = table.item(i, j)
                            row_data.append(item.text() if item else "")
                        f.write(",".join(row_data) + "\n")
                
                QMessageBox.information(None, "Success", f"Distance matrix exported to {file_path}")
            except Exception as e:
                QMessageBox.critical(None, "Error", f"Failed to export distance matrix:\n{str(e)}")
                
    def export_distance_matrix_to_xlsx(self, table, sequence_names):
        """导出距离矩阵到XLSX文件"""
        try:
            import pandas as pd
            import numpy as np
            
            file_path, _ = QFileDialog.getSaveFileName(None, "Save Distance Matrix to Excel", "", "Excel Files (*.xlsx)")
            if file_path:
                # 创建数据矩阵
                data = []
                for i in range(table.rowCount()):
                    row_data = []
                    for j in range(table.columnCount()):
                        item = table.item(i, j)
                        row_data.append(item.text() if item else "")
                    data.append(row_data)
                
                # 创建DataFrame
                df = pd.DataFrame(data, columns=sequence_names, index=sequence_names)
                
                # 保存到Excel
                df.to_excel(file_path)
                QMessageBox.information(None, "Success", f"Distance matrix exported to {file_path}")
                
        except ImportError:
            QMessageBox.warning(None, "Warning", "pandas library is required for Excel export. Please install pandas to use this feature.")
        except Exception as e:
            QMessageBox.critical(None, "Error", f"Failed to export distance matrix:\n{str(e)}")
    
    def view_phylogeny(self, phylogeny_data):
        """查看系统发育树 - 使用 IcyTree 查看器"""
        try:
            if self.parent_window:
                self.parent_window.open_icytree_wrapper(phylogeny_data=phylogeny_data)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to view phylogeny: {str(e)}")
    
    def open_icytree_wrapper(self):
        """打开IcyTree查看系统发育树"""
        from PyQt5.QtWidgets import QDialog
        try:
            dialog = QDialog()
            dialog.setWindowTitle("IcyTree - YR-MPEA")
            dialog.setWindowIcon(self.parent_window.resource_factory.get_icon("software/icytree.svg"))
            dialog.setMinimumSize(800, 600)
            dialog.setLayout(QVBoxLayout())
            
            # 使用PluginFactory获取IcyTree插件
            plugin_entry = self.parent_window.plugin_factory.get_icytree_plugin()
            icytree_wrapper = plugin_entry.run()
            
            # 从工作区获取最新的系统树数据并传递给IcyTree
            if len(self.items["phylogenies"]) > 0:
                latest_phylogeny = self.items["phylogenies"][-1]
                # 检查系统树数据格式
                if isinstance(latest_phylogeny, dict) and 'data' in latest_phylogeny and len(latest_phylogeny['data']) > 0:
                    # 如果系统树数据是字典格式，提取Newick字符串
                    tree_data = latest_phylogeny['data'][0]
                    if 'content' in tree_data:
                        icytree_wrapper.set_newick_string(tree_data['content'])
                elif isinstance(latest_phylogeny, str):
                    # 如果系统树数据直接是字符串
                    icytree_wrapper.set_newick_string(latest_phylogeny)
            
            dialog.layout().addWidget(icytree_wrapper)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open IcyTree: {str(e)}")

class YR_MPEA_entry:
    def run(self):
        return YR_MPEA_Widget()

if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication, QMainWindow
    app = QApplication(sys.argv)
    widget = YR_MPEA_Widget()
    widget.show()
    sys.exit(app.exec_())