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
from Bio import SeqIO
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout,
QMenuBar, QToolBar, QToolButton, QGroupBox, QLabel,
QAction, QMenu, QSizePolicy, QGridLayout, QFileDialog, QMessageBox, QDialog)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon

# 导入新的模块架构
from .methods import (
    WorkspaceManager, PluginManager, PluginExecutor, 
    FileOperations, UIHelpers
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
            # 这里可以添加后续处理逻辑，比如保存工作区路径或更新UI
            # 目前先简单显示一个消息框确认选择
            QMessageBox.information(
                self, 
                "Workspace Selected", 
                f"Workspace folder selected:\n{folder_path}"
            )
            # TODO: 可以在这里添加实际的工作区处理逻辑
            
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

        # 添加Select workspace folder选项
        select_workspace_action = QAction("Select workspace folder", workspace_button)
        select_workspace_action.triggered.connect(self.select_workspace_folder)
        workspace_button.addAction(select_workspace_action)
        
        # 添加分隔符
        separator = QAction(workspace_button)
        separator.setSeparator(True)
        workspace_button.addAction(separator)
        
        # Dataset子菜单 - 使用正确的Qt API
        dataset_action = QAction("Dataset", workspace_button)
        dataset_menu = QMenu(workspace_button)
        
        new_dataset_action = QAction("New", dataset_menu)
        new_dataset_action.setIcon(self.resource_factory.get_icon("new.svg"))
        new_dataset_action.triggered.connect(self.create_new_dataset)
        dataset_menu.addAction(new_dataset_action)
        
        import_nexus_action = QAction("Import from partitioned NEXUS file", dataset_menu)
        import_nexus_action.triggered.connect(self.import_dataset_from_nexus)
        dataset_menu.addAction(import_nexus_action)
        
        create_seqmatrix_action = QAction("Create by SeqMatrix", dataset_menu)
        create_seqmatrix_action.triggered.connect(self.create_dataset_by_seqmatrix)
        dataset_menu.addAction(create_seqmatrix_action)
        
        create_seqdbg_action = QAction("Create by SeqDBG", dataset_menu)
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

        # TODO: BI Programs: MrBayes

        # cons_mp_action = QAction("Cladistics - Maximum Parsimony Phylogenies (MP)", phylogeny_button)
        # cons_mp_action.setIcon(self.resource_factory.get_icon("mp.svg"))

        # TODO: MP Programs: TNT

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

        # TODO: Astral-III and other coalescent methods

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
        self.workspace = SingleGeneWorkspace(resource_factory=self.resource_factory, parent=self)
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
        
        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]
        
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

    def open_clustal_omega_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("Clustal Omega - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/clustalo.svg"))
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

        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]
        
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

        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]
        
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

        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]
        
        # use PluginFactory to get the plugin
        model_finder_entry = self.plugin_factory.get_model_finder_plugin()
        modelfinder_wrapper = model_finder_entry.run(import_from=import_from, import_data=import_data)
        modelfinder_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        modelfinder_wrapper.export_model_result_signal.connect(self.add_model_to_workspace)
        dialog.layout().addWidget(modelfinder_wrapper)
        dialog.exec_()
    
    def open_icytree_wrapper(self):
        """打开IcyTree查看系统发育树"""
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
            
            # # 如果工作区中有系统树数据，则传递给IcyTree
            # if hasattr(self, 'workspace') and len(self.workspace.items["phylogenies"]) > 0:
            #     latest_phylogeny = self.workspace.items["phylogenies"][-1]
            #     # 检查系统树数据格式
            #     if isinstance(latest_phylogeny, dict) and 'data' in latest_phylogeny and len(latest_phylogeny['data']) > 0:
            #         # 如果系统树数据是字典格式，提取Newick字符串
            #         tree_data = latest_phylogeny['data'][0]
            #         if 'content' in tree_data:
            #             icytree_wrapper.set_newick_string(tree_data['content'])
            #     elif isinstance(latest_phylogeny, str):
            #         # 如果系统树数据直接是字符串
            #         icytree_wrapper.set_newick_string(latest_phylogeny)
            
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

        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]
        
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

        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]
        
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
        dialog = QDialog()
        dialog.setWindowTitle(f"Distance Calculator [implemented from IQ-TREE] - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("dist.svg"))
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

        # Use PluginFactory to get the plugin
        ml_distance_entry = self.plugin_factory.get_ml_distance_plugin()
        plugin_wrapper = ml_distance_entry.run(import_from=import_from, import_data=import_data)
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_distance_result_signal.connect(self.add_distance_matrix_to_workspace)

        dialog.layout().addWidget(plugin_wrapper)
        dialog.exec_()
        
    def open_uncorrected_distance_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle(f"Uncorrected Distance Calculator - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("dist.svg"))
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

        # Use PluginFactory to get the plugin
        p_distance_entry = self.plugin_factory.get_p_distance_plugin()
        plugin_wrapper = p_distance_entry.run(import_from=import_from, import_data=import_data)
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_distance_result_signal.connect(self.add_distance_matrix_to_workspace)

        dialog.layout().addWidget(plugin_wrapper)
        dialog.exec_()
    def open_iqtree_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle(f"IQ-TREE 3 - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/iqtree.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        # Prepare import data
        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        best_model = ""
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]

                # 如果工作区中有模型信息，则添加到导入数据中
                if len(self.workspace.items["models"]) >= 1:
                    model_data = self.workspace.items["models"][0]
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
        plugin_wrapper = iqtree_entry.run(import_from=import_from, import_data=import_data)
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_model_result_signal.connect(self.add_model_to_workspace)
        plugin_wrapper.export_phylogeny_result_signal.connect(self.add_phylogeny_to_workspace)

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
        
        # Gamma Caterories [identify Gx]
        for item in model_entries[1:]:
            if item.startswith("G"):
                plugin_wrapper.gamma_checkbox.setChecked(True)
                plugin_wrapper.gamma_spinbox.setValue(int(item[1:]))
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
                    # TODO: 实现从NEXUS文件解析Dataset的逻辑
                    QMessageBox.information(
                        self, "Import Dataset", 
                        f"Dataset imported from NEXUS file: {nexus_file}\n"
                        "Note: Full implementation requires NEXUS parser integration."
                    )
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to import dataset from NEXUS: {str(e)}")
    
    def create_dataset_by_seqmatrix(self):
        """通过SeqMatrix创建Dataset"""
        try:
            # TODO: 实现SeqMatrix插件调用
            QMessageBox.information(
                self, "Create Dataset", 
                "Creating dataset using SeqMatrix plugin...\n"
                "Note: SeqMatrix plugin integration pending."
            )
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create dataset by SeqMatrix: {str(e)}")
    
    def create_dataset_by_seqdbg(self):
        """通过SeqDBG创建Dataset"""
        try:
            # TODO: 实现SeqDBG插件调用
            QMessageBox.information(
                self, "Create Dataset", 
                "Creating dataset using SeqDBG plugin...\n"
                "Note: SeqDBG plugin integration pending."
            )
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create dataset by SeqDBG: {str(e)}")

class SingleGeneWorkspace(QWidget):
    def __init__(self, resource_factory=None, parent=None):
        super().__init__()
        self.resource_factory = resource_factory  # 添加resource_factory引用
        self.parent_window = parent  # 添加parent引用
        self.items = {
            "alignments": [],
            "models": [],
            "distances": [],
            "phylogenies": [],
            "datasets": [],  # 添加datasets项
            # "variants": [], TODO
            # "coalescent": [], TODO
            # "clock": [], TODO
        }
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.init_ui()

    def init_ui(self):
        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)
        # background color: white
        # self.setStyleSheet("background-color: white;")
        # add a label: "Single Gene Workspace"
        self.workspace_hint = QLabel("Single Gene Workspace\nAdd an alignment or drag and drop a file here to start")
        self.workspace_hint.setAlignment(Qt.AlignCenter)
        self.workspace_hint.setStyleSheet("color: #555555;")
        self.main_layout.addWidget(self.workspace_hint)

        # TODO: drag and drop event
    
    def add_sequence(self, sequences):
        # judge if there's already an alignment
        if len(self.items["alignments"]) > 0:
            # replace the existing alignment or not
            replace_alignment = QMessageBox.question(self, "Replace Alignment", "Do you want to replace the existing alignment?", QMessageBox.Yes | QMessageBox.No)
            if replace_alignment == QMessageBox.Yes:
                self.items["alignments"] = []
            else:
                return
        else:
            # disable hint label
            self.workspace_hint.setVisible(False)
            # add a grid layout
            self.grid_widget = QWidget()
            self.grid_layout = QGridLayout()
            self.grid_widget.setLayout(self.grid_layout)
            self.main_layout.addWidget(self.grid_widget)

            # align to top-left
            self.grid_layout.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        
        # add the sequence set to items["alignments"]
        self.items["alignments"].append(sequences)
        # add an 'alignment' icon to workspace
        alignment_icon = self.resource_factory.get_icon("file/sequence.svg")
        alignment_button = QToolButton()
        alignment_button.setIcon(alignment_icon)
        alignment_button.setIconSize(QSize(45, 45))
        # alignment_button.setToolTip("Alignment")
        alignment_button.clicked.connect(lambda: self.view_alignment(sequences))
        self.grid_layout.addWidget(alignment_button, 0, len(self.items["alignments"])-1)
    
    def add_model(self, model):
        # 检查是单个模型还是模型表
        if isinstance(model, dict) and "type" in model and model["type"] == "model_table":
            # 处理模型表
            model_name = model.get("name", "Model Table")
            self.items["models"].append(model)
            
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
            
            # add a 'model' icon to workspace
            model_icon = self.resource_factory.get_icon("file/model.svg")
            model_button = QToolButton()
            model_button.setIcon(model_icon)
            model_button.setIconSize(QSize(45, 45))
            model_button.setToolTip(f"Substitution Model: {model_name}")
            model_button.clicked.connect(lambda: self.view_model_result(model))
            self.grid_layout.addWidget(model_button, 1, len(self.items["models"])-1)
        else:
            # 处理单个模型（保持向后兼容）
            model_name = getattr(model, 'model_name', 'Unknown Model')
            self.items["models"].append(model)
            
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
            
            # add a 'model' icon to workspace
            model_icon = self.resource_factory.get_icon("file/model.svg")
            model_button = QToolButton()
            model_button.setIcon(model_icon)
            model_button.setIconSize(QSize(45, 45))
            model_button.setToolTip(f"Substitution Model: {model_name}")
            model_button.clicked.connect(lambda: self.view_model_result(model))
            self.grid_layout.addWidget(model_button, 1, len(self.items["models"])-1)
    
    def add_distance(self, distance):
        # add a distance to items["distances"]
        self.items["distances"].append(distance)
        
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
        
        # add a 'distance' icon to workspace
        distance_icon = self.resource_factory.get_icon("file/distance.svg")
        distance_button = QToolButton()
        distance_button.setIcon(distance_icon)
        distance_button.setIconSize(QSize(45, 45))
        distance_button.setToolTip(f"Pairwise Distance Matrix")
        distance_button.clicked.connect(lambda: self.view_distance_matrix(distance))
        self.grid_layout.addWidget(distance_button, 2, len(self.items["distances"])-1)
    
    def add_phylogeny(self, phylogeny):
        # add a phylogeny to items["phylogenies"]
        self.items["phylogenies"].append(phylogeny)
        # add a phylogeny icon to workspace
        phylogeny_icon = self.resource_factory.get_icon("file/phylogeny.svg")
        phylogeny_button = QToolButton()
        phylogeny_button.setIcon(phylogeny_icon)
        phylogeny_button.setIconSize(QSize(45, 45))
        phylogeny_button.setToolTip(f"Phylogenetic Tree")
        phylogeny_button.clicked.connect(self.open_icytree_wrapper)
        self.grid_layout.addWidget(phylogeny_button, 3, len(self.items["phylogenies"])-1)
    
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
        
        # add a dataset icon to workspace
        dataset_icon = self.resource_factory.get_icon("file/dataset.svg")
        dataset_button = QToolButton()
        dataset_button.setIcon(dataset_icon)
        dataset_button.setIconSize(QSize(45, 45))
        dataset_button.setToolTip(f"Dataset: {getattr(dataset, 'dataset_name', 'Unnamed Dataset')}")
        # Dataset数据项应该显示在数据区域，使用合适的行号
        # 根据老版实现，Dataset应该有自己的行（比如第4行，因为0-3已被其他功能占用）
        dataset_row = 4
        
        # 实现单击选中，双击查看的交互逻辑
        dataset_button.doubleClicked = False
        dataset_button.mousePressEvent = lambda event, btn=dataset_button: self._on_dataset_mouse_press(btn, event)
        dataset_button.mouseReleaseEvent = lambda event, btn=dataset_button, ds=dataset: self._on_dataset_mouse_release(btn, ds, event)
        dataset_button.mouseDoubleClickEvent = lambda event, btn=dataset_button: self._on_dataset_double_click(btn, event)
        
        self.grid_layout.addWidget(dataset_button, dataset_row, len(self.items["datasets"])-1)
        
        # 存储dataset引用到按钮上，便于后续访问
        dataset_button.dataset_ref = dataset
    
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
        
        # 添加datasets
        for i, dataset in enumerate(self.items["datasets"]):
            dataset_icon = self.resource_factory.get_icon("file/dataset.svg")
            dataset_button = QToolButton()
            dataset_button.setIcon(dataset_icon)
            dataset_button.setIconSize(QSize(45, 45))
            dataset_button.setToolTip(f"Dataset: {getattr(dataset, 'dataset_name', 'Unnamed Dataset')}")
            dataset_row = 4
            
            # 实现单击选中，双击查看的交互逻辑
            dataset_button.doubleClicked = False
            dataset_button.mousePressEvent = lambda event, btn=dataset_button: self._on_dataset_mouse_press(btn, event)
            dataset_button.mouseReleaseEvent = lambda event, btn=dataset_button, ds=dataset: self._on_dataset_mouse_release(btn, ds, event)
            dataset_button.mouseDoubleClickEvent = lambda event, btn=dataset_button: self._on_dataset_double_click(btn, event)
            
            self.grid_layout.addWidget(dataset_button, dataset_row, i)
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
            # 使用PluginFactory获取DatasetManager插件
            dataset_manager = self.parent_window.plugin_factory.get_dataset_manager()
            
            # 保存引用防止被垃圾回收
            if not hasattr(self.parent_window, 'dataset_managers'):
                self.parent_window.dataset_managers = []
            self.parent_window.dataset_managers.append(dataset_manager)
            
            dialog = dataset_manager
            
            # 设置必要的参数
            dialog.dataset_name = getattr(dataset, 'dataset_name', 'Dataset')
            dialog.plugin_factory = self.parent_window.plugin_factory
            dialog.workspace = self
            
            # 如果dataset包含items数据，加载到dialog中
            if hasattr(dataset, 'items') and dataset.items:
                dialog.dataset_items = dataset.items
                # 重新填充表格
                dialog.table.setRowCount(0)
                for item in dataset.items:
                    dialog.add_dataset_to_table(item)
            
            # 执行对话框
            dialog.show()  # 改为非模态显示
            # 无论用户如何关闭对话框，都保存数据
            # 因为DatasetManager没有明确的"Cancel"按钮
            if hasattr(dialog, 'dataset_items'):
                dataset.items = dialog.dataset_items
                    
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
        """查看模型选择结果 - 增强版支持替换矩阵可视化"""
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
            
            # 格式化模型结果显示
            if isinstance(model_result, dict) and "type" in model_result and model_result["type"] == "model_table":
                # 显示模型表
                content = f"Model Selection Results for {model_result.get('name', 'Unknown')}\n\n"
                content += "Rank\tModel\tAIC\tBIC\tWeight\n"
                content += "-" * 50 + "\n"
                for i, model_info in enumerate(model_result.get("models", [])):
                    content += f"{i+1}\t{model_info.get('model', 'N/A')}\t{model_info.get('aic', 'N/A')}\t{model_info.get('bic', 'N/A')}\t{model_info.get('weight', 'N/A')}\n"
                
                # 模型表标签页
                model_table_tab = QWidget()
                model_table_layout = QVBoxLayout()
                text_edit = QTextEdit()
                text_edit.setReadOnly(True)
                text_edit.setText(content)
                model_table_layout.addWidget(text_edit)
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