# YR-MPEA
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

import os
from Bio import SeqIO
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, 
QMenuBar, QToolBar, QToolButton, QGroupBox, QLabel,
QAction, QMenu, QSizePolicy, QGridLayout, QFileDialog, QMessageBox)
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
        build_action = QAction("New Alignment", align_button)
        build_action.setIcon(self.resource_factory.get_icon("new.svg"))
        align_action = QAction("Align by...", align_button)
        align_action.setIcon(self.resource_factory.get_icon("alignby.svg"))
        trim_action = QAction("Trim Alignment by...", align_button)
        trim_action.setIcon(self.resource_factory.get_icon("trim.svg"))
        
        align_button.addAction(open_action)
        align_button.addAction(build_action)
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
        aligners_menu.addAction(QAction("MACSE", align_button))
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
        
        # 添加"Compute Overall Mean Distances"动作
        self.comp_mdist_action = QAction("Compute Overall Mean Distances (ML)", self)
        self.comp_mdist_action.setIcon(self.resource_factory.get_icon("mdist.svg"))
        # TODO: 实现平均距离计算功能
        # self.comp_mdist_action.triggered.connect(self.compute_mean_distances)
        distance_button.addAction(self.comp_mdist_action)

        phylogeny_button = QToolButton()
        phylogeny_button.setText("PHYLOGENY")
        phylogeny_button.setIcon(self.resource_factory.get_icon("phylogeny.svg"))
        phylogeny_button.setPopupMode(QToolButton.InstantPopup)
        phylogeny_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        phylogeny_button_menu = QMenu()
        phylogeny_button.setMenu(phylogeny_button_menu)
        # phylogeny_button.setFixedSize(60, 60)
        main_toolbar.addWidget(phylogeny_button)

        cons_ml_action = QAction("Phenetics - Maximum Likelihood Phylogenies (ML)", phylogeny_button)
        cons_ml_action.setIcon(self.resource_factory.get_icon("ml.svg"))

        ml_menu = QMenu()
        cons_ml_action.setMenu(ml_menu)
        
        # Add IQ-TREE action
        iqtree_action = QAction("IQ-Tree 3", phylogeny_button)
        iqtree_action.setIcon(self.resource_factory.get_icon("software/iqtree.svg"))
        iqtree_action.triggered.connect(self.open_iqtree_wrapper)

        ml_menu.addAction(iqtree_action)
        
        # TODO: ML Programs: IQ-TREE 3 / FastTree

        cons_bi_action = QAction("Phenetics - Bayesian Inference Phylogenies (BI)", phylogeny_button)
        cons_bi_action.setIcon(self.resource_factory.get_icon("bi.svg"))

        # TODO: BI Programs: MrBayes

        cons_mp_action = QAction("Cladistics - Maximum Parsimony Phylogenies (MP)", phylogeny_button)
        cons_mp_action.setIcon(self.resource_factory.get_icon("mp.svg"))

        # TODO: MP Programs: TNT

        cons_nj_action = QAction("Clustering - Neighbor-Joining (NJ)", phylogeny_button)
        cons_nj_action.setIcon(self.resource_factory.get_icon("nj.svg"))

        cons_upgma_action = QAction("Clustering - UPGMA", phylogeny_button)
        cons_upgma_action.setIcon(self.resource_factory.get_icon("upgma.svg"))

        cons_me_action = QAction("Clustering - Minimum Evolution (ME)", phylogeny_button)
        cons_me_action.setIcon(self.resource_factory.get_icon("me.svg"))

        cons_bionj_action = QAction("Clustering - BioNJ [IQ-TREE 3]", phylogeny_button)
        cons_bionj_action.setIcon(self.resource_factory.get_icon("bionj.svg"))

        cons_rand_action = QAction("Simulation - Random Trees (Evolver)", phylogeny_button)
        cons_rand_action.setIcon(self.resource_factory.get_icon("randtree.svg"))

        cons_seq_action = QAction("Simulation - Simulated Sequences (Evolver)", phylogeny_button)
        cons_seq_action.setIcon(self.resource_factory.get_icon("randseq.svg"))

        tree_viewer_action = QAction("Tree Viewer (IcyTree)", phylogeny_button)
        tree_viewer_action.setIcon(self.resource_factory.get_icon("software/icytree.svg"))
        tree_viewer_action.triggered.connect(self.open_icytree_wrapper)

        phylogeny_button_menu.addAction(cons_ml_action)
        phylogeny_button_menu.addAction(cons_bi_action)
        phylogeny_button_menu.addSeparator()
        phylogeny_button_menu.addAction(cons_mp_action)
        phylogeny_button_menu.addSeparator()
        phylogeny_button_menu.addAction(cons_nj_action)
        phylogeny_button_menu.addAction(cons_upgma_action)
        phylogeny_button_menu.addAction(cons_me_action)
        # phylogeny_button_menu.addSeparator()
        phylogeny_button_menu.addAction(cons_bionj_action)
        phylogeny_button_menu.addSeparator()
        phylogeny_button_menu.addAction(cons_rand_action)
        phylogeny_button_menu.addAction(cons_seq_action)
        phylogeny_button_menu.addSeparator()
        phylogeny_button_menu.addAction(tree_viewer_action)
        # phylogeny_button_menu.addSeparator()

        variants_button = QToolButton()
        variants_button.setText("VARIANTS")
        variants_button.setIcon(self.resource_factory.get_icon("variants.svg"))
        variants_button.setPopupMode(QToolButton.InstantPopup)
        variants_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        variants_button_menu = QMenu()
        variants_button.setMenu(variants_button_menu)
        # phylogeny_button.setFixedSize(60, 60)
        main_toolbar.addWidget(variants_button)

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
        self.workspace = SingleGeneWorkspace()
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
            from ..plugins.dataset_plugin import DatasetPluginEntry
            from PyQt5.QtWidgets import QDialog
            dialog = QDialog()
            dialog.setWindowTitle("Dataset Manager - YR-MPEA")
            dialog.setMinimumSize(800, 600)
            dialog.setLayout(QVBoxLayout())
            
            dataset_plugin = DatasetPluginEntry().run()
            # 连接信号以便在Dataset管理器中查看序列
            if hasattr(dataset_plugin, 'view_sequence_signal'):
                dataset_plugin.view_sequence_signal.connect(self.view_sequence_in_viewer)
            
            dialog.layout().addWidget(dataset_plugin)
            dialog.exec_()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Dataset Manager: {str(e)}")
    
    def view_sequence_in_viewer(self, sequences):
        """使用SeqViewer插件查看序列"""
        try:
            sequence_viewer = self.plugin_factory.get_sequence_viewer()
            viewer = sequence_viewer.run(sequences)
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
        from YR_MPE.icytree import IcyTreePlugin
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("IcyTree - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/icytree.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())
        icytree_wrapper = IcyTreePlugin()
        dialog.layout().addWidget(icytree_wrapper)
        dialog.exec_()

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
        dialog.setWindowIcon(self.resource_factory.get_icon("distance.svg"))
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
        """创建新的Dataset"""
        try:
            from .methods.dataset_manager import DatasetManager
            dialog = DatasetManager("New Dataset")
            dialog.exec_()
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
    def __init__(self):
        super().__init__()
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
            
            # 添加工作区主区块的图标按钮（第一行）
            self.add_workspace_main_icons()
        
        # add the sequence set to items["alignments"]
        self.items["alignments"].append(sequences)
        # add an 'alignment' icon to workspace - 这里应该是第二行及以后的内容
        alignment_icon = self.resource_factory.get_icon("file/sequence.svg")
        alignment_button = QToolButton()
        alignment_button.setIcon(alignment_icon)
        alignment_button.setIconSize(QSize(45, 45))
        # alignment_button.setToolTip("Alignment")
        alignment_button.clicked.connect(lambda: self.view_alignment(sequences))
        # 添加到第二行第一列
        self.grid_layout.addWidget(alignment_button, 1, 0)
    
    def add_workspace_main_icons(self):
        """在工作区主区块第一行添加图标按钮"""
        # Sequence (Alignments)
        seq_button = QToolButton()
        seq_button.setIcon(self.resource_factory.get_icon("file/sequence.svg"))
        seq_button.setIconSize(QSize(45, 45))
        seq_button.setToolTip("Sequence/Alignment")
        seq_button.clicked.connect(self.open_sequence_files)
        self.grid_layout.addWidget(seq_button, 0, 0)
        
        # Model
        model_button = QToolButton()
        model_button.setIcon(self.resource_factory.get_icon("file/model.svg"))
        model_button.setIconSize(QSize(45, 45))
        model_button.setToolTip("Model")
        model_button.clicked.connect(self.open_model_selection)
        self.grid_layout.addWidget(model_button, 0, 1)
        
        # Distance
        distance_button = QToolButton()
        distance_button.setIcon(self.resource_factory.get_icon("file/distance.svg"))
        distance_button.setIconSize(QSize(45, 45))
        distance_button.setToolTip("Distance")
        distance_button.clicked.connect(self.open_distance_calculation)
        self.grid_layout.addWidget(distance_button, 0, 2)
        
        # Tree (Phylogeny)
        tree_button = QToolButton()
        tree_button.setIcon(self.resource_factory.get_icon("file/phylogeny.svg"))
        tree_button.setIconSize(QSize(45, 45))
        tree_button.setToolTip("Tree")
        tree_button.clicked.connect(self.open_tree_building)
        self.grid_layout.addWidget(tree_button, 0, 3)
        
        # Dataset - NEW (无文字，只有图标)
        dataset_button = QToolButton()
        dataset_button.setIcon(self.resource_factory.get_icon("file/dataset.svg"))  # 使用专门的dataset图标
        dataset_button.setIconSize(QSize(45, 45))
        dataset_button.setToolTip("Dataset")
        dataset_button.setCheckable(True)
        dataset_button.clicked.connect(lambda: self.handle_dataset_click(dataset_button))
        dataset_button.doubleClicked = False
        
        # Override mouse events for double-click functionality
        dataset_button.mouseDoubleClickEvent_original = dataset_button.mouseDoubleClickEvent
        dataset_button.mousePressEvent_original = dataset_button.mousePressEvent
        
        def dataset_mouse_double_click(event):
            dataset_button.doubleClicked = True
            self.open_dataset_manager()
            dataset_button.mouseDoubleClickEvent_original(event)
            
        def dataset_mouse_press(event):
            dataset_button.doubleClicked = False
            dataset_button.mousePressEvent_original(event)
            
        dataset_button.mouseDoubleClickEvent = dataset_mouse_double_click
        dataset_button.mousePressEvent = dataset_mouse_press
        
        self.grid_layout.addWidget(dataset_button, 0, 4)
    
    def handle_dataset_click(self, button):
        """处理Dataset按钮点击事件"""
        if not button.doubleClicked:
            # 单击：选中Dataset功能模式
            button.setChecked(True)
    
    def open_dataset_manager(self):
        """打开Dataset管理对话框"""
        try:
            from ..plugins.dataset_plugin import DatasetPluginEntry
            from PyQt5.QtWidgets import QDialog
            dialog = QDialog()
            dialog.setWindowTitle("Dataset Manager - YR-MPEA")
            dialog.setMinimumSize(800, 600)
            dialog.setLayout(QVBoxLayout())
            
            dataset_plugin = DatasetPluginEntry().run()
            # 连接信号以便在Dataset管理器中查看序列
            dataset_plugin.view_sequence_signal = getattr(dataset_plugin, 'view_sequence_signal', None)
            if hasattr(dataset_plugin, 'view_sequence_signal'):
                dataset_plugin.view_sequence_signal.connect(self.view_sequence_in_viewer)
            
            dialog.layout().addWidget(dataset_plugin)
            dialog.exec_()
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
    
    def open_sequence_files(self):
        """打开序列文件"""
        file_dialog = QFileDialog()
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
            from Bio import SeqIO
            file_format = file_types[file.split(".")[-1].lower()]
            sequences = list(SeqIO.parse(file, file_format))
            
            if not sequences:
                QMessageBox.warning(self, "Error", "No sequences found in the file")
                return
            
            self.add_sequence(sequences)
            QMessageBox.information(self, "Success", f"Successfully loaded {len(sequences)} sequence(s) from {file}")
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error opening file: {e}")
            return
    
    def open_model_selection(self):
        """打开模型选择"""
        if len(self.items["alignments"]) == 0:
            QMessageBox.warning(self, "Warning", "Please load an alignment first.")
            return
        self.open_modelfinder_wrapper()
    
    def open_distance_calculation(self):
        """打开距离计算"""
        if len(self.items["alignments"]) == 0:
            QMessageBox.warning(self, "Warning", "Please load an alignment first.")
            return
        self.open_ml_distance_wrapper()
    
    def open_tree_building(self):
        """打开树构建"""
        if len(self.items["alignments"]) == 0:
            QMessageBox.warning(self, "Warning", "Please load an alignment first.")
            return
        self.open_iqtree_wrapper()
    
    def add_model(self, model):
        # 检查是单个模型还是模型表
        if isinstance(model, dict) and "type" in model and model["type"] == "model_table":
            # 处理完整模型表
            self.items["models"].append(model)
            # 添加模型表图标到工作区
            model_icon = QIcon(os.path.join(self.plugin_path, "icons/file/model.svg"))
            model_button = QToolButton()
            model_button.setIcon(model_icon)
            model_button.setIconSize(QSize(45, 45))
            model_button.setToolTip(f"Model Table ({len(model['data'])} models)")
            model_button.clicked.connect(lambda: self.view_model_table(model))
            self.grid_layout.addWidget(model_button, 1, 0)
        else:
            # 处理单个模型（向后兼容）
            self.items["models"].append(model)
            # add a model icon to workspace
            model_icon = QIcon(os.path.join(self.plugin_path, "icons/file/model.svg"))
            model_button = QToolButton()
            model_button.setIcon(model_icon)
            model_button.setIconSize(QSize(45, 45))
            self.grid_layout.addWidget(model_button, 1, 0)
    
    def add_distance(self, distance):
        # add a distance to items["distances"]
        self.items["distances"].append(distance)
        # add a distance icon to workspace
        distance_icon = QIcon(os.path.join(self.plugin_path, "icons/file/distance.svg"))
        distance_button = QToolButton()
        distance_button.setIcon(distance_icon)
        distance_button.setIconSize(QSize(45, 45))
        distance_button.setToolTip(f"Distance Matrix")
        distance_button.clicked.connect(lambda: self.show_distance_matrix(distance))
        self.grid_layout.addWidget(distance_button, 2, len(self.items["distances"])-1)
    
    def add_phylogeny(self, phylogeny):
        # add a phylogeny to items["phylogenies"]
        self.items["phylogenies"].append(phylogeny)
        # add a phylogeny icon to workspace
        phylogeny_icon = QIcon(os.path.join(self.plugin_path, "icons/file/phylogeny.svg"))
        phylogeny_button = QToolButton()
        phylogeny_button.setIcon(phylogeny_icon)
        phylogeny_button.setIconSize(QSize(45, 45))
        phylogeny_button.setToolTip(f"Phylogenetic Tree")
        phylogeny_button.clicked.connect(self.open_icytree_wrapper)
        self.grid_layout.addWidget(phylogeny_button, 4, len(self.items["phylogenies"])-1)
    
    def add_dataset(self, dataset):
        """添加数据集到工作区"""
        # add a dataset to items["datasets"]
        self.items["datasets"].append(dataset)
        # add a dataset icon to workspace
        dataset_icon = self.resource_factory.get_icon("file/dataset.svg")
        dataset_button = QToolButton()
        dataset_button.setIcon(dataset_icon)
        dataset_button.setIconSize(QSize(45, 45))
        dataset_button.setToolTip(f"Dataset: {dataset.dataset_name}")
        # TODO: 连接点击事件打开Dataset Manager
        self.grid_layout.addWidget(dataset_button, 3, len(self.items["datasets"])-1)
    
    def add_model_to_workspace(self, model_data):
        """添加模型结果到工作区"""
        # 检查是单个模型还是模型表
        if "type" in model_data and model_data["type"] == "model_table":
            # 处理完整模型表
            self.items["models"].append(model_data)
            # 添加模型表图标到工作区
            model_icon = QIcon(os.path.join(self.plugin_path, "icons/file/model.svg"))
            model_button = QToolButton()
            model_button.setIcon(model_icon)
            model_button.setIconSize(QSize(45, 45))
            model_button.setToolTip(f"Model Table ({len(model_data['data'])} models)")
            model_button.clicked.connect(lambda: self.view_model_table(model_data))
            self.grid_layout.addWidget(model_button, 1, 0)
        else:
            # 处理单个模型（向后兼容）
            self.items["models"].append(model_data)
            # 添加模型图标到工作区
            model_icon = QIcon(os.path.join(self.plugin_path, "icons/file/model.svg"))
            model_button = QToolButton()
            model_button.setIcon(model_icon)
            model_button.setIconSize(QSize(45, 45))
            model_button.setToolTip(f"Model: {model_data['Model']}\n"
                                   f"LogL: {model_data['LogL']}\n"
                                   f"AIC: {model_data['AIC']}\n"
                                   f"BIC: {model_data['BIC']}")
            self.grid_layout.addWidget(model_button, 1, 0)
    
    def open_icytree_wrapper(self):
        from YR_MPE.icytree import IcyTreePlugin
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("IcyTree - YR-MPEA")
        dialog.setWindowIcon(self.resource_factory.get_icon("software/icytree.svg"))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())
        icytree_wrapper = IcyTreePlugin()
        icytree_wrapper.set_newick_string(self.items["phylogenies"][-1]['data'][0]['content'])
        dialog.layout().addWidget(icytree_wrapper)
        dialog.exec_()
    def view_model_table(self, model_data):
        """查看模型表"""
        from PyQt5.QtWidgets import QDialog, QTableWidget, QTableWidgetItem, QVBoxLayout
        from PyQt5.QtCore import Qt
        import json
        
        dialog = QDialog()
        dialog.setWindowTitle("Model Selection Results")
        dialog.setMinimumSize(800, 400)
        layout = QVBoxLayout()
        dialog.setLayout(layout)
        
        # 创建菜单栏
        menubar = QMenuBar(dialog)
        layout.setMenuBar(menubar)
        
        # 创建Export菜单
        export_menu = QMenu("&Export", dialog)
        menubar.addMenu(export_menu)
        
        # 创建表格
        table = QTableWidget()
        table.setColumnCount(8)
        table.setHorizontalHeaderLabels([
            "Model", "LogL", "AIC", "w-AIC", "AICc", "w-AICc", "BIC", "w-BIC"
        ])
        table.setRowCount(len(model_data['data']))
        
        # 填充数据
        for row, model in enumerate(model_data['data']):
            for col, (key, value) in enumerate(model.items()):
                item = QTableWidgetItem(value)
                if key not in ['Model']:
                    # 数值列右对齐
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                table.setItem(row, col, item)
        
        # 调整列宽
        table.resizeColumnsToContents()
        table.setSortingEnabled(True)
        
        # 添加导出选项
        export_csv_action = QAction("&To CSV", dialog)
        export_csv_action.triggered.connect(lambda: self.export_model_table_to_csv(table, model_data))
        export_menu.addAction(export_csv_action)
        
        export_xlsx_action = QAction("To &XLSX", dialog)
        export_xlsx_action.triggered.connect(lambda: self.export_model_table_to_xlsx(table, model_data))
        export_menu.addAction(export_xlsx_action)
        
        layout.addWidget(table)
        dialog.exec_()
    
    def view_alignment(self, sequences):
        """查看序列比对结果"""
        from .sequence_editor import SequenceAlignmentViewer
        # convert sequences to diction format: {"header", "sequence"}
        sequences = [{"header": seq.name, "sequence": str(seq.seq)} for seq in sequences]
        viewer = SequenceAlignmentViewer(sequences)
        viewer.show()
        # 保存viewer引用以防被垃圾回收
        if not hasattr(self, 'viewers'):
            self.viewers = []
        self.viewers.append(viewer)

    def show_phylogeny(self):
        """显示系统发育树"""
        if not self.items["phylogenies"]:
            QMessageBox.information(self, "Info", "No phylogeny data available.")
            return
            
        # 显示系统发育树
        from .icytree import IcyTreeViewer
        viewer = IcyTreeViewer()
        print(self.items["phylogenies"][0])
        viewer.load_tree(self.items["phylogenies"][0])

    def show_distance_matrix(self, distance_data=None):
        """显示距离矩阵"""
        if not distance_data:
            QMessageBox.information(self, "Info", "No distance matrix data available.")
            return
            
        # 创建对话框显示距离矩阵
        from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, QHeaderView
        from PyQt5.QtCore import Qt
        from PyQt5.QtGui import QFont
        
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
        
        # 添加导出选项
        export_csv_action = QAction("&To CSV", dialog)
        export_csv_action.triggered.connect(lambda: self.export_distance_matrix_to_csv(table, sequence_names))
        export_menu.addAction(export_csv_action)
        
        export_xlsx_action = QAction("To &XLSX", dialog)
        export_xlsx_action.triggered.connect(lambda: self.export_distance_matrix_to_xlsx(table, sequence_names))
        export_menu.addAction(export_xlsx_action)
        
        # 创建表格控件
        table = QTableWidget()
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        
        # 获取第一个距离矩阵数据
        content = distance_data['data'][0]['content']
        
        # 解析距离矩阵内容
        lines = content.strip().split('\n')
        if not lines:
            QMessageBox.warning(self, "Warning", "Invalid distance matrix format.")
            return
            
        # 第一行是序列数量
        try:
            num_sequences = int(lines[0].strip())
        except ValueError:
            QMessageBox.warning(self, "Warning", "Invalid distance matrix format: first line should be the number of sequences.")
            return
            
        # 后续行是距离数据
        if len(lines) < num_sequences + 1:
            QMessageBox.warning(self, "Warning", "Incomplete distance matrix data.")
            return
            
        # 获取序列名称和距离数据
        sequence_names = []
        matrix_values = []
        
        for i in range(1, num_sequences + 1):
            parts = lines[i].split()
            if len(parts) < 2:
                QMessageBox.warning(self, "Warning", f"Invalid distance matrix format at line {i+1}.")
                return
                
            sequence_names.append(parts[0])
            try:
                # 剩余部分都是距离值
                row_values = [float(val) for val in parts[1:]]
                matrix_values.append(row_values)
            except ValueError:
                QMessageBox.warning(self, "Warning", f"Invalid numeric value in distance matrix at line {i+1}.")
                return
        
        # 设置表格行列数
        table.setRowCount(num_sequences)
        table.setColumnCount(num_sequences)
        
        # 设置表头
        for i in range(num_sequences):
            table.setHorizontalHeaderItem(i, QTableWidgetItem(sequence_names[i]))
            table.setVerticalHeaderItem(i, QTableWidgetItem(sequence_names[i]))
        
        # 填充数据
        for i in range(num_sequences):
            for j in range(num_sequences):
                if i == j:
                    item = QTableWidgetItem("0.0000")
                else:
                    # 注意矩阵索引，因为第一列是序列名
                    value = matrix_values[i][j]
                    item = QTableWidgetItem(f"{value:.4f}")
                    
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                table.setItem(i, j, item)
        
        # 设置表格属性
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.setSizeAdjustPolicy(QTableWidget.AdjustToContents)
        
        layout.addWidget(table)
        dialog.exec_()
        
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