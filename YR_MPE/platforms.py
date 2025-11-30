from PyQt5.QtWidgets import (QWidget, QVBoxLayout, 
QMenuBar, QToolBar, QToolButton, QGroupBox, QLabel,
QAction, QMenu, QSizePolicy, QGridLayout, QFileDialog, QMessageBox)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon
import os
from Bio import SeqIO
# from core.plugin_base import BasePlugin

class YR_MPEA_Widget(QWidget):
    def __init__(self):
        super().__init__()

        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.mode = "single_gene"
        self.init_ui()
        

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

        align_button = QToolButton()
        align_button.setText("ALIGN")
        # align_button.setStyleSheet("color: #555555;")
        align_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/align.svg")))
        align_button.setPopupMode(QToolButton.InstantPopup)
        align_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        # align_button.setFixedSize(60, 60)
        main_toolbar.addWidget(align_button)


        open_action = QAction("Open Sequence Files", align_button)
        open_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/open.svg")))
        open_action.triggered.connect(self.open_sequence_files)
        build_action = QAction("New Alignment", align_button)
        build_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/new.svg")))
        align_action = QAction("Align by...", align_button)
        align_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/alignby.svg")))
        trim_action = QAction("Trim Alignment by...", align_button)
        trim_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/trim.svg")))
        
        align_button.addAction(open_action)
        align_button.addAction(build_action)
        align_button.addAction(align_action)
        align_button.addAction(trim_action)

        # aligners submenu for align_action
        aligners_menu = QMenu()
        # aligners_menu.addAction(QAction("Clustal Omega", align_button))
        # aligners_menu.addAction(QAction("MAFFT", align_button))
        clustal_omega_action = QAction("Clustal Omega", align_button)
        clustal_omega_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/clustalo.svg")))
        clustal_omega_action.triggered.connect(self.open_clustal_omega_wrapper)
        aligners_menu.addAction(clustal_omega_action)
        mafft_action = QAction("MAFFT", align_button)
        mafft_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/mafft.svg")))
        mafft_action.triggered.connect(self.open_mafft_wrapper)
        aligners_menu.addAction(mafft_action)
        muscle5_action = QAction("Muscle 5", align_button)
        muscle5_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/muscle.svg")))
        muscle5_action.triggered.connect(self.open_muscle5_wrapper)
        aligners_menu.addAction(muscle5_action)
        aligners_menu.addAction(QAction("MACSE", align_button))
        align_action.setMenu(aligners_menu)

        trim_menu = QMenu()
        trimal_action = QAction("Trimal", trim_action)
        trimal_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/trimal.svg")))
        trimal_action.triggered.connect(self.open_trimal_wrapper)
        trim_menu.addAction(trimal_action)
        gblocks_action = QAction("GBlocks", trim_action)
        gblocks_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/gblocks.svg")))
        gblocks_action.triggered.connect(self.open_gblocks_wrapper)
        trim_menu.addAction(gblocks_action)
        trim_action.setMenu(trim_menu)

        model_button = QToolButton()
        model_button.setText("MODEL")
        model_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/model.svg")))
        model_button.setPopupMode(QToolButton.InstantPopup)
        model_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        # model_button.setFixedSize(60, 60)
        main_toolbar.addWidget(model_button)

        models_menu = QMenu()
        find_best_model_action = QAction("Find Best Substitution Model (ML)", model_button)
        find_best_model_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/find_model.svg")))
        find_best_model_action.triggered.connect(self.open_modelfinder_wrapper)
        models_menu.addAction(find_best_model_action)
        model_button.setMenu(models_menu)

        # TODO: Model Selection Algorithms: ModelFinder / ModelTest-NG

        distance_button = QToolButton()
        distance_button.setText("DISTANCE")
        distance_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/distance.svg")))
        distance_button.setPopupMode(QToolButton.InstantPopup)
        distance_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        # distance_button.setFixedSize(60, 60)
        main_toolbar.addWidget(distance_button)

        # 添加"Compute Pairwise Distances"动作
        self.comp_dist_action = QAction("Compute Pairwise Distances (ML)", self)
        self.comp_dist_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/dist.svg")))
        self.comp_dist_action.triggered.connect(self.open_ml_distance_wrapper)
        distance_button.addAction(self.comp_dist_action)
        
        # 添加"Compute Overall Mean Distances"动作
        self.comp_mdist_action = QAction("Compute Overall Mean Distances (ML)", self)
        self.comp_mdist_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/mdist.svg")))
        # TODO: 实现平均距离计算功能
        # self.comp_mdist_action.triggered.connect(self.compute_mean_distances)
        distance_button.addAction(self.comp_mdist_action)

        phylogeny_button = QToolButton()
        phylogeny_button.setText("PHYLOGENY")
        phylogeny_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/phylogeny.svg")))
        phylogeny_button.setPopupMode(QToolButton.InstantPopup)
        phylogeny_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        phylogeny_button_menu = QMenu()
        phylogeny_button.setMenu(phylogeny_button_menu)
        # phylogeny_button.setFixedSize(60, 60)
        main_toolbar.addWidget(phylogeny_button)

        cons_ml_action = QAction("Phenetics - Maximum Likelihood Phylogenies (ML)", phylogeny_button)
        cons_ml_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/ml.svg")))

        ml_menu = QMenu()
        cons_ml_action.setMenu(ml_menu)
        
        # Add IQ-TREE action
        iqtree_action = QAction("IQ-Tree 3", phylogeny_button)
        iqtree_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/iqtree.svg")))
        iqtree_action.triggered.connect(self.open_iqtree_wrapper)

        ml_menu.addAction(iqtree_action)
        
        # TODO: ML Programs: IQ-TREE 3 / FastTree

        cons_bi_action = QAction("Phenetics - Bayesian Inference Phylogenies (BI)", phylogeny_button)
        cons_bi_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/bi.svg")))

        # TODO: BI Programs: MrBayes

        cons_mp_action = QAction("Cladistics - Maximum Parsimony Phylogenies (MP)", phylogeny_button)
        cons_mp_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/mp.svg")))

        # TODO: MP Programs: TNT

        cons_nj_action = QAction("Clustering - Neighbor-Joining (NJ)", phylogeny_button)
        cons_nj_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/nj.svg")))

        cons_upgma_action = QAction("Clustering - UPGMA", phylogeny_button)
        cons_upgma_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/upgma.svg")))

        cons_me_action = QAction("Clustering - Minimum Evolution (ME)", phylogeny_button)
        cons_me_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/me.svg")))

        cons_bionj_action = QAction("Clustering - BioNJ [IQ-TREE 3]", phylogeny_button)
        cons_bionj_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/bionj.svg")))

        cons_rand_action = QAction("Simulation - Random Trees (Evolver)", phylogeny_button)
        cons_rand_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/randtree.svg")))

        cons_seq_action = QAction("Simulation - Simulated Sequences (Evolver)", phylogeny_button)
        cons_seq_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/randseq.svg")))

        tree_viewer_action = QAction("Tree Viewer (IcyTree)", phylogeny_button)
        tree_viewer_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/icytree.svg")))
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
        variants_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/variants.svg")))
        variants_button.setPopupMode(QToolButton.InstantPopup)
        variants_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        variants_button_menu = QMenu()
        variants_button.setMenu(variants_button_menu)
        # phylogeny_button.setFixedSize(60, 60)
        main_toolbar.addWidget(variants_button)

        coalescent_button = QToolButton()
        coalescent_button.setText("COALESCENT")
        coalescent_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/coalescent.svg")))
        coalescent_button.setPopupMode(QToolButton.InstantPopup)
        coalescent_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        coalescent_button_menu = QMenu()
        coalescent_button.setMenu(coalescent_button_menu)
        main_toolbar.addWidget(coalescent_button)

        # TODO: Astral-III and CASTER

        clock_button = QToolButton()
        clock_button.setText("CLOCK")
        clock_button.setIcon(QIcon(os.path.join(self.plugin_path, "icons/clock.svg")))
        clock_button.setPopupMode(QToolButton.InstantPopup)
        clock_button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        clock_button_menu = QMenu()
        clock_button.setMenu(clock_button_menu)
        main_toolbar.addWidget(clock_button)

        # TODO: Molecular Clock: LSD2 (IQ-TREE 3) / Bayesian Inference (MrBayes / PAML-mcmctree) / Penalized Likelihood (r8s)


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

        if self.mode == "single_gene":
            self.workspace = SingleGeneWorkspace()
        # elif self.mode == "multiple_genes":
        #     self.workspace = MultipleGenesWorkspace() # TODO: multiple genes workspace
        self.workspace.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(self.workspace)

    def open_muscle5_wrapper(self):
        # create a QDialog to open the muscle5_wrapper
        from PyQt5.QtWidgets import QDialog
        from .plugins.muscle5_plugin import Muscle5PluginEntry
        dialog = QDialog()
        dialog.setWindowTitle("MUSCLE5")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/muscle.svg")))
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
        muscle5_wrapper = Muscle5PluginEntry().run(import_from=import_from, import_data=import_data)
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

    def open_clustal_omega_wrapper(self):
        from .plugins.clustal_omega_plugin import ClustalOmegaPluginEntry
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("Clustal Omega - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/clustalo.svg")))
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
        
        # use QDialog to open the clustalo_wrapper
        clustalo_wrapper = ClustalOmegaPluginEntry().run(import_from=import_from, import_data=import_data)
        clustalo_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(clustalo_wrapper)
        dialog.exec_()
    
    def open_mafft_wrapper(self):
        from YR_MPE.plugins.mafft_plugin import MAFFTPluginEntry
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("MAFFT - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/mafft.svg")))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())

        import_from = None
        import_data = None
        workspace_type = type(self.workspace).__name__
        if workspace_type == "SingleGeneWorkspace":
            if len(self.workspace.items["alignments"]) >= 1:
                import_from = "YR_MPEA"
                import_data = self.workspace.items["alignments"][0]

        mafft_wrapper = MAFFTPluginEntry().run(import_from=import_from, import_data=import_data)
        mafft_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(mafft_wrapper)
        dialog.exec_()
    
    def open_modelfinder_wrapper(self):
        from .plugins.model_finder_plugin import ModelFinderPluginEntry
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("ModelFinder - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/find_model.svg")))
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
        
        # use QDialog to open the modelfinder_wrapper
        modelfinder_wrapper = ModelFinderPluginEntry().run(import_from=import_from, import_data=import_data)
        modelfinder_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        modelfinder_wrapper.export_model_result_signal.connect(self.add_model_to_workspace)
        dialog.layout().addWidget(modelfinder_wrapper)
        dialog.exec_()
    
    def open_icytree_wrapper(self):
        from YR_MPE.icytree import IcyTreePlugin
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("IcyTree - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/icytree.svg")))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())
        icytree_wrapper = IcyTreePlugin()
        dialog.layout().addWidget(icytree_wrapper)
        dialog.exec_()

    def open_trimal_wrapper(self):
        from .plugins.trimal_plugin import TrimAlPluginEntry
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("TrimAl - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/trimal.svg")))
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
        
        # use QDialog to open the trimal_wrapper
        trimal_wrapper = TrimAlPluginEntry().run(import_from=import_from, import_data=import_data)
        trimal_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        dialog.layout().addWidget(trimal_wrapper)
        dialog.exec_()
        
    def open_gblocks_wrapper(self):
        from .plugins.gblocks_plugin import GBlocksPluginEntry
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("GBlocks - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/gblocks.svg")))
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
        
        # use QDialog to open the gblocks_wrapper
        gblocks_wrapper = GBlocksPluginEntry().run(import_from=import_from, import_data=import_data)
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
        from .plugins.ml_distance_plugin import MLDistancePluginEntry
        dialog = QDialog()
        dialog.setWindowTitle(f"Distance Calculator [implemented from IQ-TREE] - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, f"icons/distance.svg")))
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
                    
                    # # 为导入数据添加模型信息
                    # if isinstance(import_data, list) and len(import_data) > 0:
                    #     # 为第一个序列添加matrix_header属性
                    #     first_seq = import_data[0]
                    #     if hasattr(first_seq, '__dict__') and best_model:
                    #         # 构造matrix_header信息
                    #         header_lines = ["# Sequence identity and model information:"]
                    #         if "Model" in best_model:
                    #             header_lines.append(f"Model of evolution: {best_model['Model']}")
                    #         # Gamma信息可能在顶层或在模型数据内
                    #         gamma_value = model_data.get("Gamma", best_model.get("Gamma"))
                    #         if gamma_value:
                    #             header_lines.append(f"Gamma categories: {gamma_value}")
                    #         first_seq.matrix_header = "\n".join(header_lines)

        
        # use QDialog to open the plugin_wrapper
        plugin_wrapper = MLDistancePluginEntry().run(import_from=import_from, import_data=import_data)
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_distance_result_signal.connect(self.add_distance_matrix_to_workspace)

        # parse model
        model_entries = best_model.split("+")
        plugin_wrapper.model_combo.setCurrentText(model_entries[0])

        # Invariable sites?
        if "I" in model_entries:
            plugin_wrapper.invar_checkbox.setChecked(True)
        
        # empirical?
        if "F" in model_entries:
            plugin_wrapper.empirical_checkbox.setChecked(True)
        
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
        
    def open_iqtree_wrapper(self):
        from PyQt5.QtWidgets import QDialog
        from .plugins.iqtree_plugin import IQTreePluginEntry
        dialog = QDialog()
        dialog.setWindowTitle(f"Phylogenetic Inference [implemented from IQ-TREE] - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, f"icons/software/iqtree.svg")))
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
        
        # use QDialog to open the plugin_wrapper
        plugin_wrapper = IQTreePluginEntry().run(import_from=import_from, import_data=import_data)
        plugin_wrapper.import_alignment_signal.connect(self.add_alignment_to_workspace)
        plugin_wrapper.export_model_result_signal.connect(self.add_model_to_workspace)
        plugin_wrapper.export_phylogeny_result_signal.connect(self.add_phylogeny_to_workspace)

        # parse model
        model_entries = best_model.split("+")
        plugin_wrapper.model_combo.setCurrentText(model_entries[0])

        # Invariable sites?
        if "I" in model_entries:
            plugin_wrapper.invar_checkbox.setChecked(True)
        
        # empirical?
        if "F" in model_entries:
            plugin_wrapper.empirical_checkbox.setChecked(True)
        
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

class SingleGeneWorkspace(QWidget):
    def __init__(self):
        super().__init__()
        self.items = {
            "alignments": [],
            "models": [],
            "distances": [],
            "phylogenies": [],
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
        alignment_icon = QIcon(os.path.join(self.plugin_path, "icons/file/sequence.svg"))
        alignment_button = QToolButton()
        alignment_button.setIcon(alignment_icon)
        alignment_button.setIconSize(QSize(45, 45))
        # alignment_button.setToolTip("Alignment")
        alignment_button.clicked.connect(lambda: self.view_alignment(sequences))
        self.grid_layout.addWidget(alignment_button, 0, 0)
    
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
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/icytree.svg")))
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
        
        dialog = QDialog()
        dialog.setWindowTitle("Model Selection Results")
        dialog.setMinimumSize(800, 400)
        layout = QVBoxLayout()
        dialog.setLayout(layout)
        
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
        
        dialog = QDialog(self)
        dialog.setWindowTitle("Distance Matrix Visualization")
        dialog.resize(800, 600)
        layout = QVBoxLayout(dialog)
        
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
