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
        # trimal_action.triggered.connect(self.open_trimal_wrapper)
        trim_menu.addAction(trimal_action)
        gblocks_action = QAction("GBlocks", trim_action)
        gblocks_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/software/gblocks.svg")))
        # gblocks_action.triggered.connect(self.open_gblocks_wrapper)
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

        comp_dist_action = QAction("Compute Pairwise Distances", distance_button)
        comp_dist_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/dist.svg")))
        comp_mdist_action = QAction("Compute Overall Mean Distances", distance_button)
        comp_mdist_action.setIcon(QIcon(os.path.join(self.plugin_path, "icons/mdist.svg")))

        distance_button.addAction(comp_dist_action)
        distance_button.addAction(comp_mdist_action)

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
        from YR_MPE.aligners import MAFFT_wrapper
        from PyQt5.QtWidgets import QDialog
        dialog = QDialog()
        dialog.setWindowTitle("MAFFT - YR-MPEA")
        dialog.setWindowIcon(QIcon(os.path.join(self.plugin_path, "icons/software/mafft.svg")))
        dialog.setMinimumSize(800, 600)
        dialog.setLayout(QVBoxLayout())
        mafft_wrapper = MAFFT_wrapper()
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
        self.grid_layout.addWidget(distance_button, 0, 0)
    
    def add_phylogeny(self, phylogeny):
        # add a phylogeny to items["phylogenies"]
        self.items["phylogenies"].append(phylogeny)
        # add a phylogeny icon to workspace
        phylogeny_icon = QIcon(os.path.join(self.plugin_path, "icons/file/phylogeny.svg"))
        phylogeny_button = QToolButton()
        phylogeny_button.setIcon(phylogeny_icon)
        phylogeny_button.setIconSize(QSize(45, 45))
        self.grid_layout.addWidget(phylogeny_button, 0, 0)

        # if clicked, open the tree by icytree
        phylogeny_button.clicked.connect(self.open_icytree_wrapper)
    
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
        icytree_wrapper.set_newick_string(self.items["phylogenies"][-1])
        dialog.layout().addWidget(icytree_wrapper)
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
