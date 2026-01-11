from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, 
                             QPushButton, QLineEdit, QTextEdit, 
                             QLabel, QComboBox, QCheckBox, QRadioButton, 
                             QSpinBox, QDoubleSpinBox, QScrollArea, 
                             QFrame, QTextEdit, QToolButton, QDialog,
                             QGroupBox, QSizePolicy, QFormLayout, QGridLayout)
from PyQt5.QtCore import Qt
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import os

class MrBayesPlugin(BasePlugin):
    def __init__(self, import_from=None, import_data=None):
        super().__init__()
        pass

    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "MrBayes-MPI-BEAGLE Phylogeny"
        self.tool_name = "MrBayes-MPI-BEAGLE"
        self.citation = [
            """Ronquist, F., Teslenko, M., van der Mark, P., Ayres, D. L., Darling, A., Höhna, S., Larget, B., Liu, L., Suchard, M. A., & Huelsenbeck, J. P. (2012). MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. Systematic biology, 61(3), 539–542. https://doi.org/10.1093/sysbio/sys029""",
            """Daniel L Ayres, Michael P Cummings, Guy Baele, Aaron E Darling, Paul O Lewis, David L Swofford, John P Huelsenbeck, Philippe Lemey, Andrew Rambaut, Marc A Suchard, BEAGLE 3: Improved Performance, Scaling, and Usability for a High-Performance Computing Library for Statistical Phylogenetics, Systematic Biology, Volume 68, Issue 6, November 2019, Pages 1052–1061, https://doi.org/10.1093/sysbio/syz020"""
        ]
        self.input_types = {"PHYLIP": ["phy"], "Chain File": [""]}
        self.output_types = {"Chain File": [".chain"], "Tree File": [".con.tre"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)

        # File input
        input_group = QGroupBox("Input")
        input_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select sequence files...")
        self.file_browse_btn = QPushButton("Browse Files")
        # self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)

        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # QCheckBox
        # If input NEXUS file contains a MrBayes Data Block, run it directly?
        self.run_data_block = QCheckBox("Run MrBayes Data Block directly if present")
        self.run_data_block.setChecked(True)
        input_layout.addRow("", self.run_data_block)

        # Phylogenetic parameters
        phy_parameters_group = QGroupBox("Phylogenetic parameters")
        phy_parameters_layout = QFormLayout()
        phy_parameters_group.setLayout(phy_parameters_layout)
        layout.addWidget(phy_parameters_group)

        # Format of MrBayes data block
        # begin mrbayes;
        # ...
        # end;

        # Rate param numbers
        self.rate_params_layout = QHBoxLayout()

        self.datatype_combo = QComboBox()
        self.datatype_combo.addItems(["DNA", "Protein"]) # default: DNA
        # self.datatype_combo.currentIndexChanged.connect(self.on_datatype_changed)
        self.datatype_combo.setCurrentIndex(0)

        self.rate_params_layout.addWidget(QLabel("Type:"))
        self.rate_params_layout.addWidget(self.datatype_combo)

        # for DNA data
        # lset nst = <dna_rate_num>, statefreq = <fixed(equal) / fixed(empirical)> // <none - estimated>;
        self.dnadata_widget = QWidget()
        self.dnadata_widget.setContentsMargins(0, 0, 0, 0)
        self.dnadata_layout = QHBoxLayout()
        self.dnadata_layout.setContentsMargins(0, 0, 0, 0)
        self.dnadata_widget.setLayout(self.dnadata_layout)

        self.dna_rate_num_spinbox = QSpinBox()
        self.dna_rate_num_spinbox.setRange(1, 6)
        self.dna_rate_num_spinbox.setValue(6)

        self.state_freq_pr_combo_dna = QComboBox()
        self.state_freq_pr_combo_dna.addItems(["estimated(dirichlet)", "fixed(equal)", "fixed(empirical)"])
        self.state_freq_pr_combo_dna.setCurrentIndex(0)

        self.dnadata_layout.addWidget(QLabel("Subst. num:"))
        self.dnadata_layout.addWidget(self.dna_rate_num_spinbox)
        self.dnadata_layout.addWidget(QLabel("State freq.:"))
        self.dnadata_layout.addWidget(self.state_freq_pr_combo_dna)

        self.rate_params_layout.addWidget(self.dnadata_widget)

        # for Protein data
        # prset aamodelpr = fixed(...) / mixed;
        self.prodata_widget = QWidget()
        self.prodata_widget.setContentsMargins(0, 0, 0, 0)
        self.prodata_widget.setVisible(False)
        self.prodata_layout = QHBoxLayout()
        self.prodata_layout.setContentsMargins(0, 0, 0, 0)
        self.prodata_widget.setLayout(self.prodata_layout)
        
        # self.prodata_layout.setVisible(False)
        self.prot_model_combo = QComboBox()
        self.prot_model_combo.addItems(["Blosum62", "Blosum", "Wag", "Lg", "gtr", "jones", "mtrev", "Poisson", "mixed"])

        self.prodata_layout.addWidget(QLabel("Protein model"))
        self.prodata_layout.addWidget(self.prot_model_combo)
        self.prodata_layout.addWidget(QLabel("State freq.:"))
        

        self.state_freq_pr_combo_prot = QComboBox()
        self.state_freq_pr_combo_prot.addItems(["estimated(dirichlet)", "fixed(equal)", "fixed(empirical)"])
        self.state_freq_pr_combo_prot.setCurrentIndex(0)

        self.prodata_layout.addWidget(self.state_freq_pr_combo_prot)

        self.rate_params_layout.addWidget(self.prodata_widget)

        phy_parameters_layout.addRow("Model Settings:", self.rate_params_layout)
        
        # Rate Heterogenity
        # lset rates = euqal/gamma/invgamma/propinv/lnorm/adgamma Ngammacat=<gamma categories>;
        rate_hetero_layout = QHBoxLayout()
        self.rate_hetero_combo = QComboBox()
        self.rate_hetero_combo.addItems(["Equal", "Gamma (+G)", "InvGamma (+G+I)", "PropInv (+I)", "Lognormal", "Adgamma"])

        self.gamma_categories_spinbox = QSpinBox()
        self.gamma_categories_spinbox.setRange(1, 100)
        self.gamma_categories_spinbox.setValue(4)

        # rate_hetero_layout.addWidget(QLabel("Rate Var.:"))
        rate_hetero_layout.addWidget(self.rate_hetero_combo)
        rate_hetero_layout.addWidget(QLabel("Gamma Categories:"))
        rate_hetero_layout.addWidget(self.gamma_categories_spinbox)

        phy_parameters_layout.addRow("Rate Heterogenity:", rate_hetero_layout)


        
        # MPI & BEAGLE settings
        mpi_beagle_group = QGroupBox("MPI && BEAGLE settings")
        mpi_beagle_layout = QFormLayout()
        mpi_beagle_group.setLayout(mpi_beagle_layout)
        layout.addWidget(mpi_beagle_group)

        # MPI parallels
        mpi_params_layout = QHBoxLayout()
        self.use_mpi = QCheckBox("Use MPI")
        self.use_mpi.setChecked(True)
        # Number of Parallels
        self.n_mpi_parallel = QSpinBox()
        n_cpu_cores = os.cpu_count()
        # Set the maximum value to the number of CPU cores
        self.n_mpi_parallel.setMaximum(n_cpu_cores)
        # default value: max(n_cpu_cores / 4, 1)
        self.n_mpi_parallel.setValue(max(int(n_cpu_cores / 4), 1))

        mpi_params_layout.addWidget(self.use_mpi)
        mpi_params_layout.addWidget(QLabel("Number of Parallels:"))
        mpi_params_layout.addWidget(self.n_mpi_parallel)
        mpi_beagle_layout.addRow("MPI settings:", mpi_params_layout) # <mpirun [UNIX] / mpiexec [WINDOWS]> -np <number>

        self.beagle_params_layout = QGridLayout()
        # use beagle?

        # set usebeagle=yes/no beagledevice=cpu/gpu beagleprecision=double/single beaglescaling=dynamic/always;
        self.use_beagle = QCheckBox("Use BEAGLE") 
        self.use_beagle.setChecked(True)

        self.beagle_device_label = QLabel("Device:")
        self.beagle_device_label.setAlignment(Qt.AlignRight)
        self.beagle_device_combo = QComboBox()
        self.beagle_device_combo.addItems(["CPU", "GPU"])
        self.beagle_device_combo.setCurrentText("CPU")

        self.beagle_precision_label = QLabel("Precision:")
        self.beagle_precision_label.setAlignment(Qt.AlignRight)
        self.beagle_precision_combo = QComboBox()
        self.beagle_precision_combo.addItems(["Double", "Single"])
        self.beagle_precision_combo.setCurrentText("Double")

        self.beagle_scaling_label = QLabel("Scaling:")
        self.beagle_scaling_label.setAlignment(Qt.AlignRight)
        self.beagle_scaling_combo = QComboBox()
        self.beagle_scaling_combo.addItems(["Dynamic", "Always"])
        self.beagle_scaling_combo.setCurrentText("Dynamic")

        self.beagle_params_layout.addWidget(self.use_beagle, 0, 0)
        self.beagle_params_layout.addWidget(self.beagle_device_label, 0, 2)
        self.beagle_params_layout.addWidget(self.beagle_device_combo, 0, 3)
        self.beagle_params_layout.addWidget(self.beagle_precision_label, 1, 0)
        self.beagle_params_layout.addWidget(self.beagle_precision_combo, 1, 1)
        self.beagle_params_layout.addWidget(self.beagle_scaling_label, 1, 2)
        self.beagle_params_layout.addWidget(self.beagle_scaling_combo, 1, 3)

        mpi_beagle_layout.addRow("BEAGLE settings", self.beagle_params_layout)

        mcmc_params_group = QGroupBox("MCMC settings")
        mcmc_params_layout = QFormLayout()
        mcmc_params_group.setLayout(mcmc_params_layout)

        layout.addWidget(mcmc_params_group)

        # generation; sampling frequency; run num; chain num;
        # mcmcp ngen=* samplefreq=* printfreq=<samplefreq> nchains=* nruns=* savebrlens=yes checkpoint=yes checkfreq=5000;
        self.generation_spinbox = QSpinBox()
        self.generation_spinbox.setMinimum(1)
        self.generation_spinbox.setValue(1000000)
        self.generation_spinbox.setMaximum(100000000000000000)

        self.sampling_frequency_spinbox = QSpinBox()
        self.sampling_frequency_spinbox.setMinimum(1)
        self.sampling_frequency_spinbox.setValue(1000)
        self.sampling_frequency_spinbox.setMaximum(100000000000000000)

        self.run_num_spinbox = QSpinBox()
        self.run_num_spinbox.setMinimum(1)
        self.run_num_spinbox.setValue(2)

        self.chain_num_spinbox = QSpinBox()
        self.chain_num_spinbox.setMinimum(1)
        self.chain_num_spinbox.setValue(4) # MC^3

        mcmcp_widget = QWidget()
        mcmcp_layout = QHBoxLayout()
        mcmcp_widget.setContentsMargins(0, 0, 0, 0)
        mcmcp_layout.setContentsMargins(0, 0, 0, 0)
        mcmcp_widget.setLayout(mcmcp_layout)

        mcmcp_layout.addWidget(QLabel("Generations:"))
        mcmcp_layout.addWidget(self.generation_spinbox)
        mcmcp_layout.addWidget(QLabel("Sampling Freq.:"))
        mcmcp_layout.addWidget(self.sampling_frequency_spinbox)
        mcmcp_layout.addWidget(QLabel("Runs:"))
        mcmcp_layout.addWidget(self.run_num_spinbox)
        mcmcp_layout.addWidget(QLabel("Chains:"))
        mcmcp_layout.addWidget(self.chain_num_spinbox)

        mcmc_params_layout.addRow("MCMC:", mcmcp_widget)

        # summary and consensus
        # sump <relburnin=yes burninfrac=*> / <burnin=*>;
        # sumt <relburnin=yes burninfrac=*> / <burnin=*> conformat=<Simple/FigTree> contype=*

        sum_widget = QWidget()
        sum_layout = QVBoxLayout()
        sum_widget.setContentsMargins(0, 0, 0, 0)
        sum_layout.setContentsMargins(0, 0, 0, 0)
        sum_widget.setLayout(sum_layout)

        self.contree_type = QComboBox()
        self.contree_type.addItems(["Majority-rule", "Halfcompat", "Allcompat"])

        self.contree_format = QComboBox()
        self.contree_format.addItems(["FigTree", "Simple"])

        sum_layout.addWidget(QLabel("Method:"))
        sum_layout.addWidget(self.contree_type)
        sum_layout.addWidget(QLabel("Output format:"))
        sum_layout.addWidget(self.contree_format)

        mcmc_params_layout.addRow("Consensus:", sum_layout)

        burnin_widget = QWidget()
        burnin_layout = QHBoxLayout()
        burnin_widget.setContentsMargins(0, 0, 0, 0)
        burnin_layout.setContentsMargins(0, 0, 0, 0)
        burnin_widget.setLayout(burnin_layout)

        self.burnin_as_fraction = QRadioButton("Fraction:")
        self.burnin_as_fraction.setChecked(True)

        self.burnin_fraction = QDoubleSpinBox()
        self.burnin_fraction.setRange(0, 1)
        self.burnin_fraction.setSingleStep(0.01)
        self.burnin_fraction.setValue(0.25)

        self.burnin_as_states = QRadioButton("States:")
        self.burnin_as_states.setChecked(False)

        self.burnin_states = QSpinBox()
        self.burnin_states.setRange(0, 1000000000000000)
        self.burnin_states.setValue(1000)

        burnin_layout.addWidget(self.burnin_as_fraction)
        burnin_layout.addWidget(self.burnin_fraction)
        burnin_layout.addWidget(self.burnin_as_states)
        burnin_layout.addWidget(self.burnin_states)

        mcmc_params_layout.addRow("Burn-in:", burnin_layout)



        mb_data_block_widget = QWidget()
        mb_data_block_widget.setContentsMargins(0,0,0,0)
        mb_data_block_layout = QHBoxLayout()
        mb_data_block_layout.setContentsMargins(0,0,0,0)
        mb_data_block_widget.setLayout(mb_data_block_layout)

        show_mb_data_block = QPushButton("Show MrBayes data block") # use QDialog to showcase MrBayes data block (Font: Consolas; same style as Console)
        copy_mb_data_block = QPushButton("Copy MrBayes data block")

        mb_data_block_layout.addWidget(show_mb_data_block)
        mb_data_block_layout.addWidget(copy_mb_data_block)

        layout.addWidget(mb_data_block_widget)

class MrBayesPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.plugin_path = plugin_path
        # self.config = config_loader()
    
    def run(self, import_from=None, import_data=None):
        return MrBayesPlugin(import_from, import_data)