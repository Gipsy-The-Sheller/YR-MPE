from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QLabel, QLineEdit,
                             QListWidget, QListWidgetItem, QPushButton, QComboBox,
                             QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
                             QSizePolicy)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_treeview import MplTreeView
from mpl_distribution import MplDistribution
import sys

class MDMRCA(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()
    
    def init_ui(self, example=False):
        self.setWindowTitle("MD-MRCA")
        # self.setGeometry(100, 100, 800, 600)

        self.main_layout = QHBoxLayout()
        self.setLayout(self.main_layout)

        self.left_layout = QVBoxLayout()
        self.left_layout.setContentsMargins(0,0,0,0)

        self.main_layout.addLayout(self.left_layout)

        self.init_taxon_set()
        self.init_tmrca_set()
        self.init_tree_view()

        if example:
            pass
    def init_taxon_set(self):
        taxon_set_layout = QHBoxLayout()
        self.left_layout.addLayout(taxon_set_layout)
        
        # part 1. OTU explorer
        otu_explor_layout = QVBoxLayout()

        taxon_set_layout.addLayout(otu_explor_layout)

        self.otu_filter = QLineEdit()
        self.otu_filter.setPlaceholderText("Filter taxa...")
        self.otu_filter.textChanged.connect(self.filter_taxa)

        self.taxon_list = QListWidget()
        self.taxon_list.setSelectionMode(QListWidget.ExtendedSelection)

        otu_explor_layout.addWidget(self.otu_filter)
        otu_explor_layout.addWidget(self.taxon_list)

        # part 2. Pushbuttons
        button_layout = QVBoxLayout()

        taxon_set_layout.addLayout(button_layout)

        self.add_taxon_button = QPushButton("Add taxa >")
        self.add_taxon_button.clicked.connect(self.add_taxa)
        button_layout.addWidget(self.add_taxon_button)

        self.remove_taxon_button = QPushButton("< Remove taxa")
        self.remove_taxon_button.clicked.connect(self.remove_taxa)
        button_layout.addWidget(self.remove_taxon_button)

        self.clear_taxon_button = QPushButton("Clear taxa")
        self.clear_taxon_button.clicked.connect(self.clear_taxa)
        button_layout.addWidget(self.clear_taxon_button)

        # part 3. selected taxa
        selected_taxa_layout = QVBoxLayout()

        taxon_set_layout.addLayout(selected_taxa_layout)

        self.taxon_set_name = QLineEdit()
        self.taxon_set_name.setPlaceholderText("Taxon set name...")
        selected_taxa_layout.addWidget(self.taxon_set_name)

        self.selected_taxa_list = QListWidget()
        selected_taxa_layout.addWidget(self.selected_taxa_list)
        
    def init_tmrca_set(self):
        tmrca_set_layout = QHBoxLayout()

        self.left_layout.addLayout(tmrca_set_layout)

        # 1. left side: tMRCA type & tMRCA params
        tmrca_left_layout = QVBoxLayout()
        tmrca_set_layout.addLayout(tmrca_left_layout)

        tmrca_type_layout = QHBoxLayout()

        tmrca_left_layout.addLayout(tmrca_type_layout)

        tmrca_type_label = QLabel("Calibration type:")
        tmrca_type_label.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Minimum)
        tmrca_type_layout.addWidget(tmrca_type_label)

        tmrca_type_combo = QComboBox()
        tmrca_type_combo.addItems(['Point', 'Uniform', 'Upper Boundary', 'Lower Boundary', 'Normal', 'Lognormal'])
        tmrca_type_layout.addWidget(tmrca_type_combo)

        tmrca_setting_group = QGroupBox("tMRCA parameters")
        self.tmrca_setting_group_layout = QFormLayout()
        tmrca_setting_group.setLayout(self.tmrca_setting_group_layout)

        tmrca_left_layout.addWidget(tmrca_setting_group)

        # matplotlib figure canvas to display tMRCA distribution
        self.tmrca_figure_canvas = MplDistribution()
        tmrca_set_layout.addWidget(self.tmrca_figure_canvas)

    def init_tree_view(self):
        # 使用封装的MplTreeView组件来绘制树状图
        self.tree_figure_canvas = MplTreeView()
        self.main_layout.addWidget(self.tree_figure_canvas)

        # 示例：绘制一个默认的树
        # newick_str = "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0,E:1.0):1.0);"
        # self.tree_figure_canvas.draw_dendrogram(newick_str=newick_str)
        
    def init_tmrca_parameters(self):
        """初始化tMRCA参数输入字段"""
        self.param_widgets = {}
        
        # Point distribution parameters
        self.point_value = QLineEdit()
        self.point_value.setText("1.0")
        self.param_widgets['Point'] = [('Point value:', self.point_value)]
        
        # Uniform distribution parameters
        self.uniform_lower = QLineEdit()
        self.uniform_lower.setText("0.5")
        self.uniform_upper = QLineEdit()
        self.uniform_upper.setText("1.5")
        self.param_widgets['Uniform'] = [('Lower bound:', self.uniform_lower), 
                                        ('Upper bound:', self.uniform_upper)]
        
        # Upper Boundary parameters
        self.upper_boundary = QLineEdit()
        self.upper_boundary.setText("2.0")
        self.param_widgets['Upper Boundary'] = [('Upper boundary:', self.upper_boundary)]
        
        # Lower Boundary parameters
        self.lower_boundary = QLineEdit()
        self.lower_boundary.setText("0.5")
        self.param_widgets['Lower Boundary'] = [('Lower boundary:', self.lower_boundary)]
        
        # Normal distribution parameters
        self.normal_mean = QLineEdit()
        self.normal_mean.setText("1.0")
        self.normal_std = QLineEdit()
        self.normal_std.setText("0.2")
        self.param_widgets['Normal'] = [('Mean:', self.normal_mean), 
                                       ('Standard deviation:', self.normal_std)]
        
        # Lognormal distribution parameters
        self.lognormal_mean = QLineEdit()
        self.lognormal_mean.setText("1.0")
        self.lognormal_std = QLineEdit()
        self.lognormal_std.setText("0.2")
        self.param_widgets['Lognormal'] = [('Mean:', self.lognormal_mean), 
                                          ('Standard deviation:', self.lognormal_std)]
        
        # 设置初始参数显示
        self.update_tmrca_parameters('Point')
        
    def update_tmrca_parameters(self, distribution_type):
        """根据分布类型更新参数输入字段"""
        # 清除现有参数
        while self.tmrca_setting_group_layout.count():
            child = self.tmrca_setting_group_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        
        # 添加新参数
        if distribution_type in self.param_widgets:
            for label_text, widget in self.param_widgets[distribution_type]:
                label = QLabel(label_text)
                self.tmrca_setting_group_layout.addRow(label, widget)
                # 连接参数变化信号
                if hasattr(widget, 'textChanged'):
                    widget.textChanged.connect(self.on_tmrca_params_changed)
        
        # 更新分布图
        self.update_tmrca_distribution()
        
    def on_tmrca_type_changed(self, distribution_type):
        """当tMRCA类型改变时的回调"""
        self.update_tmrca_parameters(distribution_type)
        
    def on_tmrca_params_changed(self):
        """当tMRCA参数改变时的回调"""
        self.update_tmrca_distribution()
        
    def update_tmrca_distribution(self):
        """更新tMRCA分布图"""
        distribution_type = self.tmrca_type_combo.currentText()
        
        try:
            if distribution_type == 'Point':
                point_val = float(self.point_value.text())
                self.tmrca_figure_canvas.plot_point_distribution(point_val)
                
            elif distribution_type == 'Uniform':
                lower = float(self.uniform_lower.text())
                upper = float(self.uniform_upper.text())
                self.tmrca_figure_canvas.plot_uniform_distribution(lower, upper)
                
            elif distribution_type == 'Upper Boundary':
                upper = float(self.upper_boundary.text())
                self.tmrca_figure_canvas.plot_upper_boundary_distribution(upper)
                
            elif distribution_type == 'Lower Boundary':
                lower = float(self.lower_boundary.text())
                self.tmrca_figure_canvas.plot_lower_boundary_distribution(lower)
                
            elif distribution_type == 'Normal':
                mean = float(self.normal_mean.text())
                std = float(self.normal_std.text())
                self.tmrca_figure_canvas.plot_normal_distribution(mean, std)
                
            elif distribution_type == 'Lognormal':
                mean = float(self.lognormal_mean.text())
                std = float(self.lognormal_std.text())
                self.tmrca_figure_canvas.plot_lognormal_distribution(mean, std)
                
        except ValueError:
            # 参数无效时清空图表
            self.tmrca_figure_canvas.clear_plot()
    
    def filter_taxa(self):
        pass
    def add_taxa(self):
        pass    
    def remove_taxa(self):
        pass
    def clear_taxa(self):
        pass



def main():
    app = QApplication(sys.argv)  # 使用sys.argv
    
    # 创建主窗口
    main_window = QMainWindow()
    main_window.setWindowTitle("MD-MRCA Application")
    
    # 创建MDMRCA部件并设置为中央部件
    md_mrca_widget = MDMRCA()
    main_window.setCentralWidget(md_mrca_widget)
    
    # 显示窗口
    main_window.show()
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()