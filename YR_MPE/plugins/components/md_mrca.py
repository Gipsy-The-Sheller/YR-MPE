from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QLabel, QLineEdit,
                             QListWidget, QListWidgetItem, QPushButton, QComboBox,
                             QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
                             QSizePolicy)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_treeview import MplTreeView
from mpl_distribution import MplDistribution
import sys

class MDMRCA(QWidget):
    def __init__(self, example=False):
        super().__init__()
        self.param_widgets = {}  # 初始化param_widgets属性
        self.init_ui()

        if example:
            self.load_example_data()
    
    def init_ui(self):
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
    
    def load_example_data(self):
        """加载示例数据：一棵系统发育树和对应的分类单元"""
        # 示例Newick格式的树
        example_newick = "((Human:0.1,Chimpanzee:0.1):0.2,(Gorilla:0.3,Orangutan:0.4):0.1,(Mouse:0.8,Rat:0.7):0.2);"
        
        # 从树中提取所有叶子节点（分类单元）
        example_taxa = ["Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse", "Rat"]
        
        # 存储所有可用的分类单元
        self.all_taxa = set(example_taxa)
        
        # 填充OTU浏览器中的分类单元列表（初始时显示所有）
        self.taxon_list.clear()
        for taxon in example_taxa:
            item = QListWidgetItem(taxon)
            self.taxon_list.addItem(item)
        
        # 绘制示例树
        self.tree_figure_canvas.draw_dendrogram(newick_str=example_newick)
        
        # 可选：预选择一些分类单元作为示例
        # 例如选择前两个作为示例分类单元集
        self.taxon_set_name.setText("Primates")
        self.selected_taxa_list.clear()
        for i in range(2):  # 选择前两个（Human, Chimpanzee）
            item = QListWidgetItem(example_taxa[i])
            self.selected_taxa_list.addItem(item)
        
        # 更新左侧显示，移除已选中的分类单元
        self._update_available_taxa_display()
    
    def _update_available_taxa_display(self):
        """更新左侧可用分类单元的显示，排除已选中的"""
        # 获取当前选中的分类单元
        selected_names = {self.selected_taxa_list.item(i).text() 
                         for i in range(self.selected_taxa_list.count())}
        
        # 获取应该显示的分类单元（所有 - 已选中）
        available_names = self.all_taxa - selected_names
        
        # 应用过滤器
        filter_text = self.otu_filter.text().lower()
        if filter_text:
            available_names = {name for name in available_names if filter_text in name.lower()}
        
        # 更新显示
        self.taxon_list.clear()
        for name in sorted(available_names):
            item = QListWidgetItem(name)
            self.taxon_list.addItem(item)
    
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

        self.tmrca_type_combo = QComboBox()
        self.tmrca_type_combo.addItems(['Point', 'Uniform', 'Upper Boundary', 'Lower Boundary', 'Normal', 'Lognormal'])
        self.tmrca_type_combo.currentTextChanged.connect(self.on_tmrca_type_changed)
        tmrca_type_layout.addWidget(self.tmrca_type_combo)

        tmrca_setting_group = QGroupBox("tMRCA parameters")
        self.tmrca_setting_group_layout = QFormLayout()
        tmrca_setting_group.setLayout(self.tmrca_setting_group_layout)

        tmrca_left_layout.addWidget(tmrca_setting_group)

        # matplotlib figure canvas to display tMRCA distribution
        self.tmrca_figure_canvas = MplDistribution()
        tmrca_set_layout.addWidget(self.tmrca_figure_canvas)

        # 初始化参数输入字段
        self.init_tmrca_parameters()
        
    def init_tree_view(self):
        # 使用封装的MplTreeView组件来绘制树状图
        self.tree_figure_canvas = MplTreeView()
        self.tree_figure_canvas.draw_dendrogram("(((1:0.1, 2:0.1):0.1,3:0.3):0.5,4:0.9);")
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
        self.point_value.setPlaceholderText("e.g., 1.5")
        self.param_widgets['Point'] = [('Age (time units):', self.point_value)]
        
        # Uniform distribution parameters
        self.uniform_lower = QLineEdit()
        self.uniform_lower.setText("0.5")
        self.uniform_lower.setPlaceholderText("e.g., 0.8")
        self.uniform_upper = QLineEdit()
        self.uniform_upper.setText("1.5")
        self.uniform_upper.setPlaceholderText("e.g., 2.0")
        self.param_widgets['Uniform'] = [('Minimum age:', self.uniform_lower), 
                                        ('Maximum age:', self.uniform_upper)]
        
        # Upper Boundary parameters
        self.upper_boundary = QLineEdit()
        self.upper_boundary.setText("2.0")
        self.upper_boundary.setPlaceholderText("e.g., 3.0")
        self.param_widgets['Upper Boundary'] = [('Maximum age:', self.upper_boundary)]
        
        # Lower Boundary parameters
        self.lower_boundary = QLineEdit()
        self.lower_boundary.setText("0.5")
        self.lower_boundary.setPlaceholderText("e.g., 0.2")
        self.param_widgets['Lower Boundary'] = [('Minimum age:', self.lower_boundary)]
        
        # Normal distribution parameters
        self.normal_mean = QLineEdit()
        self.normal_mean.setText("1.0")
        self.normal_mean.setPlaceholderText("e.g., 1.2")
        self.normal_std = QLineEdit()
        self.normal_std.setText("0.2")
        self.normal_std.setPlaceholderText("e.g., 0.3")
        self.param_widgets['Normal'] = [('Mean age:', self.normal_mean), 
                                       ('Standard deviation:', self.normal_std)]
        
        # Lognormal distribution parameters
        self.lognormal_mean = QLineEdit()
        self.lognormal_mean.setText("1.0")
        self.lognormal_mean.setPlaceholderText("e.g., 1.5")
        self.lognormal_std = QLineEdit()
        self.lognormal_std.setText("0.2")
        self.lognormal_std.setPlaceholderText("e.g., 0.4")
        self.param_widgets['Lognormal'] = [('Mean age:', self.lognormal_mean), 
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
        """根据过滤文本筛选OTU列表"""
        # 过滤逻辑已整合到 _update_available_taxa_display 中
        self._update_available_taxa_display()
    
    def add_taxa(self):
        """将选中的分类单元从OTU列表添加到选中列表"""
        selected_items = self.taxon_list.selectedItems()
        if not selected_items:
            return
            
        current_selected_names = {self.selected_taxa_list.item(i).text() 
                                for i in range(self.selected_taxa_list.count())}
        
        for item in selected_items:
            taxon_name = item.text()
            if taxon_name not in current_selected_names:
                new_item = QListWidgetItem(taxon_name)
                self.selected_taxa_list.addItem(new_item)
                current_selected_names.add(taxon_name)
        
        # 更新左侧显示
        self._update_available_taxa_display()
        
        # 清除选择
        self.taxon_list.clearSelection()
    
    def remove_taxa(self):
        """从选中列表中移除选中的分类单元"""
        selected_items = self.selected_taxa_list.selectedItems()
        if not selected_items:
            return
            
        for item in selected_items:
            row = self.selected_taxa_list.row(item)
            self.selected_taxa_list.takeItem(row)
        
        # 更新左侧显示，重新显示被移除的分类单元
        self._update_available_taxa_display()
        
        # 清除选择
        self.selected_taxa_list.clearSelection()
    
    def clear_taxa(self):
        """清空所有选中的分类单元"""
        self.selected_taxa_list.clear()
        # 重新显示所有分类单元
        self._update_available_taxa_display()

def main():
    app = QApplication(sys.argv)  # 使用sys.argv
    
    # 创建主窗口
    main_window = QMainWindow()
    main_window.setWindowTitle("MD-MRCA Application")
    
    # 创建MDMRCA部件并设置为中央部件，启用示例数据
    md_mrca_widget = MDMRCA(example=True)
    # 已经在__init__中调用了init_ui()，不需要再次调用
    main_window.setCentralWidget(md_mrca_widget)
    
    # 显示窗口
    main_window.show()
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()