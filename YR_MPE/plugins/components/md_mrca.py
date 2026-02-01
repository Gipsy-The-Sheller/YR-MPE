from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QLabel, QLineEdit,
                             QListWidget, QListWidgetItem, QPushButton, QComboBox,
                             QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
                             QSizePolicy, QMessageBox, QSplitter)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# 修复相对导入
from .mpl_treeview import MplTreeView
from .mpl_distribution import MplDistribution
# 导入streamlined_ete模块
from .methods.streamlined_ete import build_tree_block, find_mrca_position, inspect_newick_string
import sys

class MDMRCA(QWidget):
    def __init__(self, parent=None, example=False):
        super().__init__(parent)
        self.setWindowTitle("MRCA Annotation")
        self.resize(800, 600)
        
        # 初始化Newick字符串（示例数据）
        if example:
            self.newick_string = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
        else:
            self.newick_string = ""
            
        # 新增：用于存储带MRCA注释的Newick字符串
        self.annotated_newick_str = self.newick_string
        
        self.param_widgets = {}  # 初始化param_widgets属性
        self.init_ui()

        if example:
            self.load_example_data()
    
    def init_ui(self):
        self.setWindowTitle("MD-MRCA")
        # self.setGeometry(100, 100, 800, 600)

        # 使用QSplitter替代HLayout
        self.splitter = QSplitter()
        self.setLayout(QHBoxLayout())
        self.layout().addWidget(self.splitter)

        # 创建左侧控件容器
        self.left_widget = QWidget()
        self.left_layout = QVBoxLayout()
        self.left_layout.setContentsMargins(0,0,0,0)
        self.left_widget.setLayout(self.left_layout)

        self.splitter.addWidget(self.left_widget)

        self.init_taxon_set()
        self.init_tmrca_set()
        self.init_tree_view()
        
        # 设置分割器比例 6:4 (左侧60%，右侧40%)
        self.splitter.setSizes([600, 400])  # 初始大小比例
    
    def extract_all_taxa_from_newick(self, newick_string):
        """
        从Newick字符串中提取所有OTU名称
        """
        if not newick_string:
            return []
            
        taxa = []
        current_otu = ""
        in_quotes = False
        quote_char = None
        current_depth = 0
        
        i = 0
        n = len(newick_string)
        
        while i < n:
            char = newick_string[i]
            
            if char == "'" or char == '"':
                if not in_quotes:
                    in_quotes = True
                    quote_char = char
                    current_otu += char
                elif char == quote_char:
                    in_quotes = False
                    quote_char = None
                    current_otu += char
                else:
                    current_otu += char
            elif in_quotes:
                current_otu += char
            elif char == ':':
                # 遇到冒号，说明后面是分支长度，当前OTU结束
                if current_otu:
                    otu_name = current_otu.strip()
                    if otu_name:
                        taxa.append(otu_name)
                    current_otu = ""
                # 跳过分支长度部分（直到遇到分隔符）
                i += 1
                while i < n and newick_string[i] not in '(),;':
                    i += 1
                continue
            elif char == '(':
                current_depth += 1
                if current_otu:
                    otu_name = current_otu.strip()
                    if otu_name:
                        taxa.append(otu_name)
                    current_otu = ""
            elif char == ')' or char == ',' or char == ';':
                if current_otu:
                    otu_name = current_otu.strip()
                    if otu_name:
                        taxa.append(otu_name)
                    current_otu = ""
                if char == ')':
                    current_depth -= 1
                elif char == ';':
                    break
            else:
                current_otu += char
                
            i += 1
        
        return taxa
    
    def find_mrca_position(self, newick_string, target_taxa_set):
        """
        使用streamlined_ete中的鲁棒MRCA定位算法
        
        算法步骤：
        1. 将Newick字符串分割成块
        2. 使用基于栈深度的算法找到MRCA位置
        3. 返回MRCA在块列表中的索引位置
        """
        if not newick_string or not target_taxa_set:
            return None
            
        try:
            # 构建树块
            tree_blocks = build_tree_block(newick_string)
            
            # 创建目标taxa集合的副本，因为find_mrca_position会修改它
            target_taxa_list = list(target_taxa_set)
            
            # 查找MRCA位置
            mrca_index = find_mrca_position(tree_blocks, target_taxa_list)
            
            if mrca_index is None:
                return None
                
            # 验证是否所有目标taxa都被找到
            if len(target_taxa_list) > 0:
                return None
                
            # 返回MRCA位置信息
            # 注意：streamlined_ete返回的是块索引，我们需要转换为字符串位置
            # 为了兼容现有代码，我们返回一个元组 (start_pos, end_pos, depth)
            # 这里我们计算实际的字符串位置
            mrca_block = tree_blocks[mrca_index]
            if mrca_block.startswith(')'):
                # MRCA是一个右括号，找到它在原字符串中的位置
                # 重建newick字符串来定位
                reconstructed = inspect_newick_string(tree_blocks)
                # 找到对应的右括号位置
                pos = -1
                block_count = 0
                for i, char in enumerate(reconstructed):
                    if char in '()':
                        if block_count == mrca_index:
                            pos = i
                            break
                        block_count += 1
                    elif char == ',' or char == ';':
                        if block_count == mrca_index:
                            pos = i
                            break
                        block_count += 1
                
                if pos >= 0:
                    return (pos, pos, 0)  # depth信息暂时设为0，因为新算法不直接提供
                    
            return None
            
        except Exception as e:
            # 如果新算法失败，可以考虑回退到旧算法或返回None
            print(f"MRCA定位算法错误: {e}")
            return None

    def insert_mrca_label(self, newick_string, target_taxa_set, label_name):
        """
        在Newick字符串中插入MRCA标签 - 使用新的鲁棒算法
        """
        result = self.find_mrca_position(newick_string, target_taxa_set)
        if not result:
            return newick_string
            
        start_pos, end_pos, depth = result
        
        # 在MRCA节点后插入标签
        if end_pos < len(newick_string) and newick_string[end_pos] == ')':
            # 在右括号后插入
            labeled_newick = (newick_string[:end_pos + 1] + 
                            f"[&name={label_name}]" + 
                            newick_string[end_pos + 1:])
        else:
            # 单个taxa的情况或其他情况
            labeled_newick = (newick_string[:end_pos] + 
                            f"[&name={label_name}]" + 
                            newick_string[end_pos:])
        
        return labeled_newick
    
    def load_example_data(self):
        """加载示例数据：一棵系统发育树和对应的分类单元"""
        # 示例Newick格式的树
        example_newick = "(((A:0.1,B:0.1):0.1,C:0.2):0.1, (D:0.1,E:0.1):0.2);"
        
        # 存储原始Newick字符串
        self.newick_string = example_newick
        
        # 从树中提取所有叶子节点（分类单元）
        example_taxa = self.extract_all_taxa_from_newick(example_newick)
        
        # 存储所有可用的分类单元
        self.all_taxa = set(example_taxa)
        
        # 填充OTU浏览器中的分类单元列表（初始时显示所有）
        self.taxon_list.clear()
        for taxon in example_taxa:
            item = QListWidgetItem(taxon)
            self.taxon_list.addItem(item)
        
        # 绘制示例树
        self.tree_figure_canvas.draw_dendrogram(newick_str=example_newick)
        
        # 预选择一些分类单元作为示例
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
        self.selected_taxa_list.setSelectionMode(QListWidget.ExtendedSelection)
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
        # 将树视图添加到splitter而不是main_layout
        self.splitter.addWidget(self.tree_figure_canvas)

        # 示例：绘制一个默认的树
        # newick_str = "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0,E:1.0):1.0);"
        # self.tree_figure_canvas.draw_dendrogram(newick_str=newick_str)
        
    def init_tmrca_parameters(self):
        """初始化tMRCA参数输入字段 - 预创建所有控件"""
        # Point distribution parameters
        self.point_value = QLineEdit()
        self.point_value.setText("1.0")
        self.point_value.setPlaceholderText("e.g., 1.5")
        self.point_value.textChanged.connect(self.on_tmrca_params_changed)
        
        # Uniform distribution parameters
        self.uniform_lower = QLineEdit()
        self.uniform_lower.setText("0.5")
        self.uniform_lower.setPlaceholderText("e.g., 0.8")
        self.uniform_lower.textChanged.connect(self.on_tmrca_params_changed)
        self.uniform_upper = QLineEdit()
        self.uniform_upper.setText("1.5")
        self.uniform_upper.setPlaceholderText("e.g., 2.0")
        self.uniform_upper.textChanged.connect(self.on_tmrca_params_changed)
        
        # Upper Boundary parameters
        self.upper_boundary = QLineEdit()
        self.upper_boundary.setText("2.0")
        self.upper_boundary.setPlaceholderText("e.g., 3.0")
        self.upper_boundary.textChanged.connect(self.on_tmrca_params_changed)
        
        # Lower Boundary parameters
        self.lower_boundary = QLineEdit()
        self.lower_boundary.setText("0.5")
        self.lower_boundary.setPlaceholderText("e.g., 0.2")
        self.lower_boundary.textChanged.connect(self.on_tmrca_params_changed)
        
        # Normal distribution parameters
        self.normal_mean = QLineEdit()
        self.normal_mean.setText("1.0")
        self.normal_mean.setPlaceholderText("e.g., 1.2")
        self.normal_mean.textChanged.connect(self.on_tmrca_params_changed)
        self.normal_std = QLineEdit()
        self.normal_std.setText("0.2")
        self.normal_std.setPlaceholderText("e.g., 0.3")
        self.normal_std.textChanged.connect(self.on_tmrca_params_changed)
        
        # Lognormal distribution parameters
        self.lognormal_mean = QLineEdit()
        self.lognormal_mean.setText("1.0")
        self.lognormal_mean.setPlaceholderText("e.g., 1.5")
        self.lognormal_mean.textChanged.connect(self.on_tmrca_params_changed)
        self.lognormal_std = QLineEdit()
        self.lognormal_std.setText("0.2")
        self.lognormal_std.setPlaceholderText("e.g., 0.4")
        self.lognormal_std.textChanged.connect(self.on_tmrca_params_changed)
        
        # 添加所有控件到布局中，但初始时只显示Point类型的
        self.tmrca_setting_group_layout.addRow("Age (time units):", self.point_value)
        self.tmrca_setting_group_layout.addRow("Minimum age:", self.uniform_lower)
        self.tmrca_setting_group_layout.addRow("Maximum age:", self.uniform_upper)
        self.tmrca_setting_group_layout.addRow("Maximum age:", self.upper_boundary)
        self.tmrca_setting_group_layout.addRow("Minimum age:", self.lower_boundary)
        self.tmrca_setting_group_layout.addRow("Mean age:", self.normal_mean)
        self.tmrca_setting_group_layout.addRow("Standard deviation:", self.normal_std)
        self.tmrca_setting_group_layout.addRow("Mean age:", self.lognormal_mean)
        self.tmrca_setting_group_layout.addRow("Standard deviation:", self.lognormal_std)
        
        # 存储所有控件的引用
        self.all_param_widgets = {
            'Point': [self.point_value],
            'Uniform': [self.uniform_lower, self.uniform_upper],
            'Upper Boundary': [self.upper_boundary],
            'Lower Boundary': [self.lower_boundary],
            'Normal': [self.normal_mean, self.normal_std],
            'Lognormal': [self.lognormal_mean, self.lognormal_std]
        }
        
        # 设置初始参数显示
        self.update_tmrca_parameters('Point')
        
    def update_tmrca_parameters(self, distribution_type):
        """根据分布类型更新参数输入字段的可见性"""
        # 隐藏所有控件
        for widget_list in self.all_param_widgets.values():
            for widget in widget_list:
                widget.hide()
                # 隐藏对应的标签
                label = self.tmrca_setting_group_layout.labelForField(widget)
                if label:
                    label.hide()
        
        # 显示当前分布类型的控件
        if distribution_type in self.all_param_widgets:
            for widget in self.all_param_widgets[distribution_type]:
                widget.show()
                # 显示对应的标签
                label = self.tmrca_setting_group_layout.labelForField(widget)
                if label:
                    label.show()
        
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
        
        # 更新树视图高亮
        self.update_tree_highlight()
    
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
        
        # 更新树视图高亮
        self.update_tree_highlight()
    
    def clear_taxa(self):
        """清空所有选中的分类单元"""
        self.selected_taxa_list.clear()
        # 重新显示所有分类单元
        self._update_available_taxa_display()
        
        # 更新树视图高亮
        self.update_tree_highlight()
    
    def update_tree_highlight(self):
        """更新树的高亮显示"""
        # 获取当前选中的分类单元
        selected_taxa = [self.selected_taxa_list.item(i).text() 
                        for i in range(self.selected_taxa_list.count())]
        
        # 获取MRCA名称 - 优先使用输入框中的名称，如果有的话
        mrca_name = self.taxon_set_name.text().strip()
        if not mrca_name:
            # 如果输入框为空，尝试从Newick字符串中提取
            mrca_name = self.get_current_mrca_name()
        
        # 构建args字典
        args = {}
        if selected_taxa:
            args['taxon_set'] = selected_taxa
        if mrca_name:
            args['mrca_name'] = mrca_name
            
        # 如果有MRCA名称和选中的分类单元，更新带注释的Newick字符串
        if mrca_name and selected_taxa:
            self.annotated_newick_str = self.insert_mrca_label(self.newick_string, selected_taxa, mrca_name)
        else:
            # 如果没有MRCA信息，使用原始Newick字符串
            self.annotated_newick_str = self.newick_string
            
        # 更新树视图
        self.tree_figure_canvas.draw_dendrogram(newick_str=self.annotated_newick_str, args=args)

    def get_current_mrca_name(self):
        """获取当前的MRCA标注名称"""
        if not self.annotated_newick_str:
            return None
            
        # 首先尝试从输入框获取（保持向后兼容）
        taxon_set_name = self.taxon_set_name.text().strip()
        if taxon_set_name and f"[&name={taxon_set_name}]" in self.annotated_newick_str:
            return taxon_set_name
            
        # 如果输入框为空或不匹配，直接从带注释的Newick字符串中提取MRCA名称
        import re
        # 匹配 [&name=...] 格式的注释
        match = re.search(r'\[&name=([^\]]+)\]', self.annotated_newick_str)
        if match:
            extracted_name = match.group(1)
            # 验证提取的名称是否有效（非空且不包含特殊字符）
            if extracted_name.strip():
                return extracted_name.strip()
                
        return None
    
    def apply_mrca_annotation(self):
        """应用MRCA标注并更新树视图"""
        selected_taxa = [self.selected_taxa_list.item(i).text() 
                        for i in range(self.selected_taxa_list.count())]
        
        if not selected_taxa:
            QMessageBox.warning(self, "Warning", "Please select at least one taxon.")
            return
            
        taxon_set_name = self.taxon_set_name.text().strip()
        if not taxon_set_name:
            QMessageBox.warning(self, "Warning", "Please enter a taxon set name.")
            return
        
        # 检查是否已经存在相同的MRCA注释（在annotated_newick_str中）
        if taxon_set_name and f"[&name={taxon_set_name}]" in self.annotated_newick_str:
            reply = QMessageBox.question(
                self, "Confirm Overwrite", 
                f"MRCA annotation '{taxon_set_name}' already exists. Overwrite?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        # 应用MRCA标注到annotated_newick_str
        self.annotated_newick_str = self.insert_mrca_label(self.newick_string, selected_taxa, taxon_set_name)
        
        # 更新树视图（包含高亮）
        args = {
            'taxon_set': selected_taxa,
            'mrca_name': taxon_set_name
        }
        self.tree_figure_canvas.draw_dendrogram(newick_str=self.annotated_newick_str, args=args)
        
        QMessageBox.information(self, "Success", f"MRCA annotation '{taxon_set_name}' applied successfully!")
        
def main():
    from PyQt5.QtWidgets import QApplication
    import sys
    
    app = QApplication(sys.argv)
    md_mrca_widget = MDMRCA(example=True)
    md_mrca_widget.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()