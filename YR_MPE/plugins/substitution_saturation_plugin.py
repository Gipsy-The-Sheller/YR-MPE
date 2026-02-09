# substitution_saturation_plugin.py
#
# Copyright (c) 2026 YRTools Team
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
import sys
import tempfile
import warnings
warnings.filterwarnings("ignore", module="matplotlib")

import numpy as np
from typing import List, Dict, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QLabel, QComboBox, QTextEdit, QTabWidget, QFrame,
                             QDialog, QTableWidget, QTableWidgetItem, QHeaderView, QSplitter,
                             QDialogButtonBox, QScrollArea)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont

import matplotlib
matplotlib.use('Qt5Agg')  # 使用Qt5后端
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from scipy import stats

from ..templates.base_plugin_ui import BasePlugin

# 相对导入距离计算工具
distance_utils_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                                   'plugins', 'components', 'methods')
sys.path.insert(0, distance_utils_path)
try:
    import distance_utils
    from distance_utils import (
        calculate_saturation_metrics,
        validate_sequences
    )
except ImportError:
    # 备用导入方式
    from .components.methods.distance_utils import (
        calculate_saturation_metrics,
        validate_sequences
    )


class MethodDescriptionDialog(QDialog):
    """方法说明对话框"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Method Description - Substitution Saturation Analysis")
        self.setMinimumSize(600, 500)
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # 创建滚动区域
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        
        # 创建内容widget
        content_widget = QWidget()
        content_layout = QVBoxLayout()
        content_widget.setLayout(content_layout)
        
        # 方法说明文本
        method_text = QLabel("""
        <h2 style="color: #333;">Substitution Saturation Analysis (BaCoCa Method)</h2>
        
        <h3 style="color: #555;">Overview</h3>
        <p style="font-size: 20px;">This analysis calculates the substitution saturation index C based on the BaCoCa method, which evaluates the degree to which nucleotide substitutions have reached saturation in multiple sequence alignments (MSAs).</p>
        
        <h3 style="color: #555;">Calculation Formula</h3>
        <p style="font-size: 20px;"><b>C = std(Ti/Tv) / std(p)</b></p>
<p style="font-size: 20px;">where <b>std</b> represents the standard deviation.</p>
        
        <h3 style="color: #555;">Parameters</h3>
        <ul style="font-size: 20px;">
        <li><b>Ti (Transitions)</b>: Number of transition substitutions (A↔G, T↔C). Transitions are purine↔purine or pyrimidine↔pyrimidine substitutions.</li>
        <li><b>Tv (Transversions)</b>: Number of transversion substitutions (A↔T, A↔C, G↔T, G↔C). Transversions are purine↔pyrimidine substitutions.</li>
        <li><b>p (p-distance)</b>: Uncorrected genetic distance = (Ti + Tv) / total compared sites</li>
        </ul>
        
        <h3 style="color: #555;">Interpretation of C-value</h3>
        <p style="font-size: 20px;">The C value provides a quantitative measure of substitution saturation based on the ratio of standard deviations. Higher C values generally indicate greater variability in Ti/Tv ratios relative to p-distances. However, specific interpretation thresholds should be determined empirically for each dataset, as C values can vary depending on the evolutionary history and characteristics of the sequences being analyzed.</p>
        
        <h3 style="color: #555;">Plot Description</h3>
        <p style="font-size: 20px;">The plot displays substitution counts (Ti and Tv) plotted against p-distance:</p>
        <ul style="font-size: 20px;">
        <li><b>Red points (#ba3e45)</b>: Transition (Ti) counts vs. p-distance</li>
        <li><b>Blue points (#4e699a)</b>: Transversion (Tv) counts vs. p-distance</li>
        <li><b>Dashed lines</b>: Linear regression trend lines with correlation coefficients (r)</li>
        </ul>
        <p style="font-size: 20px;">Transitions (Ti) typically accumulate faster than transversions (Tv) when there is little or no saturation. Substitution saturation is indicated when the Ti and Tv trend lines converge or overlap, suggesting that the preferential accumulation of transitions has been diminished due to multiple hits at the same sites.</p>
        
        <h3 style="color: #555;">Gap Treatment Options</h3>
        <ul style="font-size: 20px;">
        <li><b>Pairwise Deletion</b>: Gaps are ignored pairwise for each sequence pair comparison. This maximizes data usage but may introduce bias.</li>
        <li><b>Complete Deletion</b>: Sites containing gaps in any sequence are completely removed from the analysis. This is more conservative but reduces the amount of usable data.</li>
        </ul>
        
        <h3 style="color: #555;">Citation</h3>
        <p style="font-size: 20px;">Kück, P., & Struck, T. H. (2014). BaCoCa--a heuristic software tool for the parallel assessment of sequence biases in hundreds of gene and taxon partitions. <i>Molecular Phylogenetics and Evolution</i>, 70, 94–98. https://doi.org/10.1016/j.ympev.2013.09.011</p>
        """)
        method_text.setWordWrap(True)
        method_text.setOpenExternalLinks(True)
        method_text.setTextFormat(Qt.RichText)
        content_layout.addWidget(method_text)
        
        scroll.setWidget(content_widget)
        layout.addWidget(scroll)
        
        # 添加关闭按钮
        button_box = QDialogButtonBox(QDialogButtonBox.Close)
        button_box.rejected.connect(self.accept)
        layout.addWidget(button_box)


class SaturationCanvas(FigureCanvas):
    """用于显示饱和度分析的Matplotlib画布"""
    
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.tight_layout()
        self.ax = self.fig.add_subplot(111)
        super(SaturationCanvas, self).__init__(self.fig)
        self.setParent(parent)
        
        # 设置中文字体支持
        plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
    
    def plot_saturation(self, p_values: List[float], ti_values: List[float], 
                       tv_values: List[float], c_value: float, 
                       gap_treatment: str = 'pairwise'):
        """
        绘制Ti vs p和Tv vs p的散点图
        
        Args:
            p_values: p值列表
            ti_values: Ti值列表
            tv_values: Tv值列表
            c_value: C值
            gap_treatment: gap处理方式
        """
        self.ax.clear()
        
        # 转换为numpy数组以便计算
        p_array = np.array(p_values)
        ti_array = np.array(ti_values)
        tv_array = np.array(tv_values)
        
        # 绘制Ti vs p散点图（红色 #ba3e45）
        self.ax.scatter(p_array, ti_array, alpha=0.6, color='#ba3e45', 
                       label='Transitions (Ti)', s=50, edgecolors='white', linewidth=0.5)
        
        # 绘制Tv vs p散点图（蓝色 #4e699a）
        self.ax.scatter(p_array, tv_array, alpha=0.6, color='#4e699a', 
                       label='Transversions (Tv)', s=50, edgecolors='white', linewidth=0.5)
        
        # 添加线性回归趋势线
        if len(p_array) > 1:
            # Ti趋势线
            ti_slope, ti_intercept, ti_r, ti_p, ti_std_err = stats.linregress(p_array, ti_array)
            ti_line_x = np.linspace(p_array.min(), p_array.max(), 100)
            ti_line_y = ti_slope * ti_line_x + ti_intercept
            self.ax.plot(ti_line_x, ti_line_y, '--', color='#ba3e45', alpha=0.8, linewidth=2,
                        label=f'Ti trend (r={ti_r:.3f})')
            
            # Tv趋势线
            tv_slope, tv_intercept, tv_r, tv_p, tv_std_err = stats.linregress(p_array, tv_array)
            tv_line_x = np.linspace(p_array.min(), p_array.max(), 100)
            tv_line_y = tv_slope * tv_line_x + tv_intercept
            self.ax.plot(tv_line_x, tv_line_y, '--', color='#4e699a', alpha=0.8, linewidth=2,
                        label=f'Tv trend (r={tv_r:.3f})')
        
        # 设置标签和标题
        self.ax.set_xlabel('p-distance (Genetic Distance)', fontsize=12, fontweight='bold')
        self.ax.set_ylabel('Substitution Count', fontsize=12, fontweight='bold')
        self.ax.set_title(f'Substitution Saturation Analysis\nC-value = {c_value:.4f} | Gap Treatment: {gap_treatment}', 
                         fontsize=14, fontweight='bold', pad=15)
        
        # 设置网格
        self.ax.grid(True, alpha=0.3, linestyle='--')
        
        # 设置图例
        self.ax.legend(loc='best', frameon=True, shadow=True, fontsize=10)
        
        # 设置背景色
        self.ax.set_facecolor('#f8f9fa')
        self.fig.patch.set_facecolor('white')
        
        # 调整布局
        self.fig.tight_layout()
        
        # 重绘画布
        self.draw()


class OverviewCanvas(FigureCanvas):
    """用于显示总体C值比较的Matplotlib画布"""
    
    def __init__(self, parent=None, width=8, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.tight_layout()
        self.ax = self.fig.add_subplot(111)
        super(OverviewCanvas, self).__init__(self.fig)
        self.setParent(parent)
        
        # 设置中文字体支持
        plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
    
    def plot_c_values(self, file_names: List[str], c_values: List[float]):
        """
        绘制总体C值柱状图
        
        Args:
            file_names: 文件名列表
            c_values: C值列表
        """
        self.ax.clear()
        
        if not c_values:
            self.ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                        transform=self.ax.transAxes, fontsize=12)
            self.draw()
            return
        
        # 创建渐变红色
        import matplotlib.colors as mcolors
        import matplotlib.cm as cm
        
        # 找到最大和最小C值
        max_c = max(c_values)
        min_c = min(c_values)
        
        # 为每个柱子创建颜色（渐变红色，最高处#ba3e45）
        colors = []
        for c in c_values:
            if max_c == min_c:
                # 所有C值相同，使用中等红色
                norm_c = 0.5
            else:
                # 归一化C值（0-1范围）
                norm_c = (c - min_c) / (max_c - min_c)
            
            # 从浅红色到#ba3e45的渐变
            # 浅红色：#ffe5e5 (255, 229, 229)
            # 目标红色：#ba3e45 (186, 62, 69)
            r = int(255 - (255 - 186) * norm_c)
            g = int(229 - (229 - 62) * norm_c)
            b = int(229 - (229 - 69) * norm_c)
            colors.append((r/255, g/255, b/255))
        
        # 绘制柱状图
        x_pos = np.arange(len(file_names))
        bars = self.ax.bar(x_pos, c_values, color=colors, edgecolor='white', linewidth=1.5, alpha=0.9)
        
        # 设置标签
        self.ax.set_xlabel('Files', fontsize=11, fontweight='bold')
        self.ax.set_ylabel('C-value', fontsize=11, fontweight='bold')
        self.ax.set_title('C-values Comparison Overview', fontsize=12, fontweight='bold', pad=10)
        
        # 设置x轴标签（旋转以避免重叠）
        short_names = [os.path.basename(name)[:15] + '...' if len(os.path.basename(name)) > 15 
                       else os.path.basename(name) for name in file_names]
        self.ax.set_xticks(x_pos)
        self.ax.set_xticklabels(short_names, rotation=45, ha='right', fontsize=9)
        
        # 在柱子上显示数值
        for bar, c in zip(bars, c_values):
            height = bar.get_height()
            self.ax.text(bar.get_x() + bar.get_width()/2., height,
                        f'{c:.3f}',
                        ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        # 设置网格
        self.ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        # 设置背景色
        self.ax.set_facecolor('#f8f9fa')
        self.fig.patch.set_facecolor('white')
        
        # 调整布局
        self.fig.tight_layout()
        
        # 重绘画布
        self.draw()


class SubstitutionSaturationPlugin(BasePlugin):
    """替代饱和度分析插件"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化替代饱和度分析插件"""
        super().__init__(import_from, import_data)
        
        # 初始化用于存储计算结果的属性
        self.saturation_results = None
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
    
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "Substitution Saturation Analysis"
        self.tool_name = "Built-in Python Implementation (BaCoCa Method)"
        self.citation = [
            """Kück, P., & Struck, T. H. (2014). BaCoCa--a heuristic software tool for the parallel assessment of sequence biases in hundreds of gene and taxon partitions. Molecular Phylogenetics and Evolution, 70, 94–98. https://doi.org/10.1016/j.ympev.2013.09.011"""
        ]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"], "NEXUS": ["nex", "nexus"]}
        self.output_types = {"Saturation Plot": ".png", "Report": ".txt"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    
    def config(self):
        """
        检查插件配置
        对于纯Python实现的插件，始终返回True
        """
        return True
    
    def setup_ui(self):
        """设置用户界面"""
        # 创建标签页
        self.tab_widget = QTabWidget()
        self.layout().addWidget(self.tab_widget)
        
        # 输入/参数标签页
        self.input_tab = QWidget()
        self.setup_input_tab()
        self.tab_widget.addTab(self.input_tab, "Input & Parameters")
        
        # 输出标签页
        self.output_tab = QWidget()
        self.setup_output_tab()
        self.tab_widget.addTab(self.output_tab, "Output")
        
        # 控制面板
        self.setup_control_panel()
        
        # 设置默认值
        self.reset_to_defaults()
    
    def setup_input_tab(self):
        """设置输入和参数标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # 输入组
        input_group = QGroupBox("Input")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # 文件输入
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select alignment files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)
        
        # 文件标签容器
        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # 文本输入
        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste alignment in FASTA format...")
        self.sequence_text.setMaximumHeight(200)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        input_layout.addRow("Sequence text:", self.sequence_text)
        
        # 处理导入的数据
        if self.import_file:
            self.file_path_edit.setText(self.import_file)
            input_group.setVisible(False)
        elif hasattr(self, 'imported_files') and self.imported_files:
            # 显示导入的文件
            for file_path in self.imported_files:
                self.add_file_tag(file_path)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            # 隐藏文本输入框
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
        
        # 参数组
        params_group = QGroupBox("Analysis Parameters")
        params_form_layout = QFormLayout()
        params_group.setLayout(params_form_layout)
        layout.addWidget(params_group)
        
        # Gap处理方式
        self.gap_treatment_combo = QComboBox()
        self.gap_treatment_combo.addItems(["Pairwise Deletion", "Complete Deletion"])
        params_form_layout.addRow("Gap Treatment:", self.gap_treatment_combo)
        
        # 方法说明按钮
        method_info_layout = QHBoxLayout()
        method_info_label = QLabel("Method: BaCoCa saturation index")
        method_info_layout.addWidget(method_info_label)
        
        self.method_info_btn = QPushButton("View Method Description")
        self.method_info_btn.clicked.connect(self.show_method_description)
        method_info_layout.addWidget(self.method_info_btn)
        method_info_layout.addStretch()
        
        params_form_layout.addRow("", method_info_layout)
        
        layout.addStretch()
        
        if not hasattr(self, 'imported_files'):
            self.imported_files = []
        if not hasattr(self, 'file_tags'):
            self.file_tags = []
    
    def setup_output_tab(self):
        """设置输出标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 文件选择区域（批量处理时显示）
        file_selection_group = QGroupBox("File Selection")
        file_selection_layout = QHBoxLayout()
        file_selection_group.setLayout(file_selection_layout)
        
        self.file_combo_label = QLabel("Select File:")
        file_selection_layout.addWidget(self.file_combo_label)
        
        self.file_combo = QComboBox()
        self.file_combo.currentIndexChanged.connect(self.on_file_selection_changed)
        file_selection_layout.addWidget(self.file_combo)
        
        self.file_combo_label.setVisible(False)
        self.file_combo.setVisible(False)
        
        layout.addWidget(file_selection_group)
        
        # 创建分割器
        splitter = QSplitter(Qt.Vertical)
        
        # 图表区域
        chart_group = QGroupBox("Substitution Saturation Plot")
        chart_layout = QVBoxLayout()
        chart_group.setLayout(chart_layout)
        
        # 创建Matplotlib画布（散点图）
        self.canvas = SaturationCanvas(self, width=8, height=6, dpi=100)
        chart_layout.addWidget(self.canvas)
        
        # 导出图表按钮
        self.export_plot_btn = QPushButton("Export Current Plot")
        self.export_plot_btn.clicked.connect(self.export_plot)
        self.export_plot_btn.setEnabled(False)
        chart_layout.addWidget(self.export_plot_btn)
        
        splitter.addWidget(chart_group)
        
        # 总体C值柱状图区域
        overview_group = QGroupBox("Overview: C-values Comparison")
        overview_layout = QVBoxLayout()
        overview_group.setLayout(overview_layout)
        
        # 创建Matplotlib画布（柱状图）
        self.overview_canvas = OverviewCanvas(self, width=8, height=4, dpi=100)
        overview_layout.addWidget(self.overview_canvas)
        
        splitter.addWidget(overview_group)
        
        # 统计信息区域
        stats_group = QGroupBox("Statistics Summary")
        stats_layout = QVBoxLayout()
        stats_group.setLayout(stats_layout)
        
        self.stats_text = QTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setMaximumHeight(200)
        stats_layout.addWidget(self.stats_text)
        
        # splitter.addWidget(stats_group)
        
        # 设置分割器比例
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        splitter.setStretchFactor(2, 1)
        
        layout.addWidget(splitter)
        
        # 初始化结果存储
        if not hasattr(self, 'saturation_results'):
            self.saturation_results = []  # 存储所有文件的结果列表
    
    def on_file_selection_changed(self, index):
        """当文件选择改变时更新图表"""
        if not self.saturation_results or index < 0 or index >= len(self.saturation_results):
            return
        
        result = self.saturation_results[index]
        
        # 更新散点图
        self.canvas.plot_saturation(
            result['p_values'],
            result['ti_values'],
            result['tv_values'],
            result['c_value'],
            result['gap_treatment']
        )
        
        # 更新统计信息
        self.display_single_result(result, result['filename'])
    
    def setup_control_panel(self):
        """设置控制面板"""
        # 先调用父类方法设置基本控件
        super().setup_control_panel()
    
    def browse_files(self):
        """浏览选择文件"""
        file_filter = "Alignment files (*.fas *.fna *.fa *.fasta *.nex *.nexus);;All files (*)"
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select alignment files", "", file_filter)
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
    
    def add_file_tag(self, file_path):
        """添加文件标签"""
        if file_path not in self.imported_files:
            self.imported_files.append(file_path)
            
            # Create file tag widget (使用IQ-TREE plugin的灰底风格)
            tag_widget = QFrame()
            tag_widget.setFrameStyle(QFrame.Box)
            tag_widget.setStyleSheet("""
                QFrame {
                    background-color: #e9ecef;
                    border-radius: 15px;
                    margin: 2px;
                }
            """)
            
            tag_layout = QHBoxLayout()
            tag_layout.setContentsMargins(8, 4, 8, 4)
            tag_widget.setLayout(tag_layout)
            
            # 文件名标签
            name_label = QLabel(os.path.basename(file_path))
            name_label.setStyleSheet("color: #495057;")
            tag_layout.addWidget(name_label)
            
            # 关闭按钮
            close_btn = QPushButton("×")
            close_btn.setFixedSize(20, 20)
            close_btn.setStyleSheet("""
                QPushButton {
                    background-color: transparent;
                    border: none;
                    color: #6c757d;
                    font-weight: bold;
                    font-size: 20px;
                }
                QPushButton:hover {
                    color: #dc3545;
                }
            """)
            close_btn.clicked.connect(lambda: self.remove_file_tag(tag_widget, file_path))
            tag_layout.addWidget(close_btn)
            
            self.file_tags_layout.addWidget(tag_widget)
            self.file_tags.append((file_path, tag_widget))
            self.file_tags_container.setVisible(True)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def remove_file_tag(self, tag_frame, file_path):
        """移除文件标签"""
        # 从列表中移除
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        
        # 从UI中移除
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        self.file_tags_layout.removeWidget(tag_frame)
        tag_frame.deleteLater()
        
        # 更新显示
        if not self.imported_files:
            self.file_tags_container.setVisible(False)
            self.file_path_edit.clear()
        else:
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def on_text_changed(self):
        """文本输入变化处理"""
        if self.sequence_text.toPlainText().strip():
            # 如果有文本输入，清空文件选择
            self.imported_files.clear()
            for tag in self.file_tags:
                self.file_tags_layout.removeWidget(tag)
                tag.deleteLater()
            self.file_tags.clear()
            self.file_tags_container.setVisible(False)
        elif not self.imported_files:
            # 如果既没有文本也没有文件，恢复初始状态
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
    
    def reset_to_defaults(self):
        """重置参数为默认值"""
        self.gap_treatment_combo.setCurrentText("Pairwise Deletion")
    
    def show_method_description(self):
        """显示方法说明对话框"""
        dialog = MethodDescriptionDialog(self)
        dialog.exec_()
    
    def prepare_input_files(self):
        """准备输入文件"""
        try:
            input_files = []
            if self.imported_files:
                for file_path in self.imported_files:
                    if os.path.exists(file_path):
                        input_files.append(file_path)
                return input_files
            elif self.sequence_text.toPlainText().strip():
                # 从文本创建临时文件
                temp_file = self.create_temp_file(suffix='.fas')
                with open(temp_file, 'w') as f:
                    f.write(self.sequence_text.toPlainText())
                self.temp_files.append(temp_file)
                return [temp_file]
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return None
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
    
    def get_parameters(self):
        """获取参数设置"""
        params = {}
        
        # Gap处理方式
        gap_treatment = self.gap_treatment_combo.currentText()
        params['gap_treatment'] = 'pairwise' if 'Pairwise' in gap_treatment else 'complete'
        
        return params
    
    def run_analysis(self):
        """运行饱和度分析"""
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip():
            QMessageBox.warning(self, "Warning", "Please provide alignment files or sequence text!")
            return
            
        # 添加控制台消息
        self.add_console_message("Starting substitution saturation analysis...", "info")
        
        # 准备输入文件
        input_files = self.prepare_input_files()
        if not input_files:
            return
        
        # 获取参数
        params = self.get_parameters()
        
        # 设置运行状态
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 未知进度
        self.export_plot_btn.setEnabled(False)
        
        try:
            # 初始化结果存储
            self.saturation_results = []
            
            # 处理每个输入文件
            for i, input_file in enumerate(input_files):
                if not self.is_running:
                    break
                
                self.progress_bar.setFormat(f"Processing file {i+1}/{len(input_files)}...")
                self.add_console_message(f"Processing file {i+1}/{len(input_files)}: {os.path.basename(input_file)}", "info")
                
                # 读取序列
                sequences = list(SeqIO.parse(input_file, "fasta"))
                
                # 验证序列
                is_valid, msg = validate_sequences(sequences)
                if not is_valid:
                    self.add_console_message(f"Validation error: {msg}", "error")
                    continue
                
                # 计算饱和度指标
                results = calculate_saturation_metrics(
                    sequences,
                    gap_treatment=params['gap_treatment']
                )
                
                # 添加文件名到结果中
                results['filename'] = os.path.basename(input_file)
                results['filepath'] = input_file
                results['gap_treatment'] = params['gap_treatment']
                
                # 保存结果
                self.saturation_results.append(results)
                
                self.add_console_message(f"Analysis completed for {len(sequences)} sequences", "info")
            
            # 所有文件处理完成后，更新UI
            if self.saturation_results:
                # 更新文件选择Combobox
                self.file_combo.clear()
                for result in self.saturation_results:
                    self.file_combo.addItem(result['filename'])
                
                # 如果有多个文件，显示文件选择器
                if len(self.saturation_results) > 1:
                    self.file_combo_label.setVisible(True)
                    self.file_combo.setVisible(True)
                    
                    # 更新总体柱状图
                    file_names = [result['filename'] for result in self.saturation_results]
                    c_values = [result['c_value'] for result in self.saturation_results]
                    self.overview_canvas.plot_c_values(file_names, c_values)
                else:
                    self.file_combo_label.setVisible(False)
                    self.file_combo.setVisible(False)
                
                # 显示第一个文件的结果
                self.on_file_selection_changed(0)
                
                # 启用导出按钮
                self.export_plot_btn.setEnabled(True)
            
            # 分析完成
            self.analysis_finished([], [])
            
        except Exception as e:
            self.analysis_error(str(e))
    
    def display_single_result(self, results: Dict, filename: str):
        """显示单个文件的结果"""
        try:
            # 绘制图表
            self.canvas.plot_saturation(
                results['p_values'],
                results['ti_values'],
                results['tv_values'],
                results['c_value'],
                results['gap_treatment']
            )
            
            # 显示统计信息
            p_values = results['p_values']
            ti_values = results['ti_values']
            tv_values = results['tv_values']
            
            # 计算统计量
            p_mean = np.mean(p_values) if p_values else 0
            p_std = np.std(p_values) if p_values else 0
            ti_mean = np.mean(ti_values) if ti_values else 0
            tv_mean = np.mean(tv_values) if tv_values else 0
            
            # 计算相关系数
            correlation_text = "N/A"
            if len(p_values) > 1:
                ti_corr, ti_p = stats.pearsonr(p_values, ti_values)
                tv_corr, tv_p = stats.pearsonr(p_values, tv_values)
                correlation_text = f"Ti-p correlation: r={ti_corr:.3f} (p={ti_p:.4f})\nTv-p correlation: r={tv_corr:.3f} (p={tv_p:.4f})"
            
            # C值作为量化指标
            c_value = results['c_value']
            
            stats_report = f"""
Substitution Saturation Analysis Results
========================================

File: {filename}

C-value: {c_value:.4f}

Statistics:
- Number of sequence pairs: {results['n_pairs']}
- Mean p-distance: {p_mean:.4f} ± {p_std:.4f}
- Mean Transitions (Ti): {ti_mean:.2f}
- Mean Transversions (Tv): {tv_mean:.2f}
- p-distance range: [{min(p_values):.4f}, {max(p_values):.4f}]
- Ti range: [{min(ti_values):.2f}, {max(ti_values):.2f}]
- Tv range: [{min(tv_values):.2f}, {max(tv_values):.2f}]

Correlation Analysis:
{correlation_text}

Method: BaCoCa saturation index
Gap Treatment: {results['gap_treatment']}

Reference:
Kück, P., & Struck, T. H. (2014). BaCoCa--a heuristic software tool for the parallel 
assessment of sequence biases in hundreds of gene and taxon partitions. 
Molecular Phylogenetics and Evolution, 70, 94–98. 
https://doi.org/10.1016/j.ympev.2013.09.011
            """
            
            self.stats_text.setPlainText(stats_report.strip())
            
        except Exception as e:
            self.add_console_message(f"Error displaying results: {str(e)}", "error")
    
    def display_results(self, results: Dict, gap_treatment: str):
        """显示分析结果（兼容旧接口）"""
        # 添加文件名如果不存在
        if 'filename' not in results:
            results['filename'] = 'Unknown'
        if 'gap_treatment' not in results:
            results['gap_treatment'] = gap_treatment
        
        # 调用新方法
        self.display_single_result(results, results['filename'])
    
    def analysis_finished(self, output_files, reports):
        """分析完成回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.progress_bar.setFormat("")
        
        # 显示完成消息
        QMessageBox.information(self, "Completed", "Substitution saturation analysis completed successfully!")
    
    def analysis_error(self, error_msg):
        """分析错误回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.progress_bar.setFormat("")
        self.export_plot_btn.setEnabled(False)
        
        QMessageBox.critical(self, "Error", f"Analysis failed: {error_msg}")
    
    def export_plot(self):
        """导出图表"""
        if not self.saturation_results:
            QMessageBox.warning(self, "Warning", "No results to export!")
            return
        
        file_filter = "PNG files (*.png);;PDF files (*.pdf);;SVG files (*.svg)"
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Plot", "saturation_plot.png", file_filter
        )
        
        if file_path:
            try:
                self.canvas.fig.savefig(file_path, dpi=300, bbox_inches='tight')
                QMessageBox.information(self, "Success", f"Plot exported to {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export plot: {e}")
    
    def handle_import_data(self, import_data):
        """处理从YR-MPEA导入的数据"""
        if isinstance(import_data, dict) and 'type' in import_data and import_data['type'] == 'alignment':
            # 处理对齐序列数据
            sequences = import_data['data']
            if isinstance(sequences, list) and len(sequences) > 0:
                # 创建临时文件来存储导入的序列数据
                temp_file = self.create_temp_file(suffix='.fas')
                with open(temp_file, 'w') as f:
                    for seq in sequences:
                        f.write(f">{seq.id}\n{seq.seq}\n")
                self.temp_files.append(temp_file)
                self.import_file = temp_file
                self.imported_files = [temp_file]
                
                # 更新UI显示导入的文件
                if hasattr(self, 'file_path_edit') and self.file_path_edit:
                    self.file_path_edit.setText(temp_file)
                    
                # 添加控制台消息
                self.add_console_message(f"Imported alignment data with {len(sequences)} sequences from YR-MPEA", "info")
        elif isinstance(import_data, list):
            # 兼容旧版本的导入数据格式
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                for seq in import_data:
                    f.write(f">{seq.id}\n{seq.seq}\n")
            self.temp_files.append(temp_file)
            self.import_file = temp_file
            self.imported_files = [temp_file]
            
            # 更新UI显示导入的文件
            if hasattr(self, 'file_path_edit') and self.file_path_edit:
                self.file_path_edit.setText(temp_file)
                
            # 添加控制台消息
            self.add_console_message(f"Imported alignment data with {len(import_data)} sequences from YR-MPEA", "info")
        else:
            self.import_file = None
            self.imported_files = []
            self.add_console_message("Invalid import data format", "error")
            
        # 确保 saturation_results 属性存在
        if not hasattr(self, 'saturation_results'):
            self.saturation_results = None


# 插件入口点
class SubstitutionSaturationPluginEntry:
    """替代饱和度分析插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return SubstitutionSaturationPlugin(import_from=import_from, import_data=import_data)
