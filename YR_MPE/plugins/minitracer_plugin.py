import os
import pandas as pd
import numpy as np
from PyQt5.QtWidgets import (QWidget, QSplitter, QVBoxLayout, QHBoxLayout, 
                             QTableWidget, QTableWidgetItem, QPushButton, 
                             QTabWidget, QFileDialog, QMessageBox, QHeaderView,
                             QSizePolicy)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QColor
from .components.traceplot import TracePlot
from .components.tracecomp import TraceComp
from .components.summary_table import SummaryTable
from .components.methods.mcmc_utils import calculate_ESS


class MiniTracerPlugin(QWidget):
    """
    MiniTracer Plugin for comprehensive MCMC trace diagnostics.
    Provides a unified interface for loading, analyzing, and comparing MCMC chains.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.mcmc_files = []  # List of loaded MCMC file paths
        self.mcmc_data = {}   # Dictionary to store parsed data for each file
        self.trace_plots = {}  # Dictionary to store pre-initialized TracePlot components: {file_path: {param_name: TracePlot}}
        self.current_trace_plot = None  # Currently displayed TracePlot component
        self.init_ui()
        
    def init_ui(self):
        """Initialize the user interface with splitter layout."""
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Create main splitter with 3:7 ratio (left:right)
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # Left panel: File management and statistics
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        left_panel.setLayout(left_layout)
        
        # Create splitter for file table/buttons and stats table with 3:7 ratio
        left_splitter = QSplitter(Qt.Vertical)
        left_layout.addWidget(left_splitter)
        
        # Upper section: File table and buttons
        upper_section = QWidget()
        upper_layout = QVBoxLayout()
        upper_layout.setContentsMargins(0, 0, 0, 0)
        upper_section.setLayout(upper_layout)
        
        # File table widget
        self.file_table = QTableWidget()
        self.file_table.setStyleSheet("QTableWidgetItem { padding: 0px;}")
        self.file_table.setColumnCount(3)
        self.file_table.setHorizontalHeaderLabels(['File', 'Generations', 'Burn-in'])
        self.file_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.file_table.setSelectionMode(QTableWidget.MultiSelection)
        
        # Apply compact styling similar to TraceComp
        font = QFont()
        font.setPointSize(9)
        self.file_table.setFont(font)
        self.file_table.horizontalHeader().setFont(font)
        self.file_table.verticalHeader().setFont(font)
        self.file_table.verticalHeader().setVisible(False)
        # Allow editing only for burn-in column (column 2)
        self.file_table.setEditTriggers(QTableWidget.DoubleClicked | QTableWidget.SelectedClicked)
        
        # Set column widths
        self.file_table.setColumnWidth(1, 120)
        self.file_table.setColumnWidth(2, 120)
        
        self.file_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.file_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.file_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        
        upper_layout.addWidget(self.file_table)
        
        # Button layout
        button_layout = QHBoxLayout()
        self.add_button = QPushButton('+')
        self.remove_button = QPushButton('-')
        self.refresh_button = QPushButton('Refresh')
        
        button_layout.addWidget(self.add_button)
        button_layout.addWidget(self.remove_button)
        button_layout.addStretch()
        button_layout.addWidget(self.refresh_button)
        
        upper_layout.addLayout(button_layout)
        
        # Lower section: Statistics table
        lower_section = QWidget()
        lower_layout = QVBoxLayout()
        lower_layout.setContentsMargins(0, 0, 0, 0)
        lower_section.setLayout(lower_layout)
        
        # Statistics table widget
        self.stats_table = QTableWidget()
        self.stats_table.setStyleSheet("QTableWidgetItem { padding: 0px;}")
        self.stats_table.setColumnCount(3)
        self.stats_table.setHorizontalHeaderLabels(['Parameter', 'Mean', 'ESS'])
        self.stats_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.stats_table.setSelectionMode(QTableWidget.SingleSelection)
        
        # Apply compact styling similar to TraceComp
        self.stats_table.setFont(font)
        self.stats_table.horizontalHeader().setFont(font)
        self.stats_table.verticalHeader().setFont(font)
        self.stats_table.verticalHeader().setVisible(False)
        self.stats_table.setEditTriggers(QTableWidget.NoEditTriggers)
        
        # Set column widths
        self.stats_table.setColumnWidth(1, 120)
        self.stats_table.setColumnWidth(2, 120)
        
        self.stats_table.setStyleSheet("QTableWidget { gridline-color: lightgray; }")
        
        self.stats_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.stats_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.stats_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        
        lower_layout.addWidget(self.stats_table)
        
        # Add sections to left splitter with 3:7 ratio
        left_splitter.addWidget(upper_section)
        left_splitter.addWidget(lower_section)
        left_splitter.setSizes([300, 700])  # 3:7 ratio
        
        # Right panel: Analysis tabs
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_panel.setLayout(right_layout)
        
        # Create the tab widget for analysis
        self.tab_widget = QTabWidget()
        # self.tab_widget.setStyleSheet("background-color: blue;")
        # 创建分析用的标签页
        # Create tabs
        self.single_chain_tab = QWidget()
        self.compare_chains_tab = QWidget()
        
        self.tab_widget.addTab(self.single_chain_tab, "Analyze Single Chain")
        self.tab_widget.addTab(self.compare_chains_tab, "Compare 2 Chains")
        
        # Initialize single chain analysis tab with a container widget
        single_layout = QVBoxLayout()
        self.single_chain_tab.setLayout(single_layout)
        
        # Create summary table for statistics
        self.summary_table = SummaryTable()
        single_layout.addWidget(self.summary_table, 1)  # Add with stretch factor
        
        # Create trace plot container
        self.trace_plot_container = QWidget()
        self.trace_plot_layout = QVBoxLayout()
        self.trace_plot_container.setLayout(self.trace_plot_layout)
        # 设置容器的大小策略，确保能正确拉伸
        self.trace_plot_container.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # 设置布局的拉伸因子
        self.trace_plot_layout.setStretch(0, 1)
        single_layout.addWidget(self.trace_plot_container, 3)  # Add with larger stretch factor
        
        # Initialize compare chains analysis tab
        compare_layout = QVBoxLayout()
        self.compare_chains_tab.setLayout(compare_layout)
        self.trace_comp = TraceComp(burnin_fraction=0.1)
        compare_layout.addWidget(self.trace_comp)
        
        # Connect tab change signal to update plot when switching tabs
        self.tab_widget.currentChanged.connect(self.on_tab_changed)
        
        # Add the tab widget to the right panel layout
        right_layout.addWidget(self.tab_widget)
        
        # Connect button signals
        self.add_button.clicked.connect(self.add_files)
        self.remove_button.clicked.connect(self.remove_selected_files)
        self.refresh_button.clicked.connect(self.refresh_selected_files)
        
        # Connect file table selection change
        self.file_table.itemSelectionChanged.connect(self.on_file_selection_changed)
        
        # Connect file table item changes for burn-in editing
        self.file_table.itemChanged.connect(self.on_file_table_item_changed)
        
        # Connect stats table selection change
        self.stats_table.itemSelectionChanged.connect(self.on_stats_selection_changed)
        
        # Set initial splitter sizes (3:7 ratio)
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([300, 700])  # Approximate 3:7 ratio
        
    def update_file_table(self):
        """Update the file table with current MCMC files and their information."""
        self.file_table.setRowCount(len(self.mcmc_files))
        
        for row, file_path in enumerate(self.mcmc_files):
            # File name column (non-editable)
            file_name = os.path.basename(file_path)
            file_item = QTableWidgetItem(file_name)
            file_item.setFlags(file_item.flags() & ~Qt.ItemIsEditable)
            self.file_table.setItem(row, 0, file_item)
            
            # Generations column (non-editable)
            generations = 0
            if file_path in self.mcmc_data:
                data = self.mcmc_data[file_path]
                if hasattr(data, 'shape') and len(data.shape) > 0:
                    generations = data.shape[0]
            
            gen_item = QTableWidgetItem(str(generations))
            gen_item.setFlags(gen_item.flags() & ~Qt.ItemIsEditable)
            self.file_table.setItem(row, 1, gen_item)
            
            # Burn-in fraction (editable) - default to 0.1
            burnin_item = QTableWidgetItem("0.100")
            # Make sure burn-in column is editable
            burnin_item.setFlags(burnin_item.flags() | Qt.ItemIsEditable)
            self.file_table.setItem(row, 2, burnin_item)
            
            # Set row height for compact layout (mimic TraceComp behavior)
            self.file_table.setRowHeight(row, 25)

    def on_file_table_item_changed(self, item):
        """Handle changes to file table items, specifically burn-in values."""
        if item.column() == 2:  # Burn-in column
            row = item.row()
            if row < len(self.mcmc_files):
                try:
                    # Parse and validate the burn-in value
                    burnin_value = float(item.text())
                    # Clamp between 0.0 and 0.999
                    burnin_value = max(0.0, min(0.999, burnin_value))
                    # Update the display with formatted value
                    item.setText(f"{burnin_value:.3f}")
                    # Trigger stats table update to reflect new burn-in value
                    self.update_stats_table()
                except ValueError:
                    # If invalid input, reset to default 0.100
                    item.setText("0.100")
                    self.update_stats_table()
    
    def update_stats_table(self):
        """Update statistics table based on selected trace file(s) and current burn-in fraction."""
        selected_items = self.file_table.selectedItems()
        if not selected_items:
            # Clear stats table
            self.stats_table.setRowCount(0)
            return
        
        selected_rows = sorted(set(item.row() for item in selected_items))
        
        # Handle single file selection (original behavior)
        if len(selected_rows) == 1:
            row = selected_rows[0]
            if row >= len(self.mcmc_files):
                self.stats_table.setRowCount(0)
                return
                
            file_path = self.mcmc_files[row]
            if file_path not in self.mcmc_data:
                self.stats_table.setRowCount(0)
                return
                
            data = self.mcmc_data[file_path]
            if data is None or data.empty:
                self.stats_table.setRowCount(0)
                return
            
            # Get burn-in fraction for this file
            burnin_frac = 0.1
            burnin_item = self.file_table.item(row, 2)
            if burnin_item:
                try:
                    burnin_frac = float(burnin_item.text())
                    burnin_frac = max(0.0, min(0.999, burnin_frac))
                except ValueError:
                    burnin_frac = 0.1
            
            # Apply burn-in to get post-burn-in data
            total_samples = len(data)
            burnin_samples = int(total_samples * burnin_frac)
            if burnin_samples >= total_samples:
                burnin_samples = total_samples - 1
            filtered_data = data.iloc[burnin_samples:]
            
            if filtered_data.empty:
                self.stats_table.setRowCount(0)
                return
            
            # Get all parameters (excluding 'iterations' column)
            params = [col for col in filtered_data.columns if col != 'iterations']
            self.stats_table.setRowCount(len(params))
            
            for i, param in enumerate(params):
                # Parameter name (non-editable)
                param_item = QTableWidgetItem(param)
                param_item.setFlags(param_item.flags() & ~Qt.ItemIsEditable)
                self.stats_table.setItem(i, 0, param_item)
                
                # Mean value calculated from post-burn-in data (non-editable)
                try:
                    numeric_data = pd.to_numeric(filtered_data[param], errors='coerce')
                    numeric_data = numeric_data.dropna()
                    if len(numeric_data) > 0:
                        mean_val = np.mean(numeric_data)
                        mean_str = f"{mean_val:.6f}"
                    else:
                        mean_str = "N/A"
                except Exception:
                    mean_str = "N/A"
                    
                mean_item = QTableWidgetItem(mean_str)
                mean_item.setFlags(mean_item.flags() & ~Qt.ItemIsEditable)
                self.stats_table.setItem(i, 1, mean_item)
                
                # Calculate ESS value
                try:
                    if len(numeric_data) > 0:
                        ess_value = calculate_ESS(numeric_data.values)
                        ess_str = f"{ess_value:.1f}"
                    else:
                        ess_str = "N/A"
                        ess_value = None
                except Exception:
                    ess_str = "N/A"
                    ess_value = None
                    
                ess_item = QTableWidgetItem(ess_str)
                ess_item.setFlags(ess_item.flags() & ~Qt.ItemIsEditable)
                self.stats_table.setItem(i, 2, ess_item)
                
                # Apply color coding based on ESS value
                if ess_value is not None:
                    if ess_value >= 200:
                        # Green: #90ee90
                        bg_color = QColor("#90ee90")
                    elif ess_value >= 100:
                        # Yellow: #eede91
                        bg_color = QColor("#eede91")
                    else:
                        # Red: #ee9191
                        bg_color = QColor("#ee9191")
                        
                    # Apply background color to all cells in the row
                    param_item.setBackground(bg_color)
                    mean_item.setBackground(bg_color)
                    ess_item.setBackground(bg_color)
                
                # Set row height for compact layout (mimic TraceComp behavior)
                self.stats_table.setRowHeight(i, 25)
                
        # Handle two file selection (show common parameters)
        elif len(selected_rows) == 2:
            row1, row2 = selected_rows
            if row1 >= len(self.mcmc_files) or row2 >= len(self.mcmc_files):
                self.stats_table.setRowCount(0)
                return
                
            file_path1 = self.mcmc_files[row1]
            file_path2 = self.mcmc_files[row2]
            
            if file_path1 not in self.mcmc_data or file_path2 not in self.mcmc_data:
                self.stats_table.setRowCount(0)
                return
                
            data1 = self.mcmc_data[file_path1]
            data2 = self.mcmc_data[file_path2]
            
            if data1 is None or data1.empty or data2 is None or data2.empty:
                self.stats_table.setRowCount(0)
                return
            
            # Get burn-in fractions for both files
            burnin_frac1 = 0.1
            burnin_item1 = self.file_table.item(row1, 2)
            if burnin_item1:
                try:
                    burnin_frac1 = float(burnin_item1.text())
                    burnin_frac1 = max(0.0, min(0.999, burnin_frac1))
                except ValueError:
                    burnin_frac1 = 0.1
                    
            burnin_frac2 = 0.1
            burnin_item2 = self.file_table.item(row2, 2)
            if burnin_item2:
                try:
                    burnin_frac2 = float(burnin_item2.text())
                    burnin_frac2 = max(0.0, min(0.999, burnin_frac2))
                except ValueError:
                    burnin_frac2 = 0.1
            
            # Apply burn-in to get post-burn-in data for both files
            total_samples1 = len(data1)
            burnin_samples1 = int(total_samples1 * burnin_frac1)
            if burnin_samples1 >= total_samples1:
                burnin_samples1 = total_samples1 - 1
            filtered_data1 = data1.iloc[burnin_samples1:]
            
            total_samples2 = len(data2)
            burnin_samples2 = int(total_samples2 * burnin_frac2)
            if burnin_samples2 >= total_samples2:
                burnin_samples2 = total_samples2 - 1
            filtered_data2 = data2.iloc[burnin_samples2:]
            
            if filtered_data1.empty or filtered_data2.empty:
                self.stats_table.setRowCount(0)
                return
            
            # Get common parameters (excluding 'iterations' column)
            params1 = set(col for col in filtered_data1.columns if col != 'iterations')
            params2 = set(col for col in filtered_data2.columns if col != 'iterations')
            common_params = sorted(list(params1 & params2))
            
            if not common_params:
                self.stats_table.setRowCount(0)
                return
                
            self.stats_table.setRowCount(len(common_params))
            
            for i, param in enumerate(common_params):
                # Parameter name (non-editable)
                param_item = QTableWidgetItem(param)
                param_item.setFlags(param_item.flags() & ~Qt.ItemIsEditable)
                self.stats_table.setItem(i, 0, param_item)
                
                # Mean values for both files
                try:
                    # File 1 mean
                    numeric_data1 = pd.to_numeric(filtered_data1[param], errors='coerce')
                    numeric_data1 = numeric_data1.dropna()
                    if len(numeric_data1) > 0:
                        mean_val1 = np.mean(numeric_data1)
                        mean_str1 = f"{mean_val1:.6f}"
                    else:
                        mean_str1 = "N/A"
                except Exception:
                    mean_str1 = "N/A"
                    
                try:
                    # File 2 mean  
                    numeric_data2 = pd.to_numeric(filtered_data2[param], errors='coerce')
                    numeric_data2 = numeric_data2.dropna()
                    if len(numeric_data2) > 0:
                        mean_val2 = np.mean(numeric_data2)
                        mean_str2 = f"{mean_val2:.6f}"
                    else:
                        mean_str2 = "N/A"
                except Exception:
                    mean_str2 = "N/A"
                
                # Combine means for display
                if mean_str1 != "N/A" and mean_str2 != "N/A":
                    mean_str = f"{mean_str1} / {mean_str2}"
                else:
                    mean_str = "N/A"
                    
                mean_item = QTableWidgetItem(mean_str)
                mean_item.setFlags(mean_item.flags() & ~Qt.ItemIsEditable)
                self.stats_table.setItem(i, 1, mean_item)
                
                # Calculate ESS values for both files
                ess_value1 = None
                ess_value2 = None
                
                try:
                    if len(numeric_data1) > 0:
                        ess_value1 = calculate_ESS(numeric_data1.values)
                        ess_str1 = f"{ess_value1:.1f}"
                    else:
                        ess_str1 = "N/A"
                except Exception:
                    ess_str1 = "N/A"
                    
                try:
                    if len(numeric_data2) > 0:
                        ess_value2 = calculate_ESS(numeric_data2.values)
                        ess_str2 = f"{ess_value2:.1f}"
                    else:
                        ess_str2 = "N/A"
                except Exception:
                    ess_str2 = "N/A"
                
                # Combine ESS values for display
                if ess_str1 != "N/A" and ess_str2 != "N/A":
                    ess_str = f"{ess_str1} / {ess_str2}"
                else:
                    ess_str = "N/A"
                    
                ess_item = QTableWidgetItem(ess_str)
                ess_item.setFlags(ess_item.flags() & ~Qt.ItemIsEditable)
                self.stats_table.setItem(i, 2, ess_item)
                
                # Apply color coding based on minimum ESS value
                if ess_value1 is not None and ess_value2 is not None:
                    min_ess = min(ess_value1, ess_value2)
                    if min_ess >= 200:
                        # Green: #90ee90
                        bg_color = QColor("#90ee90")
                    elif min_ess >= 100:
                        # Yellow: #eede91
                        bg_color = QColor("#eede91")
                    else:
                        # Red: #ee9191
                        bg_color = QColor("#ee9191")
                        
                    # Apply background color to all cells in the row
                    param_item.setBackground(bg_color)
                    mean_item.setBackground(bg_color)
                    ess_item.setBackground(bg_color)
                
                # Set row height for compact layout (mimic TraceComp behavior)
                self.stats_table.setRowHeight(i, 25)
        else:
            # More than 2 files selected - clear stats table
            self.stats_table.setRowCount(0)
    
    def add_files(self):
        """Add MCMC trace files."""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select MCMC Trace Files", "", "Trace Files (*.log *.txt *.csv);;All Files (*)"
        )
        
        if not file_paths:
            return
            
        for file_path in file_paths:
            if file_path not in self.mcmc_files:
                self.mcmc_files.append(file_path)
                try:
                    data = self.parse_trace_file(file_path)
                    self.mcmc_data[file_path] = data
                    self.update_file_table()
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"Failed to parse {os.path.basename(file_path)}: {str(e)}")
                    if file_path in self.mcmc_files:
                        self.mcmc_files.remove(file_path)
                        
    def remove_selected_files(self):
        """Remove selected files from the list."""
        selected_rows = sorted(set(item.row() for item in self.file_table.selectedItems()), reverse=True)
        for row in selected_rows:
            file_path = self.mcmc_files[row]
            del self.mcmc_data[file_path]
            del self.mcmc_files[row]
            
        self.update_file_table()
        self.update_stats_table()
        self.update_analysis_tabs()
        
    def refresh_selected_files(self):
        """Refresh selected files."""
        selected_rows = set(item.row() for item in self.file_table.selectedItems())
        if not selected_rows:
            # If no selection, refresh all files
            selected_rows = range(len(self.mcmc_files))
            
        for row in selected_rows:
            file_path = self.mcmc_files[row]
            try:
                data = self.parse_trace_file(file_path)
                self.mcmc_data[file_path] = data
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Failed to refresh {os.path.basename(file_path)}: {str(e)}")
                
        self.update_file_table()
        self.update_stats_table()
        self.update_analysis_tabs()
        
    def parse_trace_file(self, trace_file):
        """Parse trace file and return pandas DataFrame with proper numeric conversion."""
        try:
            # 读取所有行
            with open(trace_file, 'r') as f:
                lines = f.readlines()
            
            # 找到包含"iter"的header行和包含"0"的数据起始行
            header_line = None
            data_start_line = None
            
            last_col, st_num = 0, 0
            for index, line in enumerate(lines):
                line_parts = line.strip().split('\t')
                if len(line_parts) > 0 and line_parts[0] == '0':
                    st_num = index
                    break
                last_col = len(line_parts)
                
            if last_col <= 1:
                raise Exception('Trace file is empty')
            
            # 使用pandas读取，指定header行和跳过前面的行
            df = pd.read_csv(trace_file, sep='\t', header=0, skiprows=st_num-1)
            
            # 确保DataFrame不为空
            if df.empty:
                return df
                
            # 将除第一列（迭代数）外的所有列转换为数值类型
            # 第一列始终是迭代数，不管列名是什么
            numeric_columns = []
            for col in df.columns[1:]:  # 跳过第一列
                df[col] = pd.to_numeric(df[col], errors='coerce')
                numeric_columns.append(col)
            
            # 重命名第一列为'iterations'以便统一处理
            original_first_col = df.columns[0]
            df = df.rename(columns={original_first_col: 'iterations'})
            
            return df
            
        except Exception:
            return pd.DataFrame()
    
    def get_trace_info(self, trace_file):
        """Get trace file information including generations and parameters"""
        try:
            df = self.parse_trace_file(trace_file)
            if df.empty:
                return {
                    'generations': 0,
                    'burnin_fraction': 0.1,
                    'parameters': [],
                    'data': pd.DataFrame()
                }
            
            # Get generations (number of rows)
            generations = len(df)
            
            # Get parameter names (all columns except the first one which is usually iteration)
            parameters = list(df.columns)[1:] if len(df.columns) > 1 else []
            
            return {
                'generations': generations,
                'burnin_fraction': 0.1,
                'parameters': parameters,
                'data': df
            }
            
        except Exception:
            return {
                'generations': 0,
                'burnin_fraction': 0.1,
                'parameters': [],
                'data': pd.DataFrame()
            }
    
    def on_burnin_fraction_changed(self, item):
        """Handle burnin fraction changes."""
        if item.column() == 2:  # Burnin fraction column
            try:
                value = float(item.text())
                if value < 0 or value >= 1.0:
                    raise ValueError("Burnin fraction must be in [0, 1.0)")
                item.setData(Qt.UserRole, value)
            except (ValueError, TypeError):
                # Revert to previous valid value or default
                prev_value = item.data(Qt.UserRole)
                if prev_value is None:
                    prev_value = 0.1
                item.setText(str(prev_value))
                item.setData(Qt.UserRole, prev_value)
                
            # Update the corresponding trace plots and stats table if this file is selected
            row = item.row()
            selected_rows = set(item.row() for item in self.file_table.selectedItems())
            if row in selected_rows and len(selected_rows) == 1:
                # Update stats table first to get new burn-in fraction
                self.update_stats_table()
                # Then update the currently displayed trace plot if any
                stats_selected = set(item.row() for item in self.stats_table.selectedItems())
                if stats_selected:
                    stats_row = next(iter(stats_selected))
                    self.show_trace_plot_for_selection(row, stats_row)
                
    def on_file_selection_changed(self):
        """Handle file table selection changes."""
        self.update_stats_table()
        self.update_analysis_tabs()
        
        # If a file is selected and there's a parameter selected, show the corresponding trace plot
        file_selected_rows = set(item.row() for item in self.file_table.selectedItems())
        stats_selected_rows = set(item.row() for item in self.stats_table.selectedItems())
        if len(file_selected_rows) == 1 and len(stats_selected_rows) == 1:
            file_row = next(iter(file_selected_rows))
            stats_row = next(iter(stats_selected_rows))
            self.show_trace_plot_for_selection(file_row, stats_row)
        
    def on_stats_selection_changed(self):
        """Handle statistics table selection changes."""
        selected_rows = set(item.row() for item in self.stats_table.selectedItems())
        file_selected_rows = set(item.row() for item in self.file_table.selectedItems())
        
        if len(selected_rows) == 1:
            if len(file_selected_rows) == 1:
                # Single file selection - show trace plot
                file_row = next(iter(file_selected_rows))
                stats_row = next(iter(selected_rows))
                self.show_trace_plot_for_selection(file_row, stats_row)
            elif len(file_selected_rows) == 2:
                # Two files selected - update chain comparison
                self.update_analysis_tabs()
    
    def show_trace_plot_for_selection(self, file_row, param_row):
        """Show the TracePlot component for the selected file and parameter, creating it if necessary."""
        if file_row >= len(self.mcmc_files):
            return
            
        file_path = self.mcmc_files[file_row]
        if file_path not in self.mcmc_data:
            return
            
        data = self.mcmc_data[file_path]
        if data is None or data.empty:
            return
            
        # Get parameter name from stats table
        param_item = self.stats_table.item(param_row, 0)
        if not param_item:
            return
            
        param_name = param_item.text()
        
        # Initialize trace_plots dictionary if needed
        if file_path not in self.trace_plots:
            self.trace_plots[file_path] = {}
            
        # Get current burn-in fraction
        burnin_frac = 0.1
        burnin_item = self.file_table.item(file_row, 2)
        if burnin_item:
            try:
                burnin_frac = float(burnin_item.text())
                burnin_frac = max(0.0, min(0.999, burnin_frac))
            except ValueError:
                burnin_frac = 0.1
        
        # DO NOT apply burn-in here - pass complete original data to TracePlot
        # TracePlot will handle burn-in display logic internally
        
        if not data.empty:
            try:
                # Get parameter values and iterations from complete data
                param_series = data[param_name]
                numeric_mask = pd.to_numeric(param_series, errors='coerce').notna()
                valid_data = data[numeric_mask]
                
                if len(valid_data) > 0:
                    param_values = np.asarray(valid_data[param_name].values, dtype=float)
                    iterations = np.asarray(valid_data['iterations'].values, dtype=int)
                    
                    # Create TracePlot component if it doesn't exist
                    if param_name not in self.trace_plots[file_path]:
                        self.trace_plots[file_path][param_name] = TracePlot(burnin_fraction=burnin_frac)
                    
                    # ALWAYS update the data with complete original data
                    self.trace_plots[file_path][param_name].set_data(param_values, iterations)
                    self.trace_plots[file_path][param_name].update_burnin_fraction(burnin_frac)
                    
                    # Update SummaryTable with statistics for the selected parameter
                    if hasattr(self, 'summary_table'):
                        # Apply burn-in to get the actual displayed samples
                        burnin_samples = int(len(param_values) * burnin_frac)
                        displayed_samples = param_values[burnin_samples:]
                        self.summary_table.set_data(displayed_samples)
                    
            except Exception:
                return
        
        # Remove current trace plot if exists
        if self.current_trace_plot:
            self.trace_plot_layout.removeWidget(self.current_trace_plot)
            self.current_trace_plot.hide()
            
        # Show the selected trace plot
        if param_name in self.trace_plots[file_path]:
            self.current_trace_plot = self.trace_plots[file_path][param_name]
            self.trace_plot_layout.addWidget(self.current_trace_plot)
            self.current_trace_plot.show()
            # Force the container to update and resize
            self.trace_plot_container.update()
            self.trace_plot_container.adjustSize()
            
            # Ensure the trace plot stretches to fill the container
            self.current_trace_plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            # 设置布局的拉伸因子
            self.trace_plot_layout.setStretch(0, 1)
            
    def update_single_chain_analysis(self, file_row, param_row=None):
        """Update single chain analysis tab - now handled by show_trace_plot_for_selection."""
        # This method is now primarily used for compatibility, but the actual display logic
        # is handled by show_trace_plot_for_selection
        pass
        
    def update_compare_chains_analysis(self, file_row1, file_row2):
        """Update the comparison analysis for two selected files."""
        if not hasattr(self, 'trace_comp') or self.trace_comp is None:
            return
            
        # Get the selected files - use mcmc_files list instead of file table text
        # because mcmc_data dictionary uses full file paths as keys
        if file_row1 >= len(self.mcmc_files) or file_row2 >= len(self.mcmc_files):
            print(f"Debug: Invalid file row indices {file_row1}, {file_row2}")
            self.trace_comp.set_data(None, None, None, None)
            return
            
        file1 = self.mcmc_files[file_row1]
        file2 = self.mcmc_files[file_row2]
        
        data1 = self.mcmc_data.get(file1)
        data2 = self.mcmc_data.get(file2)
        
        if data1 is None or data2 is None or data1.empty or data2.empty:
            print(f"Debug: Data is None or empty for files {file1} or {file2}")
            self.trace_comp.set_data(None, None, None, None)
            return
            
        # Get burnin fractions
        burnin1 = 0.1
        burnin2 = 0.1
        
        burnin_item1 = self.file_table.item(file_row1, 2)
        if burnin_item1:
            try:
                burnin1 = float(burnin_item1.text())
                burnin1 = max(0.0, min(0.999, burnin1))
            except ValueError:
                burnin1 = 0.1
                
        burnin_item2 = self.file_table.item(file_row2, 2)
        if burnin_item2:
            try:
                burnin2 = float(burnin_item2.text())
                burnin2 = max(0.0, min(0.999, burnin2))
            except ValueError:
                burnin2 = 0.1
                
        # Apply burnin separately for each chain (don't require them to be equal)
        burnin_samples1 = int(len(data1) * burnin1)
        burnin_samples2 = int(len(data2) * burnin2)
        filtered_data1 = data1.iloc[burnin_samples1:]
        filtered_data2 = data2.iloc[burnin_samples2:]
        
        if filtered_data1.empty or filtered_data2.empty:
            print(f"Debug: Filtered data is empty after burnin")
            self.trace_comp.set_data(None, None, None, None)
            return
            
        # Get common parameters (excluding 'iterations')
        params1 = set(col for col in filtered_data1.columns if col != 'iterations')
        params2 = set(col for col in filtered_data2.columns if col != 'iterations')
        common_params = list(params1 & params2)
        
        if not common_params:
            print(f"Debug: No common parameters found between files")
            self.trace_comp.set_data(None, None, None, None)
            return
            
        # Check if a parameter is selected in stats table
        stats_selected_rows = set(item.row() for item in self.stats_table.selectedItems())
        param_name = None
        
        if len(stats_selected_rows) == 1:
            stats_row = next(iter(stats_selected_rows))
            param_item = self.stats_table.item(stats_row, 0)
            if param_item:
                selected_param = param_item.text()
                # Check if the selected parameter is common to both files
                if selected_param in common_params:
                    param_name = selected_param
                    print(f"Debug: Using selected parameter {param_name} for comparison")
        
        # If no valid parameter selected from stats table, use the first common parameter
        if param_name is None:
            param_name = common_params[0]
            print(f"Debug: Using first common parameter {param_name} for comparison")
        
        # Prepare data for tracecomp with proper numeric conversion
        try:
            # Convert parameter data to numeric, coercing errors to NaN
            param_data1 = pd.to_numeric(filtered_data1[param_name], errors='coerce')
            param_data2 = pd.to_numeric(filtered_data2[param_name], errors='coerce')
            
            # Drop NaN values and get valid data
            valid_mask1 = ~param_data1.isna()
            valid_mask2 = ~param_data2.isna()
            
            # Get valid data points
            chain1_data_valid = param_data1[valid_mask1].values
            chain2_data_valid = param_data2[valid_mask2].values
            chain1_iterations_valid = filtered_data1[valid_mask1]['iterations'].values
            chain2_iterations_valid = filtered_data2[valid_mask2]['iterations'].values
            
            # Ensure we have data to plot
            if len(chain1_data_valid) == 0 or len(chain2_data_valid) == 0:
                print(f"Debug: No valid numeric data after filtering")
                self.trace_comp.set_data(None, None, None, None)
                return
                
            print(f"Debug: Setting data with {len(chain1_data_valid)} and {len(chain2_data_valid)} points")
            self.trace_comp.set_data(chain1_data_valid, chain1_iterations_valid, chain2_data_valid, chain2_iterations_valid)
            
        except Exception as e:
            # Handle any errors during data preparation
            print(f"Error preparing comparison data: {e}")
            self.trace_comp.set_data(None, None, None, None)
            return
        
    def update_analysis_tabs(self):
        """Update analysis tabs based on current selection."""
        selected_rows = set(item.row() for item in self.file_table.selectedItems())
        current_tab = self.tab_widget.currentIndex()
        
        if current_tab == 0:  # Single chain tab
            # Single chain analysis is handled by show_trace_plot_for_selection
            pass
        elif current_tab == 1:  # Compare chains tab
            if len(selected_rows) == 2:
                rows = sorted(list(selected_rows))
                self.update_compare_chains_analysis(rows[0], rows[1])
            else:
                # Clear the comparison display when not exactly 2 files are selected
                self.trace_comp.set_data(None, None, None, None)
        
    def get_tool_name(self):
        """Return the tool name for plugin registration."""
        return "MiniTracer"
        
    def get_input_format(self):
        """Return supported input formats."""
        return ["*.log", "*.txt", "*.csv"]
        
    def get_output_format(self):
        """Return output format."""
        return []

    def on_tab_changed(self, index):
        """Handle tab change events to ensure proper layout updates."""
        if index == 0:  # Analyze Single Chain tab
            if hasattr(self, 'current_trace_plot') and self.current_trace_plot:
                # Force update the trace plot layout
                self.current_trace_plot.update_plot()
        elif index == 1:  # Compare 2 Chains tab  
            if hasattr(self, 'trace_comp') and self.trace_comp:
                # Force update the trace comparison layout
                self.trace_comp.update_plot_and_table()

class MiniTracerPluginEntry:
    def run(self):
        return MiniTracerPlugin()
