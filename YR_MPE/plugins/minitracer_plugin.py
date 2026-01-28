import os
import pandas as pd
import numpy as np
from PyQt5.QtWidgets import (QWidget, QSplitter, QVBoxLayout, QHBoxLayout, 
                             QTableWidget, QTableWidgetItem, QPushButton, 
                             QTabWidget, QFileDialog, QMessageBox, QHeaderView,
                             QSizePolicy)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
from .components.traceplot import TracePlot
from .components.tracecomp import TraceComp


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
        upper_section.setLayout(upper_layout)
        
        # File table widget
        self.file_table = QTableWidget()
        self.file_table.setColumnCount(3)
        self.file_table.setHorizontalHeaderLabels(['Trace File', 'Generations', 'Burnin Fraction'])
        self.file_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.file_table.setSelectionMode(QTableWidget.MultiSelection)
        self.file_table.itemChanged.connect(self.on_burnin_fraction_changed)
        
        # Set column widths and properties
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
        lower_section.setLayout(lower_layout)
        
        # Statistics table widget
        self.stats_table = QTableWidget()
        self.stats_table.setColumnCount(2)
        self.stats_table.setHorizontalHeaderLabels(['Parameter', 'Mean'])
        self.stats_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.stats_table.setSelectionMode(QTableWidget.SingleSelection)
        self.stats_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.stats_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        
        lower_layout.addWidget(self.stats_table)
        
        # Add sections to left splitter with 3:7 ratio
        left_splitter.addWidget(upper_section)
        left_splitter.addWidget(lower_section)
        left_splitter.setSizes([300, 700])  # 3:7 ratio
        
        # Right panel: Analysis tabs
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_panel.setLayout(right_layout)
        
        # Tab widget for different analysis modes
        self.tab_widget = QTabWidget()
        self.single_chain_tab = QWidget()
        self.compare_chains_tab = QWidget()
        
        self.tab_widget.addTab(self.single_chain_tab, "Analyze Single Chain")
        self.tab_widget.addTab(self.compare_chains_tab, "Compare 2 Chains")
        
        # Initialize single chain analysis tab with a container widget
        single_layout = QVBoxLayout()
        self.single_chain_tab.setLayout(single_layout)
        self.trace_plot_container = QWidget()
        self.trace_plot_layout = QVBoxLayout()
        self.trace_plot_container.setLayout(self.trace_plot_layout)
        single_layout.addWidget(self.trace_plot_container)
        
        # Initialize compare chains analysis tab
        compare_layout = QVBoxLayout()
        self.compare_chains_tab.setLayout(compare_layout)
        self.trace_comp = TraceComp(burnin_fraction=0.1)
        compare_layout.addWidget(self.trace_comp)
        
        # Add the tab widget to the right panel layout
        right_layout.addWidget(self.tab_widget)
        
        # Connect button signals
        self.add_button.clicked.connect(self.add_files)
        self.remove_button.clicked.connect(self.remove_selected_files)
        self.refresh_button.clicked.connect(self.refresh_selected_files)
        
        # Connect file table selection change
        self.file_table.itemSelectionChanged.connect(self.on_file_selection_changed)
        
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
            
            # Burnin fraction column (editable)
            burnin_item = QTableWidgetItem("0.1")
            burnin_item.setData(Qt.UserRole, 0.1)  # Store actual value
            self.file_table.setItem(row, 2, burnin_item)
            
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
    
    def update_stats_table(self):
        """Update statistics table based on selected trace file and current burn-in fraction."""
        selected_items = self.file_table.selectedItems()
        if not selected_items:
            # Clear stats table
            self.stats_table.setRowCount(0)
            return
        
        # Get the row of the first selected item
        row = selected_items[0].row()
        if row >= len(self.mcmc_files):
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
        if len(selected_rows) == 1:
            file_selected_rows = set(item.row() for item in self.file_table.selectedItems())
            if len(file_selected_rows) == 1:
                file_row = next(iter(file_selected_rows))
                stats_row = next(iter(selected_rows))
                self.show_trace_plot_for_selection(file_row, stats_row)
    
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
            
        # Create TracePlot component if it doesn't exist
        if param_name not in self.trace_plots[file_path]:
            # Get current burn-in fraction
            burnin_frac = 0.1
            burnin_item = self.file_table.item(file_row, 2)
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
            
            if not filtered_data.empty:
                try:
                    # Get parameter values and iterations
                    param_series = filtered_data[param_name]
                    numeric_mask = pd.to_numeric(param_series, errors='coerce').notna()
                    valid_data = filtered_data[numeric_mask]
                    
                    if len(valid_data) > 0:
                        param_values = np.asarray(valid_data[param_name].values, dtype=float)
                        iterations = np.asarray(valid_data['iterations'].values, dtype=int)
                        
                        # Create and initialize the TracePlot component
                        self.trace_plots[file_path][param_name] = TracePlot(burnin_fraction=burnin_frac)
                        self.trace_plots[file_path][param_name].set_data(param_values, iterations)
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
            
            # Update burn-in fraction for this trace plot
            burnin_item = self.file_table.item(file_row, 2)
            if burnin_item:
                try:
                    burnin_frac = float(burnin_item.text())
                    burnin_frac = max(0.0, min(0.999, burnin_frac))
                    self.current_trace_plot.burnin_fraction = burnin_frac
                except ValueError:
                    pass
    
    def update_single_chain_analysis(self, file_row, param_row=None):
        """Update single chain analysis tab - now handled by show_trace_plot_for_selection."""
        # This method is now primarily used for compatibility, but the actual display logic
        # is handled by show_trace_plot_for_selection
        pass
        
    def update_compare_chains_analysis(self, file_row1, file_row2):
        """Update compare chains analysis tab."""
        if file_row1 >= len(self.mcmc_files) or file_row2 >= len(self.mcmc_files):
            return
            
        file1 = self.mcmc_files[file_row1]
        file2 = self.mcmc_files[file_row2]
        
        if file1 not in self.mcmc_data or file2 not in self.mcmc_data:
            return
            
        data1 = self.mcmc_data[file1]
        data2 = self.mcmc_data[file2]
        
        if data1 is None or data2 is None or data1.empty or data2.empty:
            return
            
        # Check if both have same number of generations
        if len(data1) != len(data2):
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
                
        # Only compare if burnin fractions are the same
        if abs(burnin1 - burnin2) > 1e-6:
            return
            
        # Apply burnin
        burnin_samples = int(len(data1) * burnin1)
        filtered_data1 = data1.iloc[burnin_samples:]
        filtered_data2 = data2.iloc[burnin_samples:]
        
        if filtered_data1.empty or filtered_data2.empty:
            return
            
        # Get common parameters (excluding 'iterations')
        params1 = set(col for col in filtered_data1.columns if col != 'iterations')
        params2 = set(col for col in filtered_data2.columns if col != 'iterations')
        common_params = list(params1 & params2)
        
        if not common_params:
            return
            
        # Use the first common parameter
        param_name = common_params[0]
        
        # Prepare data for tracecomp
        chain1_data = filtered_data1[param_name].values
        chain2_data = filtered_data2[param_name].values
        
        # Create a simple dict structure that tracecomp can understand
        trace_data = {
            'chain1': chain1_data,
            'chain2': chain2_data,
            'param_name': param_name
        }
        
        self.trace_comp.set_data(trace_data)
        
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
        
    def get_tool_name(self):
        """Return the tool name for plugin registration."""
        return "MiniTracer"
        
    def get_input_format(self):
        """Return supported input formats."""
        return ["*.log", "*.txt", "*.csv"]
        
    def get_output_format(self):
        """Return output format."""
        return []

class MiniTracerPluginEntry:
    def run(self):
        return MiniTracerPlugin()
