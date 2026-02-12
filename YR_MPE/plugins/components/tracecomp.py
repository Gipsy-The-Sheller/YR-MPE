import pandas as pd
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QTableWidget, QTableWidgetItem, QApplication, QMainWindow
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont


class TraceComp(QWidget):
    """
    QWidget for comparing convergence between two MCMC chains.
    Displays a scatter plot of post-burnin samples from both chains vs each other,
    and a table with statistics including sample counts, means, and correlation coefficient.
    """
    
    def __init__(self, burnin_fraction=0.1, parent=None):
        super().__init__(parent)
        self.burnin_fraction = burnin_fraction
        self.chain1_data = None
        self.chain1_iterations = None
        self.chain2_data = None
        self.chain2_iterations = None
        
        self.init_ui()
    
    def init_ui(self):
        """Initialize the user interface."""
        layout = QHBoxLayout()
        
        # Create matplotlib figure for scatter plot (70% width)
        self.figure = plt.figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)
        
        # Create table widget for statistics (30% width)
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(['Metric', 'Value'])
        self.table.setRowCount(5)
        # Disable vertical header (index) to avoid display issues
        self.table.verticalHeader().setVisible(False)
        
        # Set the metric labels in the first column
        metrics = [
            'Samples (Chain 1)',
            'Mean (Chain 1)', 
            'Samples (Chain 2)',
            'Mean (Chain 2)',
            'Correlation Coefficiency (R)'
        ]
        for i, metric in enumerate(metrics):
            item = QTableWidgetItem(metric)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)  # Make non-editable
            self.table.setItem(i, 0, item)
        
        # Reduce font size and spacing for better fit
        font = QFont()
        font.setPointSize(9)  # Smaller font size
        self.table.setFont(font)
        self.table.horizontalHeader().setFont(font)
        self.table.verticalHeader().setFont(font)
        
        # Adjust row heights and column widths
        for i in range(5):
            self.table.setRowHeight(i, 25)  # Reduced row height
        self.table.setColumnWidth(0, 120)  # Adjust column widths
        self.table.setColumnWidth(1, 120)
        
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)  # Make table read-only
        
        # Set stretch factors: 70% for plot, 30% for table (7:3 ratio)
        layout.addWidget(self.canvas, 7)  # 7/(7+3) = 70%
        layout.addWidget(self.table, 3)   # 3/(7+3) = 30%
        
        self.setLayout(layout)
    
    def set_data(self, chain1_data, chain1_iterations, chain2_data, chain2_iterations):
        """
        Set the data for both chains.
        
        Parameters:
        chain1_data: array-like, data values for chain 1
        chain1_iterations: array-like, iteration numbers for chain 1
        chain2_data: array-like, data values for chain 2  
        chain2_iterations: array-like, iteration numbers for chain 2
        """
        if chain1_data is None or chain2_data is None:
            self.chain1_data = None
            self.chain2_data = None
            self.clear_plot_and_table()
            return
            
        self.chain1_data = np.array(chain1_data)
        self.chain1_iterations = np.array(chain1_iterations)
        self.chain2_data = np.array(chain2_data)
        self.chain2_iterations = np.array(chain2_iterations)
        self.update_plot_and_table()
    
    def clear_plot_and_table(self):
        """Clear the plot and table when no data is available."""
        self.ax.clear()
        self.ax.set_title('Chain Comparison (No Data)')
        self.canvas.draw()
        
        # Clear table values
        for i in range(5):
            self.table.setItem(i, 1, QTableWidgetItem(""))
    
    def update_burnin_fraction(self, burnin_fraction):
        """Update the burn-in fraction and refresh the plot and table."""
        self.burnin_fraction = burnin_fraction
        if self.chain1_data is not None:
            self.update_plot_and_table()
    
    def update_plot_and_table(self):
        """Update both the scatter plot and the statistics table."""
        if self.chain1_data is None or self.chain2_data is None:
            return
        
        # Use all data directly (data is already post-burnin from MiniTracerPlugin)
        post_burnin1 = self.chain1_data
        post_burnin2 = self.chain2_data
        
        # Clear the plot
        self.ax.clear()
        
        # Plot all points (opaque blue) - no pre-burnin since data is already filtered
        min_len = min(len(post_burnin1), len(post_burnin2))
        if min_len > 0:
            self.ax.scatter(post_burnin1[:min_len], post_burnin2[:min_len], 
                          alpha=1.0, color='#4e699a', s=20)
            
            # Calculate correlation coefficient
            if min_len > 1:
                correlation_matrix = np.corrcoef(post_burnin1[:min_len], post_burnin2[:min_len])
                r_value = correlation_matrix[0, 1]
            else:
                r_value = 0.0
                
            # Set axis limits based on all data
            x_min, x_max = post_burnin1[:min_len].min(), post_burnin1[:min_len].max()
            y_min, y_max = post_burnin2[:min_len].min(), post_burnin2[:min_len].max()
            
            # Add some padding
            x_padding = (x_max - x_min) * 0.05
            y_padding = (y_max - y_min) * 0.05
            self.ax.set_xlim(x_min - x_padding, x_max + x_padding)
            self.ax.set_ylim(y_min - y_padding, y_max + y_padding)
            
            # Add y=x reference line for convergence assessment
            line_min = min(x_min - x_padding, y_min - y_padding)
            line_max = max(x_max + x_padding, y_max + y_padding)
            self.ax.plot([line_min, line_max], [line_min, line_max], 
                        color='#ba3e45', linestyle='--', linewidth=1, alpha=0.8)
        else:
            r_value = 0.0
        
        self.ax.set_xlabel('Chain 1')
        self.ax.set_ylabel('Chain 2')
        self.ax.set_title('Chain Comparison')
        self.ax.grid(True, alpha=0.3)
        
        # Update the table with actual data counts and statistics
        count1 = len(post_burnin1)
        count2 = len(post_burnin2)
        mean1 = np.mean(post_burnin1) if count1 > 0 else 0
        mean2 = np.mean(post_burnin2) if count2 > 0 else 0
        
        self.update_statistics_table(count1, mean1, count2, mean2, r_value)
        
        self.canvas.draw()
    
    def update_statistics_table(self, count1, mean1, count2, mean2, r_value):
        """Update the statistics table with calculated values."""
        # Update values in the second column
        self.table.setItem(0, 1, QTableWidgetItem(str(count1)))
        self.table.setItem(1, 1, QTableWidgetItem(f"{mean1:.6f}"))
        self.table.setItem(2, 1, QTableWidgetItem(str(count2)))
        self.table.setItem(3, 1, QTableWidgetItem(f"{mean2:.6f}"))
        self.table.setItem(4, 1, QTableWidgetItem(f"{r_value:.6f}"))


def parse_trace_file(trace_file):
    # 读取所有行
    with open(trace_file, 'r') as f:
        lines = f.readlines()
    
    # 找到包含"iter"的header行和包含"0"的数据起始行
    header_line = None
    data_start_line = None
    
    last_col, st_num = 0,0
    for index, line in enumerate(lines):
        line = line.strip().split('\t')
        last_col = len(line)
        if line[0] == '0':
            st_num = index
            break
    if last_col <= 1:
        raise Exception('Trace file is empty')
    # for i, line in enumerate(lines):
    #     if line.strip().startswith('Gen'):
    #         header_line = i
    #     elif line.strip().startswith('0\t'):
    #         data_start_line = i
    #         break
    
    # if header_line is None or data_start_line is None:
    #     raise ValueError(f"Could not find header or data start in trace file:\n{trace_file}")
    
    # 使用pandas读取，指定header行和跳过前面的行
    df = pd.read_csv(trace_file, sep='\t', header=st_num-1, skiprows=st_num-2)
    return df


if __name__ == '__main__':
    import sys
    import os
    
    # Try to use actual files (if they exist)
    trace_file1 = r"D:\Program Files\PhyloSuite_v1.2.3_Win64_with_plugins\PhyloSuite\myWorkPlace\GenBank_File\files\MrBayes_results\2025_05_12-19_31_54\input.nex.run1.p"
    trace_file2 = r"D:\Program Files\PhyloSuite_v1.2.3_Win64_with_plugins\PhyloSuite\myWorkPlace\GenBank_File\files\MrBayes_results\2025_05_12-19_31_54\input.nex.run2.p"  # Assuming second chain file
    
    chain1_data = None
    chain1_iterations = None
    chain2_data = None
    chain2_iterations = None
    
    # Load first chain
    if os.path.exists(trace_file1):
        try:
            trace1 = parse_trace_file(trace_file1)
            print("Loaded trace file 1:")
            print(trace1.head())
            if 'Gen' in trace1.columns and 'LnL' in trace1.columns:
                chain1_iterations = trace1['Gen'].values
                chain1_data = trace1['LnL'].values
                chain1_iterations = pd.to_numeric(chain1_iterations, errors='coerce')
                chain1_data = pd.to_numeric(chain1_data, errors='coerce')
                valid_mask1 = ~(np.isnan(chain1_iterations) | np.isnan(chain1_data))
                chain1_iterations = chain1_iterations[valid_mask1]
                chain1_data = chain1_data[valid_mask1]
                print(f"Successfully loaded {len(chain1_data)} data points from chain 1")
        except Exception as e:
            print(f"Failed to load trace file 1: {e}")
            import traceback
            traceback.print_exc()
    
    # Load second chain
    if os.path.exists(trace_file2):
        try:
            trace2 = parse_trace_file(trace_file2)
            print("Loaded trace file 2:")
            print(trace2.head())
            if 'Gen' in trace2.columns and 'LnL' in trace2.columns:
                chain2_iterations = trace2['Gen'].values
                chain2_data = trace2['LnL'].values
                chain2_iterations = pd.to_numeric(chain2_iterations, errors='coerce')
                chain2_data = pd.to_numeric(chain2_data, errors='coerce')
                valid_mask2 = ~(np.isnan(chain2_iterations) | np.isnan(chain2_data))
                chain2_iterations = chain2_iterations[valid_mask2]
                chain2_data = chain2_data[valid_mask2]
                print(f"Successfully loaded {len(chain2_data)} data points from chain 2")
        except Exception as e:
            print(f"Failed to load trace file 2: {e}")
            import traceback
            traceback.print_exc()
    
    # # If files don't exist or loading failed, create simulated data
    # if chain1_data is None or len(chain1_data) == 0 or chain2_data is None or len(chain2_data) == 0:
    #     print("Using simulated data instead.")
    #     np.random.seed(42)
    #     n_samples = 1000
        
    #     # Generate correlated data for two chains
    #     true_mean = 1.0
    #     true_std = 0.2
        
    #     # Chain 1
    #     burnin_samples1 = int(n_samples * 0.2)
    #     burnin_data1 = np.linspace(2.0, true_mean, burnin_samples1) + np.random.normal(0, 0.1, burnin_samples1)
    #     post_burnin_data1 = true_mean + np.random.normal(0, true_std, n_samples - burnin_samples1)
    #     chain1_data = np.concatenate([burnin_data1, post_burnin_data1])
    #     chain1_iterations = np.arange(n_samples)
        
    #     # Chain 2 (correlated with chain 1, but ensure same post-burnin length)
    #     burnin_samples2 = int(n_samples * 0.25)
    #     burnin_data2 = np.linspace(1.8, true_mean, burnin_samples2) + np.random.normal(0, 0.12, burnin_samples2)
    #     # Make sure post-burnin data has the same length as chain1's post-burnin data
    #     post_burnin_length = n_samples - burnin_samples1
    #     # Add some correlation with chain 1's post-burnin data
    #     correlated_part = post_burnin_data1 * 0.8 + np.random.normal(true_mean * 0.2, true_std * 0.6, post_burnin_length)
    #     post_burnin_data2 = correlated_part
    #     chain2_data = np.concatenate([burnin_data2, post_burnin_data2])
    #     chain2_iterations = np.arange(len(chain2_data))
    
    # Create Qt application
    app = QApplication(sys.argv)
    
    # Create TraceComp widget
    trace_comp = TraceComp(burnin_fraction=0.1)
    trace_comp.set_data(chain1_data, chain1_iterations, chain2_data, chain2_iterations)
    
    # Create main window and set central widget
    main_window = QMainWindow()
    main_window.setCentralWidget(trace_comp)
    main_window.setWindowTitle("MCMC Chain Comparison")
    main_window.resize(1000, 600)
    main_window.show()
    
    # Example of updating burn-in fraction interactively:
    # trace_comp.update_burnin_fraction(0.2)
    
    sys.exit(app.exec_())