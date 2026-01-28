import pandas as pd
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QApplication, QMainWindow, QSizePolicy
from PyQt5.QtCore import Qt

class TracePlot(QWidget):
    def __init__(self, burnin_fraction=0.1):
        super().__init__()
        self.burnin_fraction = burnin_fraction
        self.data = None
        self.iterations = None
        self.init_ui()

    def init_ui(self):
        layout = QHBoxLayout()
        self.setLayout(layout)
        
        # 设置大小策略，使组件能自动拉伸填充父容器
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # 设置布局间距为0
        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)

        # 创建迹线图（占80%宽度）
        self.trace_figure = plt.figure(figsize=(8, 6))
        self.trace_canvas = FigureCanvas(self.trace_figure)
        layout.addWidget(self.trace_canvas, 8)  # 80%宽度
        
        # 创建密度图（占20%宽度）
        self.density_figure = plt.figure(figsize=(2, 6))
        self.density_canvas = FigureCanvas(self.density_figure)
        layout.addWidget(self.density_canvas, 2)  # 20%宽度

        self.trace_ax = self.trace_figure.add_subplot(111)
        self.density_ax = self.density_figure.add_subplot(111)

    def set_data(self, data, iterations=None):
        """设置要显示的数据
        
        Args:
            data: 时间（time）数据
            iterations: 迭代次数（iter）数据，如果为None则自动生成
        """
        self.data = np.asarray(data, dtype=float)
        if iterations is not None:
            self.iterations = np.asarray(iterations, dtype=int)
        else:
            self.iterations = np.arange(len(data))
        self.update_plot()

    def update_burnin_fraction(self, burnin_fraction):
        """更新burn-in比例并重新绘制图形"""
        self.burnin_fraction = burnin_fraction
        if self.data is not None:
            self.update_plot()

    def update_plot(self):
        """更新迹线图和密度图"""
        if self.data is None or len(self.data) == 0:
            return
            
        n_total = len(self.data)
        burnin_index = int(n_total * self.burnin_fraction)
        
        # 分离burn-in数据和有效数据
        burnin_data = self.data[:burnin_index]
        post_burnin_data = self.data[burnin_index:]
        burnin_iters = self.iterations[:burnin_index]
        post_burnin_iters = self.iterations[burnin_index:]
        
        if len(post_burnin_data) == 0:
            return
            
        # 清除之前的图形
        self.trace_ax.clear()
        self.density_ax.clear()
        
        # 设置坐标范围仅基于post-burnin数据
        y_min = np.min(post_burnin_data)
        y_max = np.max(post_burnin_data)
        y_range = y_max - y_min
        y_padding = y_range * 0.05  # 添加5%的padding
        y_lim = (y_min - y_padding, y_max + y_padding)
        
        # 绘制burn-in部分（50%透明度）
        if len(burnin_data) > 0:
            self.trace_ax.plot(burnin_iters, burnin_data, 
                             color='gray', alpha=0.5, linewidth=0.8)
            # self.trace_ax.scatter(burnin_iters, burnin_data, 
            #                     color='gray', alpha=0.5, s=10)
        
        # 绘制post-burnin部分（不透明）
        self.trace_ax.plot(post_burnin_iters, post_burnin_data, 
                         color='#4e699a', linewidth=0.8)
        # self.trace_ax.scatter(post_burnin_iters, post_burnin_data, 
        #                     color='blue', s=2)
        
        # 设置y轴范围仅基于post-burnin数据
        self.trace_ax.set_ylim(y_lim)
        self.trace_ax.set_xlabel('Iteration')
        self.trace_ax.set_ylabel('LnL')
        self.trace_ax.set_title('MCMC Trace')
        self.trace_ax.grid(True, alpha=0.3)
        
        # 绘制密度图（KDE和直方图）
        # 计算最优带宽（Silverman's rule of thumb）
        if len(post_burnin_data) > 1:
            std_val = np.std(post_burnin_data)
            h = 1.06 * std_val * (len(post_burnin_data) ** (-1/5))
            h = max(h, 1e-6)  # 避免带宽为0
            
            # KDE计算
            X, kde_vals = kernel_density_estimation_gaussian_vectorized(
                post_burnin_data, x_range=y_lim, h=h, n_points=200
            )
            
            # 绘制KDE曲线
            self.density_ax.plot(kde_vals, X, color='#ba3e45', linewidth=2)
            
            # 绘制直方图（水平方向）
            hist_bins = min(30, len(post_burnin_data) // 10)
            if hist_bins < 5:
                hist_bins = 5
            counts, bins = np.histogram(post_burnin_data, bins=hist_bins, density=True)
            bin_centers = (bins[:-1] + bins[1:]) / 2
            # 归一化直方图为概率密度
            self.density_ax.barh(bin_centers, counts, height=(bins[1] - bins[0]) * 0.8, 
                               alpha=0.6, color='lightblue', edgecolor='navy')
        
        # 设置密度图的y轴范围与迹线图一致
        self.density_ax.set_ylim(y_lim)
        self.density_ax.set_xlabel('Density')
        self.density_ax.set_ylabel('LnL')
        self.density_ax.set_title('Post-burnin Distribution')
        self.density_ax.grid(True, alpha=0.3)
        
        # 更新画布
        self.trace_canvas.draw()
        self.density_canvas.draw()


def kernel_density_estimation_gaussian_vectorized(data, x_range=None, h=1.0, n_points=100, 
                                        return_xy=True):
    """
    NOTE: This function is implemented by Deepseek-R1!
    向量化版本的核密度估计函数（效率更高）
    
    参数:
    ----------
    data : array_like
        输入的一维数据
    x_range : tuple, 可选
        计算KDE的x范围，格式为(min, max)
    h : float, 可选
        带宽(bandwidth)
    n_points : int, 可选
        在x范围内取点的数量
    kernel : str, 可选
        核函数类型，目前只实现'gaussian'（高斯核）
    return_xy : bool, 可选
        如果为True，返回x和y；如果为False，只返回y值
    
    返回:
    ----------
    如果 return_xy=True: 返回X, p
    如果 return_xy=False: 只返回p
    """
    
    data = np.asarray(data)
    n = len(data)
    
    # 如果未指定x_range，则使用数据的最小值-3h到最大值+3h
    if x_range is None:
        data_min, data_max = np.min(data), np.max(data)
        x_range = (data_min - 3*h, data_max + 3*h)
    
    # 生成x点
    X = np.linspace(x_range[0], x_range[1], n_points)
    
    # 创建网格：X (n_points,) 和 data (n,)
    # 通过广播计算所有点对的核函数值
    X_grid, data_grid = np.meshgrid(X, data, indexing='ij')
    
    p = np.exp(-((X_grid - data_grid) / h) ** 2 / 2) / (h * np.sqrt(2 * np.pi))
    
    # 对每个x点，对所有数据点的核函数值求平均
    p = np.mean(p, axis=1)
    
    if return_xy:
        return X, p
    else:
        return p

def parse_trace_file(trace_file):
    # 读取所有行
    with open(trace_file, 'r') as f:
        lines = f.readlines()
    
    # 找到包含"iter"的header行和包含"0"的数据起始行
    header_line = None
    data_start_line = None
    
    last_col, st_num = 0,0
    for index, line in enumerate(lines):
        # print(line)
        line = line.strip().split('\t')
        last_col = len(line)
        if line[0] == '0':
            st_num = index
            break
    if last_col <= 1:
        raise Exception('Trace file is empty')
    
    # 使用pandas读取，指定header行和跳过前面的行
    df = pd.read_csv(trace_file, sep='\t', header=0, skiprows=st_num-1)
    return df

if __name__ == '__main__':
    import sys
    import os
    
    # 尝试使用实际文件（如果存在）
    trace_file   = r"E:\PhyloGenetics\2026-1-9-BEAST\run3\run1.log"
    trace_file_2 = r"D:\Program Files\PhyloSuite_v1.2.3_Win64_with_plugins\PhyloSuite\myWorkPlace\GenBank_File\files\MrBayes_results\2025_05_12-19_31_54\input.nex.run1.p"
    example_data = None
    example_iterations = None
    
    if os.path.exists(trace_file):
        try:
            trace = parse_trace_file(trace_file)
            # print("Loaded trace file:")
            # print(trace.head())
            # 使用'Gen'列作为iterations，'LnL'列作为time数据
            # if 'Gen' in trace.columns and 'LnL' in trace.columns:
            example_iterations = trace[trace.columns[0]].values
            example_data = trace[trace.columns[1]].values
            print(example_data[:50])
            # 确保数据是数值类型
            example_iterations = pd.to_numeric(example_iterations, errors='coerce')
            example_data = pd.to_numeric(example_data, errors='coerce')
            # 移除任何NaN值
            valid_mask = ~(np.isnan(example_iterations) | np.isnan(example_data))
            example_iterations = example_iterations[valid_mask]
            example_data = example_data[valid_mask]
            print(f"Successfully loaded {len(example_data)} data points")
        except Exception as e:
            print(f"Failed to load trace file: {e}")
            import traceback
            traceback.print_exc()
    
    # 如果没有实际文件或加载失败，创建示例数据
    if example_data is None or len(example_data) == 0:
        print("Using simulated data instead.")
        np.random.seed(42)
        # 生成模拟的MCMC迹线数据（具有burn-in特性）
        n_samples = 1000
        # burn-in阶段：从较高值开始逐渐收敛
        burnin_samples = int(n_samples * 0.2)
        burnin_data = np.linspace(2.0, 1.0, burnin_samples) + np.random.normal(0, 0.1, burnin_samples)
        # post-burnin阶段：围绕真实值波动
        post_burnin_data = 1.0 + np.random.normal(0, 0.2, n_samples - burnin_samples)
        example_data = np.concatenate([burnin_data, post_burnin_data])
        example_iterations = np.arange(n_samples)
    
    # 创建Qt应用
    app = QApplication(sys.argv)
    
    # 创建TracePlot窗口
    trace_plot = TracePlot(burnin_fraction=0.1)
    trace_plot.set_data(example_data, example_iterations)
    
    # 创建主窗口并设置central widget
    main_window = QMainWindow()
    main_window.setCentralWidget(trace_plot)
    main_window.setWindowTitle("MCMC Trace Plot")
    main_window.resize(1000, 600)
    main_window.show()
    
    # 测试update_burnin_fraction功能
    # 可以通过以下方式交互式更新burn-in比例：
    # trace_plot.update_burnin_fraction(0.2)
    
    sys.exit(app.exec_())
