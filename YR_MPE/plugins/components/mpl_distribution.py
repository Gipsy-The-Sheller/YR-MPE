import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy import stats
import math

class MplDistribution(FigureCanvas):
    """封装的分布函数绘制器，继承自FigureCanvasQTAgg"""
    
    def __init__(self, parent=None, width=6, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        
        # 设置初始状态
        self.ax = self.fig.add_subplot(111)
        self.current_distribution = None
        self.current_params = {}
        
    def plot_point_distribution(self, point_value):
        """绘制点分布（在指定点处的垂直线）"""
        self.ax.clear()
        
        # 绘制垂直线表示点分布
        self.ax.axvline(x=point_value, color='red', linewidth=2, label=f'Point: {point_value}')
        
        # 设置图形属性
        self.ax.set_xlabel('Time')
        self.ax.set_ylabel('Probability Density')
        self.ax.set_title('Point Distribution')
        self.ax.legend()
        self.ax.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()
        self.current_distribution = 'point'
        self.current_params = {'point': point_value}
        
    def plot_uniform_distribution(self, lower_bound, upper_bound):
        """绘制均匀分布"""
        self.ax.clear()
        
        if lower_bound >= upper_bound:
            # 处理无效参数
            self.ax.text(0.5, 0.5, 'Invalid parameters: lower_bound >= upper_bound', 
                        transform=self.ax.transAxes, ha='center', va='center')
            self.ax.set_title('Uniform Distribution (Invalid Parameters)')
        else:
            # 生成x轴数据
            x_range = upper_bound - lower_bound
            x_min = lower_bound - 0.1 * x_range if x_range > 0 else lower_bound - 1
            x_max = upper_bound + 0.1 * x_range if x_range > 0 else upper_bound + 1
            x = np.linspace(x_min, x_max, 1000)
            
            # 计算均匀分布的概率密度
            y = np.where((x >= lower_bound) & (x <= upper_bound), 1.0 / (upper_bound - lower_bound), 0)
            
            # 绘制分布
            self.ax.plot(x, y, 'b-', linewidth=2, label=f'Uniform({lower_bound:.2f}, {upper_bound:.2f})')
            self.ax.fill_between(x, y, alpha=0.3, color='blue')
            
            # 标记边界
            self.ax.axvline(x=lower_bound, color='red', linestyle='--', alpha=0.7)
            self.ax.axvline(x=upper_bound, color='red', linestyle='--', alpha=0.7)
            
            # 设置图形属性
            self.ax.set_xlabel('Time')
            self.ax.set_ylabel('Probability Density')
            self.ax.set_title('Uniform Distribution')
            self.ax.legend()
            self.ax.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()
        self.current_distribution = 'uniform'
        self.current_params = {'lower': lower_bound, 'upper': upper_bound}
        
    def plot_upper_boundary_distribution(self, upper_bound):
        """绘制上边界分布（假设为指数衰减分布）"""
        self.ax.clear()
        
        # 对于上边界，我们可以使用一个从0到upper_bound的递增分布
        # 这里使用简单的线性分布作为示例
        x_max = upper_bound * 1.2 if upper_bound > 0 else 1.0
        x = np.linspace(0, x_max, 1000)
        
        # 创建一个在upper_bound处截断的分布（例如三角形分布）
        y = np.where(x <= upper_bound, (2 * x) / (upper_bound ** 2), 0)
        
        # 绘制分布
        self.ax.plot(x, y, 'g-', linewidth=2, label=f'Upper Boundary: {upper_bound:.2f}')
        self.ax.fill_between(x, y, alpha=0.3, color='green')
        self.ax.axvline(x=upper_bound, color='red', linestyle='--', alpha=0.7, label='Upper Bound')
        
        # 设置图形属性
        self.ax.set_xlabel('Time')
        self.ax.set_ylabel('Probability Density')
        self.ax.set_title('Upper Boundary Distribution')
        self.ax.legend()
        self.ax.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()
        self.current_distribution = 'upper_boundary'
        self.current_params = {'upper': upper_bound}
        
    def plot_lower_boundary_distribution(self, lower_bound):
        """绘制下边界分布（假设为从lower_bound开始的指数衰减）"""
        self.ax.clear()
        
        # 对于下边界，使用从lower_bound开始的指数衰减分布
        x_min = lower_bound * 0.8 if lower_bound > 0 else 0
        x_max = lower_bound * 3 if lower_bound > 0 else 5.0
        x = np.linspace(x_min, x_max, 1000)
        
        # 创建指数衰减分布
        if lower_bound > 0:
            decay_rate = 1.0 / lower_bound  # 衰减率基于下边界值
            y = np.where(x >= lower_bound, decay_rate * np.exp(-decay_rate * (x - lower_bound)), 0)
        else:
            # 如果下边界为0或负数，使用简单的指数分布
            y = np.exp(-x)
            y[x < lower_bound] = 0
        
        # 绘制分布
        self.ax.plot(x, y, 'orange', linewidth=2, label=f'Lower Boundary: {lower_bound:.2f}')
        self.ax.fill_between(x, y, alpha=0.3, color='orange')
        self.ax.axvline(x=lower_bound, color='red', linestyle='--', alpha=0.7, label='Lower Bound')
        
        # 设置图形属性
        self.ax.set_xlabel('Time')
        self.ax.set_ylabel('Probability Density')
        self.ax.set_title('Lower Boundary Distribution')
        self.ax.legend()
        self.ax.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()
        self.current_distribution = 'lower_boundary'
        self.current_params = {'lower': lower_bound}
        
    def plot_normal_distribution(self, mean, std_dev):
        """绘制正态分布"""
        self.ax.clear()
        
        if std_dev <= 0:
            # 处理无效参数
            self.ax.text(0.5, 0.5, 'Invalid parameters: std_dev <= 0', 
                        transform=self.ax.transAxes, ha='center', va='center')
            self.ax.set_title('Normal Distribution (Invalid Parameters)')
        else:
            # 生成x轴数据（覆盖均值±4个标准差的范围）
            x_min = mean - 4 * std_dev
            x_max = mean + 4 * std_dev
            x = np.linspace(x_min, x_max, 1000)
            
            # 计算正态分布的概率密度
            y = stats.norm.pdf(x, mean, std_dev)
            
            # 绘制分布
            self.ax.plot(x, y, 'purple', linewidth=2, label=f'Normal(μ={mean:.2f}, σ={std_dev:.2f})')
            self.ax.fill_between(x, y, alpha=0.3, color='purple')
            
            # 标记均值
            self.ax.axvline(x=mean, color='red', linestyle='--', alpha=0.7, label='Mean')
            
            # 设置图形属性
            self.ax.set_xlabel('Time')
            self.ax.set_ylabel('Probability Density')
            self.ax.set_title('Normal Distribution')
            self.ax.legend()
            self.ax.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()
        self.current_distribution = 'normal'
        self.current_params = {'mean': mean, 'std': std_dev}
        
    def plot_lognormal_distribution(self, mean, std_dev):
        """绘制对数正态分布"""
        self.ax.clear()
        
        if std_dev <= 0 or mean <= 0:
            # 处理无效参数
            self.ax.text(0.5, 0.5, 'Invalid parameters: std_dev <= 0 or mean <= 0', 
                        transform=self.ax.transAxes, ha='center', va='center')
            self.ax.set_title('Lognormal Distribution (Invalid Parameters)')
        else:
            # 对数正态分布的参数转换
            # scipy.stats.lognorm 使用形状参数 s = sigma, 位置参数 loc = 0, 尺度参数 scale = exp(mu)
            mu = math.log(mean**2 / math.sqrt(std_dev**2 + mean**2))
            sigma = math.sqrt(math.log(1 + (std_dev**2 / mean**2)))
            
            # 生成x轴数据
            x_min = 0
            x_max = mean + 4 * std_dev
            x = np.linspace(x_min, x_max, 1000)
            
            # 计算对数正态分布的概率密度
            y = stats.lognorm.pdf(x, sigma, scale=math.exp(mu))
            
            # 绘制分布
            self.ax.plot(x, y, 'brown', linewidth=2, label=f'Lognormal(mean={mean:.2f}, std={std_dev:.2f})')
            self.ax.fill_between(x, y, alpha=0.3, color='brown')
            
            # 标记均值
            self.ax.axvline(x=mean, color='red', linestyle='--', alpha=0.7, label='Mean')
            
            # 设置图形属性
            self.ax.set_xlabel('Time')
            self.ax.set_ylabel('Probability Density')
            self.ax.set_title('Lognormal Distribution')
            self.ax.legend()
            self.ax.grid(True, alpha=0.3)
        
        self.fig.tight_layout()
        self.draw()
        self.current_distribution = 'lognormal'
        self.current_params = {'mean': mean, 'std': std_dev}
        
    def plot_distribution(self, distribution_type, **params):
        """根据分布类型和参数绘制相应的分布"""
        if distribution_type == 'Point':
            self.plot_point_distribution(params.get('point', 0))
        elif distribution_type == 'Uniform':
            self.plot_uniform_distribution(params.get('lower', 0), params.get('upper', 1))
        elif distribution_type == 'Upper Boundary':
            self.plot_upper_boundary_distribution(params.get('upper', 1))
        elif distribution_type == 'Lower Boundary':
            self.plot_lower_boundary_distribution(params.get('lower', 0))
        elif distribution_type == 'Normal':
            self.plot_normal_distribution(params.get('mean', 0), params.get('std', 1))
        elif distribution_type == 'Lognormal':
            self.plot_lognormal_distribution(params.get('mean', 1), params.get('std', 1))
        else:
            self.ax.clear()
            self.ax.text(0.5, 0.5, f'Unknown distribution: {distribution_type}', 
                        transform=self.ax.transAxes, ha='center', va='center')
            self.ax.set_title('Unknown Distribution')
            self.fig.tight_layout()
            self.draw()
            
    def clear_plot(self):
        """清空当前绘图"""
        self.ax.clear()
        self.ax.set_title('No Distribution Selected')
        self.fig.tight_layout()
        self.draw()
        self.current_distribution = None
        self.current_params = {}