import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
# import numpy as np
# import re
import time  # 添加时间模块

class Node:
    def __init__(self, name="", branch_length=0.0, is_leaf=False):
        self.name = name
        self.branch_length = branch_length
        self.age = 0.0
        self.order = 0.0
        self.is_leaf = is_leaf
        self.parent = None
        self.children = []
        
    def add_child(self, child):
        child.parent = self
        self.children.append(child)

class NewickParser:
    def __init__(self, newick_str):
        self.newick_str = newick_str.strip()
        self.pos = 0
        self.node_order = 0
        
    def skip_whitespace(self):
        while self.pos < len(self.newick_str) and self.newick_str[self.pos] in ' \t\n':
            self.pos += 1
            
    def parse_label(self):
        """解析标签（节点名和分支长度）"""
        self.skip_whitespace()
        start = self.pos
        
        # 读取直到遇到分隔符
        while (self.pos < len(self.newick_str) and 
               self.newick_str[self.pos] not in ',);'):
            self.pos += 1
            
        if self.pos > start:
            label = self.newick_str[start:self.pos].strip()
            
            # 检查是否有分支长度
            if ':' in label:
                parts = label.split(':', 1)
                name = parts[0].strip() if parts[0].strip() else ""
                try:
                    branch_length = float(parts[1].strip())
                except ValueError:
                    branch_length = 0.0
            else:
                name = label
                branch_length = 0.0
                
            return name, branch_length
            
        return "", 0.0
        
    def parse(self):
        """解析Newick字符串"""
        start_time = time.time()  # 添加开始时间记录
        
        self.pos = 0
        root = Node()
        stack = [root]
        
        while self.pos < len(self.newick_str) and self.newick_str[self.pos] != ';':
            self.skip_whitespace()
            
            if self.pos >= len(self.newick_str):
                break
                
            char = self.newick_str[self.pos]
            
            if char == '(':
                # 创建新的内部节点
                new_node = Node()
                stack[-1].add_child(new_node)
                stack.append(new_node)
                self.pos += 1
                
            elif char == ')':
                # 完成当前内部节点
                completed_node = stack.pop()
                self.pos += 1  # 跳过 ')'
                
                # 解析标签
                name, branch_length = self.parse_label()
                completed_node.name = name
                completed_node.branch_length = branch_length
                
            elif char == ',':
                # 当前节点的兄弟节点
                self.pos += 1
                
            else:
                # 叶子节点
                leaf = Node(is_leaf=True)
                name, branch_length = self.parse_label()
                leaf.name = name
                leaf.branch_length = branch_length
                leaf.order = self.node_order
                self.node_order += 1
                
                stack[-1].add_child(leaf)
        
        end_time = time.time()  # 添加结束时间记录
        parse_time = end_time - start_time
        print(f"Python解析耗时: {parse_time:.4f} 秒")
        
        # 计算节点坐标
        self.calculate_coordinates(root)
        return root
        
    def calculate_coordinates(self, root):
        """计算节点坐标"""
        # 设置根节点age为0
        root.age = 0.0
        
        # 第一次遍历：计算age值
        def calculate_age(node):
            if node.parent:
                node.age = node.parent.age + node.branch_length
            for child in node.children:
                calculate_age(child)
                
        calculate_age(root)
        
        # 第二次遍历：计算内部节点的order值
        def calculate_order(node):
            # 先递归计算子节点
            for child in node.children:
                calculate_order(child)
                
            # 如果是内部节点，计算其order值（基于直接子节点的平均值）
            if not node.is_leaf and node.children:
                node.order = sum(child.order for child in node.children) / len(node.children)
            elif node.is_leaf:
                # 叶节点的order已经在解析时设置了
                pass
                
        calculate_order(root)

def collect_nodes(root):
    """收集所有节点信息"""
    nodes = []
    
    def traverse(node):
        nodes.append(node)
        for child in node.children:
            traverse(child)
            
    traverse(root)
    return nodes

class MplTreeView(FigureCanvas):
    """封装的树状图查看器，继承自FigureCanvasQTAgg"""
    
    def __init__(self, parent=None, width=12, height=8, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        
        # 设置初始状态
        self.root = None
        self.ax = self.fig.add_subplot(111)
        self.ax.axis('off')
        
    def draw_dendrogram(self, newick_str=None, root=None):
        """绘制树状图，可以传入Newick字符串或已解析的根节点"""
        start_time = time.time()  # 添加开始时间记录
        
        # 清除之前的图形
        self.ax.clear()
        
        # 如果传入的是Newick字符串，先解析
        if newick_str is not None:
            parser = NewickParser(newick_str)
            self.root = parser.parse()
        elif root is not None:
            self.root = root
        else:
            print("没有提供Newick字符串或根节点")
            return
            
        if self.root is None:
            print("没有节点可绘制")
            return
            
        nodes = collect_nodes(self.root)
        
        if not nodes:
            print("没有节点可绘制")
            return
            
        # 提取坐标范围
        max_age = max(node.age for node in nodes)
        max_order = max(node.order for node in nodes)
        min_order = min(node.order for node in nodes)
        
        # 绘制节点和连线
        for node in nodes:
            # 绘制节点
            self.ax.plot(node.age, node.order, 'o', markersize=6, color='black')
            
            # 如果节点有名称且为叶节点，添加标签
            if node.name and node.is_leaf:
                self.ax.text(node.age, node.order, f" {node.name}", 
                           verticalalignment='center', fontsize=10)
            
            # 绘制到子节点的连线
            for child in node.children:
                # 先画垂直线：从父节点向下到子节点的水平位置
                self.ax.plot([node.age, node.age], [node.order, child.order], 
                           '-', color='black', linewidth=1)
                # 再画水平线：从父节点的水平位置向右到子节点
                self.ax.plot([node.age, child.age], [child.order, child.order], 
                           '-', color='black', linewidth=1)
        
        # 设置图形属性
        self.ax.set_xlim(-0.1, max_age * 1.2 if max_age > 0 else 1.0)
        self.ax.set_ylim(min_order - 0.5, max_order + 0.5)
        self.ax.set_xlabel('Distance')
        self.ax.set_ylabel('Nodes')
        self.ax.set_title('Phylogenetic Tree')
        self.ax.grid(True, alpha=0.3)
        
        # 反转y轴，使叶节点从上到下排列
        self.ax.invert_yaxis()
        
        self.fig.tight_layout()
        self.ax.axis('off')
        self.draw()
        
        end_time = time.time()  # 添加结束时间记录
        draw_time = end_time - start_time
        # print(f"绘图耗时: {draw_time:.4f} 秒")

def draw_dendrogram(root):
    """使用matplotlib绘制树状图"""
    start_time = time.time()  # 添加开始时间记录
    
    nodes = collect_nodes(root)
    
    if not nodes:
        print("没有节点可绘制")
        return
        
    # 提取坐标范围
    max_age = max(node.age for node in nodes)
    max_order = max(node.order for node in nodes)
    min_order = min(node.order for node in nodes)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 绘制节点和连线
    for node in nodes:
        # 绘制节点
        ax.plot(node.age, node.order, 'o', markersize=6, color='black')
        
        # 如果节点有名称且为叶节点，添加标签
        if node.name and node.is_leaf:
            ax.text(node.age, node.order, f" {node.name}", 
                   verticalalignment='center', fontsize=10)
        
        # 绘制到子节点的连线
        for child in node.children:
            # 先画垂直线：从父节点向下到子节点的水平位置
            ax.plot([node.age, node.age], [node.order, child.order], 
                   '-', color='black', linewidth=1)
            # 再画水平线：从父节点的水平位置向右到子节点
            ax.plot([node.age, child.age], [child.order, child.order], 
                   '-', color='black', linewidth=1)
    
    # 设置图形属性
    ax.set_xlim(-0.1, max_age * 1.2 if max_age > 0 else 1.0)
    ax.set_ylim(min_order - 0.5, max_order + 0.5)
    ax.set_xlabel('Distance')
    ax.set_ylabel('Nodes')
    ax.set_title('Phylogenetic Tree')
    ax.grid(True, alpha=0.3)
    
    # 反转y轴，使叶节点从上到下排列
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.axis('off')
    plt.show()
    
    end_time = time.time()  # 添加结束时间记录
    draw_time = end_time - start_time
    # print(f"绘图耗时: {draw_time:.4f} 秒")

def main():
    """主函数"""
    # 默认的Newick字符串
    newick_str = "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0,E:1.0):1.0);"   
    # 如果提供了命令行参数，则使用参数中的Newick字符串
    if len(sys.argv) > 1:
        newick_str = sys.argv[1]
    
    print(f"解析Newick字符串: {newick_str}")
    
    try:
        # 解析Newick字符串
        start_total = time.time()  # 添加总时间开始记录
        parser = NewickParser(newick_str)
        root = parser.parse()
        total_time = time.time() - start_total  # 计算总时间
        
        # # 打印节点信息
        # nodes = collect_nodes(root)
        # print("\n解析得到的节点:")
        # for i, node in enumerate(nodes):
        #     node_type = "leaf" if node.is_leaf else "internal"
        #     print(f"节点 {i}: {node.name or '[unnamed]'}, "
        #           f"age={node.age:.2f}, order={node.order:.2f}, "
        #           f"type={node_type}")
        
        # 绘制树状图
        print("\n正在绘制树状图...")
        draw_dendrogram(root)
        
        # 输出总耗时
        # print(f"\n总耗时: {total_time:.4f} 秒")
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    import sys
    main()