"""
Dataset Plugin Entry
Dataset功能的插件入口点
"""
from PyQt5.QtWidgets import QToolButton
from PyQt5.QtCore import pyqtSignal


class DatasetPlugin(QToolButton):
    """Dataset插件主控件"""
    
    # 信号定义
    dataset_created = pyqtSignal(list)  # 创建的dataset列表
    view_sequence_signal = pyqtSignal(list)  # 查看序列信号
    
    def __init__(self, parent=None):
        super().__init__(parent)
        # 移除文字，只显示图标
        self.setCheckable(True)
        self.clicked.connect(self.on_clicked)
        self.double_clicked = False
        
    def mouseDoubleClickEvent(self, event):
        """双击事件处理"""
        self.double_clicked = True
        self.open_dataset_manager()
        super().mouseDoubleClickEvent(event)
        
    def mousePressEvent(self, event):
        """单击事件处理"""
        self.double_clicked = False
        super().mousePressEvent(event)
        
    def mouseReleaseEvent(self, event):
        """鼠标释放事件处理"""
        if not self.double_clicked:
            # 单击：选中Dataset功能模式
            self.setChecked(True)
        super().mouseReleaseEvent(event)
        
    def on_clicked(self):
        """点击事件处理（单击）"""
        if not self.double_clicked:
            # 单击已由mouseReleaseEvent处理
            pass
            
    def open_dataset_manager(self):
        """打开Dataset管理对话框"""
        from .dataset_manager_dialog import DatasetManagerDialog
        dialog = DatasetManagerDialog(self.view_sequence_signal)
        dialog.exec_()


class DatasetPluginEntry:
    """Dataset插件入口类"""
    
    def run(self):
        """运行Dataset插件"""
        return DatasetPlugin()