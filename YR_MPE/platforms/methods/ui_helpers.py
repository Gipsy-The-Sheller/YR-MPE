"""
UI Helpers Module
通用的UI组件创建、对话框等辅助方法
"""
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QDialog, QVBoxLayout
from PyQt5.QtCore import Qt
from typing import Optional, Tuple, List


class UIHelpers:
    """UI辅助工具类"""
    
    @staticmethod
    def select_workspace_folder(parent=None) -> Optional[str]:
        """
        选择workspace文件夹对话框
        
        Args:
            parent: 父窗口
            
        Returns:
            str: 选择的文件夹路径，如果取消则返回None
        """
        dialog = QFileDialog(parent)
        dialog.setFileMode(QFileDialog.Directory)
        dialog.setOption(QFileDialog.ShowDirsOnly, True)
        dialog.setWindowTitle("Select Workspace Folder")
        
        if dialog.exec_():
            selected_folder = dialog.selectedFiles()[0]
            return selected_folder
        return None
    
    @staticmethod
    def show_error_message(parent, title: str, message: str):
        """显示错误消息对话框"""
        QMessageBox.critical(parent, title, message)
    
    @staticmethod
    def show_info_message(parent, title: str, message: str):
        """显示信息消息对话框"""
        QMessageBox.information(parent, title, message)
    
    @staticmethod
    def create_modal_dialog(title: str, min_width: int = 800, min_height: int = 600) -> Tuple[QDialog, QVBoxLayout]:
        """
        创建模态对话框
        
        Args:
            title: 对话框标题
            min_width: 最小宽度
            min_height: 最小高度
            
        Returns:
            tuple: (dialog, layout)
        """
        dialog = QDialog()
        dialog.setWindowTitle(title)
        dialog.setMinimumSize(min_width, min_height)
        layout = QVBoxLayout()
        dialog.setLayout(layout)
        return dialog, layout