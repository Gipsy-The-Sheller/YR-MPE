"""
QDialog线程管理器 - 解决QDialog阻塞主线程的问题
"""
from PyQt5.QtWidgets import QDialog, QApplication
from PyQt5.QtCore import QThread, pyqtSignal, QObject
import threading
import sys


class DialogThread(QThread):
    """用于在独立线程中运行QDialog的线程类"""
    
    dialog_finished = pyqtSignal(int)  # 对话框关闭时发送结果信号
    
    def __init__(self, dialog_class, *args, **kwargs):
        super().__init__()
        self.dialog_class = dialog_class
        self.dialog_args = args
        self.dialog_kwargs = kwargs
        self.dialog = None
        
    def run(self):
        """在线程中创建和显示对话框"""
        try:
            # 创建对话框实例
            self.dialog = self.dialog_class(*self.dialog_args, **self.dialog_kwargs)
            
            # 显示对话框（非模态）
            self.dialog.show()
            
            # 等待对话框关闭
            while self.dialog.isVisible():
                QApplication.processEvents()
                self.msleep(10)  # 避免过度占用CPU
                
            # 发送结果信号
            result = self.dialog.result() if hasattr(self.dialog, 'result') else QDialog.Accepted
            self.dialog_finished.emit(result)
            
        except Exception as e:
            print(f"Dialog thread error: {e}")
            self.dialog_finished.emit(QDialog.Rejected)


class NonBlockingDialogManager(QObject):
    """非阻塞对话框管理器"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.active_threads = []
        
    def show_dialog(self, dialog_class, *args, **kwargs):
        """
        在独立线程中显示对话框，不阻塞主线程
        
        Args:
            dialog_class: 对话框类
            *args, **kwargs: 传递给对话框构造函数的参数
            
        Returns:
            DialogThread: 对话框线程对象
        """
        # 创建对话框线程
        dialog_thread = DialogThread(dialog_class, *args, **kwargs)
        
        # 连接信号
        dialog_thread.dialog_finished.connect(
            lambda result: self._on_dialog_finished(dialog_thread, result)
        )
        
        # 启动线程
        dialog_thread.start()
        
        # 记录活跃线程
        self.active_threads.append(dialog_thread)
        
        return dialog_thread
    
    def _on_dialog_finished(self, thread, result):
        """对话框关闭时的回调"""
        # 从活跃线程列表中移除
        if thread in self.active_threads:
            self.active_threads.remove(thread)
        
        # 清理线程
        thread.quit()
        thread.wait()
        thread.deleteLater()
    
    def close_all_dialogs(self):
        """关闭所有活跃的对话框"""
        for thread in self.active_threads[:]:  # 使用切片复制避免修改列表时出错
            if thread.dialog and thread.dialog.isVisible():
                thread.dialog.close()
    
    def get_active_dialog_count(self):
        """获取当前活跃对话框数量"""
        return len(self.active_threads)


# 全局对话框管理器实例
dialog_manager = NonBlockingDialogManager()


def show_dialog_non_blocking(dialog_class, *args, **kwargs):
    """
    便捷函数：在独立线程中显示对话框
    
    Args:
        dialog_class: 对话框类
        *args, **kwargs: 传递给对话框构造函数的参数
        
    Returns:
        DialogThread: 对话框线程对象
    """
    return dialog_manager.show_dialog(dialog_class, *args, **kwargs)


class NonBlockingDialog(QDialog):
    """非阻塞对话框基类"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setModal(False)  # 设置为非模态
        
    def show_non_blocking(self):
        """以非阻塞方式显示对话框"""
        return show_dialog_non_blocking(self.__class__, self.parent())
    
    @classmethod
    def create_and_show(cls, parent=None, *args, **kwargs):
        """创建并显示非阻塞对话框"""
        dialog = cls(parent, *args, **kwargs)
        return show_dialog_non_blocking(cls, parent, *args, **kwargs)

