from PyQt5.QtCore import QObject, pyqtSignal
from typing import Any, Dict, List, Optional
import logging

class PluginAdapter(QObject):
    """
    插件适配器 - 将现有基于PyQt信号的插件适配到新的插件系统
    """
    
    # 信号定义
    plugin_finished = pyqtSignal(str, list, list)  # (plugin_name, output_files, reports)
    plugin_error = pyqtSignal(str, str)           # (plugin_name, error_message)
    plugin_progress = pyqtSignal(str, str)        # (plugin_name, progress_message)
    plugin_console_output = pyqtSignal(str, str, str)  # (plugin_name, message, message_type)
    
    def __init__(self, plugin_instance: Any, plugin_name: str):
        super().__init__()
        self.plugin_instance = plugin_instance
        self.plugin_name = plugin_name
        self.logger = logging.getLogger(__name__)
        
        # 连接现有插件的信号到适配器的信号
        self._connect_signals()
        
    def _connect_signals(self):
        """连接现有插件的信号到适配器"""
        try:
            # 如果插件有线程，连接线程信号
            if hasattr(self.plugin_instance, 'thread') and self.plugin_instance.thread:
                thread = self.plugin_instance.thread
                if hasattr(thread, 'finished'):
                    thread.finished.connect(self._on_plugin_finished)
                if hasattr(thread, 'error'):
                    thread.error.connect(self._on_plugin_error)
                if hasattr(thread, 'progress'):
                    thread.progress.connect(self._on_plugin_progress)
                if hasattr(thread, 'console_output'):
                    thread.console_output.connect(self._on_plugin_console_output)
                    
            # 如果插件本身有信号，也连接它们
            if hasattr(self.plugin_instance, 'finished'):
                self.plugin_instance.finished.connect(self._on_plugin_finished)
            if hasattr(self.plugin_instance, 'error'):
                self.plugin_instance.error.connect(self._on_plugin_error)
            if hasattr(self.plugin_instance, 'progress'):
                self.plugin_instance.progress.connect(self._on_plugin_progress)
            if hasattr(self.plugin_instance, 'console_output'):
                self.plugin_instance.console_output.connect(self._on_plugin_console_output)
                
        except Exception as e:
            self.logger.warning(f"Could not connect all signals for plugin {self.plugin_name}: {e}")
            
    def _on_plugin_finished(self, output_files: List[str], reports: List[str]):
        """处理插件完成信号"""
        self.plugin_finished.emit(self.plugin_name, output_files, reports)
        
    def _on_plugin_error(self, error_message: str):
        """处理插件错误信号"""
        self.plugin_error.emit(self.plugin_name, error_message)
        
    def _on_plugin_progress(self, progress_message: str):
        """处理插件进度信号"""
        self.plugin_progress.emit(self.plugin_name, progress_message)
        
    def _on_plugin_console_output(self, message: str, message_type: str):
        """处理插件控制台输出信号"""
        self.plugin_console_output.emit(self.plugin_name, message, message_type)
        
    def execute(self):
        """执行插件"""
        try:
            if hasattr(self.plugin_instance, 'run_analysis'):
                self.plugin_instance.run_analysis()
            elif hasattr(self.plugin_instance, 'start_process'):
                self.plugin_instance.start_process()
            else:
                # 对于没有明确执行方法的插件，可能需要手动触发
                self.logger.warning(f"Plugin {self.plugin_name} has no standard execute method")
        except Exception as e:
            self.plugin_error.emit(self.plugin_name, f"Execution failed: {str(e)}")