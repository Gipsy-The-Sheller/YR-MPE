from PyQt5.QtCore import QObject, pyqtSignal
from typing import Dict, Any, List, Optional
from .plugin_manager import PluginManager
from .workspace_manager import WorkspaceManager
import logging

class PluginExecutor(QObject):
    """
    插件执行器 - 管理插件的执行流程、状态跟踪和结果处理
    """
    
    # 信号定义
    execution_started = pyqtSignal(str)           # (plugin_name)
    execution_finished = pyqtSignal(str, dict)    # (plugin_name, result)
    execution_error = pyqtSignal(str, str)        # (plugin_name, error_message)
    execution_progress = pyqtSignal(str, str)     # (plugin_name, progress_message)
    console_output = pyqtSignal(str, str, str)    # (plugin_name, message, message_type)
    
    def __init__(self, plugin_manager: PluginManager, workspace_manager: WorkspaceManager):
        super().__init__()
        self.plugin_manager = plugin_manager
        self.workspace_manager = workspace_manager
        self.logger = logging.getLogger(__name__)
        self.active_executions: Dict[str, Any] = {}
        
    def execute_plugin(self, plugin_name: str, input_data: Any = None, 
                      parameters: Dict[str, Any] = None) -> bool:
        """
        执行插件
        
        Args:
            plugin_name: 插件名称
            input_data: 输入数据
            parameters: 插件参数
            
        Returns:
            是否成功启动执行
        """
        try:
            self.execution_started.emit(plugin_name)
            
            # 创建并执行插件适配器
            adapter = self.plugin_manager.execute_plugin(
                plugin_name, 
                input_data=input_data, 
                parameters=parameters
            )
            
            # 连接适配器信号到执行器信号
            adapter.plugin_finished.connect(self._on_plugin_finished)
            adapter.plugin_error.connect(self._on_plugin_error)
            adapter.plugin_progress.connect(self._on_plugin_progress)
            adapter.plugin_console_output.connect(self._on_console_output)
            
            # 存储活动执行
            self.active_executions[plugin_name] = adapter
            
            # 执行插件
            adapter.execute()
            
            return True
            
        except Exception as e:
            error_msg = f"Failed to start plugin {plugin_name}: {str(e)}"
            self.logger.error(error_msg)
            self.execution_error.emit(plugin_name, error_msg)
            return False
            
    def _on_plugin_finished(self, plugin_name: str, output_files: List[str], reports: List[str]):
        """处理插件完成"""
        result = {
            "output_files": output_files,
            "reports": reports,
            "status": "success"
        }
        
        # 更新工作区状态
        if output_files:
            self.workspace_manager.add_files_to_history(output_files)
            
        self.execution_finished.emit(plugin_name, result)
        
        # 清理活动执行
        if plugin_name in self.active_executions:
            del self.active_executions[plugin_name]
            
    def _on_plugin_error(self, plugin_name: str, error_message: str):
        """处理插件错误"""
        self.execution_error.emit(plugin_name, error_message)
        
        # 清理活动执行
        if plugin_name in self.active_executions:
            del self.active_executions[plugin_name]
            
    def _on_plugin_progress(self, plugin_name: str, progress_message: str):
        """处理插件进度"""
        self.execution_progress.emit(plugin_name, progress_message)
        
    def _on_console_output(self, plugin_name: str, message: str, message_type: str):
        """处理控制台输出"""
        self.console_output.emit(plugin_name, message, message_type)
        
    def cancel_execution(self, plugin_name: str):
        """取消插件执行"""
        if plugin_name in self.active_executions:
            adapter = self.active_executions[plugin_name]
            # TODO: 实现取消逻辑
            if hasattr(adapter.plugin_instance, 'cancel'):
                adapter.plugin_instance.cancel()
            elif hasattr(adapter.plugin_instance, 'thread') and adapter.plugin_instance.thread:
                adapter.plugin_instance.thread.terminate()
                
            del self.active_executions[plugin_name]