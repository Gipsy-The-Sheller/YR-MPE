from PyQt5.QtCore import QThread, pyqtSignal
import subprocess
import os
import tempfile
import json
from typing import List, Dict, Any, Optional
import threading
import queue
import time
import platform


class BaseProcessThread(QThread):
    """
    YR-MPEA基础进程线程类
    
    这是所有需要执行外部命令的插件线程的基类，提供了通用的命令执行功能：
    - 命令执行和进度报告
    - 输出文件和报告管理
    - 错误处理和日志记录
    - 控制台消息输出
    """
    
    # 信号定义
    progress = pyqtSignal(str)                      # 进度信息
    finished = pyqtSignal(list, list)              # 完成信号 (输出文件列表, HTML报告列表)
    error = pyqtSignal(str)                        # 错误信号
    console_output = pyqtSignal(str, str)          # 控制台输出 (消息, 消息类型)
    
    def __init__(self, tool_path: str, input_files: List[str], parameters: List[str], 
                 imported_files: Optional[List[str]] = None):
        """
        初始化基础进程线程
        
        Args:
            tool_path (str): 工具可执行文件路径
            input_files (List[str]): 输入文件路径列表
            parameters (List[str]): 命令行参数列表
            imported_files (Optional[List[str]]): 导入的文件列表（用于报告生成）
        """
        super().__init__()
        self.tool_path = tool_path
        self.input_files = input_files if isinstance(input_files, list) else [input_files]
        self.parameters = parameters
        self.imported_files = imported_files or []
        self.is_running = False
        self._process = None
        
    def run(self):
        """
        执行命令的主方法
        子类应该重写此方法以实现特定的命令执行逻辑
        """
        try:
            self.is_running = True
            self.execute_commands()
        except Exception as e:
            self.error.emit(str(e))
        finally:
            self.is_running = False
    
    def execute_commands(self):
        """
        执行命令的核心方法
        子类应重写此方法来实现具体的命令执行逻辑
        """
        # 这里应该被子类重写
        self.error.emit("execute_commands method not implemented in subclass")
    
    def execute_command(self, cmd: List[str], cwd: Optional[str] = None) -> subprocess.CompletedProcess:
        """
        执行单个命令
        
        Args:
            cmd (List[str]): 命令和参数列表
            cwd (Optional[str]): 工作目录
            
        Returns:
            subprocess.CompletedProcess: 命令执行结果
        """
        self.console_output.emit(" ".join(cmd), "command")
        
        # 准备启动参数，解决Windows命令行窗口问题
        startupinfo = None
        if platform.system().lower() == "windows":
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            startupinfo.wShowWindow = subprocess.SW_HIDE

        try:
            # 使用Popen而非run，允许实时读取输出
            self._process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # 将stderr合并到stdout，以便实时捕获所有输出
                universal_newlines=True,
                bufsize=1,  # 行缓冲
                cwd=cwd,
                startupinfo=startupinfo  # Windows隐藏控制台窗口
            )

            # 实时读取输出
            stdout_lines = []
            while True:
                output = self._process.stdout.readline()
                
                if output == '' and self._process.poll() is not None:
                    # 进程结束且没有更多输出
                    break
                
                if output:
                    output = output.rstrip()  # 保留右侧空白字符可能会有用
                    stdout_lines.append(output + '\n')  # 重新加上换行符
                    # 发送输出到GUI
                    self.console_output.emit(output, "output")

            # 获取返回码
            return_code = self._process.poll()
            
            # 重构输出字符串以匹配subprocess.CompletedProcess的行为
            stdout_str = ''.join(stdout_lines)
            
            # 创建一个类似subprocess.CompletedProcess的对象
            result = subprocess.CompletedProcess(
                args=cmd,
                returncode=return_code,
                stdout=stdout_str,
                stderr=""  # 因为我们把stderr合并到了stdout中
            )
            
            return result

        except Exception as e:
            # 发生异常时，发出错误信号并抛出异常
            error_msg = f"执行命令时发生错误: {str(e)}"
            self.console_output.emit(error_msg, "error")
            raise e

    def create_temp_file(self, suffix: str = '') -> str:
        """
        创建临时文件
        
        Args:
            suffix (str): 文件后缀
            
        Returns:
            str: 临时文件路径
        """
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        temp_file.close()
        return temp_file.name
    
    def create_output_dir(self, base_path: str, dir_name: str) -> str:
        """
        创建输出目录
        
        Args:
            base_path (str): 基础路径
            dir_name (str): 目录名
            
        Returns:
            str: 输出目录路径
        """
        output_dir = os.path.join(base_path, dir_name)
        os.makedirs(output_dir, exist_ok=True)
        return output_dir
    
    def get_tool_name(self) -> str:
        """
        获取工具名称
        子类应重写此方法返回具体的工具名称
        
        Returns:
            str: 工具名称
        """
        return "Base Tool"
    
    def generate_simple_report(self, input_file: str, output_file: str, 
                              cmd: List[str], stdout: str = "", stderr: str = "") -> str:
        """
        生成简单的HTML报告
        
        Args:
            input_file (str): 输入文件路径
            output_file (str): 输出文件路径
            cmd (List[str]): 执行的命令
            stdout (str): 标准输出
            stderr (str): 标准错误输出
            
        Returns:
            str: HTML报告文件路径
        """
        try:
            base_name = os.path.splitext(os.path.basename(output_file))[0]
            html_file = f"{base_name}_report.html"
            
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>{self.get_tool_name()} Analysis Report</title>
                <style>
                    body {{ font-family: Arial, sans-serif; margin: 20px; }}
                    h1 {{ color: #2c3e50; }}
                    .info {{ background-color: #e8f4f8; padding: 10px; border-radius: 5px; }}
                    pre {{ background-color: #f8f9fa; padding: 10px; border-radius: 5px; }}
                </style>
            </head>
            <body>
                <h1>{self.get_tool_name()} Analysis Report</h1>
                <div class="info">
                    <p><strong>Input file:</strong> {input_file}</p>
                    <p><strong>Output file:</strong> {output_file}</p>
                    <p><strong>Command:</strong> {' '.join(cmd)}</p>
                </div>
                <h2>Output</h2>
                <pre>{stdout}</pre>
                {('<h2>Errors</h2><pre>' + stderr + '</pre>') if stderr else ''}
            </body>
            </html>
            """
            
            with open(html_file, 'w') as f:
                f.write(html_content)
                
            return html_file
        except Exception as e:
            self.console_output.emit(f"Failed to generate report: {e}", "error")
            return ""
    
    def stop(self):
        """
        停止线程执行
        """
        self.is_running = False
        # 如果有正在运行的进程，尝试终止它
        if self._process and self._process.poll() is None:
            try:
                # Windows下使用taskkill
                if platform.system().lower() == "windows":
                    subprocess.run(['taskkill', '/F', '/T', '/PID', str(self._process.pid)], 
                                  stdout=subprocess.DEVNULL, 
                                  stderr=subprocess.DEVNULL)
                else:
                    # Unix-like系统下发送SIGTERM信号
                    import signal
                    try:
                        os.killpg(os.getpgid(self._process.pid), signal.SIGTERM)
                    except ProcessLookupError:
                        # 进程可能已经结束
                        pass
            except Exception as e:
                print(f"终止进程时出错: {e}")