# console_plugin.py
#
# Copyright (c) 2025 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# YR-MPE program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, QDoubleSpinBox)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon
import tempfile
import os
import subprocess
import platform
import json
import shutil
from typing import List, Optional


class YRMPEConsolePlugin(QWidget):
    """YR-MPE Console插件：启动一个新终端窗口并配置工具路径"""
    
    def __init__(self):
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.init_ui()
    
    def init_ui(self):
        """初始化用户界面"""
        layout = QVBoxLayout()
        
        # 插件标题
        title_label = QLabel("YR-MPE Console")
        title_label.setStyleSheet("font-size: 18px; font-weight: bold; margin-bottom: 10px;")
        layout.addWidget(title_label)
        
        # 描述信息
        desc_label = QLabel(
            "Start up a new shell to use all YR-MPE's phylogenetic tools in a console.\n"
        )
        desc_label.setWordWrap(True)
        layout.addWidget(desc_label)
        
        # 添加间距
        layout.addSpacing(20)
        
        # 启动控制台按钮
        self.start_console_btn = QPushButton("Startup a new shell")
        self.start_console_btn.clicked.connect(self.start_console)
        layout.addWidget(self.start_console_btn)
        
        # 添加创建PhyloBOT项目按钮
        self.create_phylobot_btn = QPushButton("Create new PhyloBOT project")
        self.create_phylobot_btn.clicked.connect(self.create_phylobot_project)
        layout.addWidget(self.create_phylobot_btn)
        
        # 添加间距
        layout.addSpacing(10)
        
        # 显示将要添加的路径
        paths_group = QGroupBox("Tool PATHs")
        paths_layout = QVBoxLayout()
        
        # 读取配置文件中的路径
        self.tool_paths, self.tool_commands = self.load_tool_paths_and_commands()
        paths_text = QTextEdit()
        paths_text.setReadOnly(True)
        paths_text.setMaximumHeight(150)
        
        paths_str = "\n".join(self.tool_paths) if self.tool_paths else "Failed to load or empty..."
        paths_text.setPlainText(paths_str)
        
        paths_layout.addWidget(paths_text)
        paths_group.setLayout(paths_layout)
        layout.addWidget(paths_group)
        
        layout.addStretch()  # 添加弹性空间
        self.setLayout(layout)
    
    def load_tool_paths_and_commands(self):
        """从config.json加载工具路径和命令名称"""
        config_path = os.path.join(self.plugin_path, "config.json")
        
        if not os.path.exists(config_path):
            print(f"配置文件不存在: {config_path}")
            return [], {}
        
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
            
            # 提取所有路径并转换为绝对路径
            paths = []
            commands = {}
            
            for item in config_data:
                if "path" in item and "name" in item:
                    # 获取路径的目录部分（去掉文件名）
                    abs_path = os.path.join(self.plugin_path, '.'+item["path"])
                    dir_path = os.path.dirname(abs_path)
                    if dir_path not in paths and os.path.exists(dir_path):
                        paths.append(dir_path)
                    
                    # 提取命令名称
                    path_parts = item["path"].split('/')
                    command_name = path_parts[-1]  # 获取文件名部分
                    # 移除扩展名（如.exe, .bat等）
                    command_name = os.path.splitext(command_name)[0]
                    commands[item["name"]] = command_name

            return paths, commands
        except Exception as e:
            print(f"Failed to load config file: {e}")
            return [], {}
    
    def load_tool_paths(self):
        """从config.json加载工具路径（保持向后兼容）"""
        paths, _ = self.load_tool_paths_and_commands()
        return paths
    
    def create_phylobot_project(self):
        """创建新的PhyloBOT项目"""
        # 让用户选择项目目录
        project_dir = QFileDialog.getExistingDirectory(
            self, 
            "Select Directory for New PhyloBOT Project",
            os.path.expanduser("~")
        )
        
        if not project_dir:
            return  # 用户取消了选择
        
        try:
            # 确定模板目录
            template_dir = os.path.join(
                self.plugin_path, 
                "bin", 
                "iflow", 
                "template_phylobot_project"
            )
            
            if not os.path.exists(template_dir):
                QMessageBox.critical(
                    self, 
                    "Error", 
                    f"Template directory does not exist: {template_dir}"
                )
                return
            
            # 获取YR-MPE插件路径，用于替换模板中的路径
            yr_mpe_path = os.path.dirname(self.plugin_path)  # 获取YR-MPE根目录
            
            # 目标项目目录
            target_dir = os.path.join(project_dir, "phylobot_project")
            
            # 如果目标目录已存在，询问用户是否覆盖
            if os.path.exists(target_dir):
                reply = QMessageBox.question(
                    self, 
                    "Confirm Overwrite", 
                    f"Directory {target_dir} already exists. Do you want to overwrite it?",
                    QMessageBox.Yes | QMessageBox.No, 
                    QMessageBox.No
                )
                if reply == QMessageBox.No:
                    return
            
            # 复制模板目录
            if os.path.exists(target_dir):
                shutil.rmtree(target_dir)
            
            shutil.copytree(template_dir, target_dir)
            
            # 替换IFLOW.md中的模板路径
            iflow_md_path = os.path.join(target_dir, "IFLOW.md")
            if os.path.exists(iflow_md_path):
                with open(iflow_md_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # 替换模板变量
                content = content.replace('{YR_MPE_PATH}', yr_mpe_path.replace('\\', '/'))
                
                with open(iflow_md_path, 'w', encoding='utf-8') as f:
                    f.write(content)
            
            # 根据操作系统创建启动脚本
            system = platform.system()
            if system == "Windows":
                # 创建PowerShell启动脚本
                run_script_path = os.path.join(target_dir, "run_phylobot.ps1")
                
                # 定义PowerShell脚本内容，分别构建各部分
                node_check_part = (
                    '$nodeInstalled = Get-Command node -ErrorAction SilentlyContinue\n'
                    'if (-not $nodeInstalled) {\n'
                    '    Write-Host "Error: Node.js is not installed. Please install Node.js to continue." -ForegroundColor Red\n'
                    '    pause\n'
                    '    exit 1\n'
                    '}'
                )
                
                # 路径设置部分
                yr_mpe_path_formatted = yr_mpe_path.replace("/", "\\\\")
                setup_paths_part = (
                    f"# Set up paths\n"
                    f"$yrMpePath = \"{yr_mpe_path_formatted}\"\n"
                    f"$manualsPath = \"$yrMpePath\\\\YR_MPE\\\\bin\\\\iflow\\\\template_phylobot_project\\\\manuals\"\n"
                    f"\n"
                    f"# Copy manuals to project directory\n"
                    f"Copy-Item -Path \"$manualsPath\\\\*\" -Destination \".\\\\manuals\\\\\" -Recurse -Force"
                )
                
                # 启动iflow部分 - 现在在项目目录中启动iflow
                start_iflow_part = (
                    '# Start iflow\n'
                    f'Start-Process -FilePath "powershell" -ArgumentList "-Command", "cd \'{target_dir}\'; & \'$yrMpePath\\\\YR_MPE\\\\bin\\\\iflow\\\\iflow.ps1\'"'
                )
                
                script_content = f"  {node_check_part}\n\n{setup_paths_part}\n\n{start_iflow_part}"
                
                with open(run_script_path, 'w', encoding='utf-8') as f:
                    f.write(script_content)
            else:
                # 创建shell启动脚本
                run_script_path = os.path.join(target_dir, "run_phylobot.sh")
                
                # 定义shell脚本内容，分别构建各部分
                node_check_part = (
                    'if ! command -v node &> /dev/null; then\n'
                    '    echo "Error: Node.js is not installed. Please install Node.js to continue."\n'
                    '    read -p "Press any key to exit..." -n1 -s\n'
                    '    exit 1\n'
                    'fi'
                )
                
                # 路径设置部分
                setup_paths_part = (
                    f"# Set up paths\n"
                    f'YR_MPE_PATH="{yr_mpe_path}"\n'
                    f'MANUALS_PATH="$YR_MPE_PATH/YR_MPE/bin/iflow/template_phylobot_project/manuals"\n'
                    f"\n"
                    f"# Copy manuals to project directory\n"
                    f'cp -rf "$MANUALS_PATH/"* ./manuals/ 2>/dev/null || true'
                )
                
                # 启动iflow部分 - 现在在项目目录中启动iflow
                start_iflow_part = (
                    f'# Start iflow\n'
                    f'cd "{target_dir}"\n'
                    '"$YR_MPE_PATH/YR_MPE/bin/iflow/iflow" "$@"'
                )
                
                script_content = f"#!/bin/bash\n\n{node_check_part}\n\n{setup_paths_part}\n\n{start_iflow_part}"
                
                with open(run_script_path, 'w', encoding='utf-8') as f:
                    f.write(script_content)
                
                # Make the script executable
                os.chmod(run_script_path, 0o755)
            
            QMessageBox.information(
                self, 
                "Success", 
                f"PhyloBOT project created successfully in:\n{target_dir}\n\n" +
                "Run the project using the generated run_phylobot script."
            )
            
            # 启动PhyloBOT项目
            if system == "Windows":
                os.system(f'start powershell -ExecutionPolicy Bypass -File "{run_script_path}"')
            else:
                os.system(f'cd "{target_dir}" && ./run_phylobot.sh')
        
        except Exception as e:
            QMessageBox.critical(
                self, 
                "Error", 
                f"Failed to create PhyloBOT project: {str(e)}"
            )
    
    def start_console(self):
        """启动新的终端窗口并设置PATH"""
        try:
            # 获取当前系统
            system = platform.system()
            
            # 获取工具路径和命令
            tool_paths = self.load_tool_paths()
            _, tool_commands = self.load_tool_paths_and_commands()
            
            if not tool_paths:
                QMessageBox.warning(self, "警告", "无法加载工具路径，请检查配置文件！")
                return
            
            # 构建新的PATH环境变量
            # 首先获取当前PATH
            current_path = os.environ.get('PATH', '')
            
            # 将工具路径添加到当前PATH前面（优先级更高）
            new_paths = os.pathsep.join(tool_paths)
            final_path = new_paths + os.pathsep + current_path
            
            # 构建工具列表字符串 (使用PowerShell的转义字符)
            tools_list = "Tool name           Command`n"
            tools_list += "-" * 30 + "`n"
            for name, cmd in tool_commands.items():
                tools_list += f"{name if len(name) <= 17 else name[:18]+'...'}{' '*(20-len(name))}{cmd}`n"
            
            if system == "Windows":
                # 创建一个临时PowerShell脚本来设置环境变量并启动PowerShell
                ps_script_content = f'''
# 设置新的PATH环境变量
$env:PATH = "{final_path}"

Write-Host "Console for YR-MPE" -ForegroundColor Green
Write-Host ""
Write-Host "{tools_list}" -ForegroundColor Yellow
Write-Host ""
# 启动一个新的PowerShell会话，保持窗口开启
powershell
'''
                # 将PowerShell脚本写入临时文件
                try:
                    with open("YR_MPE_powershell_script.ps1", "w", encoding="utf-8") as temp_ps_file:
                        temp_ps_file.write(ps_script_content)
                except:
                    # ps1 is now occupied and does not need to be refreshed.
                    pass
                
                # 使用PowerShell执行临时脚本文件
                command = f'start powershell -ExecutionPolicy Bypass -File "YR_MPE_powershell_script.ps1"'
                os.system(command)
            elif system == "Darwin":  # macOS
                # 在macOS上启动终端
                # 创建一个临时脚本来设置环境变量并显示工具列表
                # 将工具列表转换为AppleScript兼容的字符串
                as_tools_list = tools_list.replace('"', '\\"')
                script = f'''
                set newPATH to "{final_path}"
                set toolsList to "{as_tools_list}"
                tell application "Terminal"
                    do script "export PATH=\\"" & quoted form of newPATH & "\\"; \\
                    echo \\"Console for YR-MPE\\"; \\
                    echo \\"\\\\"; \\
                    echo \\"" & toolsList & "\\"; \\
                    echo \\"\\\\"; \\
                    exec zsh"
                    activate
                end tell
                '''
                os.system(f"osascript -e '{script}'")
            else:  # Linux and other Unix-like systems
                # 在Linux上启动bash终端
                escaped_tools_list = tools_list.replace('"', '\\"').replace('$', '\\$')
                command = f'gnome-terminal -- bash -c "export PATH=\\"{final_path}\\":$PATH; echo \'Console for YR-MPE\'; echo \\"\\"; echo -e \\"{escaped_tools_list}\\"; echo \\"\\"; exec bash" &'
                os.system(command)
            
            QMessageBox.information(self, "提示", "YR-MPE终端已启动！\n工具路径已添加到PATH环境变量中。")
        
        except Exception as e:
            QMessageBox.critical(self, "错误", f"启动终端时发生错误：\n{str(e)}")


class ConsolePluginEntry:
    """YR-MPE Console插件入口点"""
    
    def run(self):
        return YRMPEConsolePlugin()