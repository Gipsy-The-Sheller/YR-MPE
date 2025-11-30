import sys
import os
import json
import subprocess
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                             QPushButton, QLineEdit, QFileDialog, QMessageBox, 
                             QTableWidget, QTableWidgetItem, QHeaderView, QDialog,
                             QTextEdit, QGroupBox)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont


class LicenseDialog(QDialog):
    """显示许可证信息的对话框"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("License Information")
        self.resize(600, 400)
        
        layout = QVBoxLayout()
        
        # 标题
        title_label = QLabel("GNU General Public License v3.0")
        title_font = QFont()
        title_font.setPointSize(16)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(title_label)
        
        # 许可证文本
        license_text = QTextEdit()
        license_text.setReadOnly(True)
        license_text.setText("""GNU GENERAL PUBLIC LICENSE
Version 3, 29 June 2007

Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
Everyone is permitted to copy and distribute verbatim copies
of this license document, but changing it is not allowed.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.""")
        layout.addWidget(license_text)
        
        # 关闭按钮
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)
        
        self.setLayout(layout)


class ManualConfigDialog(QDialog):
    """手动配置软件路径的对话框"""
    def __init__(self, tools_data, parent=None):
        super().__init__(parent)
        self.tools_data = tools_data
        self.tool_paths = {}
        self.setWindowTitle("Manual Configuration")
        self.resize(700, 400)
        
        layout = QVBoxLayout()
        
        # 表格
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Software Name", "Path", "Select"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Fixed)
        self.table.setColumnWidth(2, 100)
        
        # 填充表格数据
        self.table.setRowCount(len(self.tools_data))
        for i, tool in enumerate(self.tools_data):
            # 软件名
            name_item = QTableWidgetItem(tool["name"])
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(i, 0, name_item)
            
            # 路径
            path_item = QTableWidgetItem(tool.get("path", ""))
            self.table.setItem(i, 1, path_item)
            
            # 选择按钮
            select_button = QPushButton("Browse")
            select_button.clicked.connect(lambda checked, row=i: self.browse_path(row))
            self.table.setCellWidget(i, 2, select_button)
            
            # 保存路径引用
            self.tool_paths[i] = path_item
        
        layout.addWidget(self.table)
        
        # 按钮
        button_layout = QHBoxLayout()
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(ok_button)
        button_layout.addWidget(cancel_button)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def browse_path(self, row):
        """浏览选择路径"""
        path, _ = QFileDialog.getOpenFileName(self, f"Select {self.tool_paths[row].text()} executable")
        if path:
            self.tool_paths[row].setText(path)
    
    def get_configured_tools(self):
        """获取配置后的工具数据"""
        configured_tools = []
        for i, tool in enumerate(self.tools_data):
            new_tool = tool.copy()
            new_tool["path"] = self.tool_paths[i].text()
            configured_tools.append(new_tool)
        return configured_tools


class Initializer(QWidget):
    """YR-MPE初始化配置界面"""
    
    def __init__(self):
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.tools_data = []
        self.configured_tools = []
        self.init_ui()
        self.load_templates()
    
    def init_ui(self):
        """初始化用户界面"""
        self.setWindowTitle("YR-MPE Initializer")
        self.setMinimumSize(600, 500)
        
        main_layout = QVBoxLayout()
        
        # 标题
        title_label = QLabel("Welcome to YR-MPE Initialize Wizard.")
        title_font = QFont()
        # title_font.setPointSize(18)
        title_font.setBold(True)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)
        
        # 许可证查看
        license_group = QGroupBox("License")
        license_layout = QHBoxLayout()
        license_label = QLabel("View the GPL-3.0 license:")
        self.license_button = QPushButton("View License")
        self.license_button.clicked.connect(self.show_license)
        license_layout.addWidget(license_label)
        license_layout.addWidget(self.license_button)
        license_group.setLayout(license_layout)
        main_layout.addWidget(license_group)
        
        # 安装路径设置
        path_group = QGroupBox("Installation Path")
        path_layout = QHBoxLayout()
        path_layout.addWidget(QLabel("YR-MPE Path:"))
        self.path_edit = QLineEdit("YR-MPE")
        self.path_button = QPushButton("Browse")
        self.path_button.clicked.connect(self.browse_path)
        path_layout.addWidget(self.path_edit)
        path_layout.addWidget(self.path_button)
        path_group.setLayout(path_layout)
        main_layout.addWidget(path_group)
        
        # 配置工具
        config_group = QGroupBox("Tool Configuration")
        config_layout = QVBoxLayout()
        
        # 自动配置按钮
        auto_config_layout = QHBoxLayout()
        auto_config_layout.addWidget(QLabel("Automatically detect tools in PATH:"))
        self.auto_config_button = QPushButton("Auto Config")
        self.auto_config_button.clicked.connect(self.auto_config)
        auto_config_layout.addWidget(self.auto_config_button)
        config_layout.addLayout(auto_config_layout)
        
        # 手动配置按钮
        manual_config_layout = QHBoxLayout()
        manual_config_layout.addWidget(QLabel("Manually configure tool paths:"))
        self.manual_config_button = QPushButton("Manual Config")
        self.manual_config_button.clicked.connect(self.manual_config)
        manual_config_layout.addWidget(self.manual_config_button)
        config_layout.addLayout(manual_config_layout)
        
        config_group.setLayout(config_layout)
        main_layout.addWidget(config_group)
        
        # 完成安装按钮
        self.finish_button = QPushButton("Finish && Install")
        self.finish_button.clicked.connect(self.finish_install)
        main_layout.addWidget(self.finish_button)
        
        self.setLayout(main_layout)
    
    def load_templates(self):
        """加载模板文件"""
        # 加载配置模板
        config_template_path = os.path.join(self.plugin_path, "config.template.json")
        if os.path.exists(config_template_path):
            with open(config_template_path, "r", encoding='utf-8') as f:
                self.tools_data = json.load(f)
    
    def show_license(self):
        """显示许可证信息"""
        dialog = LicenseDialog(self)
        dialog.exec_()
    
    def browse_path(self):
        """浏览选择安装路径"""
        path = QFileDialog.getExistingDirectory(self, "Select Installation Path", self.plugin_path)
        if path:
            self.path_edit.setText(os.path.relpath(path, self.plugin_path))
    
    def auto_config(self):
        """自动配置工具路径"""
        if not self.tools_data:
            QMessageBox.warning(self, "Warning", "No tools data loaded.")
            return
        
        not_found_tools = []
        self.configured_tools = []
        
        for tool in self.tools_data:
            found = False
            path = ""
            
            # 检查每个可能的命令
            for cmd in tool.get("cmd", []):
                try:
                    # 尝试在系统PATH中查找命令
                    if os.name == "nt":  # Windows
                        result = subprocess.run(["where", cmd], 
                                              capture_output=True, text=True, timeout=5)
                    else:  # Unix-like systems
                        result = subprocess.run(["which", cmd], 
                                              capture_output=True, text=True, timeout=5)
                    
                    if result.returncode == 0:
                        paths = result.stdout.strip().split('\n')
                        # 获取第一个结果并去掉可能的回车符
                        path = paths[0].strip()
                        found = True
                        break
                except (subprocess.TimeoutExpired, FileNotFoundError):
                    continue
            
            # 保存配置结果
            configured_tool = tool.copy()
            configured_tool["path"] = path if found else ""
            self.configured_tools.append(configured_tool)
            
            if not found:
                not_found_tools.append(tool)
        
        # 如果有未找到的工具，显示提示
        if not_found_tools:
            message = "The following tools were not found in your PATH:\n"
            for tool in not_found_tools:
                message += f"- {tool['name']}\n"
            message += "\nYou can manually configure their paths or visit their websites to download and install them."
            
            # 显示未找到工具的网站信息
            websites = "\n\nWebsites to download tools:\n"
            for tool in not_found_tools:
                websites += f"{tool['name']}: {tool['website']}\n"
            
            QMessageBox.information(self, "Auto Config Result", message + websites)
        else:
            QMessageBox.information(self, "Auto Config Result", "All tools found successfully!")
    
    def manual_config(self):
        """手动配置工具路径"""
        if not self.tools_data:
            QMessageBox.warning(self, "Warning", "No tools data loaded.")
            return
        
        # 使用当前已配置的数据或原始模板数据
        tools_to_configure = self.configured_tools if self.configured_tools else self.tools_data
        
        dialog = ManualConfigDialog(tools_to_configure, self)
        if dialog.exec_() == QDialog.Accepted:
            self.configured_tools = dialog.get_configured_tools()
            QMessageBox.information(self, "Manual Config", "Tools configured successfully!")
    
    def finish_install(self):
        """完成安装"""
        # 检查是否进行了配置
        if not self.configured_tools:
            reply = QMessageBox.question(self, "Warning", 
                                       "Tools have not been configured. Do you want to continue anyway?",
                                       QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.No:
                return
        
        # 获取安装路径
        mpe_path = self.path_edit.text()
        
        # 1. 处理settings.template.json
        settings_template_path = os.path.join(self.plugin_path, "settings.template.json")
        settings_path = os.path.join(self.plugin_path, "settings.json")
        
        if os.path.exists(settings_template_path):
            with open(settings_template_path, "r", encoding='utf-8') as f:
                settings_content = f.read()
            
            # 替换版本号和路径
            try:
                from YR_MPE import __version__
                settings_content = settings_content.replace("%version%", __version__)
            except ImportError:
                settings_content = settings_content.replace("%version%", "0.1.0")
            
            settings_content = settings_content.replace("%MPEPATH%", mpe_path)
            
            # 保存到settings.json
            with open(settings_path, "w", encoding='utf-8') as f:
                f.write(settings_content)
        
        # 2. 创建config.json
        config_data = []
        for tool in self.configured_tools:
            config_data.append({
                "name": tool["name"],
                "path": tool["path"]
            })
        
        config_path = os.path.join(self.plugin_path, "YR_MPE", "config.json")
        with open(config_path, "w", encoding='utf-8') as f:
            json.dump(config_data, f, indent=4, ensure_ascii=False)
        
        QMessageBox.information(self, "Installation Complete", 
                              "YR-MPE has been successfully configured and installed!\nRestart YR-MPE to apply changes.")
        self.close()


class InitializerEntry:
    def run(self):
        return Initializer()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    initializer = Initializer()
    initializer.show()
    sys.exit(app.exec_())