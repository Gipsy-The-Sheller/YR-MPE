from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QLabel, QTabWidget, 
                             QTextEdit, QPushButton, QHBoxLayout, QGroupBox)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QPixmap
import os


class AboutPlugin(QWidget):
    """关于YR-MPE插件"""
    
    def __init__(self, import_from=None, import_data=None):
        """初始化关于插件"""
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.init_ui()
        
    def init_ui(self):
        """初始化用户界面"""
        self.setWindowTitle("About YR-MPE")
        self.setMinimumSize(600, 400)
        
        # 创建主布局
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # 创建标题
        header_label = QLabel("YR-MPE")
        header_font = QFont()
        header_font.setPointSize(24)
        header_font.setBold(True)
        header_label.setFont(header_font)
        header_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(header_label)

        # 副标题
        subtitle_label = QLabel("Molecular Phylogenetics & Evolution Plugins for YRTools")
        subtitle_font = QFont()
        subtitle_font.setPointSize(16)
        subtitle_label.setFont(subtitle_font)
        subtitle_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(subtitle_label)
        
        # 创建许可证信息组
        license_group = QGroupBox()
        license_layout = QVBoxLayout()
        license_group.setLayout(license_layout)
        
        # 许可证图标
        license_icon = QLabel()
        gpl_icon_path = os.path.join(self.plugin_path, "icons", "GPL3.svg")
        if os.path.exists(gpl_icon_path):
            license_icon.setText(f'<img src="{gpl_icon_path}"/>')
        else:
            license_icon.setText("GPL3")
        license_icon.setAlignment(Qt.AlignCenter)
        license_layout.addWidget(license_icon)
        
        # 简明许可证说明
        license_text = QLabel("This program is free software: you can redistribute it and/or modify "
                             "it under the terms of the GNU General Public License as published by "
                             "the Free Software Foundation, either version 3 of the License, or "
                             "(at your option) any later version.")
        license_text.setWordWrap(True)
        license_text.setAlignment(Qt.AlignLeft)
        license_layout.addWidget(license_text)
        
        main_layout.addWidget(license_group)
        


        copyright_group = QGroupBox()
        copyright_layout = QVBoxLayout()
        copyright_group.setLayout(copyright_layout)

        # 版权信息
        copyright_label = QLabel("Copyright (c) 2025 Zhi-Jie Xu")
        copyright_label.setAlignment(Qt.AlignCenter)
        copyright_layout.addWidget(copyright_label)
        
        # 联系信息
        email_label = QLabel("Email: zjxmolls@outlook.com")
        email_label.setAlignment(Qt.AlignCenter)
        copyright_layout.addWidget(email_label)

        # 存储库信息
        repo_label = QLabel("Repository: <a href='https://github.com/Gipsy-The-Sheller/YR-MPE'>https://github.com/Gipsy-The-Sheller/YR-MPE</a>")
        repo_label.setAlignment(Qt.AlignCenter)
        repo_label.setOpenExternalLinks(True)
        copyright_layout.addWidget(repo_label)

        main_layout.addWidget(copyright_group) 

        # 版本组
        version_group = QGroupBox()
        version_layout = QVBoxLayout()
        version_group.setLayout(version_layout)

        # 版本信息
        from . import __version__
        version_label = QLabel(f"Version: {__version__}")
        version_label.setAlignment(Qt.AlignCenter)
        version_layout.addWidget(version_label)  

        update_label = QLabel("Get new updates from YR-Pacman or by Git pull.")
        update_label.setAlignment(Qt.AlignCenter)
        version_layout.addWidget(update_label)

        main_layout.addWidget(version_group)      
        # 添加伸展以美化布局
        main_layout.addStretch()


class AboutPluginEntry:
    def run(self):
        return AboutPlugin()