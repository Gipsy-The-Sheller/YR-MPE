import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QToolBar, 
                             QToolButton, QMenu, QAction, QStyle)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QIcon

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle("Molecular Evolutionary Genetics Analysis")
        self.setGeometry(100, 100, 800, 600)
        
        # 创建工具栏
        self.toolbar = QToolBar("Main Toolbar")
        self.toolbar.setIconSize(QSize(32, 32))
        self.addToolBar(Qt.TopToolBarArea, self.toolbar)
        
        # 创建ALIGN按钮（带下拉菜单）
        align_button = self.create_tool_button(
            "ALIGN", 
            self.style().standardIcon(QStyle.SP_FileIcon),
            self.create_align_menu()
        )
        self.toolbar.addWidget(align_button)
        
        # 添加其他按钮...
        # 这里可以继续添加其他类似的工具栏按钮
        
    def create_tool_button(self, text, icon, menu):
        button = QToolButton()
        button.setText(text)
        button.setIcon(icon)
        button.setMenu(menu)
        button.setPopupMode(QToolButton.InstantPopup)  # 点击按钮立即显示菜单
        button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)  # 文字在图标下方
        button.setFixedSize(60, 60)  # 设置按钮大小
        return button
        
    def create_align_menu(self):
        menu = QMenu(self)
        
        # 添加菜单项
        open_action = QAction("Open a File/Session...", self)
        concatenate_action = QAction("Concatenate Sequence Alignments", self)
        # 添加更多菜单项...
        
        menu.addAction(open_action)
        menu.addAction(concatenate_action)
        
        # 连接信号和槽
        open_action.triggered.connect(self.open_file)
        concatenate_action.triggered.connect(self.concatenate_sequences)
        
        return menu
        
    def open_file(self):
        print("Open file action triggered")
        
    def concatenate_sequences(self):
        print("Concatenate sequences action triggered")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())