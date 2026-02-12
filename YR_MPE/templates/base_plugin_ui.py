from PyQt5.QtWidgets import (QWidget, QLabel, QVBoxLayout, QDialog, QTabWidget,
                             QTextEdit, QPushButton, QFileDialog, QComboBox,
                             QSpinBox, QCheckBox, QHBoxLayout, QProgressBar,
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, QSizePolicy,
                             QScrollArea, QFrame, QPlainTextEdit, QDoubleSpinBox, QApplication, QToolButton)
from PyQt5.QtCore import QThread, pyqtSignal, Qt, QUrl, QTimer
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtGui import QFont, QTextCursor, QIcon
import json
import os
import subprocess
import tempfile
import platform
import datetime
from Bio import SeqIO


class BasePlugin(QWidget):
    """
    YR-MPEA基础插件类
    
    这是所有YR-MPEA插件的基类，提供了通用的功能和接口：
    - 插件配置管理
    - 控制台日志系统
    - UI结构（标签页、控制面板等）
    - 文件导入/导出支持
    - 进度条和运行控制
    """
    
    def __init__(self, import_from=None, import_data=None):
        """
        初始化基础插件
        
        Args:
            import_from (str): 数据导入来源
            import_data (any): 导入的数据
        """
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.destroyed.connect(self.on_destroyed)
        
        # 导入相关
        self.import_from = import_from
        self.import_data = import_data
        self.import_file = None
        self.temp_files = []
        
        # 运行状态
        self.is_running = False
        self.tool_path = None
        
        # 数据存储
        self.imported_files = []  # 导入的文件路径列表
        self.file_tags = []       # 文件标签组件列表
        self.reports = []         # 报告文件列表
        self.current_report_index = 0
        self.console_output = []  # 控制台消息列表
        
        # 初始化插件信息
        self.init_plugin_info()
        
        # 检查配置
        if not self.config():
            self.show_config_guide()
    
        # 初始化UI
        self.init_ui()
        
        # 处理导入数据
        if import_from == "YR_MPEA" and isinstance(import_data, list):
            self.handle_import_data(import_data)

    def init_plugin_info(self):
        """
        初始化插件信息
        子类应重写此方法来设置插件元数据
        """
        self.plugin_name = "Base Plugin"
        self.tool_name = "base_tool"
        self.citation = ["No citation available"]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"FASTA": ".fas"}

    def info(self):
        """
        获取插件信息
        用于向用户显示插件的元数据
        """
        return {
            "name": self.plugin_name,
            "tool": self.tool_name,
            "citation": self.citation,
            "input_types": self.input_types,
            "output_types": self.output_types
        }

    def config(self):
        """
        检查插件配置
        
        Returns:
            bool: 配置是否有效
        """
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            if not os.path.exists(config_path):
                return False
            
            with open(config_path, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
            
            plugin_config = None
            for tool in config_data:
                if tool.get("name").lower() == self.tool_name.lower():
                    plugin_config = tool
                    break
            
            if not plugin_config:
                return False
            
            self.tool_path = os.path.join(self.plugin_path, plugin_config["path"].lstrip("/").lstrip("\\"))
            return os.path.exists(self.tool_path)
            
        except Exception as e:
            print(f"Config check failed: {e}")
            return False

    def show_config_guide(self):
        """
        显示配置指南
        当插件未正确配置时调用
        """
        config_guide = QDialog(self)
        config_guide.setWindowTitle("Configuration Guide")
        config_guide.setMinimumSize(600, 400)
        layout = QVBoxLayout()
        
        title = QLabel(f"Configure {self.plugin_name}")
        title.setStyleSheet("font-size: 16px; font-weight: bold; margin: 10px;")
        layout.addWidget(title)
        
        instruction = QLabel(
            f"To use this plugin, you need to configure the path to {self.tool_name}.\n\n"
            "1. Locate the {self.tool_name} executable on your system\n"
            "2. Add it to the config.json file in the plugin directory\n"
            "3. The entry should look like this:\n\n"
            f'{{\n  "name": "{self.tool_name}",\n  "path": "/path/to/executable"\n}}'
        )
        instruction.setWordWrap(True)
        layout.addWidget(instruction)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(config_guide.close)
        layout.addWidget(close_btn)
        
        config_guide.setLayout(layout)
        config_guide.exec_()

    def handle_import_data(self, import_data):
        """处理从YR-MPEA导入的数据"""
        if isinstance(import_data, list):
            # 创建临时文件来存储导入的序列数据
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                for seq in import_data:
                    f.write(f">{seq.id}\n{seq.seq}\n")
            self.temp_files.append(temp_file)
            self.import_file = temp_file
            self.imported_files = [temp_file]
            
            # 更新UI显示导入的文件
            if hasattr(self, 'file_path_edit') and self.file_path_edit:
                self.file_path_edit.setText(temp_file)
        else:
            self.import_file = None
            self.imported_files = []

    def on_destroyed(self):
        """
        清理资源
        在插件销毁时删除临时文件
        """
        for temp_file in self.temp_files:
            try:
                os.remove(temp_file)
            except:
                pass

    def init_ui(self):
        """
        初始化用户界面
        设置基本的标签页结构和控制面板
        """
        self.setWindowTitle(self.plugin_name)
        self.main_layout = QVBoxLayout()
        self.setLayout(self.main_layout)
        
        # 标签页
        self.tab_widget = QTabWidget()
        self.main_layout.addWidget(self.tab_widget)
        
        # 输入和参数标签页
        self.input_tab = QWidget()
        self.tab_widget.addTab(self.input_tab, "Input & Parameters")
        self.setup_input_tab()
        
        # 输出预览标签页
        self.output_tab = QWidget()
        self.tab_widget.addTab(self.output_tab, "Output Preview")
        self.setup_output_tab()
        
        # 控制台标签页
        self.console_tab = QWidget()
        self.tab_widget.addTab(self.console_tab, "Console")
        self.setup_console_tab()

        # 引用页
        self.citation_tab = QWidget()
        self.tab_widget.addTab(self.citation_tab, "Citation")
        self.setup_citation_tab()
        
        # 控制面板（进度条和运行按钮）
        self.setup_control_panel()
        
        # 初始化控制台
        self.auto_scroll = True
        self.add_console_message(f"{self.plugin_name} initialized", "info")

    def setup_input_tab(self):
        """
        设置输入标签页
        子类应重写此方法来创建特定的输入界面
        """
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        label = QLabel("Override setup_input_tab() to create input layout")
        layout.addWidget(label)

    def setup_output_tab(self):
        """
        设置输出标签页
        子类应重写此方法来创建特定的输出界面
        """
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 报告选择器
        report_layout = QHBoxLayout()
        report_layout.addWidget(QLabel("Select Report:"))
        
        self.report_combo = QComboBox()
        self.report_combo.setEnabled(False)
        self.report_combo.currentIndexChanged.connect(self.on_report_changed)
        report_layout.addWidget(self.report_combo)
        
        report_layout.addStretch()
        layout.addLayout(report_layout)
        
        # HTML预览
        self.web_view = QWebEngineView()
        self.web_view.setHtml("<p>Waiting for analysis results...</p>")
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(self.web_view)
        
        # 输出文件信息
        self.output_info = QLabel("Output file will be displayed here")
        layout.addWidget(self.output_info)

    def setup_console_tab(self):
        """
        设置控制台标签页
        提供命令和输出显示功能
        """
        layout = QVBoxLayout()
        self.console_tab.setLayout(layout)
        
        # 控制台控制
        console_controls = QHBoxLayout()
        
        clear_btn = QPushButton("Clear Console")
        clear_btn.clicked.connect(self.clear_console)
        console_controls.addWidget(clear_btn)
        
        auto_scroll_checkbox = QCheckBox("Auto Scroll")
        auto_scroll_checkbox.setChecked(True)
        auto_scroll_checkbox.toggled.connect(self.toggle_auto_scroll)
        console_controls.addWidget(auto_scroll_checkbox)
        
        console_controls.addStretch()
        layout.addLayout(console_controls)
        
        # 控制台文本区域
        self.console_text = QPlainTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Courier New", 10))
        
        
        # 根据平台设置样式
        if platform.system() == "Windows":
            self.console_text.setStyleSheet("""
                QPlainTextEdit {
                    background-color: #272822;
                    color: #f8f8f2;
                    border: 1px solid #3c3c3c;
                    font-family: 'Consolas', monospace;
                }
            """)  # Monokai主题
        else:
            self.console_text.setStyleSheet("""
                QPlainTextEdit {
                    background-color: #272822;
                    color: #f8f8f2;
                    border: 1px solid #3c3c3c;
                    font-family: 'Courier New', monospace;
                }
            """)  # Monokai主题
            
        layout.addWidget(self.console_text)

    def setup_citation_tab(self):
        self.citation_tab_layout = QVBoxLayout()
        self.citation_tab.setLayout(self.citation_tab_layout)
        citation_prompt = QLabel(f"If you use {self.plugin_name} in your work, please cite:")
        citation_prompt.setAlignment(Qt.AlignTop)
        citation_prompt.setWordWrap(True)
        self.citation_tab_layout.addWidget(citation_prompt)

        # find citation
        citation = self.get_citation()
        for i in citation:
            self.export_citation_layout(i)
        self.citation_tab_layout.addStretch()
    def export_citation_layout(self, citation, plugin_path=None, remarks=None):
        if plugin_path == None: 
            plugin_path = self.plugin_path
        citation_group = QGroupBox()
        citation_layout = QHBoxLayout()
        citation_group.setLayout(citation_layout)
        if remarks:
            citation_remarks = QLabel(remarks)
            citation_layout.addWidget(citation_remarks)
        citation_text = QLabel(citation)
        citation_text.setWordWrap(True)
        citation_text.setOpenExternalLinks(True)

        def copy_citation(citation):
            """
            复制引用信息到剪贴板
            """
            QApplication.clipboard().setText(citation)
            self.add_console_message("Citation copied to clipboard", "info")
        copy_btn = QPushButton()
        copy_btn.setIcon(QIcon(os.path.join(self.plugin_path, "icons/copy.svg")))
        copy_btn.clicked.connect(lambda:copy_citation(citation))
        copy_btn.setFixedSize(24, 24)

        citation_layout.addWidget(citation_text)
        citation_layout.addWidget(copy_btn)

        citation_group.setLayout(citation_layout)
        citation_group.setAlignment(Qt.AlignTop)
        citation_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.citation_tab_layout.addWidget(citation_group)

    def setup_control_panel(self):
        """
        设置控制面板
        包含进度条和运行按钮
        """
        self.control_layout = QHBoxLayout()
        
        # 进度条
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.control_layout.addWidget(self.progress_bar)
        
        # 运行按钮
        self.run_button = QPushButton(f"Run {self.tool_name}")
        self.run_button.clicked.connect(self.run_analysis)
        self.control_layout.addWidget(self.run_button)
        
        # 停止按钮
        self.stop_button = QPushButton("Stop")

        self.stop_button.clicked.connect(self.stop_analysis)
        self.stop_button.setEnabled(False)
        self.control_layout.addWidget(self.stop_button)
        
        self.main_layout.addLayout(self.control_layout)

    def clear_console(self):
        """
        清空控制台输出
        """
        self.console_text.clear()
        self.console_output.clear()
        self.add_console_message("Console cleared", "info")

    def toggle_auto_scroll(self, enabled):
        """
        切换自动滚动功能
        
        Args:
            enabled (bool): 是否启用自动滚动
        """
        self.auto_scroll = enabled

    def add_console_message(self, message, msg_type="info"):
        """
        添加消息到控制台
        
        Args:
            message (str): 消息内容
            msg_type (str): 消息类型 (info, command, output, error, warning)
        """
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        
        # 根据消息类型设置颜色
        if msg_type == "command":
            formatted_msg = f"[{timestamp}] $ {message}"
            color = "#4CAF50"  # 绿色表示命令
        elif msg_type == "output":
            formatted_msg = f"[{timestamp}] > {message}"
            color = "#2196F3"  # 蓝色表示输出
        elif msg_type == "error":
            formatted_msg = f"[{timestamp}] ERROR: {message}"
            color = "#F44336"  # 红色表示错误
        elif msg_type == "warning":
            formatted_msg = f"[{timestamp}] WARNING: {message}"
            color = "#FF9800"  # 橙色表示警告
        else:  # info
            formatted_msg = f"[{timestamp}] {message}"
            color = "#9E9E9E"  # 灰色表示信息
        
        # 添加到控制台
        self.console_text.appendPlainText(formatted_msg)
        
        # 如果启用自动滚动，则滚动到底部
        if self.auto_scroll:
            self.console_text.moveCursor(QTextCursor.End)
        
        # 存储到控制台输出列表
        self.console_output.append({
            'timestamp': timestamp,
            'message': message,
            'type': msg_type,
            'formatted': formatted_msg
        })

    def get_citation(self):
        """
        获取引用信息
        
        Returns:
            str: 插件的引用信息
        """
        return self.citation

    def run_analysis(self):
        """
        运行分析
        子类应重写此方法来实现具体的分析逻辑
        """
        self.add_console_message("Override run_analysis() to implement analysis logic", "warning")

    def stop_analysis(self):
        """
        停止分析
        子类应重写此方法来实现具体的停止逻辑
        """
        self.add_console_message("Override stop_analysis() to implement stop logic", "warning")

    def on_report_changed(self, index):
        """
        报告选择改变时的处理函数
        
        Args:
            index (int): 选中的报告索引
        """
        if 0 <= index < len(self.reports):
            self.current_report_index = index
            self.show_current_report()

    def show_current_report(self):
        """
        显示当前选中的报告
        """
        if not self.reports or self.current_report_index >= len(self.reports):
            self.web_view.setHtml("<h2>No report available</h2>")
            return
        
        report_file = self.reports[self.current_report_index]
        if report_file and os.path.exists(report_file):
            # 转换Windows路径为文件URL格式
            file_url = QUrl.fromLocalFile(report_file)
            self.web_view.load(file_url)
        else:
            self.web_view.setHtml("<h2>Report file not found</h2>")

    def update_report_combo(self):
        """
        更新报告选择下拉框
        """
        self.report_combo.clear()
        
        if not self.reports:
            self.report_combo.setEnabled(False)
            return
        
        # 添加报告选项
        for i, report_file in enumerate(self.reports):
            if os.path.exists(report_file):
                report_name = f"Report {i+1}"
                if len(self.imported_files) > 1:
                    # 如果有多个文件，显示该报告对应的文件名
                    if i < len(self.imported_files):
                        file_name = os.path.basename(self.imported_files[i])
                        report_name = f"Report {i+1} ({file_name})"
                
                self.report_combo.addItem(report_name)
        
        self.report_combo.setEnabled(len(self.reports) > 1)
        self.report_combo.setCurrentIndex(self.current_report_index)

    def create_temp_file(self, suffix=''):
        """创建临时文件"""
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        temp_file.close()
        return temp_file.name

    def on_destroyed(self):
        """清理临时文件"""
        for temp_file in self.temp_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            except Exception as e:
                print(f"Error removing temp file {temp_file}: {e}")