# IcyTree.py
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

# Copyright (C) 2025 Zhi-Jie Xu & Yi-Yang Jia
# 
# This file is part of YRTools.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from PyQt5.QtWidgets import QWidget
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtGui import QIcon
import os
from pathlib import Path
from PyQt5.QtCore import QTimer

class BasePlugin(QObject):
    # 定义标准信号
    status_changed = pyqtSignal(str)  # 状态更新
    data_ready = pyqtSignal(dict)     # 数据输出
    
    def __init__(self, config=None):
        super().__init__()
        self.config = config
        self.widget = None
        
    def create_widget(self):
        """必须由子类实现"""
        raise NotImplementedError
        
    def run(self):
        """执行插件主逻辑"""
        if not self.widget:
            self.widget = self.create_widget()
        return self.widget 

    @classmethod
    def get_icon(cls, plugin_path):
        """获取插件图标，支持多种矢量格式"""
        formats = ['eps', 'emf', 'svg']  # 按优先级排序
        icon_dir = Path(plugin_path).parent
        
        # 搜索所有支持的图标文件
        found = {}
        for f in icon_dir.glob("icon.*"):
            ext = f.suffix[1:].lower()
            if ext in formats:
                found[ext] = f
        
        # 按格式优先级选择
        for fmt in formats:
            if fmt in found:
                return cls._load_vector_icon(found[fmt])
        
        return QIcon()  # 返回空图标

    @staticmethod
    def _load_vector_icon(file_path):
        """加载矢量图标文件"""
        ext = file_path.suffix.lower()
        
        # EPS/EMF需要转换
        if ext in ('.eps', '.emf'):
            try:
                # 使用ghostscript转换EPS/EMF为临时PNG
                from subprocess import run
                temp_png = file_path.with_suffix('.png')
                
                # 添加错误检查和详细日志
                result = run(
                    ['gs', '-dSAFER', '-dBATCH', '-dNOPAUSE', 
                     '-sDEVICE=png16m', f'-sOutputFile={temp_png}',
                     '-r300', str(file_path)],
                    capture_output=True,
                    text=True
                )
                
                if result.returncode != 0:
                    print(f"Ghostscript转换失败 (代码 {result.returncode}):")
                    print(f"错误输出: {result.stderr}")
                    return QIcon()
                
                icon = QIcon(str(temp_png))
                temp_png.unlink()  # 删除临时文件
                return icon
            except Exception as e:
                print(f"EPS转换异常: {str(e)}")
                return QIcon()
        
        # 直接加载SVG
        elif ext == '.svg':
            return QIcon(str(file_path))
        
        return QIcon() 

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from PyQt5.QtCore import Qt, QUrl, pyqtSignal
from PyQt5.QtWebEngineWidgets import QWebEngineView
import os
import tempfile
# from core.plugin_base import BasePlugin

class IcyTreePlugin(QWidget, BasePlugin):
    """QWebEngine-based IcyTree plugin for phylogenetic tree visualization"""
    
    # Redefine signals
    status_changed = pyqtSignal(str)
    data_ready = pyqtSignal(dict)
    
    def __init__(self, config=None, plugin_path: str = None, **kwargs):
        self.plugin_path = plugin_path or os.path.dirname(os.path.abspath(__file__))
        QWidget.__init__(self)
        BasePlugin.__init__(self, config)
        self.nwk_string = ""
        self.icytree_path = None
        self.destroyed_flag = False  # 标记对象是否即将被销毁
        self.init_ui()
        
    def init_ui(self):
        """Initialize user interface"""
        self.setWindowTitle("IcyTree - Phylogenetic Tree Visualization")
        # self.setMinimumSize(1000, 700)
        
        # Main layout
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # Create QWebEngineView
        self.web_view = QWebEngineView()
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(self.web_view)
        
        # Initialize IcyTree
        self.init_icytree()
    
    def create_widget(self):
        """Create plugin main interface - now returns self"""
        return self
    
    def init_icytree(self):
        """Initialize IcyTree"""
        try:
            # Get icytree HTML file path
            if self.plugin_path:
                icytree_html = os.path.join(self.plugin_path, "src", "icytree", "index.html")
            else:
                # Default path
                icytree_html = os.path.join(os.path.dirname(__file__), "src", "icytree", "index.html")
            
            if os.path.exists(icytree_html):
                self.icytree_path = icytree_html
                # Load icytree page
                file_url = QUrl.fromLocalFile(icytree_html)
                
                # 断开可能存在的旧连接
                try:
                    self.web_view.loadFinished.disconnect()
                except TypeError:
                    pass  # 尚未连接
                
                # 连接页面加载完成信号
                self.web_view.loadFinished.connect(self.on_load_finished)
                
                self.web_view.load(file_url)
                self.status_changed.emit("IcyTree loaded")
            else:
                self.show_error_page("IcyTree HTML file not found")
                self.status_changed.emit("Error: IcyTree file not found")
                
        except Exception as e:
            self.show_error_page(f"Failed to initialize IcyTree: {str(e)}")
            self.status_changed.emit(f"Error: {str(e)}")
    
    def on_load_finished(self, success):
        """页面加载完成后的回调函数"""
        try:
            # 检查对象是否仍然存在或即将被销毁
            if self.destroyed_flag or not hasattr(self, 'web_view') or not self.web_view:
                return
                
            if success and self.nwk_string:
                # 页面加载成功且有Newick字符串，则注入树数据
                # 添加延时确保JavaScript完全加载
                QTimer.singleShot(1000, self._delayed_tree_injection)
        except RuntimeError:
            # 对象已被销毁，忽略
            pass
        except Exception as e:
            print(f"Error in on_load_finished: {e}")
    
    def _delayed_tree_injection(self):
        """延迟树注入，检查对象是否仍然存在"""
        try:
            # 检查对象是否仍然存在或即将被销毁
            if self.destroyed_flag or not hasattr(self, 'nwk_string') or not self.nwk_string:
                return
                
            self.send_newick_to_icytree(self.nwk_string)
        except RuntimeError:
            # 对象已被销毁，忽略
            pass
        except Exception as e:
            print(f"Error in delayed tree injection: {e}")
    
    def validate_newick(self, nwk_string):
        """Simple Newick format validation"""
        if not nwk_string:
            return False
        
        # # Basic check: should end with semicolon
        # if not nwk_string.endswith(';'):
        #     return False
        
        # Check bracket matching
        open_count = nwk_string.count('(')
        close_count = nwk_string.count(')')
        return open_count == close_count
    
    def send_newick_to_icytree(self, nwk_string):
        """Send Newick string to IcyTree"""
        try:
            # 检查对象是否仍然存在或即将被销毁
            if self.destroyed_flag or not hasattr(self, 'web_view') or not self.web_view:
                return
                
            # Escape special characters in JavaScript string
            escaped_nwk = nwk_string.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$').replace('"', '\\"')
            
            # Use JavaScript to pass Newick string to IcyTree
            js_code = f"""
            (function() {{
                try {{
                    // Set the tree data
                    window.treeData = `{escaped_nwk}`;
                    // Reload the tree data
                    if (typeof window.reloadTreeData === 'function') {{
                        window.reloadTreeData();
                    }} else {{
                        console.error("reloadTreeData function not found");
                    }}
                }} catch (e) {{
                    console.error("Error injecting tree:", e);
                }}
            }})();
            """
            
            self.web_view.page().runJavaScript(js_code)
            
        except RuntimeError:
            # 对象已被删除，忽略
            pass
        except Exception as e:
            print(f"Failed to send Newick to IcyTree: {e}")
            # Fallback: create temporary file
            self.load_tree_from_file(nwk_string)
    
    def load_tree_from_file(self, nwk_string):
        """Fallback: load tree via file"""
        try:
            # Escape special characters in JavaScript string
            escaped_nwk = nwk_string.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$').replace('"', '\\"')
            
            # Use JavaScript to simulate drag and drop operation
            js_code = f"""
            (function() {{
                try {{
                    // Set the tree data
                    window.treeData = `{escaped_nwk}`;
                    // Reload the tree data
                    if (typeof window.reloadTreeData === 'function') {{
                        window.reloadTreeData();
                    }} else {{
                        console.error("reloadTreeData function not found");
                    }}
                }} catch (e) {{
                    console.error("Error injecting tree:", e);
                }}
            }})();
            """
            
            self.web_view.page().runJavaScript(js_code)
            
        except Exception as e:
            print(f"Failed to load tree via fallback: {e}")
    
    def show_error_page(self, error_message):
        """Show error page"""
        try:
            # 检查对象是否仍然存在或即将被销毁
            if self.destroyed_flag or not hasattr(self, 'web_view') or not self.web_view:
                return
        except:
            pass        
            error_html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>IcyTree Error</title>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    margin: 50px;
                    background-color: #f5f5f5;
                }}
                .error-container {{
                    background-color: white;
                    padding: 30px;
                    border-radius: 8px;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                    text-align: center;
                }}
                .error-icon {{
                    font-size: 48px;
                    color: #e74c3c;
                    margin-bottom: 20px;
                }}
                .error-title {{
                    color: #2c3e50;
                    margin-bottom: 15px;
                }}
                .error-message {{
                    color: #7f8c8d;
                    line-height: 1.6;
                }}
            </style>
        </head>
        <body>
            <div class="error-container">
                <div class="error-icon">⚠️</div>
                <h2 class="error-title">IcyTree Loading Failed</h2>
                <p class="error-message">{error_message}</p>
                <p class="error-message">Please check if the IcyTree file path is correct.</p>
            </div>
        </body>
        </html>
        """
        
        self.web_view.setHtml(error_html)
    
    def set_newick_string(self, nwk_string):
        """External interface: set Newick string and load it to IcyTree"""
        self.nwk_string = nwk_string
        
        # Validate Newick format
        if not self.validate_newick(nwk_string):
            print(f"Invalid Newick format: {nwk_string}")
            return False
        
        # Update status
        try:
            self.status_changed.emit("Loading tree to IcyTree...")
        except RuntimeError:
            # 对象已被销毁，忽略
            return False
        
        # 如果页面已加载，则直接发送Newick字符串
        # 否则，将在on_load_finished中处理
        if self.web_view.url().isValid():
            self.send_newick_to_icytree(nwk_string)
        # else: 页面尚未加载，将在on_load_finished中处理

        # Update status after delay
        try:
            QTimer.singleShot(2000, self._delayed_status_update)
        except Exception as e:
            import logging
            logging.error(f"Failed to update status. Check whether the window is still open: {e}")
        
        return True
    
    def _delayed_status_update(self):
        """延迟状态更新，检查对象是否仍然存在"""
        try:
            # 检查对象是否仍然存在或即将被销毁
            if self.destroyed_flag or not hasattr(self, 'web_view') or not self.web_view:
                return
                
            self.status_changed.emit("Tree loaded to IcyTree")
        except RuntimeError:
            # 对象已被删除
            pass
        except Exception as e:
            # 其他异常
            import logging
            logging.error(f"Error in delayed status update: {e}")

class IcyTreePluginEntry:
    """Plugin entry point"""
    def run(self):
        plugin = IcyTreePlugin()
        return plugin.run()

    def __del__(self):
        """析构函数，标记对象即将被销毁"""
        self.destroyed_flag = True
        
        # 清理可能存在的定时器引用
        try:
            # 清理web_view引用
            if hasattr(self, 'web_view') and self.web_view:
                self.web_view.deleteLater()
                self.web_view = None
        except:
            pass
