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

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from PyQt5.QtCore import Qt, QUrl, pyqtSignal
from PyQt5.QtWebEngineWidgets import QWebEngineView
import os
import tempfile
from core.plugin_base import BasePlugin

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
        if success and self.nwk_string:
            # 页面加载成功且有Newick字符串，则注入树数据
            # 添加延时确保JavaScript完全加载
            from PyQt5.QtCore import QTimer
            QTimer.singleShot(1000, lambda: self.send_newick_to_icytree(self.nwk_string))
    
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
        self.status_changed.emit("Loading tree to IcyTree...")
        
        # 如果页面已加载，则直接发送Newick字符串
        # 否则，将在on_load_finished中处理
        if self.web_view.url().isValid():
            self.send_newick_to_icytree(nwk_string)
        # else: 页面尚未加载，将在on_load_finished中处理

        # Update status after delay
        from PyQt5.QtCore import QTimer
        import logging
        try:
            QTimer.singleShot(2000, lambda: self.status_changed.emit("Tree loaded to IcyTree"))
        except Exception as e:
            logging.error(f"Failed to update status. Check whether the window is still open: {e}")
        
        return True

class IcyTreePluginEntry:
    """Plugin entry point"""
    def run(self):
        plugin = IcyTreePlugin()
        return plugin.run()
