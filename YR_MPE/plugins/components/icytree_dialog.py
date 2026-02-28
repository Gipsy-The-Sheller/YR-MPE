# icytree_dialog.py
#
# Copyright (c) 2026 Zhi-Jie Xu
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

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, 
                            QLabel, QMessageBox, QSizePolicy, QWidget)
from PyQt5.QtCore import Qt, QTimer
import os
import tempfile


class IcyTreeDialog(QDialog):
    """IcyTree置根对话框"""
    
    def __init__(self, newick_str, plugin_path, parent=None):
        super().__init__(parent)
        self.newick_str = newick_str
        self.plugin_path = plugin_path
        self.rooted_newick = None
        self.temp_tree_file = None
        self.web_view = None
        self.js_injected = False
        self.page_loaded = False
        
        self.init_ui()
    
    def init_ui(self):
        """初始化UI"""
        self.setWindowTitle("IcyTree - Root the tree")
        self.setMinimumSize(1000, 700)
        
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # 加载状态标签
        self.loading_label = QLabel("Loading IcyTree...")
        self.loading_label.setAlignment(Qt.AlignCenter)
        self.loading_label.setStyleSheet("color: #666; font-style: italic; padding: 20px;")
        main_layout.addWidget(self.loading_label)
        
        # 说明标签
        instruction_label = QLabel(
            "Instructions: Right-click on a node in the tree and select 'Re-root' "
            "to root the tree at that node. Then click 'Apply Root' to confirm."
        )
        instruction_label.setWordWrap(True)
        instruction_label.setStyleSheet("color: #666; font-style: italic;")
        main_layout.addWidget(instruction_label)
        
        # WebEngineView 容器
        self.web_container = QWidget()
        self.web_container.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        main_layout.addWidget(self.web_container)
        
        # 按钮布局
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        
        self.apply_button = QPushButton("Apply Root")
        self.apply_button.setMinimumWidth(120)
        self.apply_button.clicked.connect(self.on_apply_root)
        self.apply_button.setEnabled(False)
        button_layout.addWidget(self.apply_button)
        
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setMinimumWidth(120)
        self.cancel_button.clicked.connect(self.reject)
        button_layout.addWidget(self.cancel_button)
        
        main_layout.addLayout(button_layout)
    
    def showEvent(self, event):
        """对话框显示时的事件处理"""
        super().showEvent(event)
        # 在对话框显示后延迟加载，避免阻塞
        QTimer.singleShot(100, self.setup_web_view)
    
    def setup_web_view(self):
        """设置WebEngineView"""
        try:
            from PyQt5.QtWebEngineWidgets import QWebEngineView
            from PyQt5.QtCore import QUrl
            
            print("Creating QWebEngineView...")
            
            # 创建WebEngineView
            self.web_view = QWebEngineView()
            
            # 获取IcyTree HTML文件路径
            icytree_html = os.path.join(self.plugin_path, "src", "icytree", "index.html")
            
            if not os.path.exists(icytree_html):
                self.loading_label.setText("IcyTree HTML file not found")
                self.loading_label.setStyleSheet("color: red; padding: 20px;")
                return
            
            # 添加到布局
            container_layout = QVBoxLayout(self.web_container)
            container_layout.setContentsMargins(0, 0, 0, 0)
            container_layout.addWidget(self.web_view)
            
            # 连接加载完成信号
            self.web_view.loadFinished.connect(self.on_page_loaded)
            
            # 加载HTML文件
            file_url = QUrl.fromLocalFile(icytree_html)
            print(f"Loading IcyTree from: {icytree_html}")
            self.web_view.load(file_url)
            
        except ImportError as e:
            print(f"Failed to import QWebEngineView: {e}")
            self.loading_label.setText("QWebEngineView not available")
            self.loading_label.setStyleSheet("color: red; padding: 20px;")
        except Exception as e:
            print(f"Failed to setup web view: {e}")
            import traceback
            traceback.print_exc()
            self.loading_label.setText(f"Failed to load IcyTree: {str(e)}")
            self.loading_label.setStyleSheet("color: red; padding: 20px;")
    
    def on_page_loaded(self, success):
        """页面加载完成回调"""
        print(f"Page loaded: {success}")
        
        if not success:
            self.loading_label.setText("Failed to load IcyTree page")
            self.loading_label.setStyleSheet("color: red; padding: 20px;")
            return
        
        self.page_loaded = True
        self.loading_label.hide()
        self.apply_button.setEnabled(True)
        
        # 延迟注入树数据和JavaScript桥接
        QTimer.singleShot(500, self.inject_tree_data)
    
    def inject_tree_data(self):
        """注入树数据和JavaScript桥接"""
        if not self.web_view:
            return
        
        print("Injecting tree data...")
        
        # 先注入树数据
        self._inject_newick()
    
    def _inject_newick(self):
        """注入Newick字符串"""
        escaped_nwk = self.newick_str.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$').replace('"', '\\"')
        
        js_code = f"""
        (function() {{
            try {{
                window.treeData = `{escaped_nwk}`;
                if (typeof window.reloadTreeData === 'function') {{
                    window.reloadTreeData();
                }}
                console.log("Tree data injected successfully");
            }} catch (e) {{
                console.error("Error injecting tree:", e);
            }}
        }})();
        """
        
        self.web_view.page().runJavaScript(js_code, lambda result: self._setup_javascript_bridge())
    
    def _setup_javascript_bridge(self):
        """设置JavaScript桥接函数"""
        js_code = """
        window.getRerootedNewick = function() {
            if (typeof currentTreeIdx !== 'undefined' && 
                currentTreeIdx >= 0 && 
                currentTreeIdx < trees.length) {
                return Write.newick(trees[currentTreeIdx]);
            }
            return null;
        }
        console.log("JavaScript bridge installed");
        """
        
        self.web_view.page().runJavaScript(js_code, lambda result: self._on_js_injected(result))
    
    def _on_js_injected(self, result):
        """JavaScript注入完成回调"""
        self.js_injected = True
        print("JavaScript bridge injected successfully")
    
    def on_apply_root(self):
        """点击Apply Root按钮的回调"""
        if not self.web_view:
            QMessageBox.warning(self, "Error", "IcyTree not available.")
            return
        
        # 显示加载提示
        self.apply_button.setEnabled(False)
        self.apply_button.setText("Getting tree...")
        
        # 调用JavaScript函数获取置根后的Newick字符串
        self.web_view.page().runJavaScript(
            "getRerootedNewick()",
            self._handle_rooted_tree
        )
    
    def _handle_rooted_tree(self, newick_str):
        """处理获取到的置根后的树"""
        self.apply_button.setEnabled(True)
        self.apply_button.setText("Apply Root")
        
        if not newick_str:
            QMessageBox.warning(self, "Error", "Failed to get the rooted tree.")
            return
        
        # 验证Newick字符串
        if not newick_str.strip().endswith(';'):
            newick_str = newick_str.strip() + ';'
        
        self.rooted_newick = newick_str
        
        # 保存到临时文件
        try:
            self.temp_tree_file = tempfile.NamedTemporaryFile(
                mode='w',
                suffix='_rooted.nwk',
                delete=False,
                encoding='utf-8'
            )
            self.temp_tree_file.write(newick_str)
            self.temp_tree_file.close()
            
            # 获取文件路径（在关闭后）
            rooted_tree_path = self.temp_tree_file.name
            
            # 确认对话框
            reply = QMessageBox.question(
                self,
                "Confirm Root",
                "Tree has been rooted. Do you want to apply this rooted tree?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes
            )
            
            if reply == QMessageBox.Yes:
                self.accept()
            else:
                # 用户取消，删除临时文件
                self.cleanup_temp_file()
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save rooted tree: {str(e)}")
            self.cleanup_temp_file()
    
    def get_rooted_tree_file(self):
        """获取置根后的树文件路径"""
        if self.temp_tree_file:
            return self.temp_tree_file.name
        return None
    
    def get_rooted_newick(self):
        """获取置根后的Newick字符串"""
        return self.rooted_newick
    
    def cleanup_temp_file(self):
        """清理临时文件"""
        if self.temp_tree_file:
            temp_path = self.temp_tree_file.name if hasattr(self.temp_tree_file, 'name') else None
            if temp_path and os.path.exists(temp_path):
                try:
                    os.unlink(temp_path)
                except Exception as e:
                    print(f"Failed to delete temporary file: {e}")
            self.temp_tree_file = None
    
    def closeEvent(self, event):
        """关闭对话框时的清理"""
        self.cleanup_temp_file()
        super().closeEvent(event)
    
    def reject(self):
        """取消时的清理"""
        self.cleanup_temp_file()
        super().reject()