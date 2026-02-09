# decenttree_plugin.py
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

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit,
                             QCheckBox, QLabel, QComboBox, QTextEdit,
                             QTabWidget, QFrame)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QIcon
import os

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class DecentTreeThread(BaseProcessThread):
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
        self.working_dir = os.path.dirname(tool_path)
    
    def get_tool_name(self):
        return "DecentTree"
        
    def execute_commands(self):
        try:
            output_files = []
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", "info")
                cmd = [self.tool_path, "-in", input_file, "-out", input_file + ".nwk", *self.parameters]
                result = self.execute_command(cmd, cwd=self.working_dir)
                if result.returncode != 0:
                    self.error.emit(f"DecentTree execution failed: {result.stderr}")
                    return
                tree_file = input_file + ".nwk"
                if os.path.exists(tree_file):
                    output_files.append(tree_file)
            self.progress.emit("Tree construction completed")
            self.finished.emit(output_files, [])
        except Exception as e:
            self.error.emit(f"Tree construction exception: {str(e)}")


class DecentTreePlugin(BasePlugin):
    export_phylogeny_result_signal = pyqtSignal(dict)
    
    def __init__(self, import_from=None, import_data=None):
        super().__init__(import_from, import_data)
        # Initialize temp_files list if not exists
        if not hasattr(self, 'temp_files'):
            self.temp_files = []
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
    
    def init_plugin_info(self):
        self.plugin_name = "DecentTree Distance Methods"
        self.tool_name = "DecentTree"
        self.citation = ["Weiwen Wang, James Barbetti, Thomas Wong, Bryan Thornlow, Russ Corbett-Detig, Yatish Turakhia, Robert Lanfear, Bui Quang Minh, DecentTree: scalable Neighbour-Joining for the genomic era, Bioinformatics, Volume 39, Issue 9, September 2023, btad536, doi: <a href='https://doi.org/10.1093/bioinformatics/btad536'>10.1093/bioinformatics/btad536</a>"]
        self.input_types = {"ML Distance Matrix": ["mldist", "dist"]}
        self.output_types = {"Newick Tree": ".nwk"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
    
    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        input_group = QGroupBox("Input")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select ML distance matrix files (.mldist)...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)
        
        # 添加文件标签容器到 input_group 内部
        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # 注意：imported_files 的处理将在 handle_import_data 中完成
        # 这里不在此处处理，避免重复操作和阻塞
        if self.import_file:
            self.file_path_edit.setText(self.import_file)
        
        params_group = QGroupBox("DecentTree Parameters")
        params_form_layout = QFormLayout()
        params_group.setLayout(params_form_layout)
        layout.addWidget(params_group)
        
        self.algorithm_combo = QComboBox()
        self.algorithm_combo.addItems([
            "BIONJ (BioNJ)", "BIONJ-R (Rapid BioNJ)", "BIONJ2009 (BioNJ with OMP)",
            "NJ (Neighbor Joining)", "NJ-R (Rapid NJ)", "NJ-R-D (Double Precision Rapid NJ)",
            "NJ-V (Vectorized NJ)", "RapidNJ (Fastest NJ)", "UNJ (Unweighted NJ)",
            "ONJ-R (Rival Rapid NJ)", "ONJ-R-V (Rival Rapid NJ Vectorized)",
            "AUCTION (Auction Joining)", "BENCHMARK (Benchmark)"
        ])
        params_form_layout.addRow("Algorithm:", self.algorithm_combo)
        
        layout.addStretch()
        
        if not hasattr(self, 'imported_files'):
            self.imported_files = []
        if not hasattr(self, 'file_tags'):
            self.file_tags = []
    
    def browse_files(self):
        file_filter = "ML Distance Matrix files (*.mldist *.dist);;All files (*)"
        file_paths, _ = QFileDialog.getOpenFileNames(self, "Select ML distance matrix files", "", file_filter)
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def add_file_tag(self, file_path):
        # 不要在这里添加到 imported_files，因为调用者已经添加过了
        # self.imported_files.append(file_path)
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.Box)
        tag_widget.setStyleSheet("QFrame { background-color: #e9ecef; border-radius: 15px; margin: 2px; }")
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)
        display_name = self.get_display_name(file_path)
        name_label = QLabel(display_name)
        name_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(name_label)
        close_btn = QPushButton("x")
        close_btn.setFixedSize(20, 20)
        close_btn.setStyleSheet("QPushButton { background-color: transparent; border: none; color: #6c757d; }")
        close_btn.clicked.connect(lambda: self.remove_file_tag(tag_widget, file_path))
        tag_layout.addWidget(close_btn)
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append(tag_widget)
        self.file_tags_container.setVisible(True)
    
    def remove_file_tag(self, tag_widget, file_path):
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        if tag_widget in self.file_tags:
            self.file_tags.remove(tag_widget)
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        if len(self.imported_files) == 0:
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
        elif len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        else:
            self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def get_display_name(self, file_path):
        base_name = os.path.basename(file_path)
        count = sum(1 for f in self.imported_files if os.path.basename(f) == base_name)
        if count > 1:
            name, ext = os.path.splitext(base_name)
            dir_name = os.path.basename(os.path.dirname(file_path))
            return f"{name} ({dir_name}){ext}"
        return base_name
    
    def setup_output_tab(self):
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # Add output info label
        self.output_info = QLabel("No tree available yet.")
        self.output_info.setAlignment(Qt.AlignCenter)
        self.output_info.setStyleSheet("color: #6c757d; padding: 20px;")
        layout.addWidget(self.output_info)
        
        self.output_preview = QTextEdit()
        self.output_preview.setReadOnly(True)
        self.output_preview.setFont(QFont("Courier", 10))
        layout.addWidget(QLabel("Phylogenetic Tree Preview (Newick Format):"))
        layout.addWidget(self.output_preview)
    
    def setup_control_panel(self):
        super().setup_control_panel()
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)
        self.control_layout.addWidget(self.import_to_platform_btn)
    
    def prepare_input_files(self):
        try:
            input_files = []
            if self.imported_files:
                for file_path in self.imported_files:
                    if os.path.exists(file_path):
                        input_files.append(file_path)
                return input_files
            elif self.import_file:
                return [self.import_file]
            QMessageBox.warning(self, "Warning", "Please provide ML distance matrix files!")
            return None
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
    
    def get_parameters(self):
        params = []
        algorithm_text = self.algorithm_combo.currentText()
        algorithm = algorithm_text.split(" (")[0]
        params.extend(["-t", algorithm])
        params.append("-no-banner")
        return params
    
    def run_analysis(self):
        if not self.imported_files and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide ML distance matrix files!")
            return
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "DecentTree executable file not found!")
            return
        self.add_console_message("Starting tree construction...", "info")
        input_files = self.prepare_input_files()
        if not input_files:
            return
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.analysis_thread = DecentTreeThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()
    
    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.current_output_files = output_files
        self.display_results(output_files)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message(f"Tree construction completed! Found {len(output_files)} result file(s)", "info")
        if self.import_from == "YR_MPEA":
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        QMessageBox.information(self, "Completed", "Tree construction completed!")
    
    def analysis_error(self, error_message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(f"Tree construction failed: {error_message}", "error")
        QMessageBox.critical(self, "Error", f"Tree construction failed: {error_message}")
    
    def stop_analysis(self):
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        QMessageBox.information(self, "Stopped", "Tree construction has been aborted.")
    
    def display_results(self, output_files):
        if not output_files:
            QMessageBox.warning(self, "Error", "No output files generated")
            return
        treefile = None
        for file in output_files:
            if file.endswith('.nwk'):
                treefile = file
                break
        if treefile:
            try:
                with open(treefile, 'r', encoding='utf-8') as f:
                    tree_content = f.read().strip()
                if not tree_content:
                    QMessageBox.warning(self, "Error", "Tree file is empty")
                    return
                from .icytree import IcyTreePlugin
                plugin_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '')
                icytree_plugin = IcyTreePlugin(plugin_path=plugin_path)
                
                # Connect to IcyTree status_changed signal to remove output_info when tree is loaded
                icytree_plugin.status_changed.connect(lambda status: self._on_icytree_status_changed(status, icytree_plugin))
                
                icytree_plugin.set_newick_string(tree_content)
                output_layout = self.output_tab.layout()
                if output_layout:
                    for i in reversed(range(output_layout.count())):
                        widget = output_layout.itemAt(i).widget()
                        if widget and widget != self.output_info:
                            widget.setParent(None)
                output_layout.addWidget(icytree_plugin)
            except ImportError:
                QMessageBox.warning(self, "Error", "IcyTree plugin not available")
            except Exception as e:
                error_msg = f"Error processing tree file: {str(e)}"
                QMessageBox.warning(self, "Error", error_msg)
                self.add_console_message(error_msg, "error")
        else:
            QMessageBox.warning(self, "Error", f"No tree file found. Generated {len(output_files)} file(s).")
    
    def _on_icytree_status_changed(self, status, icytree_plugin):
        """Handle IcyTree status changes"""
        if status == "Tree loaded to IcyTree":
            # Remove the "No tree available yet." label when tree is loaded
            if hasattr(self, 'output_info') and self.output_info:
                self.output_info.setParent(None)
                self.output_info = None
    
    def import_to_platform(self):
        if not hasattr(self, 'current_output_files') or not self.current_output_files:
            QMessageBox.warning(self, "Warning", "No phylogenetic results to import.")
            return
        try:
            phylogenies = []
            for output_file in self.current_output_files:
                if output_file.endswith('.nwk'):
                    with open(output_file, 'r', encoding='utf-8') as f:
                        content = f.read()
                        phylogenies.append({'filename': os.path.basename(output_file), 'content': content})
            if not phylogenies:
                QMessageBox.warning(self, "Warning", "No phylogenetic trees found in results.")
                return
            self.import_phylogenies_to_platform(phylogenies)
            QMessageBox.information(self, "Success", f"Successfully imported {len(phylogenies)} phylogenetic tree(s) to the platform.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import phylogenetic trees: {str(e)}")
    
    def import_phylogenies_to_platform(self, phylogenies):
        self.export_phylogeny_result_signal.emit({"type": "phylogeny", "data": phylogenies})
    
    def handle_import_data(self, import_data):
        if isinstance(import_data, list):
            if import_data:
                # 检查是否是字典格式（新格式：{'filename': '...', 'content': '...'}）
                if isinstance(import_data[0], dict):
                    # 处理带 filename 的字典
                    if 'filename' in import_data[0]:
                        if 'mldist' in import_data[0]['filename'].lower():
                            self.imported_files = [item['filename'] for item in import_data]
                            self.import_file = self.imported_files[0] if self.imported_files else None
                        else:
                            QMessageBox.warning(self, "Data Type Mismatch", "The imported data is not a distance matrix. Please use the DISTANCE menu to compute ML distances first.")
                    # 处理只有 content 没有 filename 的字典（新格式）
                    elif 'content' in import_data[0]:
                        self.imported_files = []
                        for item in import_data:
                            try:
                                temp_file = self.create_temp_file(suffix='.mldist')
                                with open(temp_file, 'w', encoding='utf-8') as f:
                                    f.write(str(item['content']))
                                self.temp_files.append(temp_file)
                                self.imported_files.append(temp_file)
                            except Exception as e:
                                continue
                        self.import_file = self.imported_files[0] if self.imported_files else None
                else:
                    QMessageBox.warning(self, "Warning", "Imported data format not supported.")
            else:
                QMessageBox.warning(self, "Warning", "No data imported.")
        elif isinstance(import_data, dict) and 'type' in import_data:
            if import_data['type'] == 'distance_matrix':
                if 'data' in import_data and isinstance(import_data['data'], list):
                    # Handle different data formats
                    self.imported_files = []
                    for item in import_data['data']:
                        try:
                            if isinstance(item, dict):
                                # 优先使用 content 创建临时文件（因为 filename 可能只是文件名，不是完整路径）
                                if 'content' in item:
                                    try:
                                        temp_file = self.create_temp_file(suffix='.mldist')
                                        with open(temp_file, 'w', encoding='utf-8') as f:
                                            f.write(str(item['content']))
                                        self.temp_files.append(temp_file)
                                        self.imported_files.append(temp_file)
                                    except Exception as e:
                                        continue
                                # 如果只有 filename 没有 content，尝试检查文件是否存在
                                elif 'filename' in item:
                                    # Check if file exists before adding
                                    if os.path.exists(item['filename']):
                                        self.imported_files.append(item['filename'])
                            elif isinstance(item, str) and os.path.exists(item):
                                self.imported_files.append(item)
                        except Exception as e:
                            continue
                    self.import_file = self.imported_files[0] if self.imported_files else None
            elif import_data['type'] == 'alignment':
                QMessageBox.warning(self, "Data Type Mismatch", "Alignment data detected. Please use the DISTANCE menu to compute ML distances first before building a tree.")
        else:
            self.import_file = None
            self.imported_files = []
        
        # 使用 QTimer.singleShot 延迟执行 UI 更新，避免阻塞
        if self.imported_files:
            from PyQt5.QtCore import QTimer
            QTimer.singleShot(0, self.update_ui_for_imported_files)
    
    def update_ui_for_imported_files(self):
        """更新 UI 以显示导入的文件"""
        if not self.imported_files:
            return
            
        if hasattr(self, 'file_path_edit'):
            # Update file path display
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            # Add file tags
            if hasattr(self, 'file_tags_layout') and hasattr(self, 'file_tags_container'):
                for file_path in self.imported_files:
                    self.add_file_tag(file_path)
                self.file_tags_container.setVisible(True)


class DecentTreePluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    def run(self, import_from=None, import_data=None):
        return DecentTreePlugin(import_from=import_from, import_data=import_data)