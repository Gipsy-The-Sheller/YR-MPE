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
from pathlib import Path

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class DecentTreeThread(BaseProcessThread):
    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters,
                 imported_files=None, workdir=None):
        super().__init__(
            tool_path, input_files, parameters, imported_files, workdir=workdir
        )
    
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
                # PyQt prepends Qt's bin directory to PATH. DecentTree links
                # its own MinGW runtime and crashes with 0xC0000005 when that
                # directory shadows the DLL search order, so run it with a
                # PATH that omits PyQt's private Qt bin.
                env = os.environ.copy()
                env["PATH"] = os.pathsep.join(
                    part for part in env.get("PATH", "").split(os.pathsep)
                    if "PyQt5" not in part.replace("\\", "/")
                )
                result = self.execute_command(
                    cmd, cwd=str(Path(input_file).parent), env=env
                )
                detail = (result.stderr or "").strip()
                if result.returncode != 0:
                    raise RuntimeError(
                        "DecentTree failed for %s%s"
                        % (os.path.basename(input_file),
                           ": " + detail if detail else "")
                    )
                tree_file = input_file + ".nwk"
                if os.path.exists(tree_file):
                    output_files.append(tree_file)
            self.progress.emit("Tree construction completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])
        except Exception as exc:
            self.error.emit("Tree construction exception: %s" % exc)


class DecentTreePlugin(BasePlugin):
    export_phylogeny_result_signal = pyqtSignal(dict)
    
    def __init__(self, import_from=None, import_data=None, **kwargs):
        self._payload = None
        super().__init__(import_from, import_data, **kwargs)
        if import_data is not None:
            self.handle_import_data(import_data)
        self.run_button.setText("Build tree")

    def init_plugin_info(self):
        self.plugin_name = "DecentTree Distance Methods"
        self.tool_name = "DecentTree"
        self.citation = []
        self.citation_keys = ["wang2023decenttree"]
        self.input_types = {"ML Distance Matrix": ["mldist", "dist"]}
        self.output_types = {"Newick Tree": ".nwk"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')

    def control_bindings(self):
        """GUI controls mapped to the DecentTree command line."""
        return [
            {"id": "input", "widget": "file_path_edit",
             "label": "Distance matrix", "kind": "path",
             "cli": ["-in", "{value}"],
             "notes": "One .mldist/.dist file per run; batch mode repeats the "
                      "command per input."},
            {"id": "algorithm", "widget": "algorithm_combo",
             "label": "Algorithm", "kind": "choice",
             "choices": {
                 "BIONJ-R (Rapid BioNJ)": "BIONJ-R",
                 "BIONJ-V (Vectorized BIONJ)": "BIONJ-V",
                 "NJ-R-D (Double Precision Rapid NJ)": "NJ-R-D",
                 "ONJ-R (Rival Rapid NJ)": "ONJ-R",
                 "ONJ-R-V (Rival Rapid NJ Vectorized)": "ONJ-R-V",
                 "AUCTION (Auction Joining)": "AUCTION",
             },
             "cli": ["-t", "{value}"]},
        ]

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
        # The current Windows build of DecentTree segfaults in the classic
        # BIONJ/NJ paths; keep the variants verified on this build.
        self.algorithm_combo.addItems([
            "BIONJ-R (Rapid BioNJ)",
            "BIONJ-V (Vectorized BIONJ)",
            "NJ-R-D (Double Precision Rapid NJ)",
            "ONJ-R (Rival Rapid NJ)",
            "ONJ-R-V (Rival Rapid NJ Vectorized)",
            "AUCTION (Auction Joining)",
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
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.output_info = QLabel("No tree available yet.")
        self.output_info.setWordWrap(True)
        self.output_info.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.output_info)
        self.output_preview = QTextEdit()
        self.output_preview.setReadOnly(True)
        self.output_preview.setFont(QFont("Courier", 10))
        layout.addWidget(self.output_preview, 1)

    def setup_control_panel(self):
        super().setup_control_panel()
        self.import_to_platform_btn = QPushButton("Send to workspace")
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
            self.tool_path, input_files, self.get_parameters(), self.imported_files,
            workdir=self.workdir,
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()
    
    def _on_result(self, payload):
        self._payload = payload

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.current_output_files = output_files
        self.display_results(output_files)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message(f"Tree construction completed! Found {len(output_files)} result file(s)", "info")
        if self.import_from == "YR_MPEA" or self.result_has_consumer():
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
        if hasattr(self, "analysis_thread") and self.analysis_thread.isRunning():
            thread = self.analysis_thread
            for signal in (
                thread.progress, thread.result_ready, thread.finished,
                thread.error, thread.console_output,
            ):
                try:
                    signal.disconnect()
                except (TypeError, RuntimeError):
                    pass
            thread.stop()
            thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Tree construction stopped by user.", "warning")

    def display_results(self, output_files):
        if not output_files:
            self.output_info.setText("No tree output.")
            self.output_preview.clear()
            return
        treefile = next((path for path in output_files if path.endswith(".nwk")), output_files[0])
        try:
            text = Path(treefile).read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            self.output_info.setText("Cannot read %s: %s" % (treefile, exc))
            self.output_preview.clear()
            return
        self.output_info.setText("Tree: %s" % Path(treefile).name)
        self.output_preview.setPlainText(text)

    def _build_result_spec(self):
        files = list(getattr(self, "current_output_files", []) or [])
        if not files:
            raise ValueError("No DecentTree output is available.")
        items = []
        for path in files:
            newick = Path(path).read_text(encoding="utf-8", errors="replace").strip()
            if not newick:
                raise ValueError("Tree file is empty: %s" % path)
            items.append(DatasetItemSpec(
                loci_name=Path(path).stem,
                item_type="phylogeny",
                file_path=str(path),
                format="newick",
                data={"newick": newick, "method": self.algorithm_combo.currentText()},
                metadata={"source": self.plugin_name, "tool": "DecentTree"},
            ))
        dataset = DatasetSpec(
            name="DecentTree tree" if len(items) == 1 else "DecentTree trees",
            items=items,
            description="%d distance-based tree(s)" % len(items),
            source="decenttree",
            summary={"trees": len(items)},
            metadata={
                "input": self.imported_files[0] if len(self.imported_files) == 1 else "",
                "parameters": self.get_parameters(),
                "algorithm": self.algorithm_combo.currentText(),
            },
        )
        return ResultSpec(source="decenttree", datasets=[dataset],
                          artifacts=[str(path) for path in files])

    def import_to_platform(self):
        if not getattr(self, "current_output_files", None):
            QMessageBox.warning(self, "No result", "Build a tree first.")
            return
        if self.result_has_consumer():
            try:
                result = self._build_result_spec()
            except Exception as exc:
                QMessageBox.critical(self, "Export failed", str(exc))
                return
            if self.export_result(result):
                self.import_to_platform_btn.setEnabled(False)
                self.import_to_platform_btn.setText("Sent to workspace")
                QMessageBox.information(
                    self, "Success", "Tree sent to the current workspace."
                )
            return
        phylogenies = []
        for path in self.current_output_files:
            if path.endswith(".nwk"):
                phylogenies.append({
                    "filename": os.path.basename(path),
                    "content": Path(path).read_text(encoding="utf-8"),
                })
        if phylogenies:
            self.export_phylogeny_result_signal.emit(
                {"type": "phylogeny", "data": phylogenies}
            )

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
    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return DecentTreePlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )