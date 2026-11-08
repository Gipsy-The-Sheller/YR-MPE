# lsd2_plugin.py
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

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
                            QMessageBox, QLabel, QFileDialog, QDialog, QTextEdit)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
import os
import tempfile

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.nexus_tree import read_tree_file
from .components.lsd2_ui import LSD2UI
from .components.icytree_dialog import IcyTreeDialog
from .components.tip_dating_dialog import TipDatingDialog
from .components.md_mrca import MDMRCA


class LSD2Plugin(BasePlugin):
    """LSD2分子钟定年插件主类"""
    
    def __init__(self, import_from=None, import_data=None, **kwargs):
        # 首先创建控制面板
        self.control_panel = LSD2UI()
        
        # 设置回调函数
        self.control_panel.open_icytree_dialog_callback = self.open_icytree_dialog
        self.control_panel.open_tip_dating_dialog_callback = self.open_tip_dating_dialog
        self.control_panel.open_md_mrca_dialog_callback = self.open_md_mrca_dialog
        
        # 调用父类初始化
        super().__init__(import_from, import_data, **kwargs)
    
    def init_plugin_info(self):
        """初始化插件信息"""
        # 必须在调用父类初始化之前设置plugin_path，指向YR_MPE目录
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
        
        self.plugin_name = "LSD2"
        self.tool_name = "LSD2"
        self.input_types = {
            "Newick Tree": ["nwk", "tree", "tre", "newick"],
            "FASTA": ["fas", "fna", "fa", "fasta"]  # 可选，用于计算序列长度
        }
        self.output_types = {
            "Newick Tree": ".nwk",
            "Result File": ".result"
        }
        self.citation = []  # citation text lives in YR_MPE/citations/*.bib
        self.citation_keys = ["to2016fast"]

    def control_bindings(self):
        """GUI controls mapped to the LSD2 command line."""
        return [
            {"id": "tree", "widget": "control_panel.tree_file_path",
             "label": "Input tree", "kind": "path",
             "cli": ["-i", "{value}"],
             "notes": "Rooted or unrooted Newick tree."},
            {"id": "sequence_file", "widget": "control_panel.seq_file_path",
             "label": "Sequence file", "kind": "path", "cli": [],
             "notes": "GUI-only: used to compute sequence length."},
            {"id": "sequence_length_mode",
             "widget": "control_panel.seq_length_combo",
             "label": "Sequence length mode", "kind": "choice", "cli": [],
             "notes": "Disabled/Manual/From sequence file; the runner always "
                      "passes -s with the resolved length (1000 if disabled)."},
            {"id": "sequence_length",
             "widget": "control_panel.sequence_length_spin",
             "label": "Manual sequence length", "kind": "int",
             "cli": ["-s", "{value}"],
             "enabled_when": {"widget": "control_panel.seq_length_combo",
                              "equals": "Manual"}},
            {"id": "root_chosen",
             "widget": "control_panel.use_chosen_root_radio",
             "label": "Use chosen root", "kind": "bool",
             "cli": ["-r", "l"],
             "notes": "Requires a rooted tree loaded from IcyTree."},
            {"id": "root_estimate",
             "widget": "control_panel.estimate_root_radio",
             "label": "Estimate the root", "kind": "bool",
             "cli": ["-r", "a"]},
            {"id": "confidence_interval",
             "widget": "control_panel.calculate_ci_checkbox",
             "label": "Confidence intervals", "kind": "bool",
             "cli": ["-f", "100"],
             "notes": "The plugin uses 100 bootstrap trees."},
            {"id": "lognormal_std",
             "widget": "control_panel.lognormal_std_edit",
             "label": "Lognormal std", "kind": "float",
             "cli": ["-q", "{value}"],
             "enabled_by": "control_panel.calculate_ci_checkbox"},
            {"id": "calibrations",
             "widget": "control_panel.calibration_table",
             "label": "Node calibrations", "kind": "table", "cli": [],
             "notes": "Written to a temporary calibration file and passed "
                      "with -d."},
            {"id": "tip_dates",
             "widget": "control_panel.use_tip_dates_checkbox",
             "label": "Use tip dates", "kind": "bool", "cli": [],
             "notes": "Tip calibrations are merged into the -d file."},
        ]

    def config(self):
        """检查插件配置"""
        # 调用父类的config()方法检查LSD2
        return super().config()
    
    def init_ui(self):
        """初始化UI"""
        # 调用父类的init_ui方法来初始化基本的UI结构
        super().init_ui()
        
        # 在input_tab中添加控制面板
        self.input_tab.layout().addWidget(self.control_panel)
    
    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        # 控制面板在init_ui中添加
    
    def setup_output_tab(self):
        """设置输出标签页"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # 添加说明标签
        self.output_info_label = QLabel("Run LSD2 analysis to see results here.")
        self.output_info_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.output_info_label)
        self.output_preview = QTextEdit()
        self.output_preview.setReadOnly(True)
        self.output_preview.setFont(QFont("Courier", 10))
        layout.addWidget(self.output_preview)
        self.import_to_platform_btn = QPushButton("Send to workspace")
        self.import_to_platform_btn.setVisible(False)
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        layout.addWidget(self.import_to_platform_btn)
    
    def open_icytree_dialog(self, newick_str):
        """打开IcyTree置根对话框"""
        dialog = IcyTreeDialog(newick_str, self.plugin_path, parent=self)
        
        if dialog.exec_() == QDialog.Accepted:
            # 获取置根后的树文件
            rooted_tree_file = dialog.get_rooted_tree_file()
            # 检查临时文件对象
            if rooted_tree_file:
                # 使用 .name 属性获取文件路径
                rooted_tree_path = rooted_tree_file.name if hasattr(rooted_tree_file, 'name') else rooted_tree_file
                if os.path.exists(rooted_tree_path):
                    # 通知UI置根完成
                    self.control_panel.on_rooted_tree_loaded(rooted_tree_path)
                    self.add_console_message("Tree has been rooted successfully.", "success")
                else:
                    QMessageBox.warning(self, "Error", "Failed to save rooted tree.")
            else:
                QMessageBox.warning(self, "Error", "No rooted tree file available.")
        
        # 对话框会自动清理临时文件
    
    def open_tip_dating_dialog(self, otu_list):
        """打开Tip Dating配置对话框"""
        # 获取已有的tip标定数据
        existing_calibrations = self.control_panel.get_tip_calibrations()
        
        dialog = TipDatingDialog(otu_list, existing_calibrations, parent=self)
        
        if dialog.exec_() == QDialog.Accepted:
            # 获取新的tip标定数据
            tip_calibrations = dialog.get_tip_calibrations()
            
            # 通知UI更新
            self.control_panel.on_tip_calibrations_updated(tip_calibrations)
            
            configured_count = len([cal for cal in tip_calibrations.values() if cal is not None])
            self.add_console_message(f"Tip dating configured: {configured_count} tips calibrated.", "success")
    
    def open_md_mrca_dialog(self, row, newick_str):
        """打开MD-MRCA对话框用于校准点配置"""
        md_mrca_dialog = MDMRCA(parent=self)
        md_mrca_dialog.newick_string = newick_str
        md_mrca_dialog.annotated_newick_str = newick_str
        
        # 加载数据
        md_mrca_dialog.load_newick_data()
        
        # 以模态对话框方式显示
        if md_mrca_dialog.exec_() == QDialog.Accepted:
            # 用户点击了Apply按钮
            selected_taxa = [md_mrca_dialog.selected_taxa_list.item(i).text() 
                           for i in range(md_mrca_dialog.selected_taxa_list.count())]
            taxon_set_name = md_mrca_dialog.taxon_set_name.text().strip()
            
            if taxon_set_name and selected_taxa:
                # 获取校准类型
                cal_type = md_mrca_dialog.tmrca_type_combo.currentText()
                
                # 获取校准参数值
                cal_values = []
                try:
                    if cal_type == 'Point':
                        cal_values = [float(md_mrca_dialog.point_value.text())]
                    elif cal_type == 'Uniform':
                        cal_values = [
                            float(md_mrca_dialog.uniform_lower.text()),
                            float(md_mrca_dialog.uniform_upper.text())
                        ]
                    elif cal_type == 'Upper Boundary':
                        cal_values = [float(md_mrca_dialog.upper_boundary.text())]
                    elif cal_type == 'Lower Boundary':
                        cal_values = [float(md_mrca_dialog.lower_boundary.text())]
                    elif cal_type == 'Normal':
                        cal_values = [
                            float(md_mrca_dialog.normal_mean.text()),
                            float(md_mrca_dialog.normal_std.text())
                        ]
                    elif cal_type == 'Lognormal':
                        cal_values = [
                            float(md_mrca_dialog.lognormal_mean.text()),
                            float(md_mrca_dialog.lognormal_std.text())
                        ]
                except ValueError as e:
                    QMessageBox.warning(self, "Invalid Parameters", 
                                      f"Please check the calibration parameters: {str(e)}")
                    return
                
                # 构建结果字典
                mrca_result = {
                    'selected_taxa': selected_taxa,
                    'taxon_set_name': taxon_set_name,
                    'cal_type': cal_type,
                    'cal_values': cal_values
                }
                
                # 通知UI更新
                self.control_panel.on_md_mrca_result(row, mrca_result)
        
        # 对话框会自动销毁
        md_mrca_dialog.deleteLater()
    
    def run_analysis(self):
        """运行LSD2分析"""
        # 获取树文件
        tree_file = None
        if self.control_panel.get_use_chosen_root():
            tree_file = self.control_panel.get_rooted_tree_file()
        else:
            tree_file = self.control_panel.get_original_tree_file()
        
        if not tree_file or not os.path.exists(tree_file):
            QMessageBox.warning(self, "No Input", "Please load a tree file first.")
            return
        
        # 获取序列长度（LSD2总是需要此参数，默认1000）
        seq_length = self.control_panel.get_sequence_length()
        if seq_length <= 0:
            seq_length = 1000  # LSD2默认值
        
        # 获取校准点数据
        node_calibrations = self.control_panel.get_calibration_data()
        tip_calibrations = self.control_panel.get_tip_calibrations()
        
        # 检查是否有已配置的校准点（过滤掉 None 值）
        has_node_calibrations = any(cal is not None for cal in node_calibrations)
        has_tip_calibrations = any(cal is not None for cal in tip_calibrations.values()) if tip_calibrations else False
        
        if not has_node_calibrations and not has_tip_calibrations:
            reply = QMessageBox.question(
                self, "No Calibrations",
                "No calibration points have been specified. LSD2 requires at least one calibration point.\n"
                "Do you want to continue anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                return
        
        # 构建LSD2参数
        parameters = []
        
        # 序列长度（总是添加）
        parameters.extend(["-s", str(seq_length)])
        
        # 根设置
        if self.control_panel.get_use_chosen_root():
            # 使用用户选择的根（置根后的树），在局部重新估计
            parameters.extend(["-r", "l"])
        else:
            # 在所有分支上估计根
            parameters.extend(["-r", "a"])
        
        # 置信区间参数
        if self.control_panel.get_calculate_ci():
            parameters.extend(["-f", "100"])  # 100个bootstrap树
            lognormal_std = self.control_panel.get_lognormal_std()
            parameters.extend(["-q", str(lognormal_std)])
        
        # 生成校准文件（使用 -d 参数）
        if node_calibrations or tip_calibrations:
            try:
                from .components.methods.lsd2_logistics import form_calibration_table_with_tips
                calibration_content = form_calibration_table_with_tips(node_calibrations, tip_calibrations)
                
                # 调试：打印校准文件内容
                self.add_console_message(f"Calibration file content:\n{calibration_content}", "info")
                
                # 根据操作系统确定行结束符
                import platform
                os_name = platform.system().lower()
                
                # LSD2 需要特定平台的行结束符
                # 在 Windows (win32) 上，LSD2 期望 CRLF (\r\n)
                # 在 Linux/Darwin 上，LSD2 可以使用 LF (\n)
                if os_name == 'windows':
                    # Windows 平台：使用 CRLF
                    line_ending = '\r\n'
                else:
                    # Linux 或 macOS 平台：使用 LF
                    line_ending = '\n'
                
                self.add_console_message(f"Platform: {platform.system()}, using line ending: {repr(line_ending)}", "info")
                
                # 保存校准文件到临时文件
                temp_calib_file = self.create_temp_file(suffix='.txt')
                with open(temp_calib_file, 'w', newline='') as f:
                    # 使用指定的行结束符
                    lines = calibration_content.split('\n')
                    content = line_ending.join(lines)
                    # 确保文件末尾有换行符
                    if not content.endswith(line_ending):
                        content += line_ending
                    f.write(content)
                
                # 调试：验证文件内容
                with open(temp_calib_file, 'r') as f:
                    first_line = f.readline().strip()
                    self.add_console_message(f"First line of calibration file: '{first_line}'", "info")
                
                # 将校准文件参数添加到命令中（使用 -d 参数）
                parameters.extend(["-d", temp_calib_file])
                
                total_calibrations = len(node_calibrations) + len([cal for cal in tip_calibrations.values() if cal is not None])
                self.add_console_message(f"Generated calibration file with {total_calibrations} calibration points", "info")
                
            except Exception as e:
                QMessageBox.warning(self, "Calibration Error", 
                                  f"Failed to generate calibration file: {str(e)}")
                import traceback
                traceback.print_exc()
                return
        
        # 启动LSD2线程
        self.lsd2_thread = LSD2Thread(
            self.tool_path, [tree_file], parameters, workdir=self.workdir
        )
        self.lsd2_thread.progress.connect(self.update_progress)
        self.lsd2_thread.result_ready.connect(self._on_result)
        self.lsd2_thread.console_output.connect(self.add_console_message)
        self.lsd2_thread.finished.connect(self.on_lsd2_finished)
        self.lsd2_thread.error.connect(self.on_lsd2_error)
        
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.lsd2_thread.start()
    
    def _on_result(self, payload):
        self._result_payload = payload

    def stop_analysis(self):
        """停止分析"""
        thread = getattr(self, 'lsd2_thread', None)
        if thread is not None and thread.isRunning():
            for signal_name in ("progress", "result_ready", "finished",
                                "error", "console_output"):
                signal = getattr(thread, signal_name, None)
                if signal is None:
                    continue
                try:
                    signal.disconnect()
                except (TypeError, RuntimeError):
                    pass
            thread.stop()
            thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
    
    def update_progress(self, message):
        """更新进度"""
        self.add_console_message(message, "info")
    
    def on_lsd2_finished(self, output_files, html_files):
        """LSD2完成回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        
        if output_files:
            self.imported_files.extend(output_files)
            self.current_output_files = output_files
            self.add_console_message(f"LSD2 analysis completed. Generated {len(output_files)} files.", "success")
            self.parse_lsd2_results(output_files)
            if (self.import_from in ["YR_MPEA", "DATASET_MANAGER"]
                    or self.result_has_consumer()):
                self.import_to_platform_btn.setVisible(True)
        else:
            self.add_console_message("LSD2 analysis completed but no output files found.", "warning")
    
    def on_lsd2_error(self, error_message):
        """LSD2错误回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.add_console_message(f"Error: {error_message}", "error")
        QMessageBox.critical(self, "LSD2 Error", error_message)
    
    def parse_lsd2_results(self, output_files):
        """解析LSD2结果文件"""
        for output_file in output_files:
            if output_file.endswith('.result'):
                self.parse_lsd2_result_file(output_file)
            elif output_file.endswith('_dated.nwk') or output_file.endswith('.nwk'):
                self.parse_lsd2_tree_file(output_file)
    
    def parse_lsd2_result_file(self, result_file):
        """解析LSD2的.result文件"""
        try:
            with open(result_file, 'r') as f:
                content = f.read()
            
            self.add_console_message(f"\n{'='*60}", "info")
            self.add_console_message(f"LSD2 Results: {os.path.basename(result_file)}", "info")
            self.add_console_message(f"{'='*60}", "info")
            
            # 解析关键信息
            lines = content.split('\n')
            for line in lines:
                if line.strip():
                    self.add_console_message(line, "output")
            
        except Exception as e:
            self.add_console_message(f"Failed to parse result file: {str(e)}", "error")
    
    def parse_lsd2_tree_file(self, tree_file):
        """解析LSD2的时间树文件"""
        try:
            with open(tree_file, 'r') as f:
                newick_str = f.read().strip()
            
            self.add_console_message(f"\n{'='*60}", "info")
            self.add_console_message(f"Time Tree: {os.path.basename(tree_file)}", "info")
            self.add_console_message(f"{'='*60}", "info")
            self.add_console_message(newick_str, "output")
            
            # 显示时间树
            self.visualize_time_tree(newick_str)
            
        except Exception as e:
            self.add_console_message(f"Failed to parse tree file: {str(e)}", "error")
    
    def visualize_time_tree(self, newick_str):
        """Show the dated tree in the native output preview."""
        self.tab_widget.setCurrentWidget(self.output_tab)
        self.output_info_label.setText("LSD2 time tree")
        self.output_preview.setPlainText(newick_str)

    def _build_result_spec(self):
        files = list(getattr(self, "current_output_files", []) or [])
        if not files:
            raise ValueError("No LSD2 output is available.")
        items = []
        artifacts = []
        for path in files:
            name = os.path.basename(path)
            if name.endswith((".nwk", ".tre", ".tree", ".newick")):
                trees = read_tree_file(path)
                if not trees:
                    raise ValueError("No tree could be parsed from %s." % name)
                items.append(DatasetItemSpec(
                    loci_name=os.path.splitext(name)[0],
                    item_type="phylogeny",
                    file_path=str(path),
                    format="newick",
                    data={"newick": trees[0], "dated": True},
                    metadata={"source": self.plugin_name, "tool": "LSD2"},
                ))
            else:
                artifacts.append(str(path))
        if not items:
            raise ValueError("LSD2 produced no dated tree item.")
        dataset = DatasetSpec(
            name="LSD2 time tree",
            items=items,
            description="%d dated tree item(s)" % len(items),
            source="lsd2",
            summary={"items": len(items)},
            metadata={"input": self.control_panel.get_original_tree_file()},
        )
        return ResultSpec(source="lsd2", datasets=[dataset], artifacts=artifacts)

    def import_to_platform(self):
        if not getattr(self, "current_output_files", None):
            QMessageBox.warning(self, "No result", "Run LSD2 first.")
            return
        if not self.result_has_consumer():
            QMessageBox.warning(
                self, "No workspace",
                "This window is not connected to a YR-MPE workspace."
            )
            return
        try:
            result = self._build_result_spec()
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        if self.export_result(result):
            self.import_to_platform_btn.setEnabled(False)
            self.import_to_platform_btn.setText("Sent to workspace")
            QMessageBox.information(
                self, "Success", "LSD2 time tree sent to the current workspace."
            )


class LSD2Thread(BaseProcessThread):
    """LSD2分子钟定年线程类"""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters, imported_files=None,
                 workdir=None):
        super().__init__(tool_path, input_files, parameters, imported_files,
                         workdir=workdir)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "LSD2 Molecular Dating"
        
    def execute_commands(self):
        """执行LSD2分子钟定年命令"""
        try:
            output_files = []
            
            if len(self.input_files) < 1:
                self.error.emit("At least one tree file is required for LSD2")
                return
            
            tree_file = self.input_files[0]
            self.progress.emit(f"Processing tree file: {os.path.basename(tree_file)}")
            self.console_output.emit(f"Processing tree file: {os.path.basename(tree_file)}", "info")
            
            # 构建命令
            base_name = os.path.splitext(tree_file)[0]
            file_dir = os.path.dirname(tree_file) or "."
            file_name = os.path.basename(tree_file)
            cmd = [
                self.tool_path,
                "-i", tree_file,
                *self.parameters
            ]

            self.console_output.emit(f"Executing LSD2 command: {' '.join(cmd)}", "command")

            # 执行命令
            result = self.execute_command(
                cmd, cwd=self.workdir if self.workdir else file_dir
            )

            if result.returncode != 0:
                detail = (result.stderr or "").strip()
                raise RuntimeError(
                    "LSD2 failed%s" % (": " + detail if detail else "")
                )

            # 查找生成的输出文件
            
            # LSD2主要输出文件
            result_file = base_name + ".result"
            if os.path.exists(result_file):
                output_files.append(result_file)
                self.console_output.emit(f"Found result file: {result_file}", "info")
            
            # 时间树文件 - 尝试多种可能的命名方式
            dated_tree_patterns = [
                tree_file + ".result.nwk",     # its1206.nwk.result.nwk
                tree_file + ".dated",          # its1206.nwk.dated
                tree_file + ".tree.dated",     # its1206.nwk.tree.dated
                base_name + "_dated.nwk",      # its1206_dated.nwk
                base_name + ".dated",          # its1206.dated
                os.path.join(file_dir, file_name + ".dated"),  # its1206.nwk.dated
            ]
            
            for dated_tree_file in dated_tree_patterns:
                if os.path.exists(dated_tree_file) and dated_tree_file not in output_files:
                    output_files.append(dated_tree_file)
                    self.console_output.emit(f"Found dated tree file: {dated_tree_file}", "info")
            
            if not output_files:
                raise RuntimeError(
                    "LSD2 produced no .result or dated tree output."
                )

            self.progress.emit("Molecular dating completed")
            self.result_ready.emit({"output_files": output_files})
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Exception during LSD2 execution: {str(e)}")


# 插件入口点
class LSD2PluginEntry:
    """LSD2插件入口点"""

    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return LSD2Plugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )
