# pdguide_plugin.py
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

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QTextEdit,
                             QTabWidget, QToolButton, QApplication, QFrame, QDoubleSpinBox,
                             QRadioButton, QButtonGroup)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QTextCursor, QIcon
import tempfile
import os
import re
import json
import datetime
from typing import List, Optional

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess
from Bio import SeqIO

from .components.pdguide_ui import PdGuideUI


class PdGuidePlugin(BasePlugin):
    """PD-Guide插件主类"""
    
    def __init__(self, import_from=None, import_data=None):
        # 首先创建控制面板
        self.control_panel = PdGuideUI()
        
        # 调用父类初始化（会调用self.init_plugin_info()、self.config()和self.init_ui()）
        super().__init__(import_from, import_data)
    
    def init_plugin_info(self):
        """初始化插件信息"""
        # 必须在调用父类初始化之前设置plugin_path，指向YR_MPE目录
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
        
        self.plugin_name = "PD-Guide"
        self.tool_name = "LSD2"  # 使用LSD2作为主要工具名称，这样BasePlugin的config()会自动查找LSD2
        self.input_types = {
            "FASTA": ["fas", "fna", "fa", "fasta"],
            "Newick Tree": ["tree", "nwk", "newick"]
        }
        self.output_types = {
            # "Newick Tree": ".nwk",
            # "JSON": ".json"
        }
        self.citation = []  # citation text lives in YR_MPE/citations/*.bib
        # PD-Guide composes LSD2 (from config) with IQ-TREE and MCMCTree.
        self.citation_keys = ["wong2026iqtree", "yang2007paml", "dosreis2011approximate"]
    
    def get_tool_path(self, tool_name):
        """获取工具路径（使用BasePlugin的标准方式）"""
        # LSD2路径已经通过父类config()方法设置在self.tool_path中
        if tool_name.lower() == "lsd2":
            return self.tool_path
        elif tool_name.lower() in ["iq-tree 3", "iqtree3", "iqtree"]:
            # IQ-TREE路径需要在配置文件中查找
            config_path = os.path.join(self.plugin_path, "config.json")
            if os.path.exists(config_path):
                try:
                    with open(config_path, 'r', encoding='utf-8') as f:
                        config_data = json.load(f)
                    for tool in config_data:
                        if tool.get("name").lower() in ["iq-tree 3", "iqtree3", "iqtree"]:
                            tool_path = os.path.join(self.plugin_path, tool["path"].lstrip("/").lstrip("\\"))
                            if os.path.exists(tool_path):
                                return tool_path
                except Exception as e:
                    print(f"Failed to load IQ-TREE config: {e}")
            return None
        else:
            return None
    
    def config(self):
        """检查插件配置"""
        # 首先调用BasePlugin的config()方法检查LSD2
        if not super().config():
            return False
        
        # 检查IQ-TREE是否可用
        iqtree_path = self.get_tool_path("iq-tree 3")
        if not iqtree_path or not os.path.exists(iqtree_path):
            return False
            
        return True
    
    def show_config_guide(self):
        """显示配置指南"""
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setWindowTitle("Configuration Required")
        msg.setText("PD-Guide requires LSD2 and IQ-TREE to be installed.")
        msg.setInformativeText(
            "Please configure the paths to LSD2 and IQ-TREE executables in the settings."
        )
        msg.exec_()
    
    def init_ui(self):
        """初始化UI"""
        # 调用父类的init_ui方法来初始化基本的UI结构（包括标签页和控制台）
        super().init_ui()
        
        # 在input_tab中添加控制面板
        self.input_tab.layout().addWidget(self.control_panel)
        
        # 连接信号
        self.control_panel.infer_tree_button.clicked.connect(self.on_infer_tree_clicked)
    
    def on_infer_tree_clicked(self):
        """处理Infer tree with seqs按钮点击"""
        if not self.imported_files:
            QMessageBox.warning(self, "No Input", "Please import sequence files first.")
            return
        
        # 获取IQ-TREE路径
        iqtree_path = self.get_tool_path("iq-tree 3")
        if not iqtree_path:
            QMessageBox.warning(self, "IQ-TREE Not Found", "Please configure IQ-TREE path in settings.")
            return
        
        # 构建参数
        parameters = []
        if self.control_panel.model_combo.currentText() != "Auto":
            parameters.extend(["-m", self.control_panel.model_combo.currentText()])
        
        if self.control_panel.bootstrap_spin.value() > 0:
            parameters.extend(["-b", str(self.control_panel.bootstrap_spin.value())])
        
        # 启动IQ-TREE线程
        self.iqtree_thread = IQTreeThread(iqtree_path, self.imported_files, parameters)
        self.iqtree_thread.progress.connect(self.update_progress)
        self.iqtree_thread.console_output.connect(self.append_console)
        self.iqtree_thread.finished.connect(self.on_iqtree_finished)
        self.iqtree_thread.error.connect(self.on_iqtree_error)
        
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.iqtree_thread.start()
    
    def on_run_lsd2_clicked(self):
        """处理Run LSD2按钮点击"""
        # 检查树文件
        tree_file = self.control_panel.tree_file_path.text().strip()
        if not tree_file:
            QMessageBox.warning(self, "No Input", "Please load a tree file first.")
            return
        
        if not os.path.exists(tree_file):
            QMessageBox.warning(self, "File Not Found", f"Tree file not found: {tree_file}")
            return
        
        # 获取LSD2路径
        lsd2_path = self.get_tool_path("lsd2")
        if not lsd2_path:
            QMessageBox.warning(self, "LSD2 Not Found", "Please configure LSD2 path in settings.")
            return
        
        # 构建LSD2参数
        parameters = self.build_lsd2_parameters()
        if not parameters:
            return
        
        # 启动LSD2线程，传入树文件作为输入
        self.lsd2_thread = LSD2Thread(lsd2_path, [tree_file], parameters)
        self.lsd2_thread.progress.connect(self.update_progress)
        self.lsd2_thread.console_output.connect(self.append_console)
        self.lsd2_thread.finished.connect(self.on_lsd2_finished)
        self.lsd2_thread.error.connect(self.on_lsd2_error)
        
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.lsd2_thread.start()
    
    def build_lsd2_parameters(self):
        """构建LSD2参数"""
        parameters = []
        
        # 检查树文件
        tree_file = self.control_panel.tree_file_path.text().strip()
        if not tree_file:
            QMessageBox.warning(self, "Tree File Required", "Please load a tree file first.")
            return None
        
        if not os.path.exists(tree_file):
            QMessageBox.warning(self, "File Not Found", f"Tree file not found: {tree_file}")
            return None
        
        # 序列长度（根据Combobox模式）
        seq_length_mode = self.control_panel.seq_length_combo.currentText()
        if seq_length_mode == 'Manual':
            seq_length = self.control_panel.sequence_length_spin.value()
            if seq_length > 0:
                parameters.extend(["-s", str(seq_length)])
        elif seq_length_mode == 'From sequence file':
            seq_file = self.control_panel.seq_file_path.text().strip()
            if seq_file and os.path.exists(seq_file):
                # 自动计算序列长度
                try:
                    from Bio import SeqIO
                    record = next(SeqIO.parse(seq_file, 'fasta'))
                    seq_length = len(record.seq)
                    if seq_length > 0:
                        parameters.extend(["-s", str(seq_length)])
                except Exception as e:
                    QMessageBox.warning(self, "Sequence File Error", 
                                      f"Failed to read sequence file: {str(e)}")
                    return None
            else:
                QMessageBox.warning(self, "Sequence File Required", 
                                  "Please select a sequence file for sequence length calculation.")
                return None
        
        # 根设置
        if self.control_panel.estimate_root_radio.isChecked():
            parameters.append("-r")  # 估计根
        else:
            # 使用指定的根（需要从选择中获取）
            # 这里暂时简化，后续需要添加根选择UI
            pass
        
        # 标记定年（tip dating）
        # 如果需要tip dating，添加 -d 参数
        # parameters.append("-d")
        
        # 生成校准文件
        calibration_data = self.control_panel.get_calibration_data()
        if calibration_data:
            try:
                from .components.methods.lsd2_logistics import form_calibration_table
                calibration_content = form_calibration_table(calibration_data)
                
                # 保存校准文件到临时文件
                temp_calib_file = self.create_temp_file(suffix='.txt')
                with open(temp_calib_file, 'w') as f:
                    f.write(calibration_content)
                
                # 将校准文件参数添加到命令中
                parameters.extend(["-c", temp_calib_file])
                
                # 存储校准文件路径以便后续清理
                self.temp_files.append(temp_calib_file)
                
                self.append_console(f"Generated calibration file with {len(calibration_data)} calibration points", "info")
                
            except Exception as e:
                QMessageBox.warning(self, "Calibration Error", 
                                  f"Failed to generate calibration file: {str(e)}")
                return None
        
        return parameters
    
    def run_analysis(self):
        """运行分析"""
        # 根据当前选择的方法决定运行哪个分析
        selected_method = self.control_panel.methods_combo.currentText()
        
        if selected_method == "LSD2":
            self.on_run_lsd2_clicked()
        elif selected_method == "PhyloBayes":
            self.add_console_message("PhyloBayes analysis not yet implemented", "warning")
        elif selected_method == "MrBayes":
            self.add_console_message("MrBayes analysis not yet implemented", "warning")
        elif selected_method == "MCMCTree":
            self.add_console_message("MCMCTree analysis not yet implemented", "warning")
        else:
            self.on_infer_tree_clicked()
    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'iqtree_thread') and self.iqtree_thread.isRunning():
            self.iqtree_thread.stop()
        if hasattr(self, 'lsd2_thread') and self.lsd2_thread.isRunning():
            self.lsd2_thread.stop()
        
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
    
    def update_progress(self, message):
        """更新进度"""
        self.add_console_message(message, "info")
    
    def append_console(self, message, level="info"):
        """添加控制台消息（兼容性方法）"""
        self.add_console_message(message, level)
    
    def on_iqtree_finished(self, output_files, html_files):
        """IQ-TREE完成回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        
        if output_files:
            self.imported_files.extend(output_files)
            self.append_console(f"IQ-TREE analysis completed. Generated {len(output_files)} files.", "success")
        else:
            self.append_console("IQ-TREE analysis completed but no output files found.", "warning")
    
    def on_iqtree_error(self, error_message):
        """IQ-TREE错误回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.append_console(f"Error: {error_message}", "error")
        QMessageBox.critical(self, "IQ-TREE Error", error_message)
    
    def on_lsd2_finished(self, output_files, html_files):
        """LSD2完成回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        
        if output_files:
            self.imported_files.extend(output_files)
            self.append_console(f"LSD2 analysis completed. Generated {len(output_files)} files.", "success")
            
            # 解析和显示结果
            self.parse_lsd2_results(output_files)
        else:
            self.append_console("LSD2 analysis completed but no output files found.", "warning")
    
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
            
            self.append_console(f"\n{'='*60}", "info")
            self.append_console(f"LSD2 Results: {os.path.basename(result_file)}", "info")
            self.append_console(f"{'='*60}", "info")
            
            # 解析关键信息
            lines = content.split('\n')
            for line in lines:
                if line.strip():
                    self.append_console(line, "output")
            
            # 提取似然值、节点年龄等信息
            # 这里可以添加更详细的解析逻辑
            
        except Exception as e:
            self.append_console(f"Failed to parse result file: {str(e)}", "error")
    
    def parse_lsd2_tree_file(self, tree_file):
        """解析LSD2的时间树文件"""
        try:
            with open(tree_file, 'r') as f:
                newick_str = f.read().strip()
            
            self.append_console(f"\n{'='*60}", "info")
            self.append_console(f"Time Tree: {os.path.basename(tree_file)}", "info")
            self.append_console(f"{'='*60}", "info")
            self.append_console(newick_str, "output")
            
            # 显示时间树
            self.visualize_time_tree(newick_str)
            
        except Exception as e:
            self.append_console(f"Failed to parse tree file: {str(e)}", "error")
    
    def visualize_time_tree(self, newick_str):
        """可视化时间树"""
        try:
            # 切换到输出标签页
            self.tabs.setCurrentWidget(self.output_tab)
            
            # 使用MplTreeView显示时间树
            from .components.mpl_treeview import MplTreeView
            
            # 创建或更新树视图
            if not hasattr(self, 'time_tree_canvas'):
                self.time_tree_canvas = MplTreeView()
                self.output_tab.layout().addWidget(self.time_tree_canvas)
            
            # 绘制时间树（带时间标尺）
            self.time_tree_canvas.draw_dendrogram(
                newick_str=newick_str,
                show_time_scale=True,
                title="Time-Scaled Phylogenetic Tree"
            )
            
            self.append_console("Time tree visualization completed.", "success")
            
        except ImportError as e:
            self.append_console(f"Could not import MplTreeView: {str(e)}", "warning")
        except Exception as e:
            self.append_console(f"Failed to visualize time tree: {str(e)}", "error")
    
    def on_lsd2_error(self, error_message):
        """LSD2错误回调"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.append_console(f"Error: {error_message}", "error")
        QMessageBox.critical(self, "LSD2 Error", error_message)


class IQTreeThread(BaseProcessThread):
    """IQ-TREE系统发育推断线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "IQ-TREE Phylogeny"
        
    def execute_commands(self):
        """执行IQ-TREE系统发育推断命令"""
        try:
            output_files = []
            
            # 分别处理每个输入文件
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                    
                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", "info")
                
                # 构建命令
                cmd = [
                    self.tool_path,
                    "-s", input_file,
                    *self.parameters
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"IQ-TREE execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 查找生成的.treefile文件
                treefile = input_file + ".treefile"
                if os.path.exists(treefile):
                    output_files.append(treefile)
                else:
                    # 尝试其他可能的命名方式
                    alternate_treefile = input_file + ".phy.treefile"
                    if os.path.exists(alternate_treefile):
                        output_files.append(alternate_treefile)
                    else:
                        self.console_output.emit(f"Warning: Could not find .treefile for {input_file}", "warning")
                        
                # 查找生成的.iqtree文件
                iqtree_file = input_file + ".iqtree"
                if os.path.exists(iqtree_file):
                    output_files.append(iqtree_file)
                
                # 查找生成的.log文件
                log_file = input_file + ".log"
                if os.path.exists(log_file):
                    output_files.append(log_file)
            
            self.progress.emit("Phylogenetic inference completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Exception during IQ-TREE execution: {str(e)}")


class LSD2Thread(BaseProcessThread):
    """LSD2分子钟定年线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "LSD2 Molecular Dating"
        
    def execute_commands(self):
        """执行LSD2分子钟定年命令"""
        try:
            output_files = []
            
            # LSD2通常需要一个树文件和一个校准文件
            if len(self.input_files) < 1:
                self.error.emit("At least one tree file is required for LSD2")
                return
            
            tree_file = self.input_files[0]
            self.progress.emit(f"Processing tree file: {os.path.basename(tree_file)}")
            self.console_output.emit(f"Processing tree file: {os.path.basename(tree_file)}", "info")
            
            # 构建命令
            cmd = [
                self.tool_path,
                "-i", tree_file,
                *self.parameters
            ]
            
            self.console_output.emit(f"Executing LSD2 command: {' '.join(cmd)}", "command")
            
            # 执行命令
            result = self.execute_command(cmd)
            
            if result.returncode != 0:
                self.error.emit(f"LSD2 execution failed with return code {result.returncode}\n{result.stderr}")
                return
            
            # 查找生成的输出文件
            base_name = os.path.splitext(tree_file)[0]
            
            # LSD2主要输出文件
            result_file = base_name + ".result"
            if os.path.exists(result_file):
                output_files.append(result_file)
                self.console_output.emit(f"Found result file: {result_file}", "info")
            
            # 时间树文件
            dated_tree_file = base_name + "_dated.nwk"
            if os.path.exists(dated_tree_file):
                output_files.append(dated_tree_file)
                self.console_output.emit(f"Found dated tree file: {dated_tree_file}", "info")
            
            # 尝试其他可能的输出文件名
            other_possible_files = [
                tree_file + ".tree",
                tree_file + ".nwk",
                tree_file + "_result.nwk"
            ]
            
            for possible_file in other_possible_files:
                if os.path.exists(possible_file) and possible_file not in output_files:
                    output_files.append(possible_file)
                    self.console_output.emit(f"Found additional output file: {possible_file}", "info")
            
            if not output_files:
                self.console_output.emit(f"Warning: Could not find LSD2 output files for {tree_file}", "warning")
                self.console_output.emit("Please check the LSD2 output in the working directory.", "warning")
            
            self.progress.emit("Molecular dating completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Exception during LSD2 execution: {str(e)}")


# 插件入口点
class PdGuidePluginEntry:
    """PD-Guide插件入口点"""
    
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path
    
    def run(self, import_from=None, import_data=None):
        return PdGuidePlugin(import_from=import_from, import_data=import_data)
