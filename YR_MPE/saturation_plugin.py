#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Nucleotide Substitution Saturation Analysis Plugin
核苷酸替代饱和度检测插件
"""

from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, 
                             QLabel, QPushButton, QFileDialog, QTextEdit, 
                             QProgressBar, QMessageBox, QGroupBox, QLineEdit,
                             QTabWidget, QScrollArea, QFrame, QPlainTextEdit,
                             QCheckBox, QSpinBox, QDoubleSpinBox, QComboBox,
                             QSizePolicy, QApplication, QToolButton)
from PyQt5.QtCore import QThread, pyqtSignal, Qt, QUrl, QTimer
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtGui import QFont, QTextCursor, QIcon, QPixmap
import os
import tempfile
import json
from pathlib import Path
import platform

try:
    from .saturation_analyzer import SaturationAnalyzer
    from .sequence_utils import SequenceUtils, SequenceItem
    from .sequence_models import SequenceType
except ImportError:
    from saturation_analyzer import SaturationAnalyzer
    from sequence_utils import SequenceUtils, SequenceItem
    from sequence_models import SequenceType


class SaturationAnalysisThread(QThread):
    """饱和度分析线程"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)  # message, msg_type
    
    def __init__(self, sequences, output_dir, iqtree_path, evolver_path, num_simulations=1000, confidence=0.95, test_type='two-sided'):
        super().__init__()
        self.sequences = sequences
        self.output_dir = output_dir
        self.iqtree_path = iqtree_path
        self.evolver_path = evolver_path
        self.analyzer = SaturationAnalyzer()
        self.num_simulations = int(num_simulations)
        self.confidence = float(confidence)
        self.test_type = str(test_type)
    
    def run(self):
        try:
            self.console_output.emit("Starting saturation analysis...", "info")
            self.progress.emit("Initializing analysis...")
            
            # 运行分析
            result = self.analyzer.analyze_saturation(
                self.sequences, 
                self.output_dir,
                self.iqtree_path,
                self.evolver_path,
                num_simulations=self.num_simulations,
                confidence=self.confidence,
                test_type=self.test_type
            )
            
            if result == -1:
                self.error.emit("IQ-TREE model selection failed")
                self.console_output.emit("IQ-TREE model selection failed", "error")
                return
            elif result == -2:
                self.error.emit("Evolver simulation failed")
                self.console_output.emit("Evolver simulation failed", "error")
                return
            
            self.console_output.emit("Analysis completed successfully!", "info")
            self.finished.emit(result)
            
        except Exception as e:
            self.error.emit(f"Analysis failed: {str(e)}")
            self.console_output.emit(f"Error: {str(e)}", "error")


class SaturationAnalysisWidget(QWidget):
    """饱和度分析主界面"""
    
    def __init__(self, importfrom=None, importdata=None):
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.imported_files = []
        self.file_tags = []
        self.analysis_result = None
        self.temp_files = []
        
        if importfrom == "YR_MPEA" and isinstance(importdata, list):
            # 从YR-MPEA导入数据
            self.imported_sequences = importdata
        else:
            self.imported_sequences = []
        
        self.init_ui()
    
    def init_ui(self):
        """初始化用户界面"""
        self.setWindowTitle("Nucleotide Substitution Saturation Analysis")
        self.setMinimumSize(1000, 700)
        
        # 创建主布局
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # 创建标签页
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # 输入标签页
        self.input_tab = QWidget()
        self.tab_widget.addTab(self.input_tab, "Input & Parameters")
        self.setup_input_tab()
        
        # 结果标签页
        self.result_tab = QWidget()
        self.tab_widget.addTab(self.result_tab, "Results")
        self.setup_result_tab()
        
        # 控制面板
        self.setup_control_panel(main_layout)
        
        # 初始化变量
        self.is_running = False
    
    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # 输入组
        input_group = QGroupBox("Input Sequences")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # 文件输入
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select FASTA files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)
        
        # 文件标签容器
        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # 文本输入
        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste FASTA format sequences...")
        self.sequence_text.setMaximumHeight(200)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        input_layout.addRow("Sequence text:", self.sequence_text)
        
        # 如果从YR-MPEA导入数据
        if self.imported_sequences:
            self.file_imported_button = QPushButton(f"Sequences imported from YR-MPEA ({len(self.imported_sequences)} sequences)")
            self.file_imported_button.setEnabled(False)
            layout.addWidget(self.file_imported_button)
            input_group.setVisible(False)
        
        # 参数组
        params_group = QGroupBox("Analysis Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # IQ-TREE路径
        iqtree_layout = QHBoxLayout()
        self.iqtree_path_edit = QLineEdit()
        self.iqtree_path_edit.setText("iqtree")  # 默认值
        self.iqtree_browse_btn = QPushButton("Browse")
        self.iqtree_browse_btn.clicked.connect(self.browse_iqtree)
        iqtree_layout.addWidget(self.iqtree_path_edit)
        iqtree_layout.addWidget(self.iqtree_browse_btn)
        params_layout.addRow("IQ-TREE path:", iqtree_layout)
        
        # Evolver路径
        evolver_layout = QHBoxLayout()
        self.evolver_path_edit = QLineEdit()
        self.evolver_path_edit.setText("evolver")  # 默认值
        self.evolver_browse_btn = QPushButton("Browse")
        self.evolver_browse_btn.clicked.connect(self.browse_evolver)
        evolver_layout.addWidget(self.evolver_path_edit)
        evolver_layout.addWidget(self.evolver_browse_btn)
        params_layout.addRow("Evolver path:", evolver_layout)
        
        # 输出目录
        output_layout = QHBoxLayout()
        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setText(os.path.join(tempfile.gettempdir(), "saturation_analysis"))
        self.output_browse_btn = QPushButton("Browse")
        self.output_browse_btn.clicked.connect(self.browse_output_dir)
        output_layout.addWidget(self.output_dir_edit)
        output_layout.addWidget(self.output_browse_btn)
        params_layout.addRow("Output directory:", output_layout)
        
        # 高级参数
        advanced_group = QGroupBox("Advanced Options")
        advanced_layout = QFormLayout()
        advanced_group.setLayout(advanced_layout)
        layout.addWidget(advanced_group)
        
        # 模拟次数
        self.simulations_spinbox = QSpinBox()
        self.simulations_spinbox.setRange(100, 10000)
        self.simulations_spinbox.setValue(1000)
        advanced_layout.addRow("Number of simulations:", self.simulations_spinbox)
        
        # 置信水平
        self.confidence_spinbox = QDoubleSpinBox()
        self.confidence_spinbox.setRange(0.5, 0.99)
        self.confidence_spinbox.setSingleStep(0.01)
        self.confidence_spinbox.setValue(0.97)
        advanced_layout.addRow("Confidence level:", self.confidence_spinbox)

        # 检验类型
        self.test_type_combo = QComboBox()
        self.test_type_combo.addItems(["two-sided", "left", "right"]) 
        self.test_type_combo.setCurrentText("two-sided")
        advanced_layout.addRow("Test type:", self.test_type_combo)
        
        # 方法说明
        method_group = QGroupBox("Method Description")
        method_layout = QVBoxLayout()
        method_group.setLayout(method_layout)
        
        method_text = QLabel("""
        <b>Index of Substitution Saturation (Iss) Analysis</b><br><br>
        
        This analysis is based on Xia (2003) method for detecting substitution saturation in nucleotide sequences.
        
        <b>Steps:</b><br>
        1. Complete deletion of gaps from alignment<br>
        2. Calculate information entropy of sequences<br>
        3. Select best substitution model using IQ-TREE with BIC<br>
        4. Generate evolver configuration file<br>
        5. Run 1000 simulations using evolver<br>
        6. Calculate 95% HPD interval for simulated Iss values<br>
        7. Compare observed Iss with simulated distribution<br><br>
        
        <b>Interpretation:</b><br>
        • If observed Iss falls within 95% HPD: sequences are NOT saturated<br>
        • If observed Iss falls outside 95% HPD: sequences are SATURATED<br><br>
        
        <b>Citation:</b><br>
        Xia, X. (2003). DAMBE: software package for data analysis in molecular biology and evolution. 
        <i>Journal of Heredity</i>, 94(4), 371-373.
        """)
        method_text.setWordWrap(True)
        method_text.setOpenExternalLinks(True)
        method_layout.addWidget(method_text)
        
        layout.addWidget(method_group)
        layout.addStretch()
    
    def setup_result_tab(self):
        """设置结果标签页"""
        layout = QVBoxLayout()
        self.result_tab.setLayout(layout)
        
        # 结果信息
        self.result_info = QLabel("Analysis results will be displayed here")
        self.result_info.setWordWrap(True)
        layout.addWidget(self.result_info)
        
        # 图片显示
        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_label.setMinimumHeight(400)
        self.image_label.setStyleSheet("border: 1px solid #ccc; background-color: #f9f9f9;")
        layout.addWidget(self.image_label)
        
        # 详细结果
        self.detailed_results = QPlainTextEdit()
        self.detailed_results.setReadOnly(True)
        self.detailed_results.setMaximumHeight(200)
        layout.addWidget(self.detailed_results)
    
    def setup_control_panel(self, main_layout):
        """设置控制面板"""
        control_layout = QHBoxLayout()
        
        # 进度条
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        # 运行按钮
        self.run_button = QPushButton("Run Analysis")
        self.run_button.clicked.connect(self.run_analysis)
        control_layout.addWidget(self.run_button)
        
        # 停止按钮
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_analysis)
        self.stop_button.setEnabled(False)
        control_layout.addWidget(self.stop_button)
        
        # 保存结果按钮
        self.save_button = QPushButton("Save Results")
        self.save_button.clicked.connect(self.save_results)
        self.save_button.setEnabled(False)
        control_layout.addWidget(self.save_button)
        
        main_layout.addLayout(control_layout)
    
    def browse_files(self):
        """浏览文件"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA files", "", 
            "FASTA files (*.fas *.fasta *.fa);;All files (*)"
        )
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
            
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
    
    def add_file_tag(self, file_path):
        """添加文件标签"""
        self.imported_files.append(file_path)
        
        tag_widget = QFrame()
        tag_widget.setFrameStyle(QFrame.Box)
        tag_widget.setStyleSheet("""
            QFrame {
                background-color: #e9ecef;
                border-radius: 15px;
                margin: 2px;
            }
        """)
        
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag_widget.setLayout(tag_layout)
        
        name_label = QLabel(os.path.basename(file_path))
        name_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(name_label)
        
        close_btn = QPushButton("×")
        close_btn.setFixedSize(20, 20)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #dc3545;
                color: white;
                border: none;
                border-radius: 10px;
                font-weight: bold;
                font-size: 12px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        close_btn.clicked.connect(lambda: self.remove_file_tag(file_path, tag_widget))
        tag_layout.addWidget(close_btn)
        
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append((file_path, tag_widget))
        self.file_tags_container.setVisible(True)
    
    def remove_file_tag(self, file_path, tag_widget):
        """移除文件标签"""
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        if not self.imported_files:
            self.file_path_edit.setText("")
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
        elif len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        else:
            self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
    
    def on_text_changed(self):
        """文本输入变化处理"""
        text = self.sequence_text.toPlainText().strip()
        if text:
            self.file_path_edit.setVisible(False)
            self.file_browse_btn.setVisible(False)
            self.file_tags_container.setVisible(False)
            self.clear_all_file_tags()
        else:
            self.file_path_edit.setVisible(True)
            self.file_browse_btn.setVisible(True)
            if self.imported_files:
                self.file_tags_container.setVisible(True)
    
    def clear_all_file_tags(self):
        """清除所有文件标签"""
        for file_path, tag_widget in self.file_tags:
            self.file_tags_layout.removeWidget(tag_widget)
            tag_widget.deleteLater()
        
        self.imported_files.clear()
        self.file_tags.clear()
        self.file_tags_container.setVisible(False)
    
    def browse_iqtree(self):
        """浏览IQ-TREE路径"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select IQ-TREE executable", "", 
            "Executable files (*.exe);;All files (*)"
        )
        if file_path:
            self.iqtree_path_edit.setText(file_path)
    
    def browse_evolver(self):
        """浏览Evolver路径"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Evolver executable", "", 
            "Executable files (*.exe);;All files (*)"
        )
        if file_path:
            self.evolver_path_edit.setText(file_path)
    
    def browse_output_dir(self):
        """浏览输出目录"""
        dir_path = QFileDialog.getExistingDirectory(self, "Select output directory")
        if dir_path:
            self.output_dir_edit.setText(dir_path)
    
    def run_analysis(self):
        """运行分析"""
        if self.is_running:
            return
        
        # 检查输入
        sequences = self.get_input_sequences()
        if not sequences:
            QMessageBox.warning(self, "Warning", "Please provide sequences!")
            return
        
        # 检查工具路径
        iqtree_path = self.iqtree_path_edit.text().strip()
        evolver_path = self.evolver_path_edit.text().strip()
        
        if not iqtree_path or not evolver_path:
            QMessageBox.warning(self, "Warning", "Please specify IQ-TREE and Evolver paths!")
            return
        
        # 更新UI状态
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        # 创建输出目录
        output_dir = self.output_dir_edit.text().strip()
        os.makedirs(output_dir, exist_ok=True)
        
        # 运行分析线程
        self.analysis_thread = SaturationAnalysisThread(
            sequences, output_dir, iqtree_path, evolver_path,
            num_simulations=int(self.simulations_spinbox.value()),
            confidence=float(self.confidence_spinbox.value()),
            test_type=str(self.test_type_combo.currentText())
        )
        self.analysis_thread.progress.connect(self.update_progress)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        # 将参数保存在实例上，供线程结束后展示用
        self._num_simulations = int(self.simulations_spinbox.value())
        self._confidence = float(self.confidence_spinbox.value())
        self._test_type = str(self.test_type_combo.currentText())
        self.analysis_thread.start()
    
    def get_input_sequences(self):
        """获取输入序列"""
        sequences = []
        
        # 优先使用导入的序列
        if self.imported_sequences:
            return self.imported_sequences
        
        # 使用文件输入
        if self.imported_files:
            for file_path in self.imported_files:
                try:
                    file_sequences = SequenceUtils.load_sequences(file_path)
                    sequences.extend(file_sequences)
                except Exception as e:
                    QMessageBox.warning(self, "Warning", f"Failed to load {file_path}: {e}")
                    return []
        
        # 使用文本输入
        elif self.sequence_text.toPlainText().strip():
            try:
                # 创建临时文件
                temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
                temp_file.write(self.sequence_text.toPlainText())
                temp_file.close()
                self.temp_files.append(temp_file.name)
                
                sequences = SequenceUtils.load_sequences(temp_file.name)
            except Exception as e:
                QMessageBox.warning(self, "Warning", f"Failed to parse sequences: {e}")
                return []
        
        return sequences
    
    def update_progress(self, message):
        """更新进度"""
        self.progress_bar.setFormat(message)
    
    def analysis_finished(self, result):
        """分析完成"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.save_button.setEnabled(True)
        
        self.analysis_result = result
        self.display_results(result)
        
        # 切换到结果标签页
        self.tab_widget.setCurrentIndex(1)
        
        QMessageBox.information(self, "Completed", "Saturation analysis completed!")
    
    def analysis_error(self, error_message):
        """分析错误"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.critical(self, "Error", f"Analysis failed: {error_message}")
    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'analysis_thread') and self.analysis_thread.isRunning():
            self.analysis_thread.terminate()
            self.analysis_thread.wait()
        
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "Analysis has been stopped.")
    
    def display_results(self, result):
        """显示结果"""
        # 更新结果信息
        info_text = f"""
        <h3>Analysis Results</h3>
        <p><b>Status:</b> {result['status']}</p>
        <p><b>Observed Iss:</b> {result['observed_iss']:.4f}</p>
        <p><b>Simulated Iss Mean:</b> {result['simulated_iss_mean']:.4f} ± {result['simulated_iss_std']:.4f}</p>
        <p><b>Test:</b> {result.get('test_type','two-sided')}, alpha={result.get('alpha', 0.05):.3f}, p={result.get('p_value', float('nan')):.4g}</p>
        <p><b>Critical values ({result.get('confidence', 0.95):.2f}):</b> [{result.get('lower_critical','N/A')}, {result.get('upper_critical','N/A')}]</p>
        <p><b>Number of sequences:</b> {result['n_sequences']}</p>
        <p><b>Sequence length:</b> {result['seq_length']}</p>
        <p><b>Model used:</b> {result['model_info'].get('model', 'Unknown')}</p>
        """
        self.result_info.setText(info_text)
        
        # 显示图片
        if os.path.exists(result['plot_file']):
            pixmap = QPixmap(result['plot_file'])
            scaled_pixmap = pixmap.scaled(self.image_label.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self.image_label.setPixmap(scaled_pixmap)
        
        # 显示详细结果
        detailed_text = f"""
Detailed Results:
================

Observed Iss: {result['observed_iss']:.6f}
Simulated Iss Statistics:
  Mean: {result['simulated_iss_mean']:.6f}
  Std:  {result['simulated_iss_std']:.6f}
  Min:  {min(result.get('simulated_iss', [0])):.6f}
  Max:  {max(result.get('simulated_iss', [0])):.6f}

Test: {result.get('test_type','two-sided')}
Alpha: {result.get('alpha', 0.05):.6f}
Empirical p-value: {result.get('p_value', float('nan')):.6f}
Critical values ({result.get('confidence', 0.95):.3f}): [{result.get('lower_critical','N/A')}, {result.get('upper_critical','N/A')}]

Model Information:
  Model: {result['model_info'].get('model', 'Unknown')}
  Alpha: {result['model_info'].get('alpha', 'Unknown')}
  Pinv:  {result['model_info'].get('pinv', 'Unknown')}

Tree Information:
  Tree used: {result.get('tree_used', 'Unknown')}
  {'Real tree from IQ-TREE' if result.get('tree_used') == 'real' else 'Random tree (IQ-TREE not available)'}

Interpretation:
{result['status']} - Decision based on empirical p-value and selected test type.

Files generated:
- Plot: {result['plot_file']}
- Summary: {result.get('summary_file', 'N/A')}
- IQ-TREE results: {result.get('iqtree_dir', 'N/A')}
- Evolver results: {result.get('evolver_dir', 'N/A')}
        """
        self.detailed_results.setPlainText(detailed_text)
    
    def save_results(self):
        """保存结果"""
        if not self.analysis_result:
            QMessageBox.warning(self, "Warning", "No results to save!")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Results", "saturation_results.json", 
            "JSON files (*.json);;All files (*)"
        )
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    json.dump(self.analysis_result, f, indent=2)
                QMessageBox.information(self, "Success", f"Results saved to {file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save results: {e}")
    
    def __del__(self):
        """清理临时文件"""
        for temp_file in self.temp_files:
            try:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
            except:
                pass


class SaturationAnalysisEntry:
    """饱和度分析插件入口"""
    
    def run(self):
        return SaturationAnalysisWidget()
