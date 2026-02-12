# macse_plugin.py
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

from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit, QToolButton, QDialog, 
                             QDoubleSpinBox, QRadioButton, QButtonGroup)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon
import tempfile
import os
import re
from typing import List, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
import subprocess


class MACSEThread(BaseProcessThread):
    """MACSE比对线程类"""
    
    def __init__(self, tool_path, input_files, parameters, mode, imported_files=None, suffix=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
        self.mode = mode  # 'align' or 'refine'
        if suffix:
            self.output_suffix = suffix
    
    def get_tool_name(self):
        """返回工具名称"""
        return "MACSE2"
        
    def execute_commands(self):
        """执行MACSE比对命令"""
        try:
            output_files = []
            html_files = []
            
            # 分别处理每个输入文件
            total_files = len(self.input_files)
            for i, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                    
                self.progress.emit(f"Processing file {i+1}/{total_files}...")
                self.console_output.emit(f"Processing file {i+1}/{total_files}: {os.path.basename(input_file)}", "info")
                
                # 创建输出文件路径
                if self.output_suffix:
                    base_name = os.path.splitext(input_file)[0] if '.' in os.path.basename(input_file) else input_file
                    output_nt = f"{base_name}{self.output_suffix}_NT.fasta"
                    output_aa = f"{base_name}{self.output_suffix}_AA.fasta"
                else:
                    output_nt = self.create_temp_file(suffix='_NT.fasta')
                    output_aa = self.create_temp_file(suffix='_AA.fasta')
                
                # 构建命令
                cmd = [
                    self.tool_path,
                    '-prog', 'alignSequences' if self.mode == 'align' else 'refineAlignment'
                ]
                
                # 根据模式添加输入参数
                if self.mode == 'align':
                    cmd.extend(['-seq', input_file])
                else:
                    cmd.extend(['-align', input_file])
                
                # 添加其他参数
                cmd.extend(self.parameters)
                
                # 添加输出参数
                cmd.extend(['-out_NT', output_nt, '-out_AA', output_aa])
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"MACSE execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 应用后处理逻辑
                self.console_output.emit(f"Applying post-processing for file {i+1}...", "info")
                self.apply_post_processing(output_nt, output_aa)
                
                # 从FASTA输出生成HTML报告
                self.console_output.emit(f"Generating HTML report for file {i+1}...", "info")
                # html_report = self.generate_html_report([output_nt, output_aa])
                # if html_report:
                #     html_files.append(html_report)
                
                output_files.extend([output_nt, output_aa])
            
            self.progress.emit("All alignments completed")
            self.finished.emit(output_files, html_files)
            
        except Exception as e:
            self.error.emit(f"MACSE alignment exception: {str(e)}")
    
    def apply_post_processing(self, nt_file, aa_file):
        """应用后处理逻辑"""
        try:
            # 读取序列
            nt_records = list(SeqIO.parse(nt_file, 'fasta'))
            aa_records = list(SeqIO.parse(aa_file, 'fasta'))
            
            # 1. Mask Exceptional Chars: 将!替换为?
            for record in nt_records:
                record.seq = Seq(str(record.seq).replace('!', '?'))
            for record in aa_records:
                record.seq = Seq(str(record.seq).replace('!', '?'))
            
            # 2. Remove Stop Codon: 处理终止密码子
            nt_records = self.remove_stop_codon(nt_records, aa_records)
            
            # 3. Remove Frameshift Gaps: 移除仅含!和-的列
            nt_records, aa_records = self.remove_frameshift_gaps(nt_records, aa_records)
            
            # 写回文件
            SeqIO.write(nt_records, nt_file, 'fasta')
            SeqIO.write(aa_records, aa_file, 'fasta')
            
        except Exception as e:
            self.console_output.emit(f"Post-processing warning: {str(e)}", "warning")
    
    def remove_stop_codon(self, nt_records, aa_records):
        """处理终止密码子"""
        processed_nt = []
        for nt_record, aa_record in zip(nt_records, aa_records):
            # 找到AA序列尾部第一个非-字符
            aa_seq_str = str(aa_record.seq)
            trimmed_aa = aa_seq_str.rstrip('-')
            if trimmed_aa and trimmed_aa[-1] == '*':
                # 找到最后一个完整密码子的位置（不包括终止密码子）
                codon_pos = len(trimmed_aa) - 1  # 终止密码子的位置
                nt_pos = codon_pos * 3  # 对应的核苷酸位置
                nt_seq_str = str(nt_record.seq)
                if nt_pos + 3 <= len(nt_seq_str):
                    # 遮蔽终止密码子为???
                    new_nt_seq = nt_seq_str[:nt_pos] + '???' + nt_seq_str[nt_pos+3:]
                    nt_record.seq = Seq(new_nt_seq)
            processed_nt.append(nt_record)
        return processed_nt
    
    def remove_frameshift_gaps(self, nt_records, aa_records):
        """移除仅包含!和-的列"""
        if not nt_records:
            return nt_records, aa_records
        
        # 转换为序列字符串列表
        nt_sequences = [str(record.seq) for record in nt_records]
        aa_sequences = [str(record.seq) for record in aa_records]
        
        # 检查所有序列长度是否一致
        seq_length = len(nt_sequences[0])
        if not all(len(seq) == seq_length for seq in nt_sequences):
            return nt_records, aa_records  # 长度不一致，跳过处理
        
        # 转换为列列表
        nt_columns = []
        for i in range(seq_length):
            column = ''.join(seq[i] for seq in nt_sequences)
            nt_columns.append(column)
        
        # 过滤列：保留不全是!和-的列
        filtered_columns = []
        for i, nt_col in enumerate(nt_columns):
            if not all(base in ['!', '-', '?'] for base in nt_col):
                filtered_columns.append(i)
        
        if not filtered_columns:
            return nt_records, aa_records  # 所有列都被过滤，保持原样
        
        # 重建序列
        new_nt_sequences = []
        new_aa_sequences = []
        for nt_seq, aa_seq in zip(nt_sequences, aa_sequences):
            new_nt_seq = ''.join(nt_seq[i] for i in filtered_columns)
            new_aa_seq = ''.join(aa_seq[i] for i in filtered_columns)
            new_nt_sequences.append(new_nt_seq)
            new_aa_sequences.append(new_aa_seq)
        
        # 更新记录
        for i, record in enumerate(nt_records):
            record.seq = Seq(new_nt_sequences[i])
        for i, record in enumerate(aa_records):
            record.seq = Seq(new_aa_sequences[i])
        
        return nt_records, aa_records


class MACSEPlugin(BasePlugin):
    """MACSE插件类"""
    import_alignment_signal = pyqtSignal(list)
    batch_import_alignment_signal = pyqtSignal(list)
    
    def __init__(self, import_from=None, import_data=None):
        """初始化MACSE插件"""
        super().__init__(import_from, import_data)
        self.batch_mode = False
        self.current_mode = 'align'  # 'align' or 'refine'
        if import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"] and import_data is not None:
            self.handle_import_data(import_data)
    
    def handle_import_data(self, import_data):
        """处理从YR-MPEA或Dataset Manager导入的数据"""
        if isinstance(import_data, list):
            # 判断是否为批量模式：[[SeqRecord, ...], [SeqRecord, ...]]
            if len(import_data) > 0 and isinstance(import_data[0], list):
                # 批量模式
                self.batch_mode = True
                self.temp_files = []
                self.imported_files = []
                
                for i, seq_list in enumerate(import_data):
                    temp_file = self.create_temp_file(suffix=f'_batch_{i}.fas')
                    with open(temp_file, 'w') as f:
                        for seq in seq_list:
                            f.write(f">{seq.id}\n{seq.seq}\n")
                    self.temp_files.append(temp_file)
                    self.imported_files.append(temp_file)
                    
                # 更新UI显示
                if hasattr(self, 'file_path_edit') and self.file_path_edit:
                    self.file_path_edit.setText(f"{len(import_data)} files selected")
                    self.file_path_edit.setEnabled(False)
                    
            else:
                # 单文件模式
                self.batch_mode = False
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
        
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "MACSE Aligner"
        self.tool_name = "MACSE2"
        self.citation = [
            """Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. <i>Molecular Biology and Evolution</i>. 2020;35(10):2582-2584. DOI: <a href="https://doi.org/10.1093/molbev/msy159">10.1093/molbev/msy159</a>""",
            """Ranwez V, Harispe S, Delsuc F, Douzery EJP. MACSE: Multiple Alignment of Coding SEquences accounting for frameshifts and stop codons. <i>PLoS One</i>. 2011;6(9):e22594. DOI: <a href="https://doi.org/10.1371/journal.pone.0022594">10.1371/journal.pone.0022594</a>"""
        ]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"FASTA": ".fasta"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    
    def setup_input_tab(self):
        """设置输入标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # 模式选择
        mode_layout = QHBoxLayout()
        mode_label = QLabel("Mode:")
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Align Sequences (alignSequences)",
            "Refine Alignment (refineAlignment)"
        ])
        self.mode_combo.currentIndexChanged.connect(self.on_mode_changed)
        mode_layout.addWidget(mode_label)
        mode_layout.addWidget(self.mode_combo)
        mode_layout.addStretch()
        layout.addLayout(mode_layout)
        
        # 输入组
        input_group = QGroupBox("Input / Output")
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
        self.sequence_text.setPlaceholderText("Or paste FASTA format sequence...")
        self.sequence_text.setMaximumHeight(200)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        input_layout.addRow("Sequence text:", self.sequence_text)
        
        # 回传选项（仅在YR-MPEA导入时显示）
        self.return_option_container = QWidget()
        return_option_layout = QHBoxLayout()
        self.return_option_container.setLayout(return_option_layout)
        return_option_label = QLabel("Return:")
        self.return_option_combo = QComboBox()
        self.return_option_combo.addItems([
            "Codon Sequences",
            "Translated AA Sequences"
        ])
        return_option_layout.addWidget(return_option_label)
        return_option_layout.addWidget(self.return_option_combo)
        return_option_layout.addStretch()
        input_layout.addRow("", self.return_option_container)
        self.return_option_container.setVisible(self.import_from == "YR_MPEA")
        
        # 处理导入的数据
        if self.import_file:
            self.file_path_edit.setText(self.import_file)
            input_group.setVisible(False)
        
        # 基本参数组
        basic_params_group = QGroupBox("Basic Parameters")
        basic_params_layout = QFormLayout()
        basic_params_group.setLayout(basic_params_layout)
        layout.addWidget(basic_params_group)
        
        # 遗传密码表
        self.genetic_code_combo = QComboBox()
        genetic_codes = [
            "1 Standard",
            "2 Vertebrate Mitochondrial", 
            "3 Yeast Mitochondrial",
            "4 Mold, Protozoan, Coelenterate Mitochondrial",
            "5 Invertebrate Mitochondrial",
            "6 Ciliate, Dasycladacean & Hexamita Nuclear",
            "9 Echinoderm & Flatworm Mitochondrial",
            "10 Euplotid Nuclear",
            "11 Bacterial Archaeal and Plant Plastid",
            "12 Alternative Yeast Nuclear",
            "13 Ascidian Mitochondrial",
            "14 Alternative Flatworm Mitochondrial",
            "15 Blepharisma Nuclear",
            "16 Chlorophycean Mitochondrial",
            "21 Trematode Mitochondrial",
            "22 Scenedesmus obliquus mitochondrial",
            "23 Thraustochytrium Mitochondrial",
            "24 Rhabdopleuridae Mitochondrial",
            "25 Candidate Division SR1 & Gracilibacteria",
            "26 Pachysolen tannophilus Nuclear",
            "27 Karyorelict Nuclear",
            "28 Condylostoma Nuclear",
            "29 Mesodinium Nuclear",
            "30 Peritrich Nuclear",
            "31 Blastocrithidia Nuclear",
            "32 Seleno Protein",
            "33 Cephalodiscidae Mitochondrial UAA-Tyr"
        ]
        self.genetic_code_combo.addItems(genetic_codes)
        basic_params_layout.addRow("Genetic Code:", self.genetic_code_combo)
        
        # 氨基酸字母表
        self.aa_alphabet_combo = QComboBox()
        aa_alphabets = [
            "SE_B_14",
            "SE_B_10", 
            "SE_V_10",
            "Li_A_10",
            "Li_B_10",
            "Solis_D_10",
            "Solis_G_10",
            "Murphy_10",
            "SE_B_8",
            "SE_B_6",
            "Dayhoff_6"
        ]
        self.aa_alphabet_combo.addItems(aa_alphabets)
        self.aa_alphabet_combo.setCurrentText("SE_B_8")
        basic_params_layout.addRow("AA Alphabet:", self.aa_alphabet_combo)
        
        # 打分矩阵
        self.score_matrix_combo = QComboBox()
        score_matrices = ["BLOSUM62", "VTML200_BIS", "VTML240"]
        self.score_matrix_combo.addItems(score_matrices)
        self.score_matrix_combo.setCurrentText("BLOSUM62")
        basic_params_layout.addRow("Score Matrix:", self.score_matrix_combo)
        
        # 优化级别
        self.optim_spinbox = QSpinBox()
        self.optim_spinbox.setRange(0, 2)
        self.optim_spinbox.setValue(2)
        basic_params_layout.addRow("Optimization Level:", self.optim_spinbox)
        
        # 后处理选项
        postprocess_group = QGroupBox("Post-processing Options (default enabled)")
        postprocess_layout = QVBoxLayout()
        postprocess_group.setLayout(postprocess_layout)
        layout.addWidget(postprocess_group)
        
        self.mask_chars_checkbox = QCheckBox("Mask Exceptional Chars (! → ?)")
        self.mask_chars_checkbox.setChecked(True)
        postprocess_layout.addWidget(self.mask_chars_checkbox)
        
        self.remove_frameshift_checkbox = QCheckBox("Remove Frameshift Gaps (columns with only ! and -)")
        self.remove_frameshift_checkbox.setChecked(True)
        postprocess_layout.addWidget(self.remove_frameshift_checkbox)
        
        self.remove_stop_codon_checkbox = QCheckBox("Remove Stop Codon (* → ??? in NT)")
        self.remove_stop_codon_checkbox.setChecked(True)
        postprocess_layout.addWidget(self.remove_stop_codon_checkbox)
        
        # 高级参数按钮
        advanced_btn = QPushButton("Advanced Parameters...")
        advanced_btn.clicked.connect(self.show_advanced_dialog)
        layout.addWidget(advanced_btn)
        
        # Output settings
        output_settings_group = QWidget()
        output_settings_group.setLayout(QVBoxLayout())
        output_settings_group.layout().setContentsMargins(0, 0, 0, 0)
        output_settings_group.setContentsMargins(0, 0, 0, 0)
        output_radio_group = QWidget()
        output_radio_group.setLayout(QHBoxLayout())
        self.save_to_cwd = QRadioButton("Current Directory")
        self.save_to_cwd.setChecked(True)
        self.save_to_tmp = QRadioButton("Temporary File")
        output_radio_group.layout().addWidget(self.save_to_cwd)
        output_radio_group.layout().addWidget(self.save_to_tmp)
        output_settings_group.layout().addWidget(output_radio_group)

        self.save_to_cwd.toggled.connect(self.on_output_radio_changed)
        self.save_to_tmp.toggled.connect(self.on_output_radio_changed)

        self.output_suffix_edit = QLineEdit('_macse')
        self.output_suffix_edit.setPlaceholderText("Output suffix")
        output_settings_group.layout().addWidget(self.output_suffix_edit)
        input_layout.addRow("Save to:", output_settings_group)
        
        layout.addStretch()
    
    def on_mode_changed(self, index):
        """处理模式切换"""
        self.current_mode = 'align' if index == 0 else 'refine'
        # 更新文件输入提示
        if self.current_mode == 'align':
            self.file_path_edit.setPlaceholderText("Select unaligned FASTA files...")
        else:
            self.file_path_edit.setPlaceholderText("Select aligned FASTA files...")
    
    def on_output_radio_changed(self):
        """处理输出选项的切换"""
        if self.save_to_cwd.isChecked():
            self.output_suffix_edit.setEnabled(True)
            self.output_suffix_edit.setVisible(True)
        elif self.save_to_tmp.isChecked():
            self.output_suffix_edit.setEnabled(False)
            self.output_suffix_edit.setVisible(False)
    
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)
    
    def show_advanced_dialog(self):
        """显示高级参数对话框"""
        dialog = MACSEAdvancedDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            # 存储高级参数
            self.advanced_params = dialog.get_parameters()
    
    def get_parameters(self):
        """获取MACSE命令参数"""
        params = []
        
        # 遗传密码表
        genetic_code_text = self.genetic_code_combo.currentText()
        genetic_code_num = genetic_code_text.split()[0]
        params.extend(['-gc_def', genetic_code_num])
        
        # 氨基酸字母表
        aa_alphabet = self.aa_alphabet_combo.currentText()
        params.extend(['-alphabet_AA', aa_alphabet])
        
        # 打分矩阵
        score_matrix = self.score_matrix_combo.currentText()
        params.extend(['-score_matrix', score_matrix])
        
        # 优化级别
        optim_level = self.optim_spinbox.value()
        params.extend(['-optim', str(optim_level)])
        
        # 后处理选项参数（这些在Python中处理，不在命令行中）
        # 但可以传递一些相关参数
        if self.current_mode == 'align':
            # alignSequences特定参数
            max_refine_iter = -1  # 默认自动
            params.extend(['-max_refine_iter', str(max_refine_iter)])
        else:
            # refineAlignment特定参数
            local_realign_init = 1.0
            local_realign_dec = 1.0
            params.extend([
                '-local_realign_init', str(local_realign_init),
                '-local_realign_dec', str(local_realign_dec)
            ])
        
        # 添加高级参数
        if hasattr(self, 'advanced_params'):
            params.extend(self.advanced_params)
        
        return params
    
    
    def on_analysis_finished(self, output_files, html_files):
        """分析完成回调"""
        # 更新UI状态
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # 保存输出文件和报告
        self.reports = html_files
        self.alignment_output_files = output_files  # 保存output_files为实例属性
        
        # 显示导入按钮（支持多种来源）
        if self.import_from in ["YR_MPEA", "seq_viewer", "DATASET_MANAGER"]:
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        
        # 设置报告文件
        if self.reports:
            self.update_report_combo()
            self.show_current_report()
        
        self.add_console_message("MACSE analysis completed successfully!", "info")
        QMessageBox.information(self, "Success", "MACSE analysis completed successfully!")
    
    def on_analysis_error(self, error_msg):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Error", f"MACSE analysis failed: {error_msg}")
        self.tab_widget.setCurrentIndex(2)  # Switch to console tab

    def import_to_platform(self):
        """导入到当前平台"""
        # if not self.reports:
        #     return
        
        # 根据回传选项决定返回哪种序列
        return_codon = self.return_option_combo.currentText() == "Codon Sequences"
        
        try:
            # 解析输出文件
            all_sequences = []
            # 使用实例属性而不是线程属性
            nt_files = [f for f in self.alignment_output_files if '_NT.fasta' in f]
            aa_files = [f for f in self.alignment_output_files if '_AA.fasta' in f]
            
            if self.batch_mode:
                # 批量模式
                batch_sequences = []
                for nt_file in nt_files:
                    if return_codon:
                        sequences = list(SeqIO.parse(nt_file, 'fasta'))
                    else:
                        # 找到对应的AA文件
                        base_name = nt_file.replace('_NT.fasta', '')
                        aa_file = base_name + '_AA.fasta'
                        if os.path.exists(aa_file):
                            sequences = list(SeqIO.parse(aa_file, 'fasta'))
                        else:
                            sequences = list(SeqIO.parse(nt_file, 'fasta'))
                    batch_sequences.append(sequences)
                self.batch_import_alignment_signal.emit(batch_sequences)
            else:
                # 单文件模式
                if return_codon and nt_files:
                    sequences = list(SeqIO.parse(nt_files[0], 'fasta'))
                elif not return_codon and aa_files:
                    sequences = list(SeqIO.parse(aa_files[0], 'fasta'))
                else:
                    sequences = list(SeqIO.parse(nt_files[0] if nt_files else aa_files[0], 'fasta'))
                self.import_alignment_signal.emit(sequences)
            
            self.close()
            
        except Exception as e:
            QMessageBox.critical(self, "Import Error", f"Failed to import sequences: {str(e)}")

    def browse_files(self):
        """浏览选择文件"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA files", "", "FASTA files (*.fas *.fasta *.fa *.fna);;All files (*)"
        )
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            # 隐藏文本输入框
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)

    def add_file_tag(self, file_path):
        """添加文件标签"""
        self.imported_files.append(file_path)
        
        # Create file tag widget
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
        
        # Get display name (handle duplicate names)
        display_name = self.get_display_name(file_path)
        
        # File name label
        name_label = QLabel(display_name)
        name_label.setStyleSheet("color: #495057;")
        tag_layout.addWidget(name_label)
        
        # Close button
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
        
        # Add to layout
        self.file_tags_layout.addWidget(tag_widget)
        self.file_tags.append((file_path, tag_widget))
        
        # Show container if it was hidden
        self.file_tags_container.setVisible(True)

    def remove_file_tag(self, file_path, tag_widget):
        """Remove a file tag"""
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        
        # Remove from tags list
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        
        # Remove widget
        self.file_tags_layout.removeWidget(tag_widget)
        tag_widget.deleteLater()
        
        # Update file path display
        if not self.imported_files:
            self.file_path_edit.setText("")
            self.file_tags_container.setVisible(False)
            # Show text input when no files
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)
        elif len(self.imported_files) == 1:
            self.file_path_edit.setText(self.imported_files[0])
        else:
            self.file_path_edit.setText(f"{len(self.imported_files)} files selected")

    def get_display_name(self, file_path):
        """Get display name for file, handling duplicates"""
        filename = os.path.basename(file_path)
        
        # Check for duplicates
        duplicate_count = 0
        for existing_path in self.imported_files:
            if os.path.basename(existing_path) == filename and existing_path != file_path:
                duplicate_count += 1
        
        if duplicate_count > 0:
            # Use last directory name for duplicates
            parent_dir = os.path.basename(os.path.dirname(file_path))
            return f"{filename} ({parent_dir})"
        else:
            return filename

    def on_text_changed(self):
        """Handle text input changes"""
        text = self.sequence_text.toPlainText().strip()
        if text:
            # Hide file input when text is present
            self.file_path_edit.setVisible(False)
            self.file_browse_btn.setVisible(False)
            self.file_tags_container.setVisible(False)
            # Clear imported files
            self.clear_all_file_tags()
        else:
            # Show file input when text is empty
            self.file_path_edit.setVisible(True)
            self.file_browse_btn.setVisible(True)
            if self.imported_files:
                self.file_tags_container.setVisible(True)

    def clear_all_file_tags(self):
        """Clear all file tags"""
        # Clear imported files list
        self.imported_files.clear()
        
        # Clear file tags list and widgets
        for file_path, tag_widget in self.file_tags:
            self.file_tags_layout.removeWidget(tag_widget)
            tag_widget.deleteLater()
        self.file_tags.clear()
        
        # Hide container
        self.file_tags_container.setVisible(False)

    def prepare_input_files(self):
        """Prepare input files from multiple sources"""
        try:
            input_files = []
            
            # If there are imported files, use them individually
            if self.imported_files:
                for file_path in self.imported_files:
                    if os.path.exists(file_path):
                        input_files.append(file_path)
                    else:
                        QMessageBox.warning(self, "Warning", f"File does not exist: {file_path}")
                return input_files
            elif self.import_file:
                return [self.import_file]
            
            # Otherwise, use text input to create a temporary file
            sequence_text = self.sequence_text.toPlainText().strip()
            if not sequence_text and not self.import_file:
                QMessageBox.warning(self, "Warning", "Please input sequence text!")
                return None
            
            if self.save_to_tmp.isChecked():
                # Create temporary file
                temp_file = self.create_temp_file(suffix='.fas')
                with open(temp_file, 'w') as f:
                    f.write(sequence_text)
                return [temp_file]
            
            elif self.save_to_cwd.isChecked():
                # directly use input paths
                temp_file = self.import_file
                return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None

    def run_analysis(self):
        """运行MACSE分析"""
        if self.is_running:
            return
            
        # 检查输入
        if not self.imported_files and not self.sequence_text.toPlainText().strip() and not self.import_file:
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
            return
            
        # 获取工具路径 - 使用config方法而不是硬编码
        if not self.config():
            QMessageBox.critical(self, "Error", f"MACSE executable configuration not found in config.json")
            return
        tool_path = self.tool_path
        
        # 添加控制台消息
        self.add_console_message("Starting MACSE alignment...", "info")
        
        # 准备输入文件
        input_files = self.prepare_input_files()
        if not input_files:
            return
            
        # 更新UI状态
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 未知进度
        
        # 如果需要添加输出文件后缀，则设置suffix属性
        suffix = None
        if self.save_to_cwd.isChecked():
            suffix = self.output_suffix_edit.text()
        
        # 在单独的线程中运行比对
        self.thread = MACSEThread(
            tool_path, 
            input_files, 
            self.get_parameters(), 
            self.current_mode,
            self.imported_files,
            suffix=suffix
        )
        self.thread.progress.connect(self.progress_bar.setFormat)
        self.thread.finished.connect(self.on_analysis_finished)
        self.thread.error.connect(self.on_analysis_error)
        self.thread.console_output.connect(self.add_console_message)
        self.thread.start()

    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'thread') and self.thread.isRunning():
            self.thread.terminate()
            self.thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "MACSE analysis has been aborted.")


class MACSEAdvancedDialog(QDialog):
    """MACSE高级参数对话框"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MACSE Advanced Parameters")
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout()
        
        # 参数输入
        self.param_edit = QTextEdit()
        self.param_edit.setPlaceholderText("Enter additional command line parameters...\nExample: -max_refine_iter 5 -local_realign_init 0.5")
        layout.addWidget(self.param_edit)
        
        # 按钮
        button_layout = QHBoxLayout()
        ok_btn = QPushButton("OK")
        cancel_btn = QPushButton("Cancel")
        ok_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(ok_btn)
        button_layout.addWidget(cancel_btn)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def get_parameters(self):
        """获取参数列表"""
        text = self.param_edit.toPlainText().strip()
        if not text:
            return []
        # 简单分割参数（需要更复杂的解析）
        return text.split()

class MACSEPluginEntry:
    def run(self, import_from=None, import_data=None):
        return MACSEPlugin(import_from, import_data)