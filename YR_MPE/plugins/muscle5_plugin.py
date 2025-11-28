from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, 
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSpinBox, QCheckBox, QLabel, QComboBox, QScrollArea,
                             QWidget, QFrame, QTextEdit)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon
import tempfile
import os
from Bio import SeqIO
import subprocess
from typing import List, Optional

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread


class Muscle5Thread(BaseProcessThread):
    """MUSCLE5比对线程类"""
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__(tool_path, input_files, parameters, imported_files)
    
    def get_tool_name(self):
        """返回工具名称"""
        return "MUSCLE5"
        
    def execute_command(self, cmd: List[str], cwd: Optional[str] = None) -> subprocess.CompletedProcess:
        """
        执行单个命令并实时更新进度
        
        Args:
            cmd (List[str]): 命令和参数列表
            cwd (Optional[str]): 工作目录
            
        Returns:
            subprocess.CompletedProcess: 命令执行结果
        """
        self.console_output.emit(" ".join(cmd), "command")
        
        # 执行命令
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
            bufsize=1,
            universal_newlines=True
        )
        
        # 实时读取输出
        stdout_lines = []
        stderr_lines = []
        
        while True:
            # 检查进程是否还在运行
            if process.poll() is not None:
                break
                
            # 非阻塞读取输出
            try:
                # 读取标准输出
                output = process.stdout.readline()
                if output:
                    stdout_lines.append(output)
                    self.console_output.emit(output.strip(), "output")
                    # 发出进度更新信号
                    self.progress.emit(f"Processing... {len(stdout_lines)} lines")
                
                # 读取错误输出
                error = process.stderr.readline()
                if error:
                    stderr_lines.append(error)
                    self.console_output.emit(error.strip(), "error")
                    # 发出进度更新信号
                    self.progress.emit(f"Processing... {len(stdout_lines)} lines")
                    
            except Exception as e:
                self.console_output.emit(f"Error reading output: {e}", "error")
        
        # 读取剩余输出
        try:
            remaining_stdout, remaining_stderr = process.communicate(timeout=1)
            if remaining_stdout:
                for line in remaining_stdout.split('\n'):
                    if line.strip():
                        stdout_lines.append(line)
                        self.console_output.emit(line.strip(), "output")
                        
            if remaining_stderr:
                for line in remaining_stderr.split('\n'):
                    if line.strip():
                        stderr_lines.append(line)
                        self.console_output.emit(line.strip(), "error")
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            if stdout:
                for line in stdout.split('\n'):
                    if line.strip():
                        stdout_lines.append(line)
                        self.console_output.emit(line.strip(), "output")
            if stderr:
                for line in stderr.split('\n'):
                    if line.strip():
                        stderr_lines.append(line)
                        self.console_output.emit(line.strip(), "error")
        
        # 创建CompletedProcess对象
        result = subprocess.CompletedProcess(
            args=cmd,
            returncode=process.returncode,
            stdout='\n'.join(stdout_lines),
            stderr='\n'.join(stderr_lines)
        )
        
        return result

    def execute_commands(self):
        """执行MUSCLE5比对命令"""
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
                
                # 创建输出文件
                output_file = self.create_temp_file(suffix='.fas')
                # html_file = self.create_temp_file(suffix='.html')
                
                # 构建命令
                cmd = [
                    self.tool_path,
                    *self.parameters,
                    input_file,
                    "-output", output_file
                ]
                
                # 执行命令
                result = self.execute_command(cmd)
                
                if result.returncode != 0:
                    self.error.emit(f"MUSCLE5 execution failed for file {i+1}: {result.stderr}")
                    return
                
                # 从FASTA输出生成HTML报告
                self.console_output.emit(f"Generating HTML report for file {i+1}...", "info")
                # self.generate_html_from_fasta(output_file, html_file)
                
                output_files.append(output_file)
                # html_files.append(html_file)
            
            self.progress.emit("All alignments completed")
            self.finished.emit(output_files, [])
            
        except Exception as e:
            self.error.emit(f"Alignment exception: {str(e)}")
    
    def generate_html_from_fasta(self, fasta_file, html_file):
        """Generate HTML visualization from FASTA alignment with interleaved display and conservation highlighting"""
        try:
            with open(fasta_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Parse FASTA sequences
            sequences = []
            current_seq = ""
            current_header = ""
            
            for line in content.split('\n'):
                line = line.strip()
                if line.startswith('>'):
                    if current_seq and current_header:
                        sequences.append((current_header, current_seq))
                    current_header = line[1:]  # Remove '>'
                    current_seq = ""
                elif line:
                    current_seq += line
            
            if current_seq and current_header:
                sequences.append((current_header, current_seq))
            
            if not sequences:
                return
            
            # Calculate conservation for each position
            conservation = self.calculate_conservation(sequences)
            
            # Pre-calculate statistics
            num_sequences = len(sequences)
            alignment_length = len(sequences[0][1]) if sequences else 0
            conserved_count = sum(1 for c in conservation if (c >= 0.8 and c != 1))
            conserved_percent = conserved_count / len(conservation) * 100 if conservation else 0
            invariable_count = sum(1 for c in conservation if c == 1)
            invariable_percent = invariable_count / len(conservation) * 100 if conservation else 0
            
            # Use list for efficient string building
            html_parts = []
            
            # HTML header and CSS
            html_parts.append(f"""<!DOCTYPE html>
<html>
<head>
    <title>MUSCLE5 Alignment Result</title>
    <style>
        body {{ 
            font-family: 'Courier New', monospace; 
            margin: 20px; 
            background-color: #f8f9fa;
        }}
        .container {{ 
            max-width: 1200px; 
            margin: 0 auto; 
            background-color: white; 
            padding: 20px; 
            border-radius: 8px; 
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .stats {{ 
            background-color: #e9ecef; 
            padding: 15px; 
            margin-bottom: 20px; 
            border-radius: 5px; 
            border-left: 4px solid #007bff;
        }}
        .alignment-container {{ 
            overflow-x: auto; 
            border: 1px solid #dee2e6; 
            border-radius: 5px; 
            background-color: white;
        }}
        .sequence-block {{ 
            margin-bottom: 20px; 
            border: 1px solid #dee2e6; 
            border-radius: 5px; 
            background-color: white;
        }}
        .sequence-row {{ 
            display: flex; 
            border-bottom: 1px solid #e9ecef; 
            padding: 3px 0;
        }}
        .sequence-row:hover {{ 
            background-color: #f8f9fa; 
        }}
        .sequence-header {{ 
            min-width: 200px; 
            padding: 3px 10px; 
            font-weight: bold; 
            color: #495057; 
            background-color: #f8f9fa;
            border-right: 2px solid #dee2e6;
            display: flex;
            align-items: center;
            font-size: 12px;
        }}
        .sequence-content {{ 
            display: flex; 
            flex-wrap: wrap; 
            font-size: 12px; 
            letter-spacing: 1px; 
            line-height: 1.1;
            font-family: 'Courier New', monospace;
        }}
        .nucleotide {{ 
            display: inline-block; 
            width: 20px; 
            text-align: center; 
            margin: 1px;
            border-radius: 2px;
        }}
        .invariable {{
            background-color: #8bb0cd;
            color: #005178;
            font-weight: bold;
        }}
        .conserved {{ 
            background-color: #d4edda; 
            color: #155724; 
            font-weight: bold;
        }}
        .semi-conserved {{ 
            background-color: #fff3cd; 
            color: #856404;
        }}
        .variable {{ 
            background-color: #f8d7da; 
            color: #721c24;
        }}
        .gap {{ 
            background-color: #e2e3e5; 
            color: #6c757d;
        }}
        .position-numbers {{ 
            display: flex; 
            padding: 5px 0; 
            background-color: #f8f9fa; 
            border-bottom: 2px solid #dee2e6;
        }}
        .pos-number {{ 
            width: 20px; 
            text-align: center; 
            font-size: 10px; 
            color: #6c757d; 
            margin: 1px;
        }}
        .legend {{ 
            margin-top: 20px; 
            padding: 15px; 
            background-color: #f8f9fa; 
            border-radius: 5px;
        }}
        .legend-item {{ 
            display: inline-block; 
            margin-right: 20px; 
            margin-bottom: 5px;
        }}
        .legend-color {{ 
            display: inline-block; 
            width: 20px; 
            height: 20px; 
            margin-right: 5px; 
            border-radius: 2px; 
            vertical-align: middle;
        }}
        h1 {{ color: #343a40; margin-bottom: 20px; }}
        h3 {{ color: #495057; margin-bottom: 10px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 MUSCLE5 Alignment Result</h1>
        
        <div class="stats">
            <h3>📊 Alignment Statistics</h3>
            <p><strong>Number of sequences:</strong> {num_sequences}</p>
            <p><strong>Alignment length:</strong> {alignment_length} bp</p>
            <p><strong>Conserved positions:</strong> {conserved_count} ({conserved_percent:.1f}%)</p>
            <p><strong>Invariable positions:</strong> {invariable_count} ({invariable_percent:.1f}%)</p>
        </div>

        <div class="legend">
            <h3>🎨 Conservation Legend</h3>
            <div class="legend-item">
                <span class="legend-color invariable"></span>
                <span>Invariable (100%)</span>
            </div>
            <div class="legend-item">
                <span class="legend-color conserved"></span>
                <span>Highly Conserved (≥80%)</span>
            </div>
            <div class="legend-item">
                <span class="legend-color semi-conserved"></span>
                <span>Moderately Conserved (50-79%)</span>
            </div>
            <div class="legend-item">
                <span class="legend-color variable"></span>
                <span>Variable (<50%)</span>
            </div>
            <div class="legend-item">
                <span class="legend-color gap"></span>
                <span>Gap</span>
            </div>
        </div>
        
        <div class="alignment-container">""")
            
            # Calculate number of blocks (each block shows 40 characters)
            chars_per_block = 40
            num_blocks = (alignment_length + chars_per_block - 1) // chars_per_block
            
            # Generate interleaved blocks efficiently
            for block_idx in range(num_blocks):
                start_pos = block_idx * chars_per_block
                end_pos = min(start_pos + chars_per_block, alignment_length)
                
                # Block header
                html_parts.append(f"""
            <div class="sequence-block">
                <div class="position-numbers">
                    <div style="min-width: 200px; padding: 5px 10px; font-weight: bold; color: #6c757d;">
                        Block {block_idx + 1} (pos {start_pos + 1}-{end_pos})
                    </div>
                    <div class="sequence-content">""")
                
                # Position numbers
                pos_numbers = []
                for i in range(start_pos, end_pos, 10):
                    pos_numbers.append(f'<span class="pos-number">{i+1}</span>')
                html_parts.append(''.join(pos_numbers))
                
                html_parts.append("""
                    </div>
                </div>""")
                
                # Sequences for this block
                for i, (header, sequence) in enumerate(sequences):
                    html_parts.append(f"""
                <div class="sequence-row">
                    <div class="sequence-header">
                        <span>{header[:35]}{'...' if len(header) > 35 else ''}</span>
                    </div>
                    <div class="sequence-content">""")
                    
                    # Nucleotides for this block
                    nucleotides = []
                    for j in range(start_pos, end_pos):
                        if j < len(sequence):
                            nucleotide = sequence[j]
                            cons_level = conservation[j] if j < len(conservation) else 0
                            
                            if nucleotide in ['-', '.']:
                                css_class = "gap"
                            elif cons_level == 1:
                                css_class = "invariable"
                            elif cons_level >= 0.8:
                                css_class = "conserved"
                            elif cons_level >= 0.5:
                                css_class = "semi-conserved"
                            else:
                                css_class = "variable"
                            
                            nucleotides.append(f'<span class="nucleotide {css_class}">{nucleotide}</span>')
                        else:
                            nucleotides.append('<span class="nucleotide gap">-</span>')
                    
                    html_parts.append(''.join(nucleotides))
                    html_parts.append("""
                    </div>
                </div>""")
                
                html_parts.append("""
            </div>""")
            
            # HTML footer
            html_parts.append("""
        </div>
    </div>
</body>
</html>""")
            
            # Join all parts efficiently
            html_content = ''.join(html_parts)
            
            # Write HTML file
            with open(html_file, 'w', encoding='utf-8') as f:
                f.write(html_content)
                
        except Exception as e:
            print(f"Error generating HTML: {e}")
            # Create a simple error HTML
            with open(html_file, 'w', encoding='utf-8') as f:
                f.write(f"""<!DOCTYPE html>
<html>
<head><title>Error</title></head>
<body>
<h2>Error generating HTML visualization</h2>
<p>Error: {e}</p>
<p>Alignment completed successfully, but HTML visualization failed.</p>
</body>
</html>""")
    
    def calculate_conservation(self, sequences):
        """Calculate conservation score for each position in the alignment"""
        if not sequences:
            return []
        
        alignment_length = len(sequences[0][1])
        conservation_scores = []
        
        for pos in range(alignment_length):
            # Get nucleotides at this position
            nucleotides = []
            for _, seq in sequences:
                if pos < len(seq):
                    nucleotide = seq[pos].upper()
                    if nucleotide not in ['-', '.', 'N']:  # Ignore gaps and N
                        nucleotides.append(nucleotide)
            
            if not nucleotides:
                conservation_scores.append(0.0)
                continue
            
            # Calculate conservation
            total_count = len(nucleotides)
            most_common_count = max(nucleotides.count(nuc) for nuc in set(nucleotides))
            conservation_score = most_common_count / total_count
            
            conservation_scores.append(conservation_score)
        
        return conservation_scores


class Muscle5Plugin(BasePlugin):
    """MUSCLE5比对插件"""
    
    # 定义信号
    import_alignment_signal = pyqtSignal(list)  # 导入比对结果信号
    
    def __init__(self, import_from=None, import_data=None):
        """初始化MUSCLE5插件"""
        super().__init__(import_from, import_data)
        # 初始化插件信息
        self.init_plugin_info()
        
        # 特别处理YR-MPEA导入的数据
        if import_from == "YR_MPEA" and import_data is not None:
            self.handle_import_data(import_data)
    
    def init_plugin_info(self):
        """初始化插件信息"""
        self.plugin_name = "MUSCLE5 Aligner"
        self.tool_name = "Muscle5"
        self.citation = ["""Edgar, R.C. Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. <i>Nat Commun</i> <b>13</b>, 6968 (2022). DOI: <a href="https://doi.org/10.1038/s41467-022-34630-w">10.1038/s41467-022-34630-w</a>"""]
        self.input_types = {"FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"FASTA": ".fas"}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
            
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

    def setup_input_tab(self):
        """设置输入参数标签页"""
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # Input by file or text
        input_group = QGroupBox("Input")
        input_layout = QFormLayout()
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        # File input
        file_layout = QHBoxLayout()
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select FASTA files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.file_browse_btn)
        input_layout.addRow("File input:", file_layout)
        
        # File tags container
        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        input_layout.addRow("", self.file_tags_container)
        
        # Text input
        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste FASTA format sequence...")
        self.sequence_text.setMaximumHeight(200)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        input_layout.addRow("Sequence text:", self.sequence_text)

        # Check if file is imported
        if self.import_file:
            # File Imported
            self.file_imported_button = QPushButton(f"File Already Imported")
            self.file_imported_button.setEnabled(False)
            layout.addWidget(self.file_imported_button)

            # set the file_path_edit to the importfile
            self.file_path_edit.setText(self.import_file)
            # invisible the input_group
            input_group.setVisible(False)
        elif hasattr(self, 'imported_files') and self.imported_files:
            # 显示导入的文件
            for file_path in self.imported_files:
                self.add_file_tag(file_path)
            
            # 更新文件路径显示
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            # 隐藏文本输入框
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
        
        # Parameters setting
        params_group = QGroupBox("Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # Alignment mode
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Standard alignment (-align)", 
            "Large dataset (-super5)"
        ])
        params_layout.addRow("Alignment mode:", self.mode_combo)
        
        layout.addStretch()
        
        # Initialize variables
        if not hasattr(self, 'imported_files'):
            self.imported_files = []  # List of imported file paths
        if not hasattr(self, 'file_tags'):
            self.file_tags = []  # List of file tag widgets

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
        for file_path, tag_widget in self.file_tags:
            self.file_tags_layout.removeWidget(tag_widget)
            tag_widget.deleteLater()
        
        self.imported_files.clear()
        self.file_tags.clear()
        self.file_tags_container.setVisible(False)
    
    def setup_control_panel(self):
        """设置控制面板"""
        super().setup_control_panel()
        
        # 添加导入到平台按钮
        self.import_to_platform_btn = QPushButton("Import to Current Platform")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)  # 初始隐藏
        
        # 确保布局存在并添加按钮
        self.control_layout.addWidget(self.import_to_platform_btn)

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
                
            # Create temporary file
            temp_file = self.create_temp_file(suffix='.fas')
            with open(temp_file, 'w') as f:
                f.write(sequence_text)
            return [temp_file]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None

    def get_parameters(self):
        """获取命令行参数"""
        parameters = []
        
        # 添加比对模式参数
        if self.mode_combo.currentText().startswith("Large dataset"):
            parameters.append("-super5")
        else:
            parameters.append("-align")

        # Go fuck yourself, Qwen! Muscle 5 DO NOT HAVE THESE PARAMETERS!   
        # # 最大迭代次数
        # parameters.extend(["-maxiters", str(self.max_iterations.value())])
        
        # # 其他选项
        # if not self.diags_checkbox.isChecked():
        #     parameters.append("-diags")
            
        # if not self.smooth_checkbox.isChecked():
        #     parameters.append("-smooth")
            
        return parameters
    
    def run_analysis(self):
        """运行比对分析"""
        # 检查输入文件
        input_files = self.prepare_input_files()
        if not input_files:
            return
        
        # 检查工具路径
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", "Tool path not configured correctly.")
            return
        
        # 更新UI状态
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        # 创建并启动比对线程
        self.alignment_thread = Muscle5Thread(
            self.tool_path, 
            input_files, 
            self.get_parameters(), 
            self.imported_files
        )
        
        # 连接信号
        self.alignment_thread.progress.connect(self.progress_bar.setFormat)
        self.alignment_thread.finished.connect(self.analysis_finished)
        self.alignment_thread.error.connect(self.analysis_error)
        self.alignment_thread.console_output.connect(self.add_console_message)
        
        # 启动线程
        self.alignment_thread.start()
    
    def analysis_finished(self, output_files, html_files):
        """分析完成处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # Save output files and reports
        self.reports = html_files
        self.update_report_combo()
        
        # Show results
        if html_files:
            self.show_current_report()
            self.tab_widget.setCurrentIndex(1)  # Switch to output preview tab
        
        self.add_console_message("Alignment completed successfully!", "info")
        QMessageBox.information(self, "Success", "Alignment completed successfully!")
        
        # 显示导入按钮（仅在从平台导入数据时显示）
        if self.import_from == "YR_MPEA":
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)
        
        # 保存输出文件路径供导入使用
        self.alignment_output_files = output_files

    def analysis_error(self, error_msg):
        """分析错误处理"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "Error", f"Alignment failed: {error_msg}")
        self.tab_widget.setCurrentIndex(2)  # Switch to console tab
    
    def stop_analysis(self):
        """停止分析"""
        if hasattr(self, 'alignment_thread') and self.alignment_thread.isRunning():
            self.alignment_thread.stop()
            self.alignment_thread.terminate()
            self.alignment_thread.wait()
        
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

    def import_to_platform(self):
        """将比对结果导入到当前平台"""
        if not hasattr(self, 'alignment_output_files') or not self.alignment_output_files:
            QMessageBox.warning(self, "Warning", "No alignment results to import.")
            return
            
        try:
            # 解析比对结果文件
            sequences = []
            for output_file in self.alignment_output_files:
                # 读取FASTA文件
                for record in SeqIO.parse(output_file, "fasta"):
                    sequences.append(record)
            
            if not sequences:
                QMessageBox.warning(self, "Warning", "No sequences found in alignment results.")
                return
                
            # 发送信号将数据导入到平台
            self.import_alignment_to_platform(sequences)
            
            # 显示成功消息
            QMessageBox.information(self, "Success", f"Successfully imported {len(sequences)} sequences to the platform.")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import alignment results: {str(e)}")

    def import_alignment_to_platform(self, sequences):
        """将比对结果导入到平台的工作区"""
        # 发送信号将数据导入到平台
        self.import_alignment_signal.emit(sequences)

    def copy_citation(self):
        """复制引用信息到剪贴板"""
        from PyQt5.QtWidgets import QApplication
        QApplication.clipboard().setText(self.get_citation())
        self.add_console_message("Citation copied to clipboard", "info")


class Muscle5PluginEntry:
    """MUSCLE5插件入口类"""
    
    def run(self, import_from=None, import_data=None):
        """运行插件"""
        return Muscle5Plugin(import_from, import_data)
