#
# sequence_concatenator_plugin.py
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

import os
from collections import defaultdict
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
                             QTableWidget, QTableWidgetItem, QFileDialog,
                             QMessageBox, QHeaderView, QGroupBox, QLabel,
                             QLineEdit, QDialog, QRadioButton, QButtonGroup, QTextEdit)
from PyQt5.QtCore import Qt, QRegExp
from PyQt5.QtGui import QFont, QSyntaxHighlighter, QTextCharFormat, QColor, QTextCursor
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np


class NexusHighlighter(QSyntaxHighlighter):
    """NEXUS 语法高亮器，使用 Monokai 主题"""

    def __init__(self, parent=None):
        super().__init__(parent)

        # 定义格式
        self.begin_end_format = QTextCharFormat()
        self.begin_end_format.setForeground(QColor("#66D9EF"))  # 青色
        self.begin_end_format.setFontWeight(QFont.Bold)

        self.command_format = QTextCharFormat()
        self.command_format.setForeground(QColor("#F92772"))  # 红橙色
        self.command_format.setFontWeight(QFont.Bold)

        self.setting_format = QTextCharFormat()
        self.setting_format.setForeground(QColor("#E6DB74"))  # 黄色

        self.value_format = QTextCharFormat()
        self.value_format.setForeground(QColor("#AE81FF"))  # 紫色

        self.comment_format = QTextCharFormat()
        self.comment_format.setForeground(QColor("#75715E"))  # 灰色

        # 编译正则表达式
        self.begin_end_regex = QRegExp(r"^\s*(begin|end)\s*;?\s*$")
        self.command_regex = QRegExp(r"\b(charset|end)\b")
        self.settings_regex = QRegExp(r"([a-zA-Z_]+)=([^\s;]+)")
        self.comment_regex = QRegExp(r"\[.*?\]")

    def highlightBlock(self, text):
        # 高亮注释
        pos = 0
        index = self.comment_regex.indexIn(text, pos)
        while index >= 0:
            length = self.comment_regex.matchedLength()
            self.setFormat(index, length, self.comment_format)
            pos = index + length
            index = self.comment_regex.indexIn(text, pos)

        # 高亮 begin/end 行
        pos = self.begin_end_regex.indexIn(text)
        if pos >= 0:
            self.setFormat(0, len(text), self.begin_end_format)
            return  # 整行已经高亮，不需要再处理

        # 高亮命令
        pos = 0
        index = self.command_regex.indexIn(text, pos)
        while index >= 0:
            length = self.command_regex.matchedLength()
            self.setFormat(index, length, self.command_format)
            pos = index + length
            index = self.command_regex.indexIn(text, pos)

        # 高亮设置项和值
        pos = 0
        index = self.settings_regex.indexIn(text, pos)
        while index >= 0:
            matched_str = self.settings_regex.cap(0)
            setting_part = self.settings_regex.cap(1)
            value_part = self.settings_regex.cap(2)

            # 高亮设置项（等号前的部分）
            setting_pos = text.find(setting_part, index)
            self.setFormat(setting_pos, len(setting_part), self.setting_format)

            # 高亮值（等号后的部分）
            value_pos = text.find(value_part, setting_pos + len(setting_part))
            self.setFormat(value_pos, len(value_part), self.value_format)

            pos = index + len(matched_str)
            index = self.settings_regex.indexIn(text, pos)


class NexusDisplayDialog(QDialog):
    """显示 NEXUS 分区文件的对话框"""

    def __init__(self, parent=None, nexus_content=""):
        super().__init__(parent)
        self.setWindowTitle("NEXUS Partition File")
        self.setMinimumSize(800, 500)
        self.nexus_content = nexus_content

        layout = QVBoxLayout()

        # 创建文本编辑器
        self.text_edit = QTextEdit()
        self.text_edit.setFont(QFont("Consolas", 10))
        self.text_edit.setStyleSheet("""
            QTextEdit {
                background-color: #272822;
                color: #f8f8f2;
                border: 1px solid #3c3c3c;
            }
        """)
        self.text_edit.setPlainText(self.nexus_content)

        # 添加语法高亮
        self.highlighter = NexusHighlighter(self.text_edit.document())

        layout.addWidget(self.text_edit)

        # 按钮布局
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        self.copy_button = QPushButton("Copy to Clipboard")
        self.copy_button.clicked.connect(self.copy_to_clipboard)
        button_layout.addWidget(self.copy_button)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        button_layout.addWidget(close_btn)

        layout.addLayout(button_layout)
        self.setLayout(layout)

    def copy_to_clipboard(self):
        """复制内容到剪贴板"""
        from PyQt5.QtWidgets import QApplication
        QApplication.clipboard().setText(self.nexus_content)
        QMessageBox.information(self, "Success", "NEXUS partition file copied to clipboard!")


class SequenceConcatenatorPlugin(QWidget):
    """
    Sequence Concatenator Plugin for merging multiple sequence alignments.
    Concatenates sequences with matching descriptions, filling missing sequences with ?.
    Outputs concatenated FASTA and NEXUS partition file.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.sequence_files = []  # List of sequence file paths
        self.dataset_info = []  # Store dataset info for overview figure
        self.all_descriptions = []  # Store all sequence descriptions for overview figure
        self.init_ui()

    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # Description
        desc_label = QLabel(
            "Concatenate multiple sequence alignments (FASTA/PHYLIP). "
            "Sequences with matching descriptions will be concatenated in the order of input files. "
            "Missing sequences will be filled with '?'."
        )
        desc_label.setWordWrap(True)
        desc_label.setStyleSheet("color: #6c757d; padding: 5px;")
        main_layout.addWidget(desc_label)

        # Input files group
        input_group = QGroupBox("Input Alignment Files")
        input_layout = QVBoxLayout()
        input_group.setLayout(input_layout)

        # File table
        self.file_table = QTableWidget()
        self.file_table.setColumnCount(3)
        self.file_table.setHorizontalHeaderLabels(['File Name', 'Sequences', 'Length'])
        self.file_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.file_table.setSelectionMode(QTableWidget.MultiSelection)
        self.file_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.file_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.file_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        input_layout.addWidget(self.file_table)

        # Button layout for file management
        file_button_layout = QHBoxLayout()
        self.add_button = QPushButton("Add Files")
        self.add_button.clicked.connect(self.add_files)
        self.remove_button = QPushButton("Remove Selected")
        self.remove_button.clicked.connect(self.remove_selected_files)
        self.clear_button = QPushButton("Clear All")
        self.clear_button.clicked.connect(self.clear_all_files)

        file_button_layout.addWidget(self.add_button)
        file_button_layout.addWidget(self.remove_button)
        file_button_layout.addWidget(self.clear_button)
        file_button_layout.addStretch()

        input_layout.addLayout(file_button_layout)
        main_layout.addWidget(input_group)

        # Output options group
        output_group = QGroupBox("Output Options")
        output_layout = QVBoxLayout()
        output_group.setLayout(output_layout)

        # FASTA output file
        fasta_layout = QHBoxLayout()
        fasta_label = QLabel("Concatenated FASTA:")
        fasta_label.setMinimumWidth(150)
        self.fasta_output_edit = QLineEdit()
        self.fasta_output_edit.setPlaceholderText("Select output FASTA file...")
        self.browse_fasta_button = QPushButton("Browse")
        self.browse_fasta_button.clicked.connect(self.browse_fasta_file)

        fasta_layout.addWidget(fasta_label)
        fasta_layout.addWidget(self.fasta_output_edit)
        fasta_layout.addWidget(self.browse_fasta_button)
        output_layout.addLayout(fasta_layout)

        # NEXUS partition file output option
        nexus_layout = QHBoxLayout()
        nexus_label = QLabel("NEXUS Partition:")
        nexus_label.setMinimumWidth(150)
        self.nexus_radio_group = QButtonGroup(self)
        self.show_text_radio = QRadioButton("Show Text (Default)")
        self.show_text_radio.setChecked(True)
        self.to_file_radio = QRadioButton("To File")
        self.nexus_radio_group.addButton(self.show_text_radio)
        self.nexus_radio_group.addButton(self.to_file_radio)

        self.nexus_output_edit = QLineEdit()
        self.nexus_output_edit.setPlaceholderText("Select output NEXUS file...")
        self.nexus_output_edit.setEnabled(False)
        self.browse_nexus_button = QPushButton("Browse")
        self.browse_nexus_button.clicked.connect(self.browse_nexus_file)
        self.browse_nexus_button.setEnabled(False)

        # 连接单选按钮信号
        self.show_text_radio.toggled.connect(self.on_nexus_option_changed)
        self.to_file_radio.toggled.connect(self.on_nexus_option_changed)

        nexus_layout.addWidget(nexus_label)
        nexus_layout.addWidget(self.show_text_radio)
        nexus_layout.addWidget(self.to_file_radio)
        output_layout.addLayout(nexus_layout)

        nexus_file_layout = QHBoxLayout()
        nexus_file_layout.addSpacing(150)
        nexus_file_layout.addWidget(self.nexus_output_edit)
        nexus_file_layout.addWidget(self.browse_nexus_button)
        output_layout.addLayout(nexus_file_layout)

        main_layout.addWidget(output_group)

        # Concatenate button
        self.concat_button = QPushButton("Concatenate Sequences")
        self.concat_button.clicked.connect(self.concatenate_sequences)
        main_layout.addWidget(self.concat_button)

        # Overview figure button
        self.overview_button = QPushButton("Display Overview Figure")
        self.overview_button.clicked.connect(self.display_overview_figure)
        self.overview_button.setEnabled(False)  # Disabled until concatenation is done
        main_layout.addWidget(self.overview_button)

        # main_layout.addStretch()

    def update_file_table(self):
        """Update the file table with current sequence files and their info."""
        self.file_table.setRowCount(len(self.sequence_files))

        for row, file_path in enumerate(self.sequence_files):
            # File name column
            file_name = os.path.basename(file_path)
            name_item = QTableWidgetItem(file_name)
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            self.file_table.setItem(row, 0, name_item)

            # Get sequence info
            try:
                sequences = list(SeqIO.parse(file_path, "fasta"))
                num_seqs = len(sequences)
                if num_seqs > 0:
                    seq_length = len(str(sequences[0].seq))
                else:
                    seq_length = 0
            except Exception as e:
                num_seqs = 0
                seq_length = 0

            # Number of sequences column
            num_item = QTableWidgetItem(str(num_seqs))
            num_item.setFlags(num_item.flags() & ~Qt.ItemIsEditable)
            self.file_table.setItem(row, 1, num_item)

            # Length column
            length_item = QTableWidgetItem(str(seq_length))
            length_item.setFlags(length_item.flags() & ~Qt.ItemIsEditable)
            self.file_table.setItem(row, 2, length_item)

    def on_nexus_option_changed(self):
        """Handle NEXUS output option change."""
        is_to_file = self.to_file_radio.isChecked()
        self.nexus_output_edit.setEnabled(is_to_file)
        self.browse_nexus_button.setEnabled(is_to_file)

    def add_files(self):
        """Add sequence files to the list."""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select Alignment Files",
            "",
            "FASTA Files (*.fasta *.fa *.fas *.fna *.ffn *.faa *.frn);;PHYLIP Files (*.phy *.phylip);;All Files (*)"
        )

        if not file_paths:
            return

        for file_path in file_paths:
            if file_path not in self.sequence_files:
                self.sequence_files.append(file_path)

        self.update_file_table()

    def remove_selected_files(self):
        """Remove selected files from the list."""
        selected_rows = sorted(set(item.row() for item in self.file_table.selectedItems()), reverse=True)

        for row in selected_rows:
            if row < len(self.sequence_files):
                del self.sequence_files[row]

        self.update_file_table()

    def clear_all_files(self):
        """Clear all files from the list."""
        self.sequence_files.clear()
        self.update_file_table()

    def browse_fasta_file(self):
        """Browse for FASTA output file location."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Output FASTA File",
            "",
            "FASTA Files (*.fasta *.fa);;All Files (*)"
        )

        if file_path:
            self.fasta_output_edit.setText(file_path)

    def browse_nexus_file(self):
        """Browse for NEXUS output file location."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Output NEXUS File",
            "",
            "NEXUS Files (*.nexus *.nex);;All Files (*)"
        )

        if file_path:
            self.nexus_output_edit.setText(file_path)

    def concatenate_sequences(self):
        """Concatenate all sequence alignments."""
        # Validate inputs
        if not self.sequence_files:
            QMessageBox.warning(self, "Warning", "Please add at least one sequence file!")
            return

        fasta_output = self.fasta_output_edit.text()
        if not fasta_output:
            QMessageBox.warning(self, "Warning", "Please specify an output FASTA file!")
            return

        # Check if all input files exist
        missing_files = [f for f in self.sequence_files if not os.path.exists(f)]
        if missing_files:
            QMessageBox.warning(
                self,
                "Warning",
                f"The following files do not exist:\n" + "\n".join(missing_files)
            )
            return

        try:
            # Parse all input files
            all_datasets = []
            dataset_info = []

            for file_path in self.sequence_files:
                try:
                    sequences = list(SeqIO.parse(file_path, "fasta"))
                    if not sequences:
                        QMessageBox.warning(
                            self,
                            "Warning",
                            f"File {os.path.basename(file_path)} is empty or invalid!"
                        )
                        return

                    # Check if all sequences in this dataset have the same length
                    first_length = len(str(sequences[0].seq))
                    for i, seq in enumerate(sequences):
                        seq_length = len(str(seq.seq))
                        if seq_length != first_length:
                            QMessageBox.warning(
                                self,
                                "Warning",
                                f"File {os.path.basename(file_path)} contains sequences of different lengths!\n"
                                f"Sequence {seq.description}: {seq_length} bp (expected {first_length} bp)"
                            )
                            return

                    # Store sequences by description
                    seq_dict = {seq.description: str(seq.seq) for seq in sequences}
                    all_datasets.append(seq_dict)

                    dataset_info.append({
                        'file': os.path.basename(file_path),
                        'length': first_length,
                        'count': len(sequences)
                    })

                except Exception as e:
                    QMessageBox.warning(
                        self,
                        "Error",
                        f"Failed to parse file {os.path.basename(file_path)}: {str(e)}"
                    )
                    return

            # Collect all unique sequence descriptions
            all_descriptions = set()
            for dataset in all_datasets:
                all_descriptions.update(dataset.keys())
            all_descriptions = sorted(list(all_descriptions))

            # Concatenate sequences
            concatenated_sequences = []
            for desc in all_descriptions:
                concatenated_seq = ""
                for i, dataset in enumerate(all_datasets):
                    if desc in dataset:
                        concatenated_seq += dataset[desc]
                    else:
                        # Fill with '?' for missing sequence (use this dataset's own length)
                        concatenated_seq += "?" * dataset_info[i]['length']
                concatenated_sequences.append((desc, concatenated_seq))

            # Write concatenated FASTA
            with open(fasta_output, 'w', encoding='utf-8') as f:
                for desc, seq in concatenated_sequences:
                    f.write(f">{desc}\n{seq}\n")

            # Generate NEXUS partition file
            nexus_content = self.generate_nexus_partition(dataset_info, fasta_output)

            # Handle NEXUS output
            if self.show_text_radio.isChecked():
                # Show in dialog
                dialog = NexusDisplayDialog(self, nexus_content)
                dialog.exec_()
            else:
                # Write to file
                nexus_output = self.nexus_output_edit.text()
                if not nexus_output:
                    QMessageBox.warning(self, "Warning", "Please specify an output NEXUS file!")
                    return

                with open(nexus_output, 'w', encoding='utf-8') as f:
                    f.write(nexus_content)

            # Success message
            total_length = sum(info['length'] for info in dataset_info)
            msg = f"Successfully concatenated {len(concatenated_sequences)} sequence(s) into:\n{fasta_output}\n\n"
            msg += f"Total concatenated length: {total_length} bp\n"
            msg += f"Number of partitions: {len(dataset_info)}"

            if self.to_file_radio.isChecked():
                msg += f"\n\nNEXUS partition file written to:\n{nexus_output}"

            # Save data for overview figure
            self.dataset_info = dataset_info
            self.all_descriptions = all_descriptions
            self.overview_button.setEnabled(True)

            QMessageBox.information(self, "Success", msg)

        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to concatenate sequences: {str(e)}"
            )

    def generate_nexus_partition(self, dataset_info, fasta_output):
        """Generate NEXUS partition file content."""
        import string

        lines = ["#nexus", "begin sets;"]

        current_pos = 1
        for info in dataset_info:
            file_name = os.path.splitext(info['file'])[0]
            # Replace all non-underscore punctuation with underscores
            for char in string.punctuation:
                if char != '_':
                    file_name = file_name.replace(char, '_')
            end_pos = current_pos + info['length'] - 1
            lines.append(f"\tcharset {file_name}={current_pos}-{end_pos};")
            current_pos = end_pos + 1

        lines.append("end;")

        return "\n".join(lines)

    def display_overview_figure(self):
        """Display overview figure showing sequence presence/absence across datasets."""
        if not self.dataset_info or not self.all_descriptions:
            QMessageBox.warning(self, "Warning", "No data available. Please concatenate sequences first!")
            return

        try:
            # Get dataset names
            dataset_names = [os.path.splitext(info['file'])[0] for info in self.dataset_info]

            # Create presence/absence matrix
            # Rows: sequence descriptions, Columns: datasets
            matrix = np.zeros((len(self.all_descriptions), len(self.dataset_info)))

            # Populate matrix with presence/absence data
            # Re-parse files to get actual presence/absence information
            for i, file_path in enumerate(self.sequence_files):
                try:
                    sequences = list(SeqIO.parse(file_path, "fasta"))
                    seq_descriptions = {seq.description for seq in sequences}
                    for j, desc in enumerate(self.all_descriptions):
                        if desc in seq_descriptions:
                            matrix[j, i] = 100
                except Exception as e:
                    continue

            # Calculate retention percentage for each dataset
            retention_percentages = []
            for i in range(len(self.dataset_info)):
                count = np.sum(matrix[:, i])
                percentage = (count / len(self.all_descriptions)) * 100
                retention_percentages.append(percentage)

            # Create extended matrix with overview row at the top
            extended_matrix = np.vstack([retention_percentages, matrix])

            # Create figure
            fig, ax = plt.subplots(figsize=(max(10, len(dataset_names) * 0.8), max(9, (len(self.all_descriptions) + 1) * 0.3)))

            # Display heatmap using imshow with viridis colormap
            im = ax.imshow(extended_matrix, cmap='viridis', aspect='auto', vmin=0, vmax=100)

            # Set axis labels
            ax.set_xticks(np.arange(len(dataset_names)))
            ax.set_yticks(np.arange(len(self.all_descriptions) + 1))

            # Create y-axis labels with overview row
            y_labels = ['Retention %'] + list(self.all_descriptions)
            ax.set_yticklabels(y_labels)

            ax.set_xticklabels(dataset_names, rotation=45, ha='right')

            # Add title
            ax.set_title('Sequence Presence/Absence Overview', fontsize=14, fontweight='bold')

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            cbar.set_label('Retention % / Presence (100) / Absence (0)', rotation=270, labelpad=15)

            # Add text annotations for percentages in the overview row
            for i, pct in enumerate(retention_percentages):
                text_color = 'white' if pct < 50 else 'black'
                ax.text(i, 0, f'{pct:.1f}%', ha='center', va='center',
                       color=text_color, fontsize=10, fontweight='bold')

            # Add grid lines to separate overview row from data
            ax.axhline(y=0.5, color='white', linewidth=2)

            # Adjust layout to prevent label cutoff
            plt.tight_layout()

            # Show the figure
            plt.show()

        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to display overview figure: {str(e)}"
            )


class SequenceConcatenatorPluginEntry:
    def run(self):
        return SequenceConcatenatorPlugin()