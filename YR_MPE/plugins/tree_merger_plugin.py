#
# tree_merger_plugin.py
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
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
                             QTableWidget, QTableWidgetItem, QFileDialog,
                             QMessageBox, QHeaderView, QGroupBox, QLabel,
                             QLineEdit)
from PyQt5.QtCore import Qt


class TreeMergerPlugin(QWidget):
    """
    TreeMerger Plugin for merging multiple phylogenetic trees into a single file.
    Useful for preparing multi-locus gene trees for species tree inference.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tree_files = []  # List of tree file paths
        self.init_ui()

    def init_ui(self):
        """Initialize the user interface."""
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # Description
        desc_label = QLabel(
            "Merge multiple Newick tree files into a single file for species tree inference. "
            "This is useful when preparing gene trees for tools like ASTRAL or PhyloNet."
        )
        desc_label.setWordWrap(True)
        desc_label.setStyleSheet("color: #6c757d; padding: 5px;")
        main_layout.addWidget(desc_label)

        # Input files group
        input_group = QGroupBox("Input Tree Files")
        input_layout = QVBoxLayout()
        input_group.setLayout(input_layout)

        # File table
        self.file_table = QTableWidget()
        self.file_table.setColumnCount(2)
        self.file_table.setHorizontalHeaderLabels(['File Name', 'File Path'])
        self.file_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.file_table.setSelectionMode(QTableWidget.MultiSelection)
        self.file_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.file_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
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

        # Output file group
        output_group = QGroupBox("Output File")
        output_layout = QVBoxLayout()
        output_group.setLayout(output_layout)

        # Output file selection
        output_file_layout = QHBoxLayout()
        self.output_file_edit = QLineEdit()
        self.output_file_edit.setPlaceholderText("Select output file...")
        self.browse_output_button = QPushButton("Browse")
        self.browse_output_button.clicked.connect(self.browse_output_file)

        output_file_layout.addWidget(self.output_file_edit)
        output_file_layout.addWidget(self.browse_output_button)
        output_layout.addLayout(output_file_layout)



        main_layout.addWidget(output_group)

        # Merge button
        self.merge_button = QPushButton("Merge Trees")
        # self.merge_button.setStyleSheet(
        #     "QPushButton { background-color: #007bff; color: white; font-weight: bold; padding: 10px; } "
        #     "QPushButton:hover { background-color: #0056b3; }"
        # )
        self.merge_button.clicked.connect(self.merge_trees)
        main_layout.addWidget(self.merge_button)        

        # main_layout.addStretch()

    def update_file_table(self):
        """Update the file table with current tree files."""
        self.file_table.setRowCount(len(self.tree_files))

        for row, file_path in enumerate(self.tree_files):
            # File name column
            file_name = os.path.basename(file_path)
            name_item = QTableWidgetItem(file_name)
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            self.file_table.setItem(row, 0, name_item)

            # File path column
            path_item = QTableWidgetItem(file_path)
            path_item.setFlags(path_item.flags() & ~Qt.ItemIsEditable)
            self.file_table.setItem(row, 1, path_item)

    def add_files(self):
        """Add tree files to the list."""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select Tree Files",
            "",
            "Newick Files (*.nwk *.newick *.tre *.tree);;All Files (*)"
        )

        if not file_paths:
            return

        for file_path in file_paths:
            if file_path not in self.tree_files:
                self.tree_files.append(file_path)

        self.update_file_table()

    def remove_selected_files(self):
        """Remove selected files from the list."""
        selected_rows = sorted(set(item.row() for item in self.file_table.selectedItems()), reverse=True)

        for row in selected_rows:
            if row < len(self.tree_files):
                del self.tree_files[row]

        self.update_file_table()

    def clear_all_files(self):
        """Clear all files from the list."""
        self.tree_files.clear()
        self.update_file_table()

    def browse_output_file(self):
        """Browse for output file location."""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Select Output File",
            "",
            "Newick Files (*.nwk *.newick *.tre *.tree);;All Files (*)"
        )

        if file_path:
            self.output_file_edit.setText(file_path)

    def merge_trees(self):
        """Merge all tree files into a single output file."""
        # Validate inputs
        if not self.tree_files:
            QMessageBox.warning(self, "Warning", "Please add at least one tree file!")
            return

        output_file = self.output_file_edit.text()
        if not output_file:
            QMessageBox.warning(self, "Warning", "Please specify an output file!")
            return

        # Check if all input files exist
        missing_files = [f for f in self.tree_files if not os.path.exists(f)]
        if missing_files:
            QMessageBox.warning(
                self,
                "Warning",
                f"The following files do not exist:\n" + "\n".join(missing_files)
            )
            return

        try:
            # Read and merge tree strings
            merged_trees = []

            for tree_file in self.tree_files:
                try:
                    with open(tree_file, 'r', encoding='utf-8') as f:
                        content = f.read().strip()
                        if content:
                            # Split by newlines in case file contains multiple trees
                            trees = [line.strip() for line in content.split('\n') if line.strip()]
                            merged_trees.extend(trees)
                except Exception as e:
                    QMessageBox.warning(
                        self,
                        "Error",
                        f"Failed to read file {os.path.basename(tree_file)}: {str(e)}"
                    )
                    return

            if not merged_trees:
                QMessageBox.warning(self, "Warning", "No valid trees found in the input files!")
                return

            # Write merged trees to output file
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write('\n'.join(merged_trees) + '\n')

            # Success message
            QMessageBox.information(
                self,
                "Success",
                f"Successfully merged {len(merged_trees)} tree(s) into:\n{output_file}"
            )

        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to merge trees: {str(e)}"
            )

class TreeMergerPluginEntry:
    def run(self):
        return TreeMergerPlugin()