# seqmatrix.py
#
# Copyright (c) 2025 Zhi-Jie Xu
#
# This file is part of YR-MPE & the re-implement of SeqMatrix (https://github.com/Gipsy-The-Sheller/SeqMartix).
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

from PyQt5.QtWidgets import (QWidget, QMainWindow, 
                             QHBoxLayout, QVBoxLayout,
                             QTableWidget, QTableWidgetItem,
                             QSplitter, QTreeWidget, QTreeWidgetItem,
                             QMenu, QMenuBar, QAction,
                             QSizePolicy, QFileDialog, QMessageBox, QInputDialog, QApplication,
                             QDialog, QLineEdit, QTextEdit, QPushButton, QLabel, QHeaderView,
                             QCheckBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QMimeData, QByteArray, QDataStream, QIODevice
from PyQt5.QtGui import QIcon, QDrag, QColor, QKeyEvent

import os
import pandas as pd
import uuid

from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class NCBIdownloader(QThread):
    # Define signals for communication with the main thread
    log_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(dict)
    
    def __init__(self, accessions):
        super().__init__()
        self.accessions = accessions
        self.sequences = {}
        
    def run(self):
        Entrez.email = "example@example.com"  # Should be set to actual email
        for accession in self.accessions:
            try:
                with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
                    record = SeqIO.read(handle, "genbank")
                    full_id = record.id + " " + record.description[len(record.id):].strip()
                    self.sequences[full_id] = str(record.seq)
                    self.log_signal.emit(f"[LOG] Downloaded sequence: {full_id}")
            except Exception as e:
                self.log_signal.emit(f"[ERROR] Failed to download {accession}: {e}")
        self.finished_signal.emit(self.sequences)

class SequenceTreeWidgetItem(QTreeWidgetItem):
    """Custom TreeWidgetItem that supports dragging sequences"""
    def __init__(self, sequence_uuid, sequence_name, sequence_data, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sequence_uuid = sequence_uuid
        self.sequence_name = sequence_name
        self.sequence_data = sequence_data
        self.setFlags(self.flags() | Qt.ItemIsDragEnabled)

class DraggableTableWidget(QTableWidget):
    """Custom TableWidget that supports dropping sequences and keyboard operations"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.parent_window = parent
        self.horizontalHeader().sectionDoubleClicked.connect(self.on_header_double_clicked)
        
    def dragEnterEvent(self, event):
        if event.mimeData().hasFormat('application/x-seqmatrix-sequence'):
            event.acceptProposedAction()
        else:
            event.ignore()
            
    def dragMoveEvent(self, event):
        if event.mimeData().hasFormat('application/x-seqmatrix-sequence'):
            event.acceptProposedAction()
        else:
            event.ignore()
            
    def dropEvent(self, event):
        if event.mimeData().hasFormat('application/x-seqmatrix-sequence'):
            # Retrieve the data
            data = event.mimeData().data('application/x-seqmatrix-sequence')
            stream = QDataStream(data, QIODevice.ReadOnly)
            sequence_uuid = bytes(stream.readString()).decode('utf-8')
            sequence_name = bytes(stream.readString()).decode('utf-8')
            
            # Get drop position
            pos = event.pos()
            row = self.rowAt(pos.y())
            col = self.columnAt(pos.x())
            
            # Only allow dropping in partition columns (not sequence name column)
            if row >= 0 and col >= 0 and col > 0 and col < self.columnCount():  # col > 0 to skip sequence name column
                # Set the item with UUID as data and sequence name as display text
                item = QTableWidgetItem(sequence_name)
                item.setData(Qt.UserRole, sequence_uuid)  # Store UUID in UserRole
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)  # Make non-editable
                self.setItem(row, col, item)
                event.acceptProposedAction()
                
                # Update highlighting after drop
                if hasattr(self.parent_window, 'update_highlighting'):
                    self.parent_window.update_highlighting()
            else:
                event.ignore()
        else:
            event.ignore()
            
    def keyPressEvent(self, event):
        """Handle key press events for the table"""
        # Check if Backspace or Delete key is pressed
        if event.key() in (Qt.Key_Backspace, Qt.Key_Delete):
            # Get selected items
            selected_items = self.selectedItems()
            cleared_items = False
            
            for item in selected_items:
                row = self.row(item)
                col = self.column(item)
                
                # Only process items in partition columns (col > 0) that have UUID data
                if col > 0 and item.data(Qt.UserRole):
                    # Clear the item content and make it editable
                    item.setText("")
                    item.setData(Qt.UserRole, None)  # Remove UUID data
                    item.setFlags(item.flags() | Qt.ItemIsEditable)  # Make editable
                    cleared_items = True
            
            # If we cleared any items, update highlighting
            if cleared_items and hasattr(self.parent_window, 'update_highlighting'):
                self.parent_window.update_highlighting()
        else:
            # For all other keys, use the default behavior
            super().keyPressEvent(event)
            
    def on_header_double_clicked(self, logical_index):
        """Handle double click on horizontal header"""
        # Process clicks on all columns
        if hasattr(self.parent_window, 'on_header_double_clicked'):
            self.parent_window.on_header_double_clicked(logical_index)

class SequenceNameEditDialog(QDialog):
    """Dialog for editing sequence name column"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle("Edit Sequence Names")
        self.setModal(True)
        self.resize(400, 300)
        
        layout = QVBoxLayout()
        
        # Info label
        layout.addWidget(QLabel("Sequence Names (one per line):"))
        self.sequence_list_edit = QTextEdit()
        layout.addWidget(self.sequence_list_edit)
        
        # Checkbox for table header
        self.include_header_checkbox = QCheckBox("Include table header")
        self.include_header_checkbox.setChecked(False)
        layout.addWidget(self.include_header_checkbox)
        
        # OK button
        self.ok_button = QPushButton("OK")
        self.ok_button.clicked.connect(self.accept)
        layout.addWidget(self.ok_button)
        
        self.setLayout(layout)
        
    def get_data(self):
        """Return the entered data"""
        return self.sequence_list_edit.toPlainText(), self.include_header_checkbox.isChecked()

class PartitionEditDialog(QDialog):
    """Dialog for editing partition columns"""
    def __init__(self, current_name, parent=None):
        super().__init__(parent)
        self.current_name = current_name
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle("Edit Partition")
        self.setModal(True)
        self.resize(400, 300)
        
        layout = QVBoxLayout()
        
        # Name input
        layout.addWidget(QLabel("Partition Name:"))
        self.name_edit = QLineEdit(self.current_name)
        layout.addWidget(self.name_edit)
        
        # Sequence list input
        layout.addWidget(QLabel("Sequence Accessions (one per line):"))
        self.sequence_list_edit = QTextEdit()
        layout.addWidget(self.sequence_list_edit)
        
        # Checkbox for table header
        self.include_header_checkbox = QCheckBox("Include table header")
        self.include_header_checkbox.setChecked(False)
        layout.addWidget(self.include_header_checkbox)
        
        # OK button
        self.ok_button = QPushButton("OK")
        self.ok_button.clicked.connect(self.accept)
        layout.addWidget(self.ok_button)
        
        self.setLayout(layout)
        
    def get_data(self):
        """Return the entered data"""
        return self.name_edit.text(), self.sequence_list_edit.toPlainText(), self.include_header_checkbox.isChecked()

class SeqMatrixUI(QMainWindow):
    def __init__(self):
        super().__init__()
        
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))

        self.init_ui()
        
    def init_ui(self):
        # 3:7 QSplitter
        self.splitter = QSplitter(Qt.Horizontal)
        self.setCentralWidget(self.splitter)

        # Left side: Sequence Explorer (Tree)
        self.tree_explorer = QTreeWidget()
        self.tree_explorer.setHeaderLabel("Sequence Explorer")
        self.tree_explorer.setDragEnabled(True)
        self.tree_explorer.setDragDropMode(QTreeWidget.DragOnly)
        self.splitter.addWidget(self.tree_explorer)

        # Right Side: Sequence Matrix
        self.right_side_widget = QWidget()
        self.right_side_widget.setContentsMargins(0, 0, 0, 0)
        self.right_side_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.right_layout = QVBoxLayout()
        self.right_layout.setContentsMargins(0, 0, 0, 0)
        self.right_side_widget.setLayout(self.right_layout)
        self.init_right_side()
        self.splitter.addWidget(self.right_side_widget)

        self.splitter.setSizes([300, 700])

        # Window properties
        self.setWindowTitle("SeqMatrix")
        self.resize(1000, 600)

    def init_right_side(self):
        
        # Top: menu bar with QAction

        menu_bar = QMenuBar()
        self.right_layout.addWidget(menu_bar)

        # QAction group

        self.import_sequence_action = QAction("Import from Sequences", self)
        menu_bar.addAction(self.import_sequence_action)
        import_sequence_action_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "import_sequence.svg"))
        self.import_sequence_action.setIcon(import_sequence_action_icon)

        self.import_from_nexus_action = QAction("Import from NEXUS", self)
        menu_bar.addAction(self.import_from_nexus_action)
        import_from_nexus_action_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "import_nexus.svg"))
        self.import_from_nexus_action.setIcon(import_from_nexus_action_icon)

        self.import_from_excel_action = QAction("Import from Excel", self)
        menu_bar.addAction(self.import_from_excel_action)
        import_from_excel_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "import_excel.svg"))
        self.import_from_excel_action.setIcon(import_from_excel_icon)
        
        self.append_new_row_action = QAction("Append New Row", self)
        menu_bar.addAction(self.append_new_row_action)
        append_new_row_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "append_row.svg"))
        self.append_new_row_action.setIcon(append_new_row_icon)

        self.append_new_column_action = QAction("Append New Column", self)
        menu_bar.addAction(self.append_new_column_action)
        append_new_column_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "append_column.svg"))
        self.append_new_column_action.setIcon(append_new_column_icon)

        self.delete_row_action = QAction("Delete Row", self)
        menu_bar.addAction(self.delete_row_action)
        delete_row_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "delete_row.svg"))
        self.delete_row_action.setIcon(delete_row_icon)

        self.delete_column_action = QAction("Delete Column", self)
        menu_bar.addAction(self.delete_column_action)
        delete_column_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "delete_column.svg"))
        self.delete_column_action.setIcon(delete_column_icon)

        self.export_action = QAction("Export", self)
        menu_bar.addAction(self.export_action)
        export_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "export.svg"))
        self.export_action.setIcon(export_icon)

        self.download_from_NCBI_action = QAction("Download from NCBI", self)
        menu_bar.addAction(self.download_from_NCBI_action)
        download_from_NCBI_icon = QIcon(os.path.join(self.plugin_path, "..", "icons", "seqmatrix", "ncbi.svg"))
        self.download_from_NCBI_action.setIcon(download_from_NCBI_icon)

        # QTableWidget
        self.seq_matrix = DraggableTableWidget(self)
        self.right_layout.addWidget(self.seq_matrix)
        
        # essential column: sequence name
        self.seq_matrix.setColumnCount(2)
        self.seq_matrix.setRowCount(4)
        self.seq_matrix.setHorizontalHeaderItem(0, QTableWidgetItem("Sequence Name"))
        self.seq_matrix.setHorizontalHeaderItem(1, QTableWidgetItem("Partition 1"))
        
        # Enable selection
        self.seq_matrix.setSelectionMode(QTableWidget.ContiguousSelection)
        self.seq_matrix.setEditTriggers(QTableWidget.DoubleClicked | QTableWidget.EditKeyPressed)
        
        # Make sequence name column editable
        # Other columns are handled in the drop event to make them non-editable

class Lg_SeqMatrix(SeqMatrixUI):
    
    # Logics of SeqMatrix
    def __init__(self, import_from = None):
        super().__init__()

        # seqmatrix_data now stores mapping from UUID to sequence data
        self.seqmatrix_data = {}  # uuid -> {name: sequence_name, sequence: sequence_data}
        # sequence_data structure:
        # {
        #   "/path/to/directory": {
        #     "filename.fasta": [{"uuid": uuid, "name": sequence_name, "sequence": sequence_data}, ...],
        #     "filename.nex": [{"uuid": uuid, "name": sequence_name, "sequence": sequence_data}, ...]
        #   },
        #   "Downloaded": {
        #     "NCBI Downloads": [{"uuid": uuid, "name": sequence_name, "sequence": sequence_data}, ...]
        #   }
        # }
        self.sequence_data = {}
        
        # Keep track of tree items by UUID for highlighting
        self.tree_items_by_uuid = {}  # uuid -> [SequenceTreeWidgetItem, ...]
        
        # Store accession items for NCBI download mapping
        self.accession_items = {}  # accession -> [(row, col), ...]
        
        # Connect actions to methods
        self.setup_connections()
    
    def setup_connections(self):
        """Setup signal-slot connections"""
        self.import_sequence_action.triggered.connect(self.import_sequences)
        self.import_from_nexus_action.triggered.connect(self.import_nexus)
        self.import_from_excel_action.triggered.connect(self.import_excel)
        self.append_new_row_action.triggered.connect(self.append_row)
        self.append_new_column_action.triggered.connect(self.append_column)
        self.delete_row_action.triggered.connect(self.delete_row)
        self.delete_column_action.triggered.connect(self.delete_column)
        self.export_action.triggered.connect(self.export_data)
        self.download_from_NCBI_action.triggered.connect(self.download_ncbi)

        # Setup custom mime data for tree items
        self.tree_explorer.startDrag = self.tree_start_drag
    
    def tree_start_drag(self, supported_actions):
        """Custom drag start method for the tree widget"""
        item = self.tree_explorer.currentItem()
        if isinstance(item, SequenceTreeWidgetItem):
            drag = QDrag(self.tree_explorer)
            mime_data = QMimeData()
            
            # Create byte array to store UUID and name
            data = QByteArray()
            stream = QDataStream(data, QIODevice.WriteOnly)
            stream.writeString(item.sequence_uuid.encode('utf-8'))
            stream.writeString(item.sequence_name.encode('utf-8'))
            
            mime_data.setData('application/x-seqmatrix-sequence', data)
            drag.setMimeData(mime_data)
            
            drag.exec_(Qt.CopyAction)
    
    def import_sequences(self):
        """Import sequences from FASTA files"""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA files", "", "FASTA files (*.fas *.fasta *.fa)")
        
        if not files:
            return
            
        for file in files:
            file_name = os.path.basename(file)
            dir_path = os.path.dirname(file)
            
            # Group files by directory path
            if dir_path not in self.sequence_data:
                self.sequence_data[dir_path] = {}
                
            if file_name not in self.sequence_data[dir_path]:
                self.sequence_data[dir_path][file_name] = []
            
            try:
                records = list(SeqIO.parse(file, "fasta"))
                for record in records:
                    # Generate unique UUID for each sequence
                    sequence_uuid = str(uuid.uuid4())
                    full_id = record.id + " " + record.description[len(record.id):].strip()
                    sequence = str(record.seq)
                    
                    # Store in both structures
                    self.sequence_data[dir_path][file_name].append({
                        "uuid": sequence_uuid,
                        "name": full_id,
                        "sequence": sequence
                    })
                    self.seqmatrix_data[sequence_uuid] = {
                        "name": full_id,
                        "sequence": sequence
                    }
                    
                QMessageBox.information(self, "Success", f"Imported {len(records)} sequences from {file_name}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to import {file_name}: {str(e)}")
        
        self.update_tree_view()
        self.update_table()
    
    def import_nexus(self):
        """Import sequences from NEXUS files"""
        file, _ = QFileDialog.getOpenFileName(
            self, "Select NEXUS file", "", "NEXUS files (*.nex *.nexus)")
        
        if not file:
            return
            
        try:
            file_name = os.path.basename(file)
            dir_path = os.path.dirname(file)
            
            # Group files by directory path
            if dir_path not in self.sequence_data:
                self.sequence_data[dir_path] = {}
                
            if file_name not in self.sequence_data[dir_path]:
                self.sequence_data[dir_path][file_name] = []
            
            # Parse NEXUS file
            records = list(SeqIO.parse(file, "nexus"))
            for record in records:
                # Generate unique UUID for each sequence
                sequence_uuid = str(uuid.uuid4())
                full_id = record.id + " " + record.description[len(record.id):].strip()
                sequence = str(record.seq)
                
                # Store in both structures
                self.sequence_data[dir_path][file_name].append({
                    "uuid": sequence_uuid,
                    "name": full_id,
                    "sequence": sequence
                })
                self.seqmatrix_data[sequence_uuid] = {
                    "name": full_id,
                    "sequence": sequence
                }
                
            QMessageBox.information(self, "Success", f"Imported {len(records)} sequences from {file_name}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import {file_name}: {str(e)}")
        
        self.update_tree_view()
        self.update_table()
    
    def import_excel(self):
        """Import data from Excel files"""
        file, _ = QFileDialog.getOpenFileName(
            self, "Select Excel file", "", "Excel files (*.xlsx *.xls)")
        
        if not file:
            return
            
        try:
            file_name = os.path.basename(file)
            dir_path = os.path.dirname(file)
            
            # Group files by directory path
            if dir_path not in self.sequence_data:
                self.sequence_data[dir_path] = {}
                
            if file_name not in self.sequence_data[dir_path]:
                self.sequence_data[dir_path][file_name] = []
            
            df = pd.read_excel(file)
            
            # Assuming the first column contains sequence names and other columns contain sequences
            if len(df.columns) < 2:
                raise ValueError("Excel file must have at least 2 columns (sequence names and sequences)")
            
            for _, row in df.iterrows():
                # Generate unique UUID for each sequence
                sequence_uuid = str(uuid.uuid4())
                seq_name = str(row.iloc[0])
                sequence = str(row.iloc[1])  # Taking the second column as sequence
                
                # Store in both structures
                self.sequence_data[dir_path][file_name].append({
                    "uuid": sequence_uuid,
                    "name": seq_name,
                    "sequence": sequence
                })
                self.seqmatrix_data[sequence_uuid] = {
                    "name": seq_name,
                    "sequence": sequence
                }
            
            QMessageBox.information(self, "Success", f"Imported data from {file_name}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import {file_name}: {str(e)}")
        
        self.update_tree_view()
        self.update_table()
    
    def append_row(self):
        """Add a new row to the table"""
        current_row_count = self.seq_matrix.rowCount()
        self.seq_matrix.setRowCount(current_row_count + 1)
        
    def append_column(self):
        """Add a new column to the table"""
        current_col_count = self.seq_matrix.columnCount()
        self.seq_matrix.setColumnCount(current_col_count + 1)
        
        # Set header for new column
        partition_num = current_col_count  # Since first column is "Sequence Name"
        self.seq_matrix.setHorizontalHeaderItem(current_col_count, QTableWidgetItem(f"Partition {partition_num}"))
    
    def delete_row(self):
        """Delete selected rows from the table"""
        selected_ranges = self.seq_matrix.selectedRanges()
        if not selected_ranges:
            QMessageBox.warning(self, "Warning", "Please select rows to delete")
            return
            
        # Delete rows in reverse order to maintain indices
        rows_to_delete = []
        for selected_range in selected_ranges:
            for row in range(selected_range.topRow(), selected_range.bottomRow() + 1):
                rows_to_delete.append(row)
                
        rows_to_delete.sort(reverse=True)
        for row in rows_to_delete:
            self.seq_matrix.removeRow(row)
            
        # Update highlighting after deletion
        self.update_highlighting()
    
    def delete_column(self):
        """Delete selected columns from the table"""
        selected_ranges = self.seq_matrix.selectedRanges()
        if not selected_ranges:
            QMessageBox.warning(self, "Warning", "Please select columns to delete")
            return
            
        # Prevent deletion of the first column (Sequence Name)
        cols_to_delete = []
        for selected_range in selected_ranges:
            for col in range(selected_range.leftColumn(), selected_range.rightColumn() + 1):
                if col > 0:  # Don't allow deletion of Sequence Name column
                    cols_to_delete.append(col)
                    
        if not cols_to_delete:
            QMessageBox.warning(self, "Warning", "Cannot delete Sequence Name column")
            return
            
        cols_to_delete.sort(reverse=True)
        for col in cols_to_delete:
            self.seq_matrix.removeColumn(col)
            
        # Update highlighting after deletion
        self.update_highlighting()
    
    def export_data(self):
        """Export partitions as FASTA files"""
        directory = QFileDialog.getExistingDirectory(self, "Select Export Directory")
        if not directory:
            return
            
        try:
            # Process each partition column (skip sequence name column)
            for col in range(1, self.seq_matrix.columnCount()):
                # Get partition name
                header_item = self.seq_matrix.horizontalHeaderItem(col)
                partition_name = header_item.text() if header_item else f"Partition_{col}"
                
                # Collect sequences for this partition
                sequences = []
                for row in range(self.seq_matrix.rowCount()):
                    # Get sequence name from first column
                    name_item = self.seq_matrix.item(row, 0)
                    seq_name = name_item.text() if name_item and name_item.text() else f"Sequence_{row+1}"
                    
                    # Get UUID from partition column
                    partition_item = self.seq_matrix.item(row, col)
                    if partition_item and partition_item.data(Qt.UserRole):
                        uuid = partition_item.data(Qt.UserRole)
                        if uuid in self.seqmatrix_data:
                            sequence_data = self.seqmatrix_data[uuid]
                            # Create SeqRecord
                            seq_record = SeqRecord(
                                Seq(sequence_data["sequence"]),
                                id=seq_name,
                                description=""
                            )
                            sequences.append(seq_record)
                
                # Write FASTA file if we have sequences
                if sequences:
                    fasta_file = os.path.join(directory, f"{partition_name}.fasta")
                    with open(fasta_file, "w") as f:
                        SeqIO.write(sequences, f, "fasta")
            
            QMessageBox.information(self, "Success", f"Partitions exported as FASTA files to {directory}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export data: {str(e)}")
    
    def download_ncbi(self):
        """Download sequences from NCBI"""
        # Get accession numbers from selected cells
        selected_items = self.seq_matrix.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "Warning", "Please select cells containing accession numbers")
            return
            
        # Filter out items from the first column (Sequence Name column)
        filtered_items = [item for item in selected_items if self.seq_matrix.column(item) > 0]
        
        # If no items left after filtering, show warning
        if not filtered_items:
            QMessageBox.warning(self, "Warning", "Please select cells containing accession numbers from partition columns (not Sequence Name column)")
            return
            
        # Clear previous accession items tracking
        self.accession_items.clear()
        
        accessions = []
        for item in filtered_items:
            if item.text().strip():
                accession = item.text().strip()
                accessions.append(accession)
                # Store item position for later update
                row = self.seq_matrix.row(item)
                col = self.seq_matrix.column(item)
                if accession not in self.accession_items:
                    self.accession_items[accession] = []
                self.accession_items[accession].append((row, col))
                
        if not accessions:
            QMessageBox.warning(self, "Warning", "No valid accession numbers found in selected cells")
            return
            
        # Start download thread
        self.download_thread = NCBIdownloader(accessions)
        self.download_thread.log_signal.connect(self.handle_download_log)
        self.download_thread.finished_signal.connect(self.handle_download_finished)
        self.download_thread.start()
        
        QMessageBox.information(self, "Info", f"Started downloading {len(accessions)} sequences from NCBI")
    
    def handle_download_log(self, message):
        """Handle log messages from download thread"""
        print(message)  # In a real application, you might want to show this in a log window
    
    def handle_download_finished(self, sequences):
        """Handle completion of NCBI download"""
        if sequences:
            # Add downloaded sequences to a special "Downloaded" group
            if "Downloaded" not in self.sequence_data:
                self.sequence_data["Downloaded"] = {}
                
            if "NCBI Downloads" not in self.sequence_data["Downloaded"]:
                self.sequence_data["Downloaded"]["NCBI Downloads"] = []
                
            # Process downloaded sequences
            for name, sequence in sequences.items():
                # Generate unique UUID for each sequence
                sequence_uuid = str(uuid.uuid4())
                
                # Store in both structures
                self.sequence_data["Downloaded"]["NCBI Downloads"].append({
                    "uuid": sequence_uuid,
                    "name": name,
                    "sequence": sequence
                })
                self.seqmatrix_data[sequence_uuid] = {
                    "name": name,
                    "sequence": sequence
                }
                
                # Update the corresponding table cells with the downloaded sequence info
                # Find cells that requested this accession
                for accession in self.accession_items:
                    # Match by accession (first part of the name)
                    if name.startswith(accession):
                        positions = self.accession_items[accession]
                        for row, col in positions:
                            if row < self.seq_matrix.rowCount() and col < self.seq_matrix.columnCount():
                                # Update the cell with sequence name and UUID
                                item = QTableWidgetItem(accession)
                                item.setData(Qt.UserRole, sequence_uuid)
                                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                                self.seq_matrix.setItem(row, col, item)
                
            QMessageBox.information(self, "Success", f"Downloaded {len(sequences)} sequences")
            self.update_tree_view()
            # 注意：这里不再调用 update_table()，因为它会清空整个表格
            self.update_highlighting()  # Update highlighting for newly downloaded sequences
        else:
            QMessageBox.warning(self, "Warning", "No sequences were downloaded")
    
    def on_header_double_clicked(self, logical_index):
        """Handle double click on any column header"""
        # Handle Sequence Name column (index 0)
        if logical_index == 0:
            self.on_sequence_name_header_double_clicked()
        # Handle Partition columns (index > 0)
        elif logical_index > 0:
            self.on_partition_header_double_clicked(logical_index)
    
    def on_sequence_name_header_double_clicked(self):
        """Handle double click on sequence name column header"""
        # Show dialog for sequence names
        dialog = SequenceNameEditDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            sequence_list, include_header = dialog.get_data()
            
            # Check if sequence list contains tabs
            if '\t' in sequence_list:
                QMessageBox.critical(self, "Error", "Sequence list cannot contain tabs. Please use one name per line.")
                return
                
            # If sequence list is not empty, warn user about overwriting
            if sequence_list.strip():
                reply = QMessageBox.question(
                    self, 
                    "Warning", 
                    "Previous changes in this column will be overwritten. Continue?",
                    QMessageBox.Yes | QMessageBox.No,
                    QMessageBox.No
                )
                
                if reply == QMessageBox.No:
                    return
                    
                # Process sequence list
                names = [line for line in sequence_list.split('\n')]  # Keep empty lines
                
                # If include_header is checked, remove the first line
                if include_header and names:
                    names = names[1:]
                    
                self.populate_sequence_name_column(names)
    
    def on_partition_header_double_clicked(self, logical_index):
        """Handle double click on partition column header"""
        # Get current header text
        current_header = self.seq_matrix.horizontalHeaderItem(logical_index)
        current_name = current_header.text() if current_header else f"Partition {logical_index}"
        
        # Show dialog
        dialog = PartitionEditDialog(current_name, self)
        if dialog.exec_() == QDialog.Accepted:
            new_name, sequence_list, include_header = dialog.get_data()
            
            # Check if sequence list contains tabs
            if '\t' in sequence_list:
                QMessageBox.critical(self, "Error", "Sequence list cannot contain tabs. Please use one accession per line.")
                return
                
            # If sequence list is not empty, warn user about overwriting
            if sequence_list.strip():
                reply = QMessageBox.question(
                    self, 
                    "Warning", 
                    "Previous changes in this column will be overwritten. Continue?",
                    QMessageBox.Yes | QMessageBox.No,
                    QMessageBox.No
                )
                
                if reply == QMessageBox.No:
                    return
                    
                # Process sequence list
                accessions = [line for line in sequence_list.split('\n')]  # Keep empty lines
                
                # If include_header is checked, remove the first line
                if include_header and accessions:
                    accessions = accessions[1:]
                    
                self.populate_column_with_accessions(logical_index, accessions)
            
            # Update header name
            self.seq_matrix.setHorizontalHeaderItem(logical_index, QTableWidgetItem(new_name))
    
    def populate_sequence_name_column(self, names):
        """Populate the sequence name column with names, preserving empty lines"""
        # Ensure we have enough rows
        while self.seq_matrix.rowCount() < len(names):
            self.seq_matrix.setRowCount(self.seq_matrix.rowCount() + 1)
            
        # Fill the sequence name column with names (including empty lines)
        for row, name in enumerate(names):
            item = QTableWidgetItem(name)  # This will be an empty string for empty lines
            item.setFlags(item.flags() | Qt.ItemIsEditable)  # Make editable
            self.seq_matrix.setItem(row, 0, item)
    
    def populate_column_with_accessions(self, column_index, accessions):
        """Populate a column with accession numbers, preserving empty lines for missing genes"""
        # Ensure we have enough rows
        while self.seq_matrix.rowCount() < len(accessions):
            self.seq_matrix.setRowCount(self.seq_matrix.rowCount() + 1)
            
        # Fill the column with accession numbers (including empty lines)
        for row, accession in enumerate(accessions):
            item = QTableWidgetItem(accession)  # This will be an empty string for empty lines
            item.setFlags(item.flags() | Qt.ItemIsEditable)  # Make editable
            self.seq_matrix.setItem(row, column_index, item)
    
    def update_tree_view(self):
        """Update the tree view with current sequence data"""
        # Clear the tree
        self.tree_explorer.clear()
        self.tree_items_by_uuid.clear()
        
        # Add items to tree
        for dir_path, files in self.sequence_data.items():
            # Create parent item for directory or special group
            if dir_path == "Downloaded":
                parent_item = QTreeWidgetItem(self.tree_explorer, ["Downloaded"])
            else:
                parent_item = QTreeWidgetItem(self.tree_explorer, [dir_path])
                
            # Add child items for each file in the directory
            for file_name, sequences in files.items():
                file_item = QTreeWidgetItem(parent_item, [file_name])
                
                # Add sequence names as children of the file
                for seq_data in sequences:
                    sequence_item = SequenceTreeWidgetItem(
                        seq_data["uuid"],
                        seq_data["name"],
                        seq_data["sequence"],
                        file_item,
                        [seq_data["name"]]
                    )
                    
                    # Keep track of tree items by UUID for highlighting
                    if seq_data["uuid"] not in self.tree_items_by_uuid:
                        self.tree_items_by_uuid[seq_data["uuid"]] = []
                    self.tree_items_by_uuid[seq_data["uuid"]].append(sequence_item)
        
        self.tree_explorer.expandAll()
    
    def update_table(self):
        """Update the table with current sequence data"""
        # 注意：这个方法现在只在初始化时调用，不会清空已有数据
        # 确保有足够的行数
        while self.seq_matrix.rowCount() < 4:
            self.seq_matrix.insertRow(self.seq_matrix.rowCount())
            
        # Update highlighting
        self.update_highlighting()
    
    def update_highlighting(self):
        """Update highlighting of sequences based on occurrence count and source"""
        # Reset all highlights first
        self._reset_highlights()
        
        # Count occurrences of each UUID in the table
        uuid_counts = {}
        uuid_sources = {}  # Track if UUID comes from NCBI downloads
        table_items_by_uuid = {}  # Store table items by UUID for highlighting
        
        # Check all table items
        for row in range(self.seq_matrix.rowCount()):
            for col in range(self.seq_matrix.columnCount()):
                item = self.seq_matrix.item(row, col)
                if item and item.data(Qt.UserRole):  # Has UUID data
                    uuid = item.data(Qt.UserRole)
                    if uuid not in uuid_counts:
                        uuid_counts[uuid] = 0
                        table_items_by_uuid[uuid] = []
                    
                    uuid_counts[uuid] += 1
                    table_items_by_uuid[uuid].append(item)
                    
                    # Check if this UUID is from NCBI downloads
                    if uuid not in uuid_sources:
                        uuid_sources[uuid] = self._is_ncbi_sequence(uuid)
        
        # Apply highlighting based on counts and sources
        for uuid, count in uuid_counts.items():
            # Determine color based on count and source
            if count == 1:
                # Single occurrence - green for normal, blue for NCBI
                color = QColor(173, 216, 230) if uuid_sources.get(uuid, False) else QColor(144, 238, 144)  # Blue or Green
            else:
                # Multiple occurrences - red for duplicates
                color = QColor(255, 182, 193)  # Light red
            
            # Highlight table items
            for item in table_items_by_uuid.get(uuid, []):
                item.setBackground(color)
            
            # Highlight tree items
            for tree_item in self.tree_items_by_uuid.get(uuid, []):
                tree_item.setBackground(0, color)
    
    def _reset_highlights(self):
        """Reset all highlights in both tree and table"""
        # Reset tree highlights
        for uuid, items in self.tree_items_by_uuid.items():
            for item in items:
                item.setBackground(0, QColor(255, 255, 255))  # White background
        
        # Reset table highlights
        for row in range(self.seq_matrix.rowCount()):
            for col in range(self.seq_matrix.columnCount()):
                item = self.seq_matrix.item(row, col)
                if item:
                    item.setBackground(QColor(255, 255, 255))  # White background
    
    def _is_ncbi_sequence(self, uuid):
        """Check if a sequence UUID comes from NCBI downloads"""
        # Check if the UUID exists in the NCBI downloads section
        if "Downloaded" in self.sequence_data:
            if "NCBI Downloads" in self.sequence_data["Downloaded"]:
                for seq_data in self.sequence_data["Downloaded"]["NCBI Downloads"]:
                    if seq_data.get("uuid") == uuid:
                        return True
        return False
    
    def Lg_MenuBar(self):
        pass
        
    def ShowTreeView(self):
        pass
        
    def Lg_Drag_and_Drop(self):
        pass
        
    def load_from_file(self):
        pass
    
    def parse_tree_view(self, tree_dict):
        paths = []
        for key, value in tree_dict.items():
            if isinstance(value, dict):
                for file_name, sequences in value.items():
                    for seq_data in sequences:
                        paths.append(f"{key}/{file_name}/{seq_data['name']}")
            elif isinstance(value, str):
                paths.append(key)  
        return paths
        

class SeqMatrixPluginEntry:
    def __init__(self):
        self.window = None
        
    def run(self, import_from=None):
        return Lg_SeqMatrix(import_from=import_from)

if __name__ == "__main__":
    import sys
    
    app = QApplication(sys.argv)
    window = Lg_SeqMatrix()
    window.show()
    sys.exit(app.exec_())

def parse_tree_view(tree_dict):
    paths = []
    for key, value in tree_dict.items():
        if isinstance(value, list):
            for item in value:
                for sub_path in parse_tree_view(item):
                    paths.append(f"{key}/{sub_path}")
        elif isinstance(value, str):
            paths.append(key)  
    return paths