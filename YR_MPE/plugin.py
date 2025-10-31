from PyQt5.QtWidgets import (QWidget, QLabel, QVBoxLayout, QDialog, QTabWidget, 
                             QTextEdit, QPushButton, QFileDialog, QComboBox, 
                             QSpinBox, QCheckBox, QHBoxLayout, QProgressBar,
                             QMessageBox, QGroupBox, QFormLayout, QLineEdit, QSizePolicy,
                             QScrollArea, QFrame, QPlainTextEdit, QDoubleSpinBox, QApplication, QToolButton)
from PyQt5.QtCore import QThread, pyqtSignal, Qt, QUrl, QTimer
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtGui import QFont, QTextCursor, QIcon
# from core.plugin_base import 
import json
import os
import subprocess
import tempfile
import platform
import threading
from .template import Wrapper, CommandThread
from .sequence_editor import SequenceEditor
import re
from Bio import SeqIO

class Muscle5_wrapper(QWidget):
    def __init__(self, importfrom = None, importdata = None):
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.destroyed.connect(self.on_destroyed)
        self.importfile = None
        self.temp_files = []

        if not self.config():
            self.show_config_guide()
    
        if importfrom == None:
            self.init_ui()
        elif importfrom == "YR_MPEA":
            # check if importdata is a list of sequences
            # format. list(SeqIO.parse(file, "fasta"))
            if isinstance(importdata, list):
                # create a temporary file
                temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False)
                self.temp_files.append(temp_file.name)
                with open(temp_file.name, 'w') as f:
                    for seq in importdata:
                        f.write(f">{seq.id}\n{seq.seq}\n")
                    # SeqIO.write(importdata, f, "fasta")
                temp_file.close()
                self.importfile = temp_file.name
                self.init_ui()
            
            else:
                self.init_ui()
    
    def on_destroyed(self):
        # delete the temporary file
        for temp_file in self.temp_files:
            os.remove(temp_file)

    def get_citation(self):
        return """Edgar, R.C. Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. <i>Nat Commun</i> <b>13</b>, 6968 (2022). DOI: <a href="https://doi.org/10.1038/s41467-022-34630-w">10.1038/s41467-022-34630-w</a>"""

    def copy_citation(self):
        QApplication.clipboard().setText(self.get_citation())
        self.add_console_message("Citation copied to clipboard", "info")

    def config(self):
        
        # load config.json and find "muscle5"
        # then check if "path" is valid
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            if not os.path.exists(config_path):
                return False
            
            with open(config_path, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
            
            muscle5_config = None
            for tool in config_data:
                if tool.get("name") == "Muscle5":
                    muscle5_config = tool
                    break
            
            if not muscle5_config:
                return False
            
            muscle5_path = os.path.join(self.plugin_path, muscle5_config["path"].lstrip("/"))
            return os.path.exists(muscle5_path)
            
        except Exception as e:
            print(f"Config check failed: {e}")
            return False
        
    def show_config_guide(self):
        config_guide = QDialog()
        config_guide.setWindowTitle("Config Guide")
        config_guide.setMinimumSize(800, 600)
        config_guide.setLayout(QVBoxLayout())
        config_guide.show()

    def init_ui(self):
        self.setWindowTitle("MUSCLE5 Aligner Wrapper")
        # self.setMinimumSize(1000, 700)
        
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)

        # options:
        # Tab 1. Input
        # 1. File input / Seq input (autofill with importdata)
        # 2. Parameters
        # mode (-align / -super5[for big dataset])
        # 3. run button and progress bar (isolated pyqtthread)

        # command: %muscle% {-align / -super5} input.fas -html %aln.html% (tempfile) -output aln.fas

        # Tab 2. Output preview (html), using QWebEngineView
        
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Tab 1: Input and parameters
        self.input_tab = QWidget()
        self.tab_widget.addTab(self.input_tab, "Input & parameters")
        self.setup_input_tab()
        
        # Tab 2: Output preview
        self.output_tab = QWidget()
        self.tab_widget.addTab(self.output_tab, "Output preview")
        self.setup_output_tab()
        
        # Tab 3: Console
        self.console_tab = QWidget()
        self.tab_widget.addTab(self.console_tab, "Console")
        self.setup_console_tab()
        
        # Progress bar and run button
        self.setup_control_panel(main_layout)
        
        # Initialize variables
        self.muscle5_path = self.get_muscle5_path()
        self.is_running = False
        self.imported_files = []  # List of imported file paths
        self.file_tags = []  # List of file tag widgets
        self.html_reports = []  # List of HTML report files
        self.current_report_index = 0
        self.console_output = []  # List of console messages
        
    def setup_input_tab(self):
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

        if self.importfile:
            # File Imported
            self.file_imported_button = QPushButton(f"File Already Imported from: {self.importfrom}")
            self.file_imported_button.setEnabled(False)
            layout.addWidget(self.file_imported_button)

            # set the file_path_edit to the importfile
            self.file_path_edit.setText(self.importfile)
            # invisible the input_group
            input_group.setVisible(False)
        
        # Parameters setting
        params_group = QGroupBox("Parameters")
        params_layout = QFormLayout()
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # Alignment mode
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Standard alignment (-align)", "Large dataset (-super5)"])
        params_layout.addRow("Alignment mode:", self.mode_combo)
        
        # Citation Group
        citation_group = QGroupBox("Citation")
        # citation_layout = QFormLayout()
        citation_layout = QVBoxLayout()
        citation_group.setLayout(citation_layout)
        citation_prompt = QLabel("If you use MUSCLE5 in your work, please cite:")
        citation_layout.addWidget(citation_prompt)
        citation_text = self.get_citation()
        # remove html tags
        citation_text = re.sub(r'<[^>]*>', '', citation_text)
        citation_text = QLabel(citation_text)
        # enable word wrap
        citation_text.setWordWrap(True)
        citation_text.setOpenExternalLinks(True)

        # copy button (icon only)
        copy_btn = QToolButton()
        copy_btn.setIcon(QIcon(os.path.join(self.plugin_path, "icons/copy.svg")))
        copy_btn.clicked.connect(self.copy_citation)

        citation_widget = QWidget()
        citation_widget_layout = QHBoxLayout()
        citation_widget.setLayout(citation_widget_layout)
        citation_widget_layout.addWidget(citation_text)
        citation_widget_layout.addWidget(copy_btn)
        citation_layout.addWidget(citation_widget)

        layout.addWidget(citation_group)
        
        layout.addStretch()
        
    def setup_output_tab(self):
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # Report selector
        report_layout = QHBoxLayout()
        report_layout.addWidget(QLabel("Select Report:"))
        
        self.report_combo = QComboBox()
        self.report_combo.setEnabled(False)
        self.report_combo.currentIndexChanged.connect(self.on_report_changed)
        report_layout.addWidget(self.report_combo)
        
        report_layout.addStretch()
        layout.addLayout(report_layout)
        
        # HTML preview
        self.web_view = QWebEngineView()
        self.web_view.setHtml("<p>Waiting for alignment result...</p>")
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(self.web_view)
        
        # Output file information
        self.output_info = QLabel("Output file will be displayed here")
        layout.addWidget(self.output_info)
        
    def setup_console_tab(self):
        """Setup console tab for command and output display"""
        layout = QVBoxLayout()
        self.console_tab.setLayout(layout)
        
        # Console controls
        console_controls = QHBoxLayout()
        
        clear_btn = QPushButton("Clear Console")
        clear_btn.clicked.connect(self.clear_console)
        console_controls.addWidget(clear_btn)
        
        auto_scroll_checkbox = QCheckBox("Auto Scroll")
        auto_scroll_checkbox.setChecked(True)
        auto_scroll_checkbox.toggled.connect(self.toggle_auto_scroll)
        console_controls.addWidget(auto_scroll_checkbox)
        
        console_controls.addStretch()
        layout.addLayout(console_controls)
        
        # Console text area
        self.console_text = QPlainTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setFont(QFont("Courier New", 10))
        # If the platform is Windows, use Consolas font
        if platform.system() == "Windows":
            self.console_text.setStyleSheet("""
                QPlainTextEdit {
                    background-color: #272822;
                    color: #f8f8f2;
                    border: 1px solid #3c3c3c;
                    font-family: 'Consolas', monospace;
                }
            """) # Monokai theme
        else:
            self.console_text.setStyleSheet("""
                QPlainTextEdit {
                    background-color: #272822;
                    color: #f8f8f2;
                    border: 1px solid #3c3c3c;
                    font-family: 'Courier New', monospace;
                }
            """) # Monokai theme            
        layout.addWidget(self.console_text)
        
        # Initialize console
        self.auto_scroll = True
        self.add_console_message("MUSCLE5 Console initialized", "info")
        
    def clear_console(self):
        """Clear console output"""
        self.console_text.clear()
        self.console_output.clear()
        self.add_console_message("Console cleared", "info")
        
    def toggle_auto_scroll(self, enabled):
        """Toggle auto scroll functionality"""
        self.auto_scroll = enabled
        
    def add_console_message(self, message, msg_type="info"):
        """Add message to console"""
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        
        # Color coding for different message types
        if msg_type == "command":
            formatted_msg = f"[{timestamp}] $ {message}"
            color = "#4CAF50"  # Green for commands
        elif msg_type == "output":
            formatted_msg = f"[{timestamp}] > {message}"
            color = "#2196F3"  # Blue for output
        elif msg_type == "error":
            formatted_msg = f"[{timestamp}] ERROR: {message}"
            color = "#F44336"  # Red for errors
        elif msg_type == "warning":
            formatted_msg = f"[{timestamp}] WARNING: {message}"
            color = "#FF9800"  # Orange for warnings
        else:  # info
            formatted_msg = f"[{timestamp}] {message}"
            color = "#9E9E9E"  # Gray for info
        
        # Add to console
        self.console_text.appendPlainText(formatted_msg)
        
        # Auto scroll to bottom if enabled
        if self.auto_scroll:
            self.console_text.moveCursor(QTextCursor.End)
        
        # Store in console output list
        if not hasattr(self, 'console_output'):
            self.console_output = []
        self.console_output.append({
            'timestamp': timestamp,
            'message': message,
            'type': msg_type,
            'formatted': formatted_msg
        })
        
    def setup_control_panel(self, main_layout):
        control_layout = QHBoxLayout()
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        control_layout.addWidget(self.progress_bar)
        
        # Run button
        self.run_button = QPushButton("Run alignment")
        self.run_button.clicked.connect(self.run_alignment)
        control_layout.addWidget(self.run_button)
        
        # Abort button
        self.stop_button = QPushButton("Abort")
        self.stop_button.clicked.connect(self.stop_alignment)
        self.stop_button.setEnabled(False)
        control_layout.addWidget(self.stop_button)
        
        main_layout.addLayout(control_layout)
        
    def get_muscle5_path(self):
        """Get muscle5 executable file path"""
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            with open(config_path, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
            
            for tool in config_data:
                if tool.get("name") == "Muscle5":
                    return os.path.join(self.plugin_path, tool["path"].lstrip("/"))
            return None
        except Exception as e:
            print(f"Get muscle5 path failed: {e}\nPlease re-configure")
            return None
            
    def browse_files(self):
        """Browse multiple files"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTA files", "", "FASTA file (*.fas *.fasta *.fa);;All files (*)"
        )
        if file_paths:
            for file_path in file_paths:
                if file_path not in self.imported_files:
                    self.add_file_tag(file_path)
            
            # Update file path display
            if len(self.imported_files) == 1:
                self.file_path_edit.setText(self.imported_files[0])
            else:
                self.file_path_edit.setText(f"{len(self.imported_files)} files selected")
            
            # Hide text input when files are imported
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)
    
    def add_file_tag(self, file_path):
        """Add a file tag widget"""
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
            
    def run_alignment(self):
        """Run MUSCLE5 alignment"""
        if self.is_running:
            return
            
        # Check input
        if (not self.imported_files and not self.sequence_text.toPlainText().strip()) and not self.importfile:
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
            return
            
        if not self.muscle5_path or not os.path.exists(self.muscle5_path):
            QMessageBox.critical(self, "Error", "MUSCLE5 executable file not found!")
            return
            
        # Add console message
        self.add_console_message("Starting MUSCLE5 alignment...", "info")
        
        # Prepare input files
        input_files = self.prepare_input_files()
        if not input_files:
            return
            
        # Update UI state
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Unknown progress
        
        # Run alignment in a separate thread
        self.alignment_thread = AlignmentThread(
            self.muscle5_path, input_files, self.get_parameters(), self.imported_files
        )
        self.alignment_thread.progress.connect(self.update_progress)
        self.alignment_thread.finished.connect(self.alignment_finished)
        self.alignment_thread.error.connect(self.alignment_error)
        self.alignment_thread.console_output.connect(self.add_console_message)
        self.alignment_thread.start()
        
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
            elif self.importfile:
                return [self.importfile]
            
            # Otherwise, use text input to create a temporary file
            sequence_text = self.sequence_text.toPlainText().strip()
            if not sequence_text and not self.importfile:
                QMessageBox.warning(self, "Warning", "Please input sequence text!")
                return None
                
            # Create temporary file
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False)
            temp_file.write(sequence_text)
            temp_file.close()
            return [temp_file.name]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
            
    def get_parameters(self):
        """Get alignment parameters"""
        params = []
        
        # Alignment mode
        if self.mode_combo.currentText().startswith("Large dataset"):
            params.append("-super5")
        else:
            params.append("-align")
            
        # # 最大迭代次数
        # params.extend(["-maxiters", str(self.max_iterations.value())])
        
        # # 其他选项
        # if not self.diags_checkbox.isChecked():
        #     params.append("-diags")
            
        # if not self.smooth_checkbox.isChecked():
        #     params.append("-smooth")
            
        return params
        
    def update_progress(self, message):
        self.progress_bar.setFormat(message)
        
    def alignment_finished(self, output_files, html_files):
        """Handle alignment completion"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # Add console message
        self.add_console_message(f"Alignment completed successfully! Generated {len(html_files)} report(s)", "info")
        
        # Store HTML reports
        self.html_reports = html_files if isinstance(html_files, list) else [html_files]
        self.current_report_index = 0
        
        # Update report combo box
        self.update_report_combo()
        
        # Switch to output tab
        self.tab_widget.setCurrentIndex(1)
        
        # Show first HTML result
        self.show_current_report()
            
        # Show output file information
        if output_files:
            if isinstance(output_files, list):
                self.output_info.setText(f"Output files: {len(output_files)} files generated")
            else:
                self.output_info.setText(f"Output file: {output_files}")
            
        QMessageBox.information(self, "Completed", "Sequence alignment completed!")
    
    def update_report_combo(self):
        """Update report combo box with available reports"""
        self.report_combo.clear()
        
        if not self.html_reports:
            self.report_combo.setEnabled(False)
            return
        
        # Add report options
        for i, html_file in enumerate(self.html_reports):
            if os.path.exists(html_file):
                report_name = f"Report {i+1}"
                if len(self.imported_files) > 1:
                    # If multiple files, show which file this report corresponds to
                    if i < len(self.imported_files):
                        file_name = os.path.basename(self.imported_files[i])
                        report_name = f"Report {i+1} ({file_name})"
                
                self.report_combo.addItem(report_name)
        
        self.report_combo.setEnabled(len(self.html_reports) > 1)
        self.report_combo.setCurrentIndex(self.current_report_index)
    
    def on_report_changed(self, index):
        """Handle report selection change"""
        if 0 <= index < len(self.html_reports):
            self.current_report_index = index
            self.show_current_report()
    
    def show_current_report(self):
        """Show the currently selected report"""
        if not self.html_reports or self.current_report_index >= len(self.html_reports):
            self.web_view.setHtml("<h2>No report available</h2>")
            return
        
        html_file = self.html_reports[self.current_report_index]
        if html_file and os.path.exists(html_file):
            # Convert Windows path to file URL format
            file_url = QUrl.fromLocalFile(html_file)
            self.web_view.load(file_url)
        else:
            self.web_view.setHtml("<h2>Report file not found</h2>")
        
    def alignment_error(self, error_message):
        """Handle alignment error"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # Add console message
        self.add_console_message(f"Alignment failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"Alignment failed: {error_message}")
        
    def stop_alignment(self):
        """Stop alignment"""
        if hasattr(self, 'alignment_thread') and self.alignment_thread.isRunning():
            self.alignment_thread.terminate()
            self.alignment_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", "Alignment has been aborted.")

class Muscle5_wrapper_entry:
    def run(self):
        return Muscle5_wrapper()

class AlignmentThread(QThread):
    """Alignment thread"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, list)  # output_files, html_files
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)  # message, msg_type
    
    def __init__(self, muscle5_path, input_files, parameters, imported_files=None):
        super().__init__()
        self.muscle5_path = muscle5_path
        self.input_files = input_files if isinstance(input_files, list) else [input_files]
        self.parameters = parameters
        self.imported_files = imported_files or []
    
    def get_tool_name(self):
        """Return the tool name for MUSCLE5"""
        return "MUSCLE5"
        
    def run(self):
        try:
            output_files = []
            html_files = []
            
            # Process each input file separately
            for i, input_file in enumerate(self.input_files):
                self.progress.emit(f"Processing file {i+1}/{len(self.input_files)}...")
                self.console_output.emit(f"Processing file {i+1}/{len(self.input_files)}: {os.path.basename(input_file)}", "info")
                
                # Create output files for this input
                output_file = tempfile.NamedTemporaryFile(suffix='.fas', delete=False).name
                html_file = tempfile.NamedTemporaryFile(suffix='.html', delete=False).name
                
                # Build command
                cmd = [
                    self.muscle5_path,
                    *self.parameters,
                    input_file,
                    "-output", output_file
                ]
                
                # Send command to console
                cmd_str = ' '.join(cmd)
                self.console_output.emit(cmd_str, "command")
                
                # Execute command
                result = subprocess.run(
                    cmd, 
                    capture_output=True, 
                    text=True, 
                    timeout=300
                )
                
                # Send results to console
                self.console_output.emit(f"Return code: {result.returncode}", "info")
                if result.stdout:
                    # Split stdout into lines and send each line
                    for line in result.stdout.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "output")
                if result.stderr:
                    # Split stderr into lines and send each line
                    for line in result.stderr.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "error")
                
                if result.returncode != 0:
                    self.error.emit(f"MUSCLE5 execution failed for file {i+1}: {result.stderr}")
                    return
                
                # Generate HTML from FASTA output
                self.console_output.emit(f"Generating HTML report for file {i+1}...", "info")
                self.generate_html_from_fasta(output_file, html_file)
                
                output_files.append(output_file)
                html_files.append(html_file)
            
            self.progress.emit("All alignments completed")
            self.finished.emit(output_files, html_files)
            
        except subprocess.TimeoutExpired:
            self.error.emit("Alignment timed out")
        except Exception as e:
            self.error.emit(f"Alignment exception: {str(e)}")
        finally:
            # Clean up temporary input files (if created from text input)
            for input_file in self.input_files:
                if input_file.startswith(tempfile.gettempdir()):
                    try:
                        os.unlink(input_file)
                    except:
                        pass
    
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


if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    
    # Create application
    app = QApplication(sys.argv)
    
    # Create MUSCLE5 wrapper instance
    muscle5_wrapper = Muscle5_wrapper()
    
    # Check configuration
    if not muscle5_wrapper.config():
        print("MUSCLE5 configuration check failed!")
        print("Please ensure the correct Muscle5 path configuration exists in config.json")
        sys.exit(1)
    
    # Show window
    muscle5_wrapper.show()
    
    # Run application
    sys.exit(app.exec_())


class MAFFT_wrapper(Wrapper):
    """MAFFT wrapper class inheriting from the generic Wrapper"""
    
    def __init__(self, importfrom=None, importdata=None):
        super().__init__(importfrom, importdata)
        
    def get_citation(self):
        return """Kazutaka Katoh, Kazuharu Misawa, Kei‐ichi Kuma, Takashi Miyata, MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform, <i>Nucleic Acids Research</i>, Volume 30, Issue 14, 15 July 2002, Pages 3059–3066. DOI: <a href="https://doi.org/10.1093/nar/gkf436">10.1093/nar/gkf436</a>"""
    
    def copy_citation(self):
        QApplication.clipboard().setText(self.get_citation())
        self.add_console_message("Citation copied to clipboard", "info")
        
    def get_tool_name(self):
        """Return the tool name for MAFFT"""
        return "MAFFT"
    
    def get_input_format(self):
        """Return input format for MAFFT"""
        return {"type": "FASTA", "suffix": [".fas", ".fasta", ".fa", ".fna"]}
    
    def get_output_format(self):
        """Return output format for MAFFT"""
        return {"type": "FASTA", "suffix": ".fas"}
    
    def setup_input_tab(self):
        """Setup input tab with MAFFT-specific parameters"""
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
        
        # Basic Parameters setting
        basic_params_group = QGroupBox("Basic Parameters")
        basic_params_layout = QFormLayout()
        basic_params_group.setLayout(basic_params_layout)
        layout.addWidget(basic_params_group)
        
        # Alignment mode
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Auto mode (--auto)",
            "High speed (default)",
            "High accuracy - localpair (--localpair)",
            "High accuracy - genafpair (--genafpair)", 
            "High accuracy - globalpair (--globalpair)",
            "Fast (--retree 1)"
        ])
        basic_params_layout.addRow("Alignment mode:", self.mode_combo)
        
        # Thread count
        self.thread_spinbox = QSpinBox()
        self.thread_spinbox.setRange(-1, 32)
        self.thread_spinbox.setValue(-1)  # Default: auto
        self.thread_spinbox.setToolTip("-1 for automatic thread detection")
        basic_params_layout.addRow("Threads:", self.thread_spinbox)
        
        # Advanced parameters button
        advanced_btn = QPushButton("Advanced Parameters...")
        advanced_btn.clicked.connect(self.show_advanced_dialog)
        basic_params_layout.addRow("", advanced_btn)
        
        # Citation Group
        citation_group = QGroupBox("Citation")
        citation_layout = QVBoxLayout()
        citation_group.setLayout(citation_layout)
        citation_prompt = QLabel("If you use MAFFT in your work, please cite:")
        citation_layout.addWidget(citation_prompt)
        citation_text = QLabel(self.get_citation())
        citation_text.setWordWrap(True)
        citation_text.setOpenExternalLinks(True)
        copy_btn = QToolButton()
        copy_btn.setIcon(QIcon(os.path.join(self.plugin_path, "icons/copy.svg")))
        copy_btn.clicked.connect(self.copy_citation)
        citation_widget = QWidget()
        citation_widget_layout = QHBoxLayout()
        citation_widget.setLayout(citation_widget_layout)
        citation_widget_layout.addWidget(citation_text)
        citation_widget_layout.addWidget(copy_btn)
        citation_layout.addWidget(citation_widget)
        layout.addWidget(citation_group)
        
        layout.addStretch()
        
    def setup_output_tab(self):
        """Setup output tab with report selector and HTML preview"""
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        # Report selector
        report_layout = QHBoxLayout()
        report_layout.addWidget(QLabel("Select Report:"))
        
        self.report_combo = QComboBox()
        self.report_combo.setEnabled(False)
        self.report_combo.currentIndexChanged.connect(self.on_report_changed)
        report_layout.addWidget(self.report_combo)
        
        report_layout.addStretch()
        layout.addLayout(report_layout)
        
        # HTML preview
        self.web_view = QWebEngineView()
        self.web_view.setHtml("<p>Waiting for alignment result...</p>")
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(self.web_view)
        
        # Output file information
        self.output_info = QLabel("Output file will be displayed here")
        layout.addWidget(self.output_info)
        
    def show_advanced_dialog(self):
        """Show advanced parameters dialog"""
        dialog = MAFFTAdvancedDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            # Store advanced parameters
            self.advanced_params = dialog.get_parameters()
    
    def get_parameters(self):
        """Get MAFFT command parameters"""
        params = []
        
        # Basic parameters based on mode selection
        mode_text = self.mode_combo.currentText()
        if "Auto mode" in mode_text:
            params.append("--auto")
        elif "High accuracy - localpair" in mode_text:
            params.extend(["--maxiterate", "1000", "--localpair"])
        elif "High accuracy - genafpair" in mode_text:
            params.extend(["--maxiterate", "1000", "--genafpair"])
        elif "High accuracy - globalpair" in mode_text:
            params.extend(["--maxiterate", "1000", "--globalpair"])
        elif "Fast" in mode_text:
            params.append("--retree")
            params.append("1")
        # High speed (default) - no additional parameters
        
        # Thread count
        thread_count = self.thread_spinbox.value()
        if thread_count != -1:
            params.extend(["--thread", str(thread_count)])
        
        # Advanced parameters
        if hasattr(self, 'advanced_params'):
            params.extend(self.advanced_params)
        
        return params
    
    def run_command(self):
        """Override run_command to use MAFFTCommandThread"""
        if self.is_running:
            return
            
        # Check input
        if not self.imported_files and not self.sequence_text.toPlainText().strip():
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
            return
            
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", f"{self.get_tool_name()} executable file not found!")
            return
            
        # Add console message
        self.add_console_message(f"Starting {self.get_tool_name()} command...", "info")
        
        # Prepare input files
        input_files = self.prepare_input_files()
        if not input_files:
            return
            
        # Update UI state
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Unknown progress
        
        # Run alignment in a separate thread using MAFFTCommandThread
        self.command_thread = MAFFTCommandThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
        )
        self.command_thread.progress.connect(self.update_progress)
        self.command_thread.finished.connect(self.command_finished)
        self.command_thread.error.connect(self.command_error)
        self.command_thread.console_output.connect(self.add_console_message)
        self.command_thread.start()


class MAFFTAdvancedDialog(QDialog):
    """Advanced parameters dialog for MAFFT"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MAFFT Advanced Parameters")
        self.setModal(True)
        self.setMinimumSize(500, 400)
        
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # Create scroll area for parameters
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout()
        scroll_widget.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # Gap penalties
        gap_group = QGroupBox("Gap Penalties")
        gap_layout = QFormLayout()
        gap_group.setLayout(gap_layout)
        scroll_layout.addWidget(gap_group)
        
        self.op_spinbox = QDoubleSpinBox()
        self.op_spinbox.setRange(0.0, 10.0)
        self.op_spinbox.setSingleStep(0.1)
        self.op_spinbox.setValue(1.53)
        self.op_spinbox.setDecimals(2)
        self.op_spinbox.setToolTip("Gap opening penalty (default: 1.53)")
        gap_layout.addRow("Gap opening penalty (--op):", self.op_spinbox)
        
        self.ep_spinbox = QDoubleSpinBox()
        self.ep_spinbox.setRange(0.0, 5.0)
        self.ep_spinbox.setSingleStep(0.1)
        self.ep_spinbox.setValue(0.0)
        self.ep_spinbox.setDecimals(2)
        self.ep_spinbox.setToolTip("Offset/gap extension penalty (default: 0.0)")
        gap_layout.addRow("Gap extension penalty (--ep):", self.ep_spinbox)
        
        # Iteration settings
        iter_group = QGroupBox("Iteration Settings")
        iter_layout = QFormLayout()
        iter_group.setLayout(iter_layout)
        scroll_layout.addWidget(iter_group)
        
        self.maxiterate_spinbox = QSpinBox()
        self.maxiterate_spinbox.setRange(0, 10000)
        self.maxiterate_spinbox.setValue(0)
        self.maxiterate_spinbox.setToolTip("Maximum number of iterative refinement (default: 0)")
        iter_layout.addRow("Max iterations (--maxiterate):", self.maxiterate_spinbox)
        
        # Output options
        output_group = QGroupBox("Output Options")
        output_layout = QFormLayout()
        output_group.setLayout(output_layout)
        scroll_layout.addWidget(output_group)
        
        self.clustalout_checkbox = QCheckBox("Output in Clustal format")
        self.clustalout_checkbox.setToolTip("Output alignment in Clustal format (--clustalout)")
        output_layout.addRow("", self.clustalout_checkbox)
        
        self.reorder_checkbox = QCheckBox("Reorder sequences")
        self.reorder_checkbox.setToolTip("Reorder sequences in aligned order (--reorder)")
        output_layout.addRow("", self.reorder_checkbox)
        
        self.quiet_checkbox = QCheckBox("Quiet mode")
        self.quiet_checkbox.setToolTip("Do not report progress (--quiet)")
        output_layout.addRow("", self.quiet_checkbox)
        
        # Special options
        special_group = QGroupBox("Special Options")
        special_layout = QFormLayout()
        special_group.setLayout(special_layout)
        scroll_layout.addWidget(special_group)
        
        self.dash_checkbox = QCheckBox("Add structural information")
        self.dash_checkbox.setToolTip("Add structural information (--dash)")
        special_layout.addRow("", self.dash_checkbox)
        
        scroll_layout.addStretch()
        
        # Dialog buttons
        button_layout = QHBoxLayout()
        
        reset_btn = QPushButton("Reset to Defaults")
        reset_btn.clicked.connect(self.reset_to_defaults)
        button_layout.addWidget(reset_btn)
        
        button_layout.addStretch()
        
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)
        
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self.accept)
        button_layout.addWidget(ok_btn)
        
        layout.addLayout(button_layout)
    
    def reset_to_defaults(self):
        """Reset all parameters to default values"""
        self.op_spinbox.setValue(1.53)
        self.ep_spinbox.setValue(0.0)
        self.maxiterate_spinbox.setValue(0)
        self.clustalout_checkbox.setChecked(False)
        self.reorder_checkbox.setChecked(False)
        self.quiet_checkbox.setChecked(False)
        self.dash_checkbox.setChecked(False)
    
    def get_parameters(self):
        """Get the advanced parameters as a list"""
        params = []
        
        # Gap penalties
        if self.op_spinbox.value() != 1.53:  # Only add if not default
            params.extend(["--op", str(self.op_spinbox.value())])
        
        if self.ep_spinbox.value() != 0.0:  # Only add if not default
            params.extend(["--ep", str(self.ep_spinbox.value())])
        
        # Iteration settings
        if self.maxiterate_spinbox.value() != 0:  # Only add if not default
            params.extend(["--maxiterate", str(self.maxiterate_spinbox.value())])
        
        # Output options
        if self.clustalout_checkbox.isChecked():
            params.append("--clustalout")
        
        if self.reorder_checkbox.isChecked():
            params.append("--reorder")
        
        if self.quiet_checkbox.isChecked():
            params.append("--quiet")
        
        # Special options
        if self.dash_checkbox.isChecked():
            params.append("--dash")
        
        return params


class MAFFTCommandThread(CommandThread):
    """MAFFT-specific command thread"""
    
    def get_tool_name(self):
        """Return the tool name for MAFFT"""
        return "MAFFT"
    
    def run(self):
        try:
            output_files = []
            html_files = []
            
            # Process each input file separately
            for i, input_file in enumerate(self.input_files):
                self.progress.emit(f"Processing file {i+1}/{len(self.input_files)}...")
                self.console_output.emit(f"Processing file {i+1}/{len(self.input_files)}: {os.path.basename(input_file)}", "info")
                
                # Create output files for this input
                output_file = tempfile.NamedTemporaryFile(suffix='.fas', delete=False).name
                html_file = tempfile.NamedTemporaryFile(suffix='.html', delete=False).name
                
                # Build MAFFT command (MAFFT outputs to stdout, not -output)
                cmd = [
                    self.tool_path,
                    *self.parameters,
                    input_file
                ]
                
                # Send command to console
                cmd_str = ' '.join(cmd)
                self.console_output.emit(cmd_str, "command")
                
                # Execute command and capture output
                result = subprocess.run(
                    cmd, 
                    capture_output=True, 
                    text=True, 
                    timeout=300
                )
                
                # Send results to console
                self.console_output.emit(f"Return code: {result.returncode}", "info")
                if result.stdout:
                    # Split stdout into lines and send each line
                    for line in result.stdout.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "output")
                if result.stderr:
                    # Split stderr into lines and send each line
                    for line in result.stderr.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "error")
                
                if result.returncode != 0:
                    self.error.emit(f"MAFFT execution failed for file {i+1}: {result.stderr}")
                    return
                
                # Write MAFFT output to file
                with open(output_file, 'w', encoding='utf-8') as f:
                    f.write(result.stdout)
                
                # Generate HTML from FASTA output
                self.console_output.emit(f"Generating HTML report for file {i+1}...", "info")
                self.generate_html_from_fasta(output_file, html_file)
                
                output_files.append(output_file)
                html_files.append(html_file)
            
            self.progress.emit("All alignments completed")
            self.finished.emit(output_files, html_files)
            
        except subprocess.TimeoutExpired:
            self.error.emit("MAFFT timed out")
        except Exception as e:
            self.error.emit(f"MAFFT exception: {str(e)}")
        finally:
            # Clean up temporary input files (if created from text input)
            for input_file in self.input_files:
                if input_file.startswith(tempfile.gettempdir()):
                    try:
                        os.unlink(input_file)
                    except:
                        pass
    
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
    <title>MAFFT Alignment Result</title>
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
            border-left: 4px solid #28a745;
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
        <h1>🧬 MAFFT Alignment Result</h1>
        
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


class MAFFT_wrapper_entry:
    def run(self):
        return MAFFT_wrapper()


class ClustalOmega_wrapper(Wrapper):
    """Clustal Omega wrapper using the generic Wrapper base."""
    
    def __init__(self, importfrom=None, importdata=None):
        super().__init__(importfrom, importdata)
        
    def get_citation(self):
        return """Sievers F, Wilm A, Dineen D, et al. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. <i>Mol Syst Biol.</i> 2011;7:539. Published 2011 Oct 11. DOI: <a href=\"https://doi.org/10.1038/msb.2011.75\">10.1038/msb.2011.75</a>"""
    
    def copy_citation(self):
        QApplication.clipboard().setText(self.get_citation())
        self.add_console_message("Citation copied to clipboard", "info")
        
    def get_tool_name(self):
        return "Clustal Omega"
    
    def get_input_format(self):
        return {"type": "FASTA", "suffix": [".fas", ".fasta", ".fa", ".fna", ".faa"]}
    
    def get_output_format(self):
        return {"type": "FASTA", "suffix": ".fas"}
    
    def setup_input_tab(self):
        layout = QVBoxLayout()
        self.input_tab.setLayout(layout)
        
        # Input group
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
        
        # Basic params
        basic_group = QGroupBox("Basic Parameters")
        basic_layout = QFormLayout()
        basic_group.setLayout(basic_layout)
        layout.addWidget(basic_group)
        
        # Auto mode
        self.auto_checkbox = QCheckBox("Auto options (--auto)")
        self.auto_checkbox.setChecked(True)
        basic_layout.addRow("", self.auto_checkbox)
        
        # Sequence type
        self.seqtype_combo = QComboBox()
        self.seqtype_combo.addItems(["Auto", "Protein", "RNA", "DNA"])
        basic_layout.addRow("Sequence type:", self.seqtype_combo)
        
        # Threads
        self.thread_spinbox = QSpinBox()
        self.thread_spinbox.setRange(1, 64)
        self.thread_spinbox.setValue(4)
        basic_layout.addRow("Threads:", self.thread_spinbox)
        
        # Output order
        self.order_combo = QComboBox()
        self.order_combo.addItems(["input-order", "tree-order"]) 
        basic_layout.addRow("Output order:", self.order_combo)
        
        # Advanced button
        adv_btn = QPushButton("Advanced Parameters...")
        adv_btn.clicked.connect(self.show_advanced_dialog)
        basic_layout.addRow("", adv_btn)
        
        # Citation Group
        citation_group = QGroupBox("Citation")
        citation_layout = QVBoxLayout()
        citation_group.setLayout(citation_layout)
        citation_prompt = QLabel("If you use Clustal Omega in your work, please cite:")
        citation_layout.addWidget(citation_prompt)
        citation_text = QLabel(self.get_citation())
        citation_text.setWordWrap(True)
        citation_text.setOpenExternalLinks(True)
        copy_btn = QToolButton()
        copy_btn.setIcon(QIcon(os.path.join(self.plugin_path, "icons/copy.svg")))
        copy_btn.clicked.connect(self.copy_citation)
        citation_widget = QWidget()
        citation_widget_layout = QHBoxLayout()
        citation_widget.setLayout(citation_widget_layout)
        citation_widget_layout.addWidget(citation_text)
        citation_widget_layout.addWidget(copy_btn)
        citation_layout.addWidget(citation_widget)
        layout.addWidget(citation_group)
        
        layout.addStretch()
    
    def setup_output_tab(self):
        layout = QVBoxLayout()
        self.output_tab.setLayout(layout)
        
        report_layout = QHBoxLayout()
        report_layout.addWidget(QLabel("Select Report:"))
        
        self.report_combo = QComboBox()
        self.report_combo.setEnabled(False)
        self.report_combo.currentIndexChanged.connect(self.on_report_changed)
        report_layout.addWidget(self.report_combo)
        report_layout.addStretch()
        layout.addLayout(report_layout)
        
        self.web_view = QWebEngineView()
        self.web_view.setHtml("<p>Waiting for alignment result...</p>")
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout.addWidget(self.web_view)
        
        self.output_info = QLabel("Output file will be displayed here")
        layout.addWidget(self.output_info)
    
    def show_advanced_dialog(self):
        dialog = ClustalAdvancedDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            self.advanced_params = dialog.get_parameters()
    
    def get_parameters(self):
        params = []
        
        # Auto
        if self.auto_checkbox.isChecked():
            params.append("--auto")
        
        # Seqtype
        seqtype = self.seqtype_combo.currentText()
        if seqtype != "Auto":
            params.extend(["--seqtype", seqtype])
        
        # Threads
        params.extend(["--threads", str(self.thread_spinbox.value())])
        
        # Output order
        params.extend(["--output-order", self.order_combo.currentText()])
        
        # Advanced
        if hasattr(self, 'advanced_params'):
            params.extend(self.advanced_params)
        
        return params
    
    def run_command(self):
        if self.is_running:
            return
        
        if not self.imported_files and not self.sequence_text.toPlainText().strip():
            QMessageBox.warning(self, "Warning", "Please provide sequence files or sequence text!")
            return
        
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(self, "Error", f"{self.get_tool_name()} executable file not found!")
            return
        
        self.add_console_message(f"Starting {self.get_tool_name()} command...", "info")
        
        input_files = self.prepare_input_files()
        if not input_files:
            return
        
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        self.command_thread = ClustalOmegaCommandThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
        )
        self.command_thread.progress.connect(self.update_progress)
        self.command_thread.finished.connect(self.command_finished)
        self.command_thread.error.connect(self.command_error)
        self.command_thread.console_output.connect(self.add_console_message)
        self.command_thread.start()


class ClustalAdvancedDialog(QDialog):
    """Advanced parameters dialog for Clustal Omega"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Clustal Omega Advanced Parameters")
        self.setModal(True)
        self.setMinimumSize(520, 460)
        
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        scroll = QScrollArea()
        content = QWidget()
        form = QFormLayout()
        content.setLayout(form)
        scroll.setWidget(content)
        scroll.setWidgetResizable(True)
        layout.addWidget(scroll)
        
        # Formats
        self.outfmt_combo = QComboBox()
        self.outfmt_combo.addItems(["fasta", "clustal", "msf", "phylip", "selex", "stockholm", "vienna"])
        form.addRow("Output format (--outfmt):", self.outfmt_combo)
        
        # Wrap
        self.wrap_spin = QSpinBox()
        self.wrap_spin.setRange(0, 1000)
        self.wrap_spin.setValue(0)
        form.addRow("Wrap (--wrap):", self.wrap_spin)
        
        # Iterations
        self.iter_spin = QSpinBox()
        self.iter_spin.setRange(0, 10000)
        self.iter_spin.setValue(0)
        form.addRow("Iterations (--iterations):", self.iter_spin)
        
        self.max_gt_iter = QSpinBox()
        self.max_gt_iter.setRange(0, 10000)
        form.addRow("Max guide-tree iter (--max-guidetree-iterations):", self.max_gt_iter)
        
        self.max_hmm_iter = QSpinBox()
        self.max_hmm_iter.setRange(0, 10000)
        form.addRow("Max HMM iter (--max-hmm-iterations):", self.max_hmm_iter)
        
        # Clustering toggles
        self.full_checkbox = QCheckBox("Use full distance matrix (--full)")
        form.addRow("", self.full_checkbox)
        self.full_iter_checkbox = QCheckBox("Full matrix during iteration (--full-iter)")
        form.addRow("", self.full_iter_checkbox)
        self.pileup_checkbox = QCheckBox("Sequentially align (--pileup)")
        form.addRow("", self.pileup_checkbox)
        
        # Distance options
        self.use_kimura = QCheckBox("Use Kimura correction (--use-kimura)")
        form.addRow("", self.use_kimura)
        self.percent_id = QCheckBox("Report percent identity (--percent-id)")
        form.addRow("", self.percent_id)
        
        # Limits
        self.maxnumseq = QSpinBox()
        self.maxnumseq.setRange(0, 10000000)
        form.addRow("Max sequences (--maxnumseq):", self.maxnumseq)
        
        self.maxseqlen = QSpinBox()
        self.maxseqlen.setRange(0, 10000000)
        form.addRow("Max sequence length (--maxseqlen):", self.maxseqlen)
        
        # Output tweaks
        self.resno_checkbox = QCheckBox("Clustal format residue numbers (--resno)")
        form.addRow("", self.resno_checkbox)
        
        # Buttons
        buttons = QHBoxLayout()
        reset_btn = QPushButton("Reset")
        reset_btn.clicked.connect(self.reset_to_defaults)
        buttons.addWidget(reset_btn)
        buttons.addStretch()
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        buttons.addWidget(cancel_btn)
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(self.accept)
        buttons.addWidget(ok_btn)
        layout.addLayout(buttons)
    
    def reset_to_defaults(self):
        self.outfmt_combo.setCurrentText("fasta")
        self.wrap_spin.setValue(0)
        self.iter_spin.setValue(0)
        self.max_gt_iter.setValue(0)
        self.max_hmm_iter.setValue(0)
        self.full_checkbox.setChecked(False)
        self.full_iter_checkbox.setChecked(False)
        self.pileup_checkbox.setChecked(False)
        self.use_kimura.setChecked(False)
        self.percent_id.setChecked(False)
        self.maxnumseq.setValue(0)
        self.maxseqlen.setValue(0)
        self.resno_checkbox.setChecked(False)
    
    def get_parameters(self):
        params = []
        # Output format
        fmt = self.outfmt_combo.currentText()
        if fmt and fmt != "fasta":
            params.extend(["--outfmt", fmt])
        
        # Wrap
        if self.wrap_spin.value() > 0:
            params.extend(["--wrap", str(self.wrap_spin.value())])
        
        # Iterations
        if self.iter_spin.value() > 0:
            params.extend(["--iterations", str(self.iter_spin.value())])
        if self.max_gt_iter.value() > 0:
            params.extend(["--max-guidetree-iterations", str(self.max_gt_iter.value())])
        if self.max_hmm_iter.value() > 0:
            params.extend(["--max-hmm-iterations", str(self.max_hmm_iter.value())])
        
        # Clustering toggles
        if self.full_checkbox.isChecked():
            params.append("--full")
        if self.full_iter_checkbox.isChecked():
            params.append("--full-iter")
        if self.pileup_checkbox.isChecked():
            params.append("--pileup")
        
        # Distance options
        if self.use_kimura.isChecked():
            params.append("--use-kimura")
        if self.percent_id.isChecked():
            params.append("--percent-id")
        
        # Limits
        if self.maxnumseq.value() > 0:
            params.extend(["--maxnumseq", str(self.maxnumseq.value())])
        if self.maxseqlen.value() > 0:
            params.extend(["--maxseqlen", str(self.maxseqlen.value())])
        
        # Output tweaks
        if self.resno_checkbox.isChecked():
            params.append("--resno")
        
        return params


class ClustalOmegaCommandThread(CommandThread):
    """Clustal Omega execution thread; capture stdout to FASTA, then render HTML."""
    
    def get_tool_name(self):
        """Return the tool name for Clustal Omega"""
        return "Clustal Omega"
    
    def run(self):
        try:
            output_files = []
            html_files = []
            for i, input_file in enumerate(self.input_files):
                self.progress.emit(f"Processing file {i+1}/{len(self.input_files)}...")
                self.console_output.emit(f"Processing file {i+1}/{len(self.input_files)}: {os.path.basename(input_file)}", "info")
                
                output_file = tempfile.NamedTemporaryFile(suffix='.fas', delete=False).name
                html_file = tempfile.NamedTemporaryFile(suffix='.html', delete=False).name
                
                cmd = [
                    self.tool_path,
                    *self.parameters,
                    "-i", input_file,
                ]
                
                # Prefer stdout to capture; could also use -o output_file
                cmd_str = ' '.join(cmd)
                self.console_output.emit(cmd_str, "command")
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                self.console_output.emit(f"Return code: {result.returncode}", "info")
                if result.stdout:
                    for line in result.stdout.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "output")
                if result.stderr:
                    for line in result.stderr.strip().split('\n'):
                        if line.strip():
                            self.console_output.emit(line.strip(), "error")
                
                if result.returncode != 0:
                    self.error.emit(f"Clustal Omega execution failed for file {i+1}: {result.stderr}")
                    return
                
                # Write stdout to output file
                with open(output_file, 'w', encoding='utf-8') as f:
                    f.write(result.stdout)
                
                # Render HTML (reuse base implementation)
                self.console_output.emit(f"Generating HTML report for file {i+1}...", "info")
                self.generate_html_from_output(output_file, html_file)
                output_files.append(output_file)
                html_files.append(html_file)
            
            self.progress.emit("All alignments completed")
            self.finished.emit(output_files, html_files)
        except subprocess.TimeoutExpired:
            self.error.emit("Clustal Omega timed out")
        except Exception as e:
            self.error.emit(f"Clustal Omega exception: {str(e)}")
        finally:
            for input_file in self.input_files:
                if input_file.startswith(tempfile.gettempdir()):
                    try:
                        os.unlink(input_file)
                    except:
                        pass

class ClustalOmega_wrapper_entry:
    def run(self):
        return ClustalOmega_wrapper()




class SequenceEditor_wrapper(QWidget):
    """序列编辑器包装器"""
    
    def __init__(self, importfrom=None, importdata=None):
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        self.sequence_editor = None
        self.init_ui()
        
    def init_ui(self):
        """初始化用户界面"""
        self.setWindowTitle("YR-MPE Sequence Editor")
        self.setMinimumSize(1200, 800)
        
        # 创建主布局
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
        # 创建工具栏
        toolbar = QHBoxLayout()
        
        # 新建按钮
        new_btn = QPushButton("New Sequence")
        new_btn.clicked.connect(self.new_sequence)
        toolbar.addWidget(new_btn)
        
        # 打开按钮
        open_btn = QPushButton("Open File")
        open_btn.clicked.connect(self.open_file)
        toolbar.addWidget(open_btn)
        
        # 保存按钮
        save_btn = QPushButton("Save")
        save_btn.clicked.connect(self.save_file)
        toolbar.addWidget(save_btn)
        
        toolbar.addStretch()
        
        # 帮助按钮
        help_btn = QPushButton("Help")
        help_btn.clicked.connect(self.show_help)
        toolbar.addWidget(help_btn)
        
        main_layout.addLayout(toolbar)
        
        # 创建序列编辑器
        self.sequence_editor = SequenceEditor(self)
        main_layout.addWidget(self.sequence_editor)
        
        # 连接信号
        self.sequence_editor.file_opened.connect(self.on_file_opened)
        self.sequence_editor.file_saved.connect(self.on_file_saved)
        
    def new_sequence(self):
        """新建序列"""
        self.sequence_editor.new_file()
        
    def open_file(self):
        """打开文件"""
        self.sequence_editor.open_file()
        
    def save_file(self):
        """保存文件"""
        self.sequence_editor.save_file()
        
    def show_help(self):
        """显示帮助"""
        help_text = """
YR-MPE Sequence Editor Help

Features:
• Support for DNA, RNA, and protein sequences
• Syntax highlighting for different sequence types
• Undo/Redo functionality
• Find and replace operations
• Sequence statistics and quality analysis
• Multiple file format support (FASTA, Phylip, Nexus, GenBank)
• Sequence manipulation (reverse, complement, translate)
• Open Reading Frame (ORF) detection
• Consensus sequence generation

Keyboard Shortcuts:
• Ctrl+N: New file
• Ctrl+O: Open file
• Ctrl+S: Save file
• Ctrl+Z: Undo
• Ctrl+Y: Redo
• Ctrl+F: Find
• Ctrl+H: Replace

For more information, visit the YR-MPE documentation.
        """
        
        QMessageBox.information(self, "Help", help_text)
        
    def on_file_opened(self, file_path):
        """文件打开事件"""
        self.setWindowTitle(f"YR-MPE Sequence Editor - {os.path.basename(file_path)}")
        
    def on_file_saved(self, file_path):
        """文件保存事件"""
        self.setWindowTitle(f"YR-MPE Sequence Editor - {os.path.basename(file_path)}")


class SequenceEditor_wrapper_entry:
    def run(self):
        return SequenceEditor_wrapper()