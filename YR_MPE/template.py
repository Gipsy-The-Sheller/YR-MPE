from PyQt5.QtWidgets import (QWidget, QDialog, QVBoxLayout, QTabWidget, QHBoxLayout, 
                             QPushButton, QCheckBox, QPlainTextEdit, QProgressBar, 
                             QMessageBox, QFileDialog, QFrame, QLabel, QTextEdit, 
                             QComboBox, QSpinBox, QGroupBox, QFormLayout, QLineEdit, 
                             QSizePolicy, QScrollArea)
from PyQt5.QtCore import QThread, pyqtSignal, QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtGui import QFont, QTextCursor
import os
import json
import tempfile
import platform
import subprocess

class Wrapper(QWidget):
    def __init__(self, importfrom = None, importdata = None):
        super().__init__()
        self.plugin_path = os.path.dirname(os.path.abspath(__file__))
        if not self.config():
            self.show_config_guide()
    
        if importfrom == None:
            self.init_ui()
        elif importfrom == "YR_MPEA":
            pass

    def get_tool_name(self):
        # This function can be overridden by the subclass to provide the tool name
        return
    

    def get_input_format(self):
        # This function can be overridden by the subclass to provide the input type
        return {"type": "FASTA", "suffix": [".fas", ".fasta", ".fa", ".fna"]}
    
    def get_output_format(self):
        # This function can be overridden by the subclass to provide the output type
        return {"type": "FASTA", "suffix": ".fas"}


    def config(self):
        
        tool_name = self.get_tool_name()
        # load config.json and find the corresponding tool
        # then check if "path" is valid
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            if not os.path.exists(config_path):
                return False
            
            with open(config_path, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
            
            tool_config = None
            for tool in config_data:
                if tool.get("name").lower() == tool_name.lower():
                    tool_config = tool
                    break
            
            if not tool_config:
                return False
            
            muscle5_path = os.path.join(self.plugin_path, tool_config["path"].lstrip("/"))
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
        tool_name = self.get_tool_name()
        self.setWindowTitle(f"{tool_name} Wrapper")
        # self.setMinimumSize(1000, 700)
        
        main_layout = QVBoxLayout()
        self.setLayout(main_layout)
        
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
        self.tool_path = self.get_tool_path()
        self.is_running = False
        self.imported_files = []  # List of imported file paths
        self.file_tags = []  # List of file tag widgets
        self.html_reports = []  # List of HTML report files
        self.current_report_index = 0
        self.console_output = []  # List of console messages
        
    def setup_input_tab(self):
        # This function can be overridden by the subclass to provide the input type
        pass
        
    def setup_output_tab(self):
        # This function can be overridden by the subclass to provide the output type
        pass
        
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
        self.add_console_message(f"{self.get_tool_name()} Console initialized", "info")
        
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
        self.run_button = QPushButton(f"Run {self.get_tool_name()}")
        self.run_button.clicked.connect(self.run_command)
        control_layout.addWidget(self.run_button)
        
        # Abort button
        self.stop_button = QPushButton("Abort")
        self.stop_button.clicked.connect(self.stop_command)
        self.stop_button.setEnabled(False)
        control_layout.addWidget(self.stop_button)
        
        main_layout.addLayout(control_layout)
        
    def get_tool_path(self):
        """Get tool executable file path"""
        tool_name = self.get_tool_name()
        try:
            config_path = os.path.join(self.plugin_path, "config.json")
            with open(config_path, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
            
            for tool in config_data:
                if tool.get("name").lower() == tool_name.lower():
                    return os.path.join(self.plugin_path, tool["path"].lstrip("/"))
            return None
        except Exception as e:
            print(f"Get tool path failed: {e}\nPlease re-configure")
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
            
    def run_command(self):

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
        
        # Run alignment in a separate thread
        self.command_thread = CommandThread(
            self.tool_path, input_files, self.get_parameters(), self.imported_files
        )
        self.command_thread.progress.connect(self.update_progress)
        self.command_thread.finished.connect(self.command_finished)
        self.command_thread.error.connect(self.command_error)
        self.command_thread.console_output.connect(self.add_console_message)
        self.command_thread.start()
        
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
            
            # Otherwise, use text input to create a temporary file
            sequence_text = self.sequence_text.toPlainText().strip()
            if not sequence_text:
                QMessageBox.warning(self, "Warning", "Please input sequence text!")
                return None
                
            # Create temporary file
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix=self.get_output_format()['suffix'], delete=False)
            temp_file.write(sequence_text)
            temp_file.close()
            return [temp_file.name]
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Prepare input files failed: {e}")
            return None
            
    def get_parameters(self):
        # This function can be overridden by the subclass to provide the command parameters
        params = []
            
        return params
        
    def update_progress(self, message):
        self.progress_bar.setFormat(message)
        
    def command_finished(self, output_files, html_files):
        """Handle command completion"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # Add console message
        self.add_console_message(f"{self.get_tool_name()} completed successfully! Generated {len(html_files)} report(s)", "info")
        
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
            
        QMessageBox.information(self, "Completed", f"{self.get_tool_name()} completed!")
    
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
        
    def command_error(self, error_message):
        """Handle alignment error"""
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        # Add console message
        self.add_console_message(f"{self.get_tool_name()} failed: {error_message}", "error")
        
        QMessageBox.critical(self, "Error", f"{self.get_tool_name()} failed: {error_message}")
        
    def stop_command(self):
        """Stop command"""
        if hasattr(self, 'command_thread') and self.command_thread.isRunning():
            self.command_thread.terminate()
            self.command_thread.wait()
            
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        
        QMessageBox.information(self, "Stopped", f"{self.get_tool_name()} has been aborted.")

class CommandThread(QThread):
    """Command thread"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, list)  # output_files, html_files
    error = pyqtSignal(str)
    console_output = pyqtSignal(str, str)  # message, msg_type
    
    def __init__(self, tool_path, input_files, parameters, imported_files=None):
        super().__init__()
        self.tool_path = tool_path
        self.input_files = input_files if isinstance(input_files, list) else [input_files]
        self.parameters = parameters
        self.imported_files = imported_files or []
        
    def run(self):
        try:
            output_files = []
            html_files = []
            
            # Process each input file separately
            for i, input_file in enumerate(self.input_files):
                self.progress.emit(f"Processing file {i+1}/{len(self.input_files)}...")
                self.console_output.emit(f"Processing file {i+1}/{len(self.input_files)}: {os.path.basename(input_file)}", "info")
                
                # Create output files for this input
                output_file = tempfile.NamedTemporaryFile(suffix=self.get_output_format()['suffix'], delete=False).name
                html_file = tempfile.NamedTemporaryFile(suffix='.html', delete=False).name
                
                # Build command
                cmd = [
                    self.tool_path,
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
                    self.error.emit(f"{self.get_tool_name()} execution failed for file {i+1}: {result.stderr}")
                    return
                
                # Generate HTML from FASTA output
                self.console_output.emit(f"Generating HTML report for file {i+1}...", "info")
                self.generate_html_from_output(output_file, html_file)
                
                output_files.append(output_file)
                html_files.append(html_file)
            
            self.progress.emit("All alignments completed")
            self.finished.emit(output_files, html_files)
            
        except subprocess.TimeoutExpired:
            self.error.emit(f"{self.get_tool_name()} timed out")
        except Exception as e:
            self.error.emit(f"{self.get_tool_name()} exception: {str(e)}")
        finally:
            # Clean up temporary input files (if created from text input)
            for input_file in self.input_files:
                if input_file.startswith(tempfile.gettempdir()):
                    try:
                        os.unlink(input_file)
                    except:
                        pass
    
    def generate_html_from_output(self, output_file, html_file):
        """Generate HTML visualization from output with interleaved display and conservation highlighting"""
        try:
            with open(output_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Parse output sequences
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
    <title>{self.get_tool_name()} Result</title>
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
        <h1>🧬 {self.get_tool_name()} Result</h1>
        
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