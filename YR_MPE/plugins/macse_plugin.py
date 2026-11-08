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

"""MACSE2 codon-aware aligner, standard ResultSpec plugin."""

from __future__ import annotations

import os
import shlex
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QFileDialog, QFormLayout, QFrame,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton,
    QRadioButton, QScrollArea, QSpinBox, QTextEdit, QVBoxLayout, QWidget,
)

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.alignment_common import (
    prepare_import, read_alignment, records_temp, setup_preview, show_preview,
)

GENETIC_CODES = [
    (1, "1 Standard"),
    (2, "2 Vertebrate Mitochondrial"),
    (3, "3 Yeast Mitochondrial"),
    (4, "4 Mold/Protozoan/Coelenterate Mitochondrial"),
    (5, "5 Invertebrate Mitochondrial"),
    (6, "6 Ciliate/Dasycladacean/Hexamita Nuclear"),
    (9, "9 Echinoderm/Flatworm Mitochondrial"),
    (10, "10 Euplotid Nuclear"),
    (11, "11 Bacterial/Archaeal/Plant Plastid"),
    (12, "12 Alternative Yeast Nuclear"),
    (13, "13 Ascidian Mitochondrial"),
    (14, "14 Alternative Flatworm Mitochondrial"),
    (15, "15 Blepharisma Nuclear"),
    (16, "16 Chlorophycean Mitochondrial"),
    (21, "21 Trematode Mitochondrial"),
    (22, "22 Scenedesmus obliquus Mitochondrial"),
    (23, "23 Thraustochytrium Mitochondrial"),
    (24, "24 Rhabdopleuridae Mitochondrial"),
    (25, "25 Candidate Division SR1/Gracilibacteria"),
    (26, "26 Pachysolen tannophilus Nuclear"),
    (27, "27 Karyorelict Nuclear"),
    (28, "28 Condylostoma Nuclear"),
    (29, "29 Mesodinium Nuclear"),
    (30, "30 Peritrich Nuclear"),
    (31, "31 Blastocrithidia Nuclear"),
    (32, "32 Seleno Protein"),
    (33, "33 Cephalodiscidae Mitochondrial UAA-Tyr"),
]

AA_ALPHABETS = [
    "SE_B_14", "SE_B_10", "SE_V_10", "Li_A_10", "Li_B_10",
    "Solis_D_10", "Solis_G_10", "Murphy_10", "SE_B_8", "SE_B_6", "Dayhoff_6",
]


class MACSEThread(BaseProcessThread):
    """Run MACSE alignSequences/refineAlignment for one or more files."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, input_files, parameters, mode,
                 imported_files=None, suffix=None, workdir=None,
                 mask_chars=True, remove_stop=True, remove_frameshift=True):
        super().__init__(
            tool_path, input_files, parameters, imported_files, workdir=workdir
        )
        self.mode = mode
        self.output_suffix = suffix or ""
        self.mask_chars = bool(mask_chars)
        self.remove_stop = bool(remove_stop)
        self.remove_frameshift = bool(remove_frameshift)

    def get_tool_name(self):
        return "MACSE2"

    def _output_paths(self, input_file):
        if self.output_suffix:
            base = os.path.splitext(input_file)[0]
            return (
                "%s%s_NT.fasta" % (base, self.output_suffix),
                "%s%s_AA.fasta" % (base, self.output_suffix),
            )
        return (
            self.create_temp_file(suffix="_NT.fasta"),
            self.create_temp_file(suffix="_AA.fasta"),
        )

    def execute_commands(self):
        try:
            output_files = []
            pairs = []
            total = len(self.input_files)
            for index, input_file in enumerate(self.input_files):
                if not self.is_running:
                    break
                self.progress.emit("Processing file %d/%d..." % (index + 1, total))
                self.console_output.emit(
                    "Processing file %d/%d: %s"
                    % (index + 1, total, os.path.basename(input_file)),
                    "info",
                )
                output_nt, output_aa = self._output_paths(input_file)
                cmd = [
                    self.tool_path,
                    "-prog", "alignSequences" if self.mode == "align" else "refineAlignment",
                    "-seq" if self.mode == "align" else "-align", input_file,
                    *self.parameters,
                    "-out_NT", output_nt,
                    "-out_AA", output_aa,
                ]
                result = self.execute_command(
                    cmd, cwd=str(Path(output_nt).parent)
                )
                detail = (result.stderr or "").strip()
                if result.returncode != 0:
                    raise RuntimeError(
                        "MACSE2 failed for %s%s"
                        % (os.path.basename(input_file),
                           ": " + detail if detail else "")
                    )
                if not Path(output_nt).is_file() or not Path(output_aa).is_file():
                    raise RuntimeError(
                        "MACSE2 did not write both output files for %s"
                        % os.path.basename(input_file)
                    )
                self.apply_post_processing(output_nt, output_aa)
                pairs.append((output_nt, output_aa))
                output_files.extend([output_nt, output_aa])

            self.progress.emit("All alignments completed")
            self.result_ready.emit({"output_files": output_files, "pairs": pairs})
            self.finished.emit(output_files, [])
        except Exception as exc:
            self.error.emit("MACSE2 alignment exception: %s" % exc)

    def apply_post_processing(self, nt_file, aa_file):
        try:
            nt_records = list(SeqIO.parse(nt_file, "fasta"))
            aa_records = list(SeqIO.parse(aa_file, "fasta"))
            if self.mask_chars:
                for record in nt_records:
                    record.seq = Seq(str(record.seq).replace("!", "?"))
                for record in aa_records:
                    record.seq = Seq(str(record.seq).replace("!", "?"))
            if self.remove_stop:
                nt_records = self.remove_stop_codon(nt_records, aa_records)
            if self.remove_frameshift:
                nt_records, aa_records = self.remove_frameshift_gaps(
                    nt_records, aa_records
                )
            SeqIO.write(nt_records, nt_file, "fasta")
            SeqIO.write(aa_records, aa_file, "fasta")
        except Exception as exc:
            self.console_output.emit(
                "Post-processing warning: %s" % exc, "warning"
            )

    @staticmethod
    def remove_stop_codon(nt_records, aa_records):
        for nt_record, aa_record in zip(nt_records, aa_records):
            aa_seq = str(aa_record.seq).rstrip("-")
            if aa_seq.endswith("*"):
                nt_pos = (len(aa_seq) - 1) * 3
                nt_seq = str(nt_record.seq)
                if nt_pos + 3 <= len(nt_seq):
                    nt_record.seq = Seq(nt_seq[:nt_pos] + "???" + nt_seq[nt_pos + 3:])
        return nt_records

    @staticmethod
    def remove_frameshift_gaps(nt_records, aa_records):
        """Drop codon columns where every NT character is !/-/?.

        The AA alignment has one column per codon, so the filter must be
        applied at codon granularity, not at raw nucleotide index.
        """
        if not nt_records or not aa_records:
            return nt_records, aa_records
        nt_sequences = [str(record.seq) for record in nt_records]
        aa_sequences = [str(record.seq) for record in aa_records]
        nt_length = len(nt_sequences[0])
        if nt_length % 3 or any(len(seq) != nt_length for seq in nt_sequences):
            return nt_records, aa_records
        if any(len(seq) != nt_length // 3 for seq in aa_sequences):
            return nt_records, aa_records

        keep_codons = []
        for codon_index in range(nt_length // 3):
            start = codon_index * 3
            columns = [seq[start:start + 3] for seq in nt_sequences]
            if any(not all(base in "!-?" for base in codon) for codon in columns):
                keep_codons.append(codon_index)

        if len(keep_codons) == nt_length // 3:
            return nt_records, aa_records

        for index, record in enumerate(nt_records):
            seq = nt_sequences[index]
            record.seq = Seq("".join(
                seq[codon * 3:codon * 3 + 3] for codon in keep_codons
            ))
        for index, record in enumerate(aa_records):
            seq = aa_sequences[index]
            record.seq = Seq("".join(seq[codon] for codon in keep_codons))
        return nt_records, aa_records


class MACSEPlugin(BasePlugin):
    """MACSE2 codon-aware alignment panel."""

    # Deprecated compatibility channels; the host uses result_ready.
    import_alignment_signal = pyqtSignal(list)
    batch_import_alignment_signal = pyqtSignal(list)

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self.batch_mode = False
        self.current_mode = "align"
        self._payload = None
        super().__init__(import_from, import_data, **kwargs)
        if import_data is not None and (
            not isinstance(import_data, list) or import_from != "YR_MPEA"
        ):
            self.handle_import_data(import_data)
        self.run_button.setText("Run MACSE2")

    def init_plugin_info(self):
        self.plugin_name = "MACSE2 Codon Aligner"
        self.tool_name = "MACSE2"
        self.citation = []
        self.citation_keys = ["ranwez2018macse", "ranwez2011macse"]
        self.input_types = {"Nucleotide FASTA": ["fas", "fna", "fa", "fasta"]}
        self.output_types = {"Codon/AA alignment": ["fasta"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    def control_bindings(self):
        """GUI controls mapped to the MACSE2 command line."""
        genetic_codes = {
            self.genetic_code_combo.itemText(i):
                self.genetic_code_combo.itemData(i)
            for i in range(self.genetic_code_combo.count())
        }
        return [
            {"id": "mode", "widget": "mode_combo",
             "label": "Mode", "kind": "choice",
             "choices": {
                 "Align sequences (alignSequences)":
                     ["-prog", "alignSequences", "-seq"],
                 "Refine an existing alignment (refineAlignment)":
                     ["-prog", "refineAlignment", "-align"],
             },
             "cli": ["{value}"],
             "notes": "The mode selects both -prog and -seq/-align."},
            {"id": "input", "widget": "file_path_edit",
             "label": "Input nucleotide FASTA", "kind": "path", "cli": [],
             "notes": "Passed with -seq in align mode or -align in refine mode."},
            {"id": "genetic_code", "widget": "genetic_code_combo",
             "label": "Genetic code", "kind": "choice",
             "choices": genetic_codes,
             "cli": ["-gc_def", "{value}"]},
            {"id": "aa_alphabet", "widget": "aa_alphabet_combo",
             "label": "AA alphabet", "kind": "choice",
             "cli": ["-alphabet_AA", "{value}"]},
            {"id": "score_matrix", "widget": "score_matrix_combo",
             "label": "Score matrix", "kind": "choice",
             "cli": ["-score_matrix", "{value}"]},
            {"id": "optimization", "widget": "optim_spinbox",
             "label": "Optimization cycles", "kind": "int",
             "cli": ["-optim", "{value}"]},
            {"id": "import_molecule", "widget": "return_option_combo",
             "label": "Preferred result", "kind": "choice", "cli": [],
             "notes": "GUI-only: chooses which of the NT/AA outputs is "
                      "selected for import."},
            {"id": "mask_chars", "widget": "mask_chars_checkbox",
             "label": "Mask exceptional characters", "kind": "bool", "cli": [],
             "notes": "GUI post-processing; not a MACSE option."},
            {"id": "remove_stop", "widget": "remove_stop_codon_checkbox",
             "label": "Mask terminal stop codon", "kind": "bool", "cli": [],
             "notes": "GUI post-processing; not a MACSE option."},
            {"id": "remove_frameshift", "widget": "remove_frameshift_checkbox",
             "label": "Drop all-gap codon columns", "kind": "bool", "cli": [],
             "notes": "GUI post-processing; not a MACSE option."},
        ]

    # ------------------------------------------------------------------
    # input
    # ------------------------------------------------------------------
    def setup_input_tab(self):
        outer = QVBoxLayout(self.input_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        outer.addWidget(scroll)
        content = QWidget()
        scroll.setWidget(content)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        title = QLabel("MACSE2")
        title.setStyleSheet("font-size: 18px; font-weight: 600;")
        layout.addWidget(title)
        subtitle = QLabel(
            "Codon-aware alignment that tolerates frameshifts and stop codons. "
            "Produces paired nucleotide and amino-acid alignments."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        layout.addWidget(subtitle)

        mode_group = QGroupBox("Mode")
        mode_form = QFormLayout(mode_group)
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "Align sequences (alignSequences)",
            "Refine an existing alignment (refineAlignment)",
        ])
        self.mode_combo.currentIndexChanged.connect(self.on_mode_changed)
        mode_form.addRow("Mode:", self.mode_combo)
        layout.addWidget(mode_group)

        input_group = QGroupBox("Input")
        form = QFormLayout(input_group)
        self.file_path_edit = QLineEdit()
        self.file_path_edit.setPlaceholderText("Select nucleotide FASTA files...")
        self.file_browse_btn = QPushButton("Browse Files")
        self.file_browse_btn.clicked.connect(self.browse_files)
        row = QHBoxLayout()
        row.addWidget(self.file_path_edit, 1)
        row.addWidget(self.file_browse_btn)
        form.addRow("Files:", row)

        self.file_tags_container = QFrame()
        self.file_tags_layout = QVBoxLayout()
        self.file_tags_container.setLayout(self.file_tags_layout)
        self.file_tags_container.setVisible(False)
        form.addRow("", self.file_tags_container)

        self.sequence_text = QTextEdit()
        self.sequence_text.setPlaceholderText("Or paste nucleotide FASTA...")
        self.sequence_text.setMaximumHeight(160)
        self.sequence_text.textChanged.connect(self.on_text_changed)
        form.addRow("Paste:", self.sequence_text)

        self.return_option_combo = QComboBox()
        self.return_option_combo.addItems(["Codon sequences", "Translated AA sequences"])
        form.addRow("Preferred result:", self.return_option_combo)
        layout.addWidget(input_group)

        params = QGroupBox("Parameters")
        pform = QFormLayout(params)
        self.genetic_code_combo = QComboBox()
        for code, label in GENETIC_CODES:
            self.genetic_code_combo.addItem(label, code)
        pform.addRow("Genetic code:", self.genetic_code_combo)
        self.aa_alphabet_combo = QComboBox()
        self.aa_alphabet_combo.addItems(AA_ALPHABETS)
        self.aa_alphabet_combo.setCurrentText("SE_B_8")
        pform.addRow("AA alphabet:", self.aa_alphabet_combo)
        self.score_matrix_combo = QComboBox()
        self.score_matrix_combo.addItems(["BLOSUM62", "VTML200_BIS", "VTML240"])
        pform.addRow("Score matrix:", self.score_matrix_combo)
        self.optim_spinbox = QSpinBox()
        self.optim_spinbox.setRange(0, 2)
        self.optim_spinbox.setValue(2)
        pform.addRow("Optimization:", self.optim_spinbox)
        layout.addWidget(params)

        post = QGroupBox("Post-processing")
        post_layout = QVBoxLayout(post)
        self.mask_chars_checkbox = QCheckBox("Mask exceptional characters (! → ?)")
        self.mask_chars_checkbox.setChecked(True)
        post_layout.addWidget(self.mask_chars_checkbox)
        self.remove_frameshift_checkbox = QCheckBox("Drop all-gap codon columns")
        self.remove_frameshift_checkbox.setChecked(True)
        post_layout.addWidget(self.remove_frameshift_checkbox)
        self.remove_stop_codon_checkbox = QCheckBox("Mask terminal stop codon")
        self.remove_stop_codon_checkbox.setChecked(True)
        post_layout.addWidget(self.remove_stop_codon_checkbox)
        advanced_btn = QPushButton("Advanced parameters...")
        advanced_btn.clicked.connect(self.show_advanced_dialog)
        post_layout.addWidget(advanced_btn)
        layout.addWidget(post)

        output_group = QGroupBox("Output")
        oform = QFormLayout(output_group)
        mode_row = QHBoxLayout()
        self.save_to_cwd = QRadioButton("Beside input")
        self.save_to_cwd.setChecked(True)
        self.save_to_tmp = QRadioButton("Temporary file")
        mode_row.addWidget(self.save_to_cwd)
        mode_row.addWidget(self.save_to_tmp)
        oform.addRow("Mode:", mode_row)
        self.output_suffix_edit = QLineEdit("_macse")
        oform.addRow("Suffix:", self.output_suffix_edit)
        self.save_to_cwd.toggled.connect(self.on_output_radio_changed)
        layout.addWidget(output_group)
        layout.addStretch()
        self.file_tags = []

    def on_mode_changed(self, index):
        self.current_mode = "align" if index == 0 else "refine"
        self.file_path_edit.setPlaceholderText(
            "Select unaligned nucleotide FASTA files..."
            if self.current_mode == "align"
            else "Select aligned nucleotide FASTA files..."
        )

    def on_output_radio_changed(self):
        enabled = self.save_to_cwd.isChecked()
        self.output_suffix_edit.setEnabled(enabled)
        self.output_suffix_edit.setVisible(enabled)

    def setup_output_tab(self):
        setup_preview(self, "The MACSE2 nucleotide/AA output will appear here.")

    def setup_control_panel(self):
        super().setup_control_panel()
        self.import_to_platform_btn = QPushButton("Send to workspace")
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        self.import_to_platform_btn.setVisible(False)
        self.control_layout.addWidget(self.import_to_platform_btn)

    def handle_import_data(self, import_data):
        files, batch = prepare_import(self, import_data)
        self.batch_mode = batch
        self.imported_files = files
        self.import_file = files[0] if len(files) == 1 else None
        if hasattr(self, "file_path_edit"):
            if not files:
                self.file_path_edit.clear()
            elif len(files) == 1:
                self.file_path_edit.setText(files[0])
            else:
                self.file_path_edit.setText("%d files selected" % len(files))
                self.file_path_edit.setEnabled(False)

    # ------------------------------------------------------------------
    # input helpers
    # ------------------------------------------------------------------
    def browse_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select nucleotide FASTA files", "",
            "FASTA files (*.fas *.fasta *.fa *.fna);;All files (*)",
        )
        for path in paths:
            if path not in self.imported_files:
                self.add_file_tag(path)
        if self.imported_files:
            self.file_path_edit.setText(
                self.imported_files[0] if len(self.imported_files) == 1
                else "%d files selected" % len(self.imported_files)
            )
            self.sequence_text.setVisible(False)
            self.sequence_text.setEnabled(False)

    def add_file_tag(self, file_path):
        self.imported_files.append(file_path)
        self.import_file = self.imported_files[0] if len(self.imported_files) == 1 else None
        self.batch_mode = len(self.imported_files) > 1
        tag = QFrame()
        tag.setStyleSheet(
            "QFrame { background-color: palette(base); border: 1px solid palette(mid);"
            " border-radius: 12px; margin: 2px; }"
        )
        tag_layout = QHBoxLayout()
        tag_layout.setContentsMargins(8, 4, 8, 4)
        tag.setLayout(tag_layout)
        tag_layout.addWidget(QLabel(self.get_display_name(file_path)))
        close = QPushButton("×")
        close.setFlat(True)
        close.setFixedSize(20, 20)
        close.clicked.connect(lambda: self.remove_file_tag(file_path, tag))
        tag_layout.addWidget(close)
        self.file_tags_layout.addWidget(tag)
        self.file_tags.append((file_path, tag))
        self.file_tags_container.setVisible(True)

    def remove_file_tag(self, file_path, tag_widget):
        if file_path in self.imported_files:
            self.imported_files.remove(file_path)
        self.file_tags = [(fp, tw) for fp, tw in self.file_tags if fp != file_path]
        tag_widget.setParent(None)
        tag_widget.deleteLater()
        self.batch_mode = len(self.imported_files) > 1
        self.import_file = self.imported_files[0] if len(self.imported_files) == 1 else None
        if self.imported_files:
            self.file_path_edit.setText(
                self.imported_files[0] if len(self.imported_files) == 1
                else "%d files selected" % len(self.imported_files)
            )
        else:
            self.file_path_edit.clear()
            self.file_tags_container.setVisible(False)
            self.sequence_text.setVisible(True)
            self.sequence_text.setEnabled(True)

    def clear_all_file_tags(self):
        for _, tag in self.file_tags:
            tag.setParent(None)
            tag.deleteLater()
        self.file_tags = []
        self.imported_files = []
        self.import_file = None
        self.batch_mode = False
        self.file_tags_container.setVisible(False)

    def get_display_name(self, file_path):
        filename = os.path.basename(file_path)
        duplicates = [
            path for path in self.imported_files
            if os.path.basename(path) == filename and path != file_path
        ]
        if duplicates:
            return "%s (%s)" % (filename, os.path.basename(os.path.dirname(file_path)))
        return filename

    def on_text_changed(self):
        text = self.sequence_text.toPlainText().strip()
        if text:
            self.file_path_edit.setVisible(False)
            self.file_browse_btn.setVisible(False)
            self.file_tags_container.setVisible(False)
            self.save_to_tmp.setChecked(True)
            self.save_to_cwd.setEnabled(False)
            self.clear_all_file_tags()
        else:
            self.file_path_edit.setVisible(True)
            self.file_browse_btn.setVisible(True)
            self.save_to_cwd.setEnabled(True)
            if self.imported_files:
                self.file_tags_container.setVisible(True)

    def prepare_input_files(self):
        if self.imported_files:
            return [path for path in self.imported_files if os.path.exists(path)]
        text = self.sequence_text.toPlainText().strip()
        if not text:
            QMessageBox.warning(self, "Input missing", "Select files or paste nucleotide FASTA.")
            return None
        try:
            records = list(SeqIO.parse(StringIO(text), "fasta"))
        except Exception as exc:
            QMessageBox.warning(self, "Invalid FASTA", str(exc))
            return None
        if len(records) < 2:
            QMessageBox.warning(self, "Need more sequences", "Paste at least two FASTA records.")
            return None
        return [records_temp(self, records, "pasted")]

    def show_advanced_dialog(self):
        dialog = MACSEAdvancedDialog(self)
        if dialog.exec_() == QDialog.Accepted:
            self.advanced_params = dialog.get_parameters()

    def get_parameters(self):
        params = [
            "-gc_def", str(self.genetic_code_combo.currentData()),
            "-alphabet_AA", self.aa_alphabet_combo.currentText(),
            "-score_matrix", self.score_matrix_combo.currentText(),
            "-optim", str(self.optim_spinbox.value()),
        ]
        params.extend(getattr(self, "advanced_params", []))
        return params

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    def run_analysis(self):
        if self.is_running:
            return
        if not self.tool_path or not Path(self.tool_path).is_file():
            QMessageBox.critical(self, "Tool missing", "MACSE2 executable not found.")
            return
        input_files = self.prepare_input_files()
        if not input_files:
            return
        suffix = None
        if self.save_to_cwd.isChecked():
            suffix = self.output_suffix_edit.text()
            if any(sep in suffix for sep in ("/", "\\")):
                QMessageBox.warning(self, "Invalid suffix", "Suffix must not contain a path separator.")
                return

        self._payload = None
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running MACSE2...")
        self.import_to_platform_btn.setEnabled(True)
        self.import_to_platform_btn.setText("Send to workspace")
        self.tab_widget.setCurrentIndex(2)

        self.thread = MACSEThread(
            self.tool_path, input_files, self.get_parameters(),
            self.current_mode, self.imported_files, suffix=suffix,
            workdir=self.workdir,
            mask_chars=self.mask_chars_checkbox.isChecked(),
            remove_stop=self.remove_stop_codon_checkbox.isChecked(),
            remove_frameshift=self.remove_frameshift_checkbox.isChecked(),
        )
        self.thread.progress.connect(self.progress_bar.setFormat)
        self.thread.result_ready.connect(self._on_result)
        self.thread.finished.connect(self.on_analysis_finished)
        self.thread.error.connect(self.on_analysis_error)
        self.thread.console_output.connect(self.add_console_message)
        self.thread.start()

    def _on_result(self, payload):
        self._payload = payload

    def on_analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.alignment_output_files = output_files
        nt_files = [path for path in output_files if "_NT.fasta" in path]
        show_preview(self, nt_files or output_files)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("MACSE2 alignment completed.", "info")
        if self.result_has_consumer() or self.import_from in [
            "YR_MPEA", "seq_viewer", "DATASET_MANAGER"
        ]:
            self.import_to_platform_btn.setVisible(True)
        self.import_to_platform_btn.setEnabled(True)

    def on_analysis_error(self, error_msg):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(error_msg, "error")
        QMessageBox.critical(self, "MACSE2 failed", error_msg)
        self.tab_widget.setCurrentIndex(2)

    def stop_analysis(self):
        if hasattr(self, "thread") and self.thread.isRunning():
            thread = self.thread
            for signal in (
                thread.progress, thread.result_ready, thread.finished,
                thread.error, thread.console_output,
            ):
                try:
                    signal.disconnect()
                except (TypeError, RuntimeError):
                    pass
            thread.stop()
            thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("MACSE2 alignment stopped by user.", "warning")

    # ------------------------------------------------------------------
    # result
    # ------------------------------------------------------------------
    def _build_result_spec(self):
        if not getattr(self, "alignment_output_files", None):
            raise ValueError("No MACSE2 output is available.")
        preferred = "NT" if self.return_option_combo.currentIndex() == 0 else "AA"
        items = []
        nt_alignments = 0
        total_sequences = 0
        for path in self.alignment_output_files:
            records, file_format = read_alignment(path)
            if len(records) < 2:
                raise ValueError("%s contains fewer than two sequences." % path)
            lengths = {len(str(record.seq)) for record in records}
            if len(lengths) > 1:
                raise ValueError("%s is not rectangular." % path)
            molecule = "AA" if "_AA" in Path(path).stem else "NT"
            if molecule == "NT":
                nt_alignments += 1
                total_sequences += len(records)
            items.append(DatasetItemSpec(
                loci_name=Path(path).stem,
                item_type="alignment",
                sequences=records,
                file_path=str(path),
                format=file_format,
                is_aligned=True,
                length=max(lengths),
                sequence_count=len(records),
                selected=(molecule == preferred),
                data={
                    "tool": "MACSE2",
                    "molecule": molecule,
                    "mode": self.current_mode,
                    "parameters": self.get_parameters(),
                },
                metadata={"source": self.plugin_name},
            ))
        inputs = list(self.imported_files or [])
        dataset = DatasetSpec(
            name="MACSE2 alignment" if nt_alignments == 1 else "MACSE2 alignments",
            items=items,
            description="%d codon alignment(s), %d sequence records"
                        % (nt_alignments, total_sequences),
            source="macse2",
            summary={
                "alignments": nt_alignments,
                "sequences": total_sequences,
                "files": len(items),
            },
            metadata={
                "input": inputs[0] if len(inputs) == 1 else "",
                "input_files": inputs,
                "parameters": self.get_parameters(),
                "mode": self.current_mode,
                "mask_chars": self.mask_chars_checkbox.isChecked(),
                "remove_stop": self.remove_stop_codon_checkbox.isChecked(),
                "remove_frameshift": self.remove_frameshift_checkbox.isChecked(),
            },
        )
        return ResultSpec(source="macse2", datasets=[dataset])

    def import_to_platform(self):
        if not getattr(self, "alignment_output_files", None):
            QMessageBox.warning(self, "No result", "Run MACSE2 first.")
            return
        if self.result_has_consumer():
            try:
                result = self._build_result_spec()
            except Exception as exc:
                QMessageBox.critical(self, "Result failed", str(exc))
                return
            if self.export_result(result):
                self.import_to_platform_btn.setEnabled(False)
                self.import_to_platform_btn.setText("Sent to workspace")
                QMessageBox.information(
                    self, "Success", "MACSE2 result sent to the current workspace."
                )
            return

        # Legacy fallback for external callers.
        try:
            preferred_nt = self.return_option_combo.currentIndex() == 0
            nt_files = [f for f in self.alignment_output_files if "_NT.fasta" in f]
            aa_files = [f for f in self.alignment_output_files if "_AA.fasta" in f]
            selected = nt_files if preferred_nt else aa_files
            if self.batch_mode:
                self.batch_import_alignment_signal.emit(
                    [list(SeqIO.parse(path, "fasta")) for path in selected]
                )
            else:
                sequences = []
                for path in selected:
                    sequences.extend(list(SeqIO.parse(path, "fasta")))
                if sequences:
                    self.import_alignment_signal.emit(sequences)
        except Exception as exc:
            QMessageBox.critical(self, "Import failed", str(exc))


class MACSEAdvancedDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MACSE2 Advanced Parameters")
        self.setMinimumSize(520, 220)
        layout = QVBoxLayout(self)
        hint = QLabel(
            "Additional MACSE2 command-line parameters. Quoted values are supported."
        )
        hint.setWordWrap(True)
        layout.addWidget(hint)
        self.param_edit = QTextEdit()
        self.param_edit.setPlaceholderText(
            "Example: -max_refine_iter 5 -local_realign_init 0.5"
        )
        layout.addWidget(self.param_edit)
        buttons = QHBoxLayout()
        reset = QPushButton("Reset")
        reset.clicked.connect(self.param_edit.clear)
        buttons.addWidget(reset)
        buttons.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        buttons.addWidget(cancel)
        ok = QPushButton("OK")
        ok.clicked.connect(self.accept)
        buttons.addWidget(ok)
        layout.addLayout(buttons)

    def get_parameters(self):
        text = self.param_edit.toPlainText().strip()
        if not text:
            return []
        try:
            return shlex.split(text)
        except ValueError:
            return text.split()


class MACSEPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return MACSEPlugin(
            import_from=import_from, import_data=import_data,
            workdir=workdir, **kwargs,
        )
