# orthofinder_plugin.py
#
# Copyright (c) 2025 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Genomic Dataset Formation (OrthoFinder).

GBFF/GB/GBK -> CDS/AA -> OrthoFinder -> single-copy orthogroups -> CDS.

This plugin is a thin GUI around ``gbff_utils`` and ``orthogroup_utils``.
It is intentionally non-modal: OrthoFinder runs in a BaseProcessThread and
streams its log to the Console tab.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import shutil
import tempfile
import time
from pathlib import Path

from Bio import SeqIO
from PyQt5.QtCore import Qt, QUrl, pyqtSignal
from PyQt5.QtGui import QDesktopServices, QFont
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from .components.methods.gbff_utils import parse_folder, write_dataset
from .components.methods.orthogroup_utils import (
    derive_single_copy_ids,
    find_orthogroups_dir,
    read_cds_map,
    read_orthogroups,
    read_single_copy_ids,
    write_single_copy_ogs,
)


SOURCE_TYPES = ("gbff", "gbff_records", "proteome")
RUN_MODES = ("full", "parse_only", "existing")
INTERNAL_STOPS = ("skip", "warn", "keep")
GENETIC_CODES = (None, 1, 2, 4, 5, 11)


def _timestamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def _ui_font(point_size: int = 10, bold: bool = False) -> QFont:
    """Derive a font from the host application font (YRTools: MiSans Medium)."""
    font = QFont(QApplication.font())
    font.setPointSize(point_size)
    if bold:
        font.setBold(True)
    return font


def _path_row(edit: QLineEdit, button: QPushButton) -> QWidget:
    row = QWidget()
    layout = QHBoxLayout(row)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(8)
    layout.addWidget(edit, 1)
    layout.addWidget(button)
    return row


def _section(title: str) -> tuple[QGroupBox, QVBoxLayout]:
    group = QGroupBox(title)
    layout = QVBoxLayout(group)
    layout.setSpacing(8)
    return group, layout


class OrthoFinderThread(BaseProcessThread):
    """Parse input, run OrthoFinder, and post-process single-copy OGs."""

    result_summary = pyqtSignal(dict)

    def __init__(
        self,
        input_dir: str,
        output_dir: str,
        source_mode: str,
        default_table: int,
        skip_pseudo: bool,
        internal_stop: str,
        keep_stop: bool,
        run_mode: str,
        orthofinder_path: str,
        blast_dir: str,
        threads: int,
        existing_root: str | None = None,
        **kwargs,
    ):
        super().__init__(orthofinder_path or "", [input_dir], [], workdir=output_dir, **kwargs)
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.source_mode = source_mode
        self.default_table = default_table
        self.skip_pseudo = skip_pseudo
        self.internal_stop = internal_stop
        self.keep_stop = keep_stop
        self.run_mode = run_mode
        self.orthofinder_path = Path(orthofinder_path) if orthofinder_path else None
        self.blast_dir = Path(blast_dir) if blast_dir else None
        self.threads = threads
        self.existing_root = Path(existing_root) if existing_root else None

    def get_tool_name(self) -> str:
        return "OrthoFinder"

    # ------------------------------------------------------------------
    # validation helpers
    # ------------------------------------------------------------------
    def _validate_runtime(self) -> None:
        if not self.orthofinder_path or not self.orthofinder_path.exists():
            raise ValueError(f"OrthoFinder executable not found: {self.orthofinder_path}")
        missing = []
        for name in ("blastp.exe", "makeblastdb.exe", "mcl.exe"):
            candidates = []
            if self.blast_dir:
                candidates.append(self.blast_dir / name)
            candidates.append(self.orthofinder_path.parent / name)
            if not any(p.exists() for p in candidates):
                missing.append(name)
        if missing:
            raise ValueError(
                "Missing runtime tools: " + ", ".join(missing)
                + f"\nLooked in: {self.blast_dir} and {self.orthofinder_path.parent}"
            )

    def _copy_proteomes(self) -> dict[str, int]:
        proteomes = self.output_dir / "proteomes"
        proteomes.mkdir(parents=True, exist_ok=True)
        suffixes = {".fa", ".faa", ".fasta", ".fas"}
        summary: dict[str, int] = {}
        for path in sorted(self.input_dir.iterdir()):
            if not path.is_file() or path.suffix.lower() not in suffixes:
                continue
            target = proteomes / f"{path.stem}.faa"
            shutil.copy2(path, target)
            count = 0
            with target.open("r", encoding="utf-8", errors="replace") as handle:
                count = sum(1 for line in handle if line.startswith(">"))
            summary[path.stem] = count
        if not summary:
            raise ValueError(f"No protein FASTA files found in {self.input_dir}")
        return summary

    def _run_orthofinder(self) -> None:
        self._validate_runtime()
        results_root = (self.output_dir / "orthofinder").resolve()
        if results_root.exists():
            shutil.rmtree(results_root)
        # OrthoFinder 2.5.5 uses a relative '..\temp\WorkingDirectory' path.
        (self.output_dir.parent / "temp" / "WorkingDirectory").mkdir(parents=True, exist_ok=True)

        env = os.environ.copy()
        path_parts = []
        if self.blast_dir:
            path_parts.append(str(self.blast_dir))
        path_parts.append(str(self.orthofinder_path.parent))
        env["PATH"] = os.pathsep.join(path_parts + [env.get("PATH", "")])

        cmd = [
            str(self.orthofinder_path),
            "-f", str((self.output_dir / "proteomes").resolve()),
            "-o", str(results_root),
            "-t", str(self.threads),
            "-S", "blast",
            "-og",
            "-X",
        ]
        old_path = os.environ.get("PATH")
        os.environ["PATH"] = env["PATH"]
        try:
            result = self.execute_command(cmd, cwd=str(self.output_dir), env=env)
        finally:
            if old_path is None:
                os.environ.pop("PATH", None)
            else:
                os.environ["PATH"] = old_path
        if result.returncode != 0:
            raise RuntimeError(
                f"OrthoFinder failed with exit code {result.returncode}.\n"
                f"Last output:\n{(result.stderr or result.stdout)[-2000:]}"
            )

    # ------------------------------------------------------------------
    # main
    # ------------------------------------------------------------------
    def execute_commands(self) -> None:
        started = time.time()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if self.source_mode == "proteome":
            self.progress.emit("Preparing proteomes...")
            parsed = self._copy_proteomes()
            cds_count = 0
            self.console_output.emit(
                f"Copied {len(parsed)} proteomes; protein-only mode (no CDS mapping)", "info"
            )
        else:
            self.progress.emit("Parsing GBFF/GenBank files...")
            records = parse_folder(
                self.input_dir,
                default_table=self.default_table,
                skip_pseudo=self.skip_pseudo,
                internal_stop=self.internal_stop,
                keep_stop=self.keep_stop,
                record_as_species=(self.source_mode == "gbff_records"),
            )
            if not records:
                raise ValueError(
                    "No CDS records parsed. Check the input folder and genetic code."
                )
            parsed = write_dataset(records, self.output_dir)
            cds_count = sum(parsed.values())
            self.console_output.emit(
                f"Parsed {len(parsed)} species / {cds_count} CDS", "info"
            )

        if self.run_mode == "parse_only":
            summary = {
                "species": len(parsed),
                "cds": cds_count,
                "orthogroups": 0,
                "single_copy_ogs": 0,
                "output_dir": str(self.output_dir),
                "elapsed": round(time.time() - started, 1),
                "run_mode": "parse_only",
                "cds_mapping": self.source_mode in ("gbff", "gbff_records"),
            }
            self.result_summary.emit(summary)
            self.finished.emit(self._output_files(), [])
            return

        if self.run_mode == "full":
            self.progress.emit("Running OrthoFinder (BLAST + MCL)...")
            self._run_orthofinder()

        # Post-process either the new run or an existing results folder.
        root = self.existing_root if self.run_mode == "existing" else self.output_dir
        if root is None:
            root = self.output_dir
        og_dir = find_orthogroups_dir(Path(root))
        cds_map_path = self.output_dir / "cds_map.tsv"
        cds_map = read_cds_map(cds_map_path) if cds_map_path.exists() else {}
        orthogroups = read_orthogroups(og_dir / "Orthogroups.tsv")
        single_path = og_dir / "Orthogroups_SingleCopyOrthologues.txt"
        if single_path.exists():
            single_ids = read_single_copy_ids(single_path)
        else:
            single_ids = derive_single_copy_ids(orthogroups)

        if cds_map:
            post = write_single_copy_ogs(
                orthogroups, single_ids, cds_map, self.output_dir / "single_copy_ogs"
            )
            single_count = post["orthogroups_written"]
        else:
            single_count = len(single_ids)
            self.console_output.emit(
                "No cds_map.tsv; single-copy OGs are counted but CDS back-mapping is skipped "
                "(proteome-only input).",
                "warning",
            )

        summary = {
            "species": len(parsed),
            "cds": cds_count,
            "orthogroups": len(orthogroups),
            "single_copy_ogs": single_count,
            "output_dir": str(self.output_dir),
            "elapsed": round(time.time() - started, 1),
            "run_mode": self.run_mode,
            "cds_mapping": bool(cds_map),
        }
        self.result_summary.emit(summary)
        self.finished.emit(self._output_files(), [])

    def _output_files(self) -> list[str]:
        files = []
        for name in ("cds_map.tsv",):
            path = self.output_dir / name
            if path.exists():
                files.append(str(path))
        summary = self.output_dir / "single_copy_ogs" / "single_copy_og_summary.tsv"
        if summary.exists():
            files.append(str(summary))
        return files


class OrthoFinderPlugin(BasePlugin):
    """Three-step panel for GBFF -> OrthoFinder -> single-copy OG -> CDS."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self._last_summary = None
        self.analysis_thread = None
        self._opened_from_yrmpea = import_from == "YR_MPEA"
        self._status_state = "idle"
        self._started_at = None
        self._completed_at = None
        self._last_error = ""
        self._progress_text = ""
        self._show_run_dialogs = True
        super().__init__(import_from, import_data, **kwargs)
        self.run_button.setText("Run pipeline")

    def init_plugin_info(self):
        self.plugin_name = "One-Step OrthoFinder (for Organelle Genomes)"
        self.tool_name = "OrthoFinder"
        self.citation = []  # citation text lives in YR_MPE/citations/*.bib
        self.input_types = {
            "GenBank": ["gb", "gbff", "gbk"],
            "Proteome FASTA": ["fa", "faa", "fasta"],
        }
        self.output_types = {"OrthoFinder dataset": ["faa", "fna", "tsv"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    # ------------------------------------------------------------------
    # input tab / wizard
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

        title = QLabel("One-Step OrthoFinder (for Organelle Genomes)")
        title.setFont(_ui_font(15, bold=True))
        layout.addWidget(title)

        subtitle = QLabel("GBFF/GB → CDS/AA → OrthoFinder → single-copy orthogroups → CDS")
        subtitle.setStyleSheet("color: palette(mid);")
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        # step indicator
        self.step_labels = []
        step_row = QHBoxLayout()
        step_row.setSpacing(16)
        for index, text in enumerate(("1 Input", "2 OrthoFinder", "3 Review")):
            label = QLabel(text)
            label.setMinimumHeight(26)
            self.step_labels.append(label)
            step_row.addWidget(label)
        step_row.addStretch()
        layout.addLayout(step_row)

        self.step_stack = QStackedWidget()
        self.step_stack.addWidget(self._build_input_page())
        self.step_stack.addWidget(self._build_orthofinder_page())
        self.step_stack.addWidget(self._build_review_page())
        layout.addWidget(self.step_stack, 1)

        nav = QHBoxLayout()
        self.back_button = QPushButton("← Back")
        self.back_button.clicked.connect(self._go_back)
        self.next_button = QPushButton("Next: OrthoFinder →")
        self.next_button.clicked.connect(self._go_next)
        nav.addWidget(self.back_button)
        nav.addStretch()
        nav.addWidget(self.next_button)
        layout.addLayout(nav)

        self._set_step(0)
        self._update_review()

    def _build_input_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setSpacing(12)
        layout.setContentsMargins(0, 0, 0, 0)

        source_group, source_layout = _section("Input source")
        source_form = QFormLayout()
        source_form.setSpacing(8)

        self.source_combo = QComboBox()
        self.source_combo.addItem("GenBank folder (one file per species)", "gbff")
        self.source_combo.addItem("GenBank multi-record file (each record = taxon)", "gbff_records")
        self.source_combo.addItem("Proteome FASTA (protein-only)", "proteome")
        self.source_combo.currentIndexChanged.connect(self._on_source_changed)
        source_form.addRow("Source type:", self.source_combo)

        self.input_edit = QLineEdit()
        self.input_edit.setPlaceholderText("Folder with .gb/.gbff/.gbk files")
        self.input_browse_button = QPushButton("Browse")
        self.input_browse_button.clicked.connect(self._browse_input)
        self.input_label = QLabel("Input folder:")
        source_form.addRow(self.input_label, _path_row(self.input_edit, self.input_browse_button))
        source_layout.addLayout(source_form)

        self.scan_label = QLabel("No folder selected.")
        self.scan_label.setWordWrap(True)
        self.scan_label.setStyleSheet("color: palette(mid);")
        scan_row = QHBoxLayout()
        scan_row.addWidget(self.scan_label, 1)
        self.scan_button = QPushButton("Scan folder")
        self.scan_button.clicked.connect(self._scan_input)
        scan_row.addWidget(self.scan_button)
        source_layout.addLayout(scan_row)

        layout.addWidget(source_group)

        parse_group, parse_layout = _section("CDS parsing")
        parse_form = QFormLayout()
        parse_form.setSpacing(8)
        self.genetic_code_combo = QComboBox()
        self.genetic_code_combo.addItem("Auto (from /transl_table)", None)
        self.genetic_code_combo.addItem("1 — Standard", 1)
        self.genetic_code_combo.addItem("2 — Vertebrate mitochondrial", 2)
        self.genetic_code_combo.addItem("4 — Mold mitochondrial", 4)
        self.genetic_code_combo.addItem("5 — Invertebrate mitochondrial", 5)
        self.genetic_code_combo.addItem("11 — Bacterial/Plastid", 11)
        parse_form.addRow("Genetic code:", self.genetic_code_combo)

        self.skip_pseudo_check = QCheckBox("Skip pseudogenes")
        self.skip_pseudo_check.setChecked(True)
        parse_form.addRow("", self.skip_pseudo_check)

        self.internal_stop_combo = QComboBox()
        self.internal_stop_combo.addItem("Skip", "skip")
        self.internal_stop_combo.addItem("Warn and keep", "warn")
        self.internal_stop_combo.addItem("Keep", "keep")
        parse_form.addRow("Internal stop codons:", self.internal_stop_combo)

        self.keep_stop_check = QCheckBox("Keep terminal stop codon in CDS")
        parse_form.addRow("", self.keep_stop_check)
        parse_layout.addLayout(parse_form)
        layout.addWidget(parse_group)

        output_group, output_layout = _section("Output")
        output_form = QFormLayout()
        output_form.setSpacing(8)
        self.output_edit = QLineEdit()
        self.output_edit.setPlaceholderText("Output folder")
        self.output_browse_button = QPushButton("Browse")
        self.output_browse_button.clicked.connect(self._browse_output)
        output_form.addRow("Output folder:", _path_row(self.output_edit, self.output_browse_button))
        output_layout.addLayout(output_form)
        layout.addWidget(output_group)

        self.tip_label = QLabel(
            "Tip: one file per species is recommended. /translation is used when "
            "present; otherwise the CDS is translated with the selected genetic code."
        )
        self.tip_label.setWordWrap(True)
        self.tip_label.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.tip_label)
        layout.addStretch()
        return page

    def _default_blast_dir(self) -> str:
        if not self.tool_path:
            return ""
        base = os.path.dirname(self.tool_path)
        # Development layout: orthofinder.exe lives next to tools/blast_tb.
        candidate = os.path.join(base, "tools", "blast_tb")
        if os.path.isdir(candidate):
            return candidate
        return base

    def _build_orthofinder_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setSpacing(12)
        layout.setContentsMargins(0, 0, 0, 0)

        runtime_group, runtime_layout = _section("OrthoFinder runtime")
        runtime_form = QFormLayout()
        runtime_form.setSpacing(8)

        self.orthofinder_edit = QLineEdit()
        self.orthofinder_edit.setPlaceholderText("Path to orthofinder.exe")
        self.orthofinder_edit.setText(self.tool_path or "")
        self.orthofinder_browse_button = QPushButton("Browse")
        self.orthofinder_browse_button.clicked.connect(self._browse_orthofinder)
        runtime_form.addRow("Executable:", _path_row(self.orthofinder_edit, self.orthofinder_browse_button))

        self.blast_edit = QLineEdit()
        self.blast_edit.setPlaceholderText("Directory with blastp/makeblastdb/mcl")
        self.blast_edit.setText(self._default_blast_dir())
        self.blast_browse_button = QPushButton("Browse")
        self.blast_browse_button.clicked.connect(self._browse_blast)
        runtime_form.addRow("BLAST+ / MCL dir:", _path_row(self.blast_edit, self.blast_browse_button))

        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 16)
        self.threads_spin.setValue(2)
        runtime_form.addRow("Threads:", self.threads_spin)

        runtime_layout.addLayout(runtime_form)
        status_row = QHBoxLayout()
        self.runtime_status = QLabel("Runtime not checked.")
        self.runtime_status.setWordWrap(True)
        self.check_runtime_button = QPushButton("Check runtime")
        self.check_runtime_button.clicked.connect(lambda: self._check_runtime(show=True))
        status_row.addWidget(self.runtime_status, 1)
        status_row.addWidget(self.check_runtime_button)
        runtime_layout.addLayout(status_row)
        layout.addWidget(runtime_group)

        mode_group, mode_layout = _section("Run mode")
        self.run_mode_combo = QComboBox()
        self.run_mode_combo.addItem("Parse + OrthoFinder + single-copy OGs", "full")
        self.run_mode_combo.addItem("Parse only", "parse_only")
        self.run_mode_combo.addItem("Use existing OrthoFinder results", "existing")
        self.run_mode_combo.currentIndexChanged.connect(self._on_run_mode_changed)
        mode_form = QFormLayout()
        mode_form.setSpacing(8)
        mode_form.addRow("Mode:", self.run_mode_combo)

        self.existing_edit = QLineEdit()
        self.existing_edit.setPlaceholderText("Folder containing Orthogroups.tsv")
        self.existing_browse_button = QPushButton("Browse")
        self.existing_browse_button.clicked.connect(self._browse_existing)
        self.existing_row = _path_row(self.existing_edit, self.existing_browse_button)
        self.existing_label = QLabel("Existing results:")
        self.existing_container = QWidget()
        existing_layout = QHBoxLayout(self.existing_container)
        existing_layout.setContentsMargins(0, 0, 0, 0)
        existing_layout.setSpacing(8)
        existing_layout.addWidget(self.existing_label)
        existing_layout.addWidget(self.existing_row, 1)
        mode_layout.addLayout(mode_form)
        mode_layout.addWidget(self.existing_container)

        self.mode_hint = QLabel(
            "Parse only writes proteomes + CDS + cds_map and does not require OrthoFinder."
        )
        self.mode_hint.setWordWrap(True)
        self.mode_hint.setStyleSheet("color: palette(mid);")
        mode_layout.addWidget(self.mode_hint)
        layout.addWidget(mode_group)
        layout.addStretch()
        self._on_run_mode_changed()
        return page

    def _build_review_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setSpacing(12)
        layout.setContentsMargins(0, 0, 0, 0)

        group, group_layout = _section("Review")
        self.review_label = QLabel()
        self.review_label.setWordWrap(True)
        self.review_label.setTextFormat(Qt.RichText)
        group_layout.addWidget(self.review_label)
        layout.addWidget(group)

        note = QLabel(
            "Click <b>Run pipeline</b> below to parse the input and run OrthoFinder. "
            "Long runs are non-modal; watch the Console tab."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color: palette(mid);")
        layout.addWidget(note)
        layout.addStretch()
        return page

    # ------------------------------------------------------------------
    # wizard helpers
    # ------------------------------------------------------------------
    def _set_step(self, index: int) -> None:
        self.step_stack.setCurrentIndex(index)
        for i, label in enumerate(self.step_labels):
            if i < index:
                label.setText(f"✓ {i + 1} " + ("Input", "OrthoFinder", "Review")[i])
                label.setStyleSheet("color: #2e7d32; font-weight: 600;")
            elif i == index:
                label.setText(f"● {i + 1} " + ("Input", "OrthoFinder", "Review")[i])
                label.setStyleSheet("color: #1565c0; font-weight: 600; border-bottom: 2px solid #1565c0;")
            else:
                label.setText(f"○ {i + 1} " + ("Input", "OrthoFinder", "Review")[i])
                label.setStyleSheet("color: palette(mid);")
        self.back_button.setEnabled(index > 0)
        if index == 0:
            self.next_button.setText("Next: OrthoFinder →")
            self.next_button.setEnabled(True)
        elif index == 1:
            self.next_button.setText("Next: Review →")
            self.next_button.setEnabled(True)
        else:
            self.next_button.setText("Review complete")
            self.next_button.setEnabled(False)
        self._update_review()

    def _go_next(self) -> None:
        index = self.step_stack.currentIndex()
        if index == 0 and not self._validate_inputs(show=True):
            return
        if index == 1 and not self._validate_runtime_inputs(show=True):
            return
        self._set_step(min(index + 1, 2))

    def _go_back(self) -> None:
        self._set_step(max(self.step_stack.currentIndex() - 1, 0))

    def _on_source_changed(self) -> None:
        mode = self.source_combo.currentData()
        is_genbank = mode in ("gbff", "gbff_records")
        self.genetic_code_combo.setEnabled(is_genbank)
        self.skip_pseudo_check.setEnabled(is_genbank)
        self.internal_stop_combo.setEnabled(is_genbank)
        self.keep_stop_check.setEnabled(is_genbank)
        if mode == "gbff":
            placeholder = "Folder with .gb/.gbff/.gbk files"
            input_label = "Input folder:"
            scan_label = "Scan folder"
            tip = ("Tip: one file per species is recommended. /translation is used "
                   "when present; otherwise the CDS is translated with the selected "
                   "genetic code.")
        elif mode == "gbff_records":
            placeholder = "Single GenBank file with multiple records"
            input_label = "Input file:"
            scan_label = "Scan file"
            tip = ("Tip: a single GenBank file can contain multiple records; each "
                   "record is treated as a separate taxon. /translation is used when "
                   "present; otherwise the CDS is translated with the selected "
                   "genetic code.")
        else:
            placeholder = "Folder with .fa/.faa/.fasta files"
            input_label = "Input folder:"
            scan_label = "Scan folder"
            tip = "Tip: one protein FASTA per species."
        self.input_edit.setPlaceholderText(placeholder)
        if hasattr(self, "input_label"):
            self.input_label.setText(input_label)
        if hasattr(self, "scan_button"):
            self.scan_button.setText(scan_label)
        if hasattr(self, "tip_label"):
            self.tip_label.setText(tip)
        if hasattr(self, "scan_label"):
            self.scan_label.setText(
                "No input selected." if not self.input_edit.text() else self.scan_label.text()
            )

    def _on_run_mode_changed(self) -> None:
        mode = self.run_mode_combo.currentData()
        self.existing_container.setVisible(mode == "existing")
        if mode == "parse_only":
            self.mode_hint.setText("Parse only writes proteomes + CDS + cds_map; no OrthoFinder needed.")
        elif mode == "existing":
            self.mode_hint.setText("Reuse an existing OrthoFinder result folder; CDS mapping needs cds_map.tsv in the output folder.")
        else:
            self.mode_hint.setText("Runs OrthoFinder with BLAST+ and writes single-copy orthogroups.")

    # ------------------------------------------------------------------
    # browse / validation
    # ------------------------------------------------------------------
    def _browse_input(self) -> None:
        if self.source_combo.currentData() == "gbff_records":
            path, _ = QFileDialog.getOpenFileName(
                self,
                "Select multi-record GenBank file",
                self.input_edit.text() or os.path.expanduser("~"),
                "GenBank (*.gb *.gbff *.gbk);;All files (*)",
            )
        else:
            path = QFileDialog.getExistingDirectory(
                self, "Select input folder", self.input_edit.text() or os.path.expanduser("~")
            )
        if path:
            self.input_edit.setText(path)
            if not self.output_edit.text():
                self.output_edit.setText(
                    os.path.join(self.workdir or tempfile.gettempdir(), f"orthofinder_{_timestamp()}")
                )
            self._scan_input()

    def _browse_output(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output folder", self.output_edit.text() or os.path.expanduser("~"))
        if path:
            self.output_edit.setText(path)

    def _browse_orthofinder(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select orthofinder executable", self.orthofinder_edit.text(), "Executable (*.exe);;All files (*)")
        if path:
            self.orthofinder_edit.setText(path)

    def _browse_blast(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select BLAST+/MCL directory", self.blast_edit.text() or os.path.expanduser("~"))
        if path:
            self.blast_edit.setText(path)

    def _browse_existing(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select existing OrthoFinder results", self.existing_edit.text() or os.path.expanduser("~"))
        if path:
            self.existing_edit.setText(path)

    def _input_suffixes(self) -> set[str]:
        if self.source_combo.currentData() == "proteome":
            return {".fa", ".faa", ".fasta", ".fas"}
        return {".gb", ".gbff", ".gbk", ".genbank"}

    def _scan_input(self) -> None:
        input_path = Path(self.input_edit.text())
        if not input_path.exists():
            self.scan_label.setText("⚠ Input path not found.")
            self.scan_label.setStyleSheet("color: #b00020;")
            return

        if self.source_combo.currentData() == "gbff_records":
            if not input_path.is_file() or input_path.suffix.lower() not in self._input_suffixes():
                self.scan_label.setText("⚠ Select a .gb/.gbff/.gbk file.")
                self.scan_label.setStyleSheet("color: #b00020;")
                return
            try:
                count = sum(1 for _ in SeqIO.parse(str(input_path), "genbank"))
            except Exception as exc:
                self.scan_label.setText(f"⚠ Could not read records: {exc}")
                self.scan_label.setStyleSheet("color: #b00020;")
                return
            self.scan_label.setText(
                f"✓ {input_path.name} · {count} record(s) (each record = taxon) · CDS counted on run."
            )
            self.scan_label.setStyleSheet("color: #2e7d32;")
            return

        if not input_path.is_dir():
            self.scan_label.setText("⚠ Select an input folder.")
            self.scan_label.setStyleSheet("color: #b00020;")
            return
        suffixes = self._input_suffixes()
        files = [p for p in input_path.iterdir() if p.is_file() and p.suffix.lower() in suffixes]
        species = {p.stem for p in files}
        if not files:
            self.scan_label.setText("⚠ No matching files in this folder.")
            self.scan_label.setStyleSheet("color: #b00020;")
        else:
            self.scan_label.setText(
                f"✓ {len(files)} file(s) · {len(species)} species · CDS counted on run."
            )
            self.scan_label.setStyleSheet("color: #2e7d32;")

    def _validate_inputs(self, show: bool = True) -> bool:
        input_path = Path(self.input_edit.text())
        mode = self.source_combo.currentData()
        if mode == "gbff_records":
            if not input_path.is_file() or input_path.suffix.lower() not in self._input_suffixes():
                if show:
                    QMessageBox.warning(self, "Input missing", "Please select a multi-record .gb/.gbff/.gbk file.")
                return False
        else:
            if not input_path.is_dir():
                if show:
                    QMessageBox.warning(self, "Input missing", "Please select an existing input folder.")
                return False
            if not any(p.is_file() and p.suffix.lower() in self._input_suffixes() for p in input_path.iterdir()):
                if show:
                    QMessageBox.warning(self, "No input files", "No matching input files were found in the folder.")
                return False
        if not self.output_edit.text().strip():
            if show:
                QMessageBox.warning(self, "Output missing", "Please choose an output folder.")
            return False
        return True

    def _check_runtime(self, show: bool = True) -> bool:
        orthofinder = Path(self.orthofinder_edit.text())
        blast_dir = Path(self.blast_edit.text()) if self.blast_edit.text() else None
        missing = []
        if not orthofinder.is_file():
            missing.append("orthofinder")
        for name in ("blastp.exe", "makeblastdb.exe", "mcl.exe"):
            candidates = []
            if blast_dir:
                candidates.append(blast_dir / name)
            candidates.append(orthofinder.parent / name)
            if not any(p.is_file() for p in candidates):
                missing.append(name)
        if missing:
            self.runtime_status.setText("⚠ Missing: " + ", ".join(missing))
            self.runtime_status.setStyleSheet("color: #b00020;")
            ok = False
        else:
            self.runtime_status.setText("✓ OrthoFinder, blastp, makeblastdb, mcl found.")
            self.runtime_status.setStyleSheet("color: #2e7d32;")
            ok = True
        if show and not ok:
            QMessageBox.warning(self, "Runtime incomplete", "Missing: " + ", ".join(missing))
        return ok

    def _validate_runtime_inputs(self, show: bool = True) -> bool:
        mode = self.run_mode_combo.currentData()
        if mode == "parse_only":
            return True
        if mode == "existing":
            root = Path(self.existing_edit.text()) if self.existing_edit.text() else None
            if not root or not root.exists():
                if show:
                    QMessageBox.warning(self, "Existing results missing", "Please select a folder containing Orthogroups.tsv.")
                return False
            try:
                find_orthogroups_dir(root)
            except FileNotFoundError:
                if show:
                    QMessageBox.warning(self, "Orthogroups.tsv missing", "No Orthogroups.tsv was found under the selected folder.")
                return False
            return True
        return self._check_runtime(show=show)

    def _update_review(self) -> None:
        if not hasattr(self, "review_label"):
            return
        mode = self.run_mode_combo.currentData() if hasattr(self, "run_mode_combo") else "full"
        source_mode = self.source_combo.currentData()
        if source_mode == "gbff":
            source = "GenBank folder (one file per species)"
        elif source_mode == "gbff_records":
            source = "GenBank multi-record file (each record = taxon)"
        else:
            source = "Proteome FASTA"
        self.review_label.setText(
            f"<b>Input</b>: {source} — {self.input_edit.text() or '<i>not selected</i>'}<br>"
            f"<b>Output</b>: {self.output_edit.text() or '<i>not selected</i>'}<br>"
            f"<b>Run mode</b>: {self.run_mode_combo.currentText() if hasattr(self, 'run_mode_combo') else mode}<br>"
            f"<b>Runtime</b>: {self.orthofinder_edit.text() or '<i>not set</i>'} · "
            f"{self.threads_spin.value() if hasattr(self, 'threads_spin') else 2} thread(s)"
        )

    # ------------------------------------------------------------------
    # output tab
    # ------------------------------------------------------------------
    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        self.output_status = QLabel("No run yet.")
        self.output_status.setWordWrap(True)
        self.output_status.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.output_status)

        cards = QHBoxLayout()
        cards.setSpacing(12)
        self.card_labels = {}
        for key, caption in (
            ("species", "Species"),
            ("cds", "CDS"),
            ("orthogroups", "Orthogroups"),
            ("single_copy_ogs", "Single-copy OGs"),
        ):
            card = QFrame()
            card.setObjectName("statCard")
            card.setFrameShape(QFrame.StyledPanel)
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(12, 8, 12, 8)
            number = QLabel("—")
            number.setAlignment(Qt.AlignCenter)
            number.setFont(_ui_font(18, bold=True))
            caption_label = QLabel(caption)
            caption_label.setAlignment(Qt.AlignCenter)
            caption_label.setStyleSheet("color: palette(mid);")
            card_layout.addWidget(number)
            card_layout.addWidget(caption_label)
            cards.addWidget(card, 1)
            self.card_labels[key] = number
        layout.addLayout(cards)

        table_label = QLabel("Single-copy orthogroups")
        table_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(table_label)
        self.og_table = QTableWidget(0, 4)
        self.og_table.setHorizontalHeaderLabels(["Orthogroup", "Species", "Proteins", "Status"])
        self.og_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.og_table.horizontalHeader().setDefaultAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.og_table.verticalHeader().setVisible(False)
        self.og_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.og_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        layout.addWidget(self.og_table, 1)

        self.output_info = QLabel("Output files will appear here.")
        self.output_info.setWordWrap(True)
        self.output_info.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.output_info)

        actions = QHBoxLayout()
        self.open_output_button = QPushButton("Open output folder")
        self.open_output_button.clicked.connect(self._open_output_folder)
        self.export_button = QPushButton("Export summary CSV")
        self.export_button.clicked.connect(self._export_summary)
        self.workspace_button = QPushButton("Export to workspace")
        self.workspace_button.clicked.connect(self._export_to_workspace)
        self.workspace_button.setEnabled(False)
        self.workspace_button.setToolTip(
            "Create a YR-MPEA dataset with one item per single-copy orthogroup. "
            "Available when this plugin is opened from YR-MPEA."
        )
        actions.addWidget(self.open_output_button)
        actions.addWidget(self.export_button)
        actions.addWidget(self.workspace_button)
        actions.addStretch()
        layout.addLayout(actions)

    def setup_control_panel(self):
        super().setup_control_panel()
        self.run_button.setText("Run pipeline")
        self.run_button.setMinimumWidth(160)
        self.stop_button.setMinimumWidth(100)
        self.run_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.stop_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.control_layout.insertStretch(0, 1)

    # ------------------------------------------------------------------
    # run lifecycle
    # ------------------------------------------------------------------
    def run_analysis(self):
        return self.start_pipeline(show_dialogs=True)

    def start_pipeline(self, show_dialogs=True):
        """Start the pipeline. Agents call this with show_dialogs=False."""
        self._show_run_dialogs = bool(show_dialogs)
        if not self._validate_inputs(show=show_dialogs):
            return False
        mode = self.run_mode_combo.currentData()
        if mode != "parse_only" and not self._validate_runtime_inputs(show=show_dialogs):
            return False

        output_dir = Path(self.output_edit.text())
        output_dir.mkdir(parents=True, exist_ok=True)
        source_mode = self.source_combo.currentData()
        default_table = self.genetic_code_combo.currentData() or 1

        self.is_running = True
        self._status_state = "running"
        self._started_at = time.time()
        self._completed_at = None
        self._last_error = ""
        self._last_summary = None
        self._progress_text = "Starting pipeline..."
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running pipeline...")
        self.tab_widget.setCurrentIndex(2)  # Console
        self.add_console_message("Starting OrthoFinder pipeline...", "info")

        try:
            self.analysis_thread = OrthoFinderThread(
                input_dir=self.input_edit.text(),
                output_dir=str(output_dir),
                source_mode=source_mode,
                default_table=default_table,
                skip_pseudo=self.skip_pseudo_check.isChecked(),
                internal_stop=self.internal_stop_combo.currentData(),
                keep_stop=self.keep_stop_check.isChecked(),
                run_mode=mode,
                orthofinder_path=self.orthofinder_edit.text(),
                blast_dir=self.blast_edit.text(),
                threads=self.threads_spin.value(),
                existing_root=self.existing_edit.text() or None,
            )
            self.analysis_thread.progress.connect(self._on_progress)
            self.analysis_thread.console_output.connect(self.add_console_message)
            self.analysis_thread.result_summary.connect(self._on_summary)
            self.analysis_thread.finished.connect(self.analysis_finished)
            self.analysis_thread.error.connect(self.analysis_error)
            self.analysis_thread.start()
            return True
        except Exception as exc:
            self.is_running = False
            self._status_state = "failed"
            self._last_error = str(exc)
            self._progress_text = "Failed to start"
            self.run_button.setEnabled(True)
            self.stop_button.setEnabled(False)
            self.progress_bar.setVisible(False)
            self.add_console_message("Pipeline failed to start: %s" % exc, "error")
            if show_dialogs:
                QMessageBox.critical(self, "Pipeline failed", str(exc))
            return False

    def _on_progress(self, message: str) -> None:
        self._progress_text = message
        self.progress_bar.setFormat(message)
        self.add_console_message(message, "info")

    def _on_summary(self, summary: dict) -> None:
        self._last_summary = summary

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self._status_state = "completed"
        self._completed_at = time.time()
        self._progress_text = "Completed"
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        if self._last_summary:
            self.update_results(self._last_summary)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("Pipeline completed.", "info")

    def analysis_error(self, error_message: str):
        self.is_running = False
        self._status_state = "failed"
        self._completed_at = time.time()
        self._last_error = error_message
        self._progress_text = "Failed"
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(f"Pipeline failed: {error_message}", "error")
        if self._show_run_dialogs:
            QMessageBox.critical(self, "Pipeline failed", error_message)

    def stop_analysis(self):
        if self.analysis_thread and self.analysis_thread.isRunning():
            self.analysis_thread.stop()
            self.analysis_thread.wait()
        self.is_running = False
        self._status_state = "stopped"
        self._completed_at = time.time()
        self._progress_text = "Stopped"
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Stopped by user.", "warning")

    def display_results(self, output_files):
        if self._last_summary:
            self.update_results(self._last_summary)

    def update_results(self, summary: dict):
        self._last_summary = summary
        self.output_status.setText(
            f"Pipeline complete — {summary.get('species', 0)} species · "
            f"{summary.get('cds', 0)} CDS · {summary.get('orthogroups', 0)} orthogroups · "
            f"{summary.get('single_copy_ogs', 0)} single-copy OGs "
            f"({summary.get('elapsed', 0)} s)"
        )
        for key, label in self.card_labels.items():
            label.setText(str(summary.get(key, "—")))
        output_dir = Path(summary.get("output_dir", ""))
        self.output_info.setText(f"Output: {output_dir}")
        summary_tsv = output_dir / "single_copy_ogs" / "single_copy_og_summary.tsv"
        self.og_table.setRowCount(0)
        if summary_tsv.exists():
            import csv

            with summary_tsv.open("r", encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.og_table.setRowCount(len(rows))
            for row_index, row in enumerate(rows):
                values = [row.get("orthogroup", ""), row.get("species_count", ""),
                          row.get("protein_count", ""), row.get("status", "")]
                for col, value in enumerate(values):
                    self.og_table.setItem(row_index, col, QTableWidgetItem(str(value)))
        elif summary.get("single_copy_ogs", 0) == 0:
            self.output_info.setText(
                f"Output: {output_dir}\nNo single-copy orthogroups were written; check the OrthoFinder log."
            )
        if hasattr(self, "workspace_button"):
            self.workspace_button.setEnabled(
                self._opened_from_yrmpea and summary.get("single_copy_ogs", 0) > 0
            )

    # ------------------------------------------------------------------
    # semantic API for agents / scripts
    # ------------------------------------------------------------------
    def _collect_options(self):
        return {
            "source_type": self.source_combo.currentData(),
            "input_path": self.input_edit.text(),
            "output_dir": self.output_edit.text(),
            "genetic_code": self.genetic_code_combo.currentData(),
            "skip_pseudogenes": self.skip_pseudo_check.isChecked(),
            "internal_stop": self.internal_stop_combo.currentData(),
            "keep_terminal_stop": self.keep_stop_check.isChecked(),
            "run_mode": self.run_mode_combo.currentData(),
            "existing_results": self.existing_edit.text(),
            "orthofinder_path": self.orthofinder_edit.text(),
            "blast_dir": self.blast_edit.text(),
            "threads": self.threads_spin.value(),
        }

    def inspect_options(self, as_json=False):
        payload = {
            "schema": "yrmpe.orthofinder.options.v1",
            "tool": "OrthoFinder",
            "options": self._collect_options(),
            "constraints": {
                "source_type": list(SOURCE_TYPES),
                "genetic_code": list(GENETIC_CODES),
                "internal_stop": list(INTERNAL_STOPS),
                "run_mode": list(RUN_MODES),
                "threads": {"min": self.threads_spin.minimum(),
                            "max": self.threads_spin.maximum()},
            },
        }
        if as_json:
            return json.dumps(payload, indent=2, ensure_ascii=False)
        return payload

    def _set_combo_data(self, combo, value):
        for index in range(combo.count()):
            if combo.itemData(index) == value:
                combo.setCurrentIndex(index)
                return
        raise ValueError("no combo entry for %r" % (value,))

    @staticmethod
    def _as_bool(value):
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)) and value in (0, 1):
            return bool(value)
        text = str(value).strip().lower()
        if text in ("true", "yes", "1", "on"):
            return True
        if text in ("false", "no", "0", "off"):
            return False
        raise ValueError("expected a boolean, got %r" % (value,))

    def modify_options(self, options):
        """Apply a dict/JSON option set atomically.

        Returns {ok, applied, rejected, errors, options}. A single invalid key
        rejects the whole patch so callers never get a half-applied state.
        """
        if isinstance(options, str):
            try:
                options = json.loads(options)
            except Exception as exc:
                return {"ok": False, "applied": [], "rejected": {},
                        "errors": ["invalid JSON: %s" % exc],
                        "options": self._collect_options()}
        if not isinstance(options, dict):
            return {"ok": False, "applied": [], "rejected": {},
                    "errors": ["options must be an object"],
                    "options": self._collect_options()}
        if isinstance(options.get("options"), dict):
            options = options["options"]

        if self.is_running:
            return {"ok": False, "applied": [], "rejected": options,
                    "errors": ["pipeline is running; stop it before changing options"],
                    "options": self._collect_options()}

        aliases = {
            "input": "input_path",
            "output": "output_dir",
            "transl_table": "genetic_code",
            "existing": "existing_results",
            "orthofinder_executable": "orthofinder_path",
            "blast": "blast_dir",
        }
        source_aliases = {
            "genbank": "gbff", "gb": "gbff", "genbank_folder": "gbff",
            "multi_record": "gbff_records", "multi-record": "gbff_records",
            "fasta": "proteome", "proteomes": "proteome",
        }
        run_aliases = {
            "parse-only": "parse_only", "parse": "parse_only",
            "existing_results": "existing", "existing_results_only": "existing",
        }

        normalized = {}
        errors = []
        rejected = {}
        for raw_key, raw_value in options.items():
            key = aliases.get(raw_key, raw_key)
            try:
                if key == "source_type":
                    value = source_aliases.get(str(raw_value).strip().lower(), str(raw_value).strip().lower())
                    if value not in SOURCE_TYPES:
                        raise ValueError("source_type must be one of %s" % (SOURCE_TYPES,))
                elif key in ("input_path", "output_dir", "existing_results",
                             "orthofinder_path", "blast_dir"):
                    value = str(raw_value)
                elif key == "genetic_code":
                    if raw_value in (None, "", "auto", "Auto", "None", "none"):
                        value = None
                    else:
                        value = int(raw_value)
                    if value not in GENETIC_CODES:
                        raise ValueError("genetic_code must be one of %s" % (GENETIC_CODES,))
                elif key in ("skip_pseudogenes", "keep_terminal_stop"):
                    value = self._as_bool(raw_value)
                elif key == "internal_stop":
                    value = str(raw_value).strip().lower()
                    if value not in INTERNAL_STOPS:
                        raise ValueError("internal_stop must be one of %s" % (INTERNAL_STOPS,))
                elif key == "run_mode":
                    value = run_aliases.get(str(raw_value).strip().lower(), str(raw_value).strip().lower())
                    if value not in RUN_MODES:
                        raise ValueError("run_mode must be one of %s" % (RUN_MODES,))
                elif key == "threads":
                    value = int(raw_value)
                    if not self.threads_spin.minimum() <= value <= self.threads_spin.maximum():
                        raise ValueError("threads must be between %d and %d"
                                         % (self.threads_spin.minimum(), self.threads_spin.maximum()))
                else:
                    raise KeyError("unknown option %r" % raw_key)
                normalized[key] = value
            except Exception as exc:
                errors.append(str(exc))
                rejected[raw_key] = str(exc)

        if errors:
            return {"ok": False, "applied": [], "rejected": rejected,
                    "errors": errors, "options": self._collect_options()}

        setters = {
            "source_type": lambda v: self._set_combo_data(self.source_combo, v),
            "input_path": lambda v: self.input_edit.setText(v),
            "output_dir": lambda v: self.output_edit.setText(v),
            "genetic_code": lambda v: self._set_combo_data(self.genetic_code_combo, v),
            "skip_pseudogenes": lambda v: self.skip_pseudo_check.setChecked(v),
            "internal_stop": lambda v: self._set_combo_data(self.internal_stop_combo, v),
            "keep_terminal_stop": lambda v: self.keep_stop_check.setChecked(v),
            "run_mode": lambda v: self._set_combo_data(self.run_mode_combo, v),
            "existing_results": lambda v: self.existing_edit.setText(v),
            "orthofinder_path": lambda v: self.orthofinder_edit.setText(v),
            "blast_dir": lambda v: self.blast_edit.setText(v),
            "threads": lambda v: self.threads_spin.setValue(v),
        }
        for key, value in normalized.items():
            setters[key](value)
        self._update_review()
        return {"ok": True, "applied": list(normalized), "rejected": {},
                "errors": [], "options": self._collect_options()}

    def inspect_status(self):
        running = bool(self.is_running or (self.analysis_thread and self.analysis_thread.isRunning()))
        state = "running" if running else self._status_state
        elapsed = None
        if running and self._started_at:
            elapsed = round(time.time() - self._started_at, 2)
        elif self._last_summary and self._last_summary.get("elapsed") is not None:
            elapsed = self._last_summary.get("elapsed")
        elif self._started_at and self._completed_at:
            elapsed = round(self._completed_at - self._started_at, 2)

        max_value = self.progress_bar.maximum()
        can_export = bool(
            self._opened_from_yrmpea
            and self.result_has_consumer()
            and not running
            and self._last_summary
            and self._last_summary.get("single_copy_ogs", 0) > 0
        )
        return {
            "schema": "yrmpe.orthofinder.status.v1",
            "state": state,
            "running": running,
            "wizard_step": self.step_stack.currentIndex(),
            "progress": {
                "text": self._progress_text,
                "visible": self.progress_bar.isVisible(),
                "indeterminate": max_value == 0 and self.progress_bar.minimum() == 0,
                "value": None if self.progress_bar.value() < 0 else self.progress_bar.value(),
            },
            "elapsed": elapsed,
            "summary": dict(self._last_summary or {}),
            "last_error": self._last_error,
            "console_lines": len(self.console_output),
            "can_run": not running,
            "can_stop": running,
            "can_export": can_export,
            "can_close": not running,
        }

    def inspect_console(self, since=0, limit=200, level=None):
        lines = self.console_output
        start = max(0, int(since))
        limit = max(1, min(int(limit), 1000))
        window = list(enumerate(lines[start:start + limit], start=start))
        selected = [
            (index, line) for index, line in window
            if not level or line.get("type") == level
        ]
        return {
            "schema": "yrmpe.orthofinder.console.v1",
            "total": len(lines),
            "since": start,
            "next": start + len(window),
            "lines": [
                {
                    "index": index,
                    "timestamp": line.get("timestamp", ""),
                    "type": line.get("type", ""),
                    "message": line.get("message", ""),
                }
                for index, line in selected
            ],
        }

    def preview_results(self):
        summary = dict(self._last_summary or {})
        output_dir = Path(summary.get("output_dir") or self.output_edit.text())
        files = []

        def add_file(role, path):
            path = Path(path)
            item = {"role": role, "path": str(path), "exists": path.exists()}
            if path.is_file() and path.exists():
                item["size"] = path.stat().st_size
            files.append(item)

        add_file("cds_map", output_dir / "cds_map.tsv")
        add_file("single_copy_summary", output_dir / "single_copy_ogs" / "single_copy_og_summary.tsv")
        for path in sorted((output_dir / "proteomes").glob("*.faa")):
            add_file("proteome", path)
        for path in sorted((output_dir / "cds").glob("*.fna")):
            add_file("cds", path)

        orthogroups = None
        orthofinder_root = output_dir / "orthofinder"
        if self.run_mode_combo.currentData() == "existing" and self.existing_edit.text():
            orthofinder_root = Path(self.existing_edit.text())
        try:
            orthogroups = str(find_orthogroups_dir(orthofinder_root) / "Orthogroups.tsv")
        except Exception:
            pass

        single_copy = []
        for path in sorted((output_dir / "single_copy_ogs").glob("*.cds.fna")):
            og = path.name[: -len(".cds.fna")]
            single_copy.append({
                "orthogroup": og,
                "cds_fasta": str(path),
                "protein_fasta": str(path.with_name(og + ".faa")),
                "map_tsv": str(path.with_name(og + ".map.tsv")),
            })

        return {
            "schema": "yrmpe.orthofinder.results.v1",
            "ready": bool(summary),
            "output_dir": str(output_dir),
            "summary": summary,
            "orthogroups_tsv": orthogroups,
            "single_copy_ogs": single_copy,
            "files": files,
        }

    def validate_configuration(self):
        input_ok = self._validate_inputs(show=False)
        runtime_ok = self._validate_runtime_inputs(show=False)
        return {
            "ok": bool(input_ok and runtime_ok),
            "input": bool(input_ok),
            "runtime": bool(runtime_ok),
            "status": self.inspect_status(),
        }

    def can_close(self):
        running = bool(self.is_running or (self.analysis_thread and self.analysis_thread.isRunning()))
        return {
            "can_close": not running,
            "running": running,
            "reason": "pipeline is running; stop it before closing" if running else "",
        }

    def inspect_actions(self):
        status = self.inspect_status()
        can_export = status["can_export"]
        output_ready = bool(self.output_edit.text()) and Path(self.output_edit.text()).is_dir()
        actions = [
            {"id": "orthofinder.run", "label": "Run pipeline", "enabled": status["can_run"],
             "reason": "" if status["can_run"] else "pipeline already running",
             "parameters": {"show_dialogs": {"type": "boolean", "required": False}}},
            {"id": "orthofinder.stop", "label": "Stop pipeline", "enabled": status["can_stop"],
             "reason": "" if status["can_stop"] else "pipeline is not running"},
            {"id": "orthofinder.validate", "label": "Validate inputs/runtime",
             "enabled": status["can_run"]},
            {"id": "orthofinder.modify_options", "label": "Modify options",
             "enabled": status["can_run"],
             "parameters": {"options": {"type": "object", "required": True}}},
            {"id": "orthofinder.preview_results", "label": "Preview result files",
             "enabled": True},
            {"id": "orthofinder.export_to_workspace", "label": "Import to YR-MPEA dataset",
             "enabled": can_export,
             "reason": "" if can_export else "run must complete and the plugin must be opened from YR-MPEA"},
            {"id": "orthofinder.open_output", "label": "Open output folder",
             "enabled": output_ready,
             "reason": "" if output_ready else "output folder does not exist"},
            {"id": "orthofinder.set_step", "label": "Switch wizard step",
             "enabled": status["can_run"],
             "parameters": {"step": {"type": "integer", "minimum": 0, "maximum": 2}}},
            {"id": "orthofinder.close", "label": "Close window",
             "enabled": status["can_close"],
             "reason": "" if status["can_close"] else "stop the pipeline first"},
        ]
        return {
            "schema": "yrmpe.orthofinder.actions.v1",
            "actions": actions,
            "next_import": {
                "action": "orthofinder.export_to_workspace",
                "target": "YR-MPEA dataset",
                "enabled": can_export,
            },
        }

    def agent_contract(self):
        return {
            "schema": "yrmpe.ui.orthofinder.v1",
            "id": "one_step_orthofinder",
            "options": self.inspect_options(),
            "status": self.inspect_status(),
            "actions": self.inspect_actions(),
        }

    def agent_invoke(self, action_id, params=None):
        params = params or {}
        name = action_id[len("orthofinder."):] if action_id.startswith("orthofinder.") else action_id
        try:
            if name == "inspect_options":
                return {"ok": True, "result": self.inspect_options()}
            if name == "inspect_status":
                return {"ok": True, "result": self.inspect_status()}
            if name == "inspect_console":
                return {"ok": True, "result": self.inspect_console(
                    since=params.get("since", 0),
                    limit=params.get("limit", 200),
                    level=params.get("level"),
                )}
            if name == "inspect_actions":
                return {"ok": True, "result": self.inspect_actions()}
            if name == "preview_results":
                return {"ok": True, "result": self.preview_results()}
            if name == "modify_options":
                return self.modify_options(params.get("options", params))
            if name == "validate":
                return {"ok": True, "result": self.validate_configuration()}
            if name == "run":
                started = self.start_pipeline(show_dialogs=bool(params.get("show_dialogs", False)))
                outcome = {"ok": started, "status": self.inspect_status()}
                if not started:
                    outcome["error"] = "input/runtime validation failed; call orthofinder.validate for details"
                return outcome
            if name == "stop":
                self.stop_analysis()
                return {"ok": True, "status": self.inspect_status()}
            if name == "export_to_workspace":
                outcome = self.export_to_workspace()
                outcome["status"] = self.inspect_status()
                return outcome
            if name == "open_output":
                folder = self.output_edit.text()
                if not folder or not Path(folder).is_dir():
                    return {"ok": False, "error": "output folder does not exist"}
                self._open_output_folder()
                return {"ok": True, "path": folder}
            if name == "set_step":
                self._set_step(int(params["step"]))
                return {"ok": True, "status": self.inspect_status()}
            if name == "close":
                state = self.can_close()
                if not state["can_close"]:
                    return {"ok": False, "error": state["reason"]}
                self.close()
                return {"ok": True, "closed": True}
            raise KeyError("unknown action %r" % action_id)
        except Exception as exc:
            return {"ok": False, "error": str(exc)}

    def closeEvent(self, event):
        state = self.can_close()
        if not state["can_close"]:
            reply = QMessageBox.question(
                self,
                "Pipeline running",
                "The OrthoFinder pipeline is still running. Stop it and close?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                event.ignore()
                return
            self.stop_analysis()
        super().closeEvent(event)

    def _build_result_spec(self):
        """Build the dataset exported to YR-MPEA (one item per OG)."""
        summary = self._last_summary
        if not summary:
            raise ValueError("Run the pipeline before exporting to the workspace.")

        output_dir = Path(summary.get("output_dir", ""))
        single_dir = output_dir / "single_copy_ogs"
        if not single_dir.is_dir():
            raise ValueError("No single_copy_ogs folder found under %s" % output_dir)

        items = []
        for cds_path in sorted(single_dir.glob("*.cds.fna")):
            og = cds_path.name[: -len(".cds.fna")]
            records = list(SeqIO.parse(str(cds_path), "fasta"))
            if not records:
                continue
            map_path = single_dir / (og + ".map.tsv")
            faa_path = single_dir / (og + ".faa")
            items.append(
                DatasetItemSpec(
                    loci_name=og,
                    item_type="sequence",
                    sequences=records,
                    file_path=str(cds_path),
                    format="fasta",
                    is_aligned=False,
                    data={
                        "orthogroup": og,
                        "protein_fasta": str(faa_path),
                        "map_tsv": str(map_path),
                    },
                    metadata={
                        "source": self.plugin_name,
                        "output_dir": str(output_dir),
                    },
                )
            )

        if not items:
            raise ValueError(
                "No single-copy CDS FASTA files were found under %s." % single_dir
            )

        genome_label = output_dir.name or "dataset"
        dataset = DatasetSpec(
            name="OrthoFinder_%s" % genome_label,
            items=items,
            description=(
                "%d species, %d CDS, %d orthogroups, %d single-copy orthogroups"
                % (
                    summary.get("species", 0),
                    summary.get("cds", 0),
                    summary.get("orthogroups", 0),
                    summary.get("single_copy_ogs", 0),
                )
            ),
            source="one_step_orthofinder",
            summary={
                key: summary.get(key)
                for key in ("species", "cds", "orthogroups", "single_copy_ogs", "output_dir")
            },
            metadata={
                "plugin": self.plugin_name,
                "input": self.input_edit.text(),
                "source_type": self.source_combo.currentText(),
                "genetic_code": self.genetic_code_combo.currentText(),
                "output": str(output_dir),
            },
        )
        artifacts = [str(single_dir / "single_copy_og_summary.tsv")]
        return ResultSpec(
            source="one_step_orthofinder",
            datasets=[dataset],
            artifacts=artifacts,
            summary=dict(summary),
        )

    def export_to_workspace(self):
        """Emit the dataset result. Returns a JSON-safe status dict."""
        if not self._opened_from_yrmpea or not self.result_has_consumer():
            return {
                "ok": False,
                "error": "Not connected to YR-MPEA; open this plugin from the workspace.",
            }
        try:
            result = self._build_result_spec()
        except Exception as exc:
            return {"ok": False, "error": str(exc)}

        dataset = result.datasets[0]
        if not self.export_result(result):
            return {"ok": False, "error": "Result validation failed; see the console."}

        self.add_console_message(
            "Sent dataset %s to the YR-MPEA workspace (%d loci)."
            % (dataset.name, len(dataset.items)),
            "info",
        )
        return {"ok": True, "dataset": dataset.name, "loci": len(dataset.items)}

    def _export_to_workspace(self):
        """UI wrapper for export_to_workspace()."""
        outcome = self.export_to_workspace()
        if not outcome.get("ok"):
            QMessageBox.warning(self, "Export failed", outcome.get("error", "Unknown error"))
            return
        QMessageBox.information(
            self,
            "Dataset exported",
            "Dataset: %s\nLoci: %d\n\nThe dataset is now in the YR-MPEA workspace."
            % (outcome["dataset"], outcome["loci"]),
        )

    def _open_output_folder(self):
        folder = self.output_edit.text()
        if folder and Path(folder).exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(folder))

    def _export_summary(self):
        source = Path(self.output_edit.text()) / "single_copy_ogs" / "single_copy_og_summary.tsv"
        if not source.exists():
            QMessageBox.warning(self, "No summary", "Run the pipeline first; no summary file found.")
            return
        target, _ = QFileDialog.getSaveFileName(self, "Export summary CSV", str(source.with_suffix(".csv")), "CSV (*.csv)")
        if target:
            shutil.copy2(source, target)
            self.add_console_message(f"Summary exported to {target}", "info")


class OrthoFinderPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return OrthoFinderPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
