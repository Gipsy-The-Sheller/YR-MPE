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
import os
import shutil
import tempfile
import time
from pathlib import Path

from PyQt5.QtCore import Qt, QUrl, pyqtSignal
from PyQt5.QtGui import QDesktopServices, QFont
from PyQt5.QtWidgets import (
    QAbstractItemView,
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
from .components.methods.gbff_utils import parse_folder, write_dataset
from .components.methods.orthogroup_utils import (
    derive_single_copy_ids,
    find_orthogroups_dir,
    read_cds_map,
    read_orthogroups,
    read_single_copy_ids,
    write_single_copy_ogs,
)


def _timestamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


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
                "cds_mapping": self.source_mode == "gbff",
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
        super().__init__(import_from, import_data, **kwargs)
        self.run_button.setText("Run pipeline")

    def init_plugin_info(self):
        self.plugin_name = "Genomic Dataset Formation (OrthoFinder)"
        self.tool_name = "OrthoFinder"
        self.citation = [
            "Emms, D. M. & Kelly, S. (2019) OrthoFinder: phylogenetic orthology "
            "inference for comparative genomics. Genome Biology 20:238. "
            "doi: <a href='https://doi.org/10.1186/s13059-019-1816-6'>"
            "10.1186/s13059-019-1816-6</a>",
            "Camacho, C. et al. (2009) BLAST+: architecture and applications. "
            "BMC Bioinformatics 10:421. "
            "doi: <a href='https://doi.org/10.1186/1471-2105-10-421'>"
            "10.1186/1471-2105-10-421</a>",
            "van Dongen, S. (2000) Graph clustering by flow simulation. "
            "PhD thesis, University of Utrecht (MCL).",
            "Cock, P. J. A. et al. (2009) Biopython: freely available Python tools "
            "for computational molecular biology and bioinformatics. "
            "Bioinformatics 25:1422-1423. "
            "doi: <a href='https://doi.org/10.1093/bioinformatics/btp163'>"
            "10.1093/bioinformatics/btp163</a>",
        ]
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

        title = QLabel("Genomic Dataset Formation (OrthoFinder)")
        title_font = QFont()
        title_font.setPointSize(15)
        title_font.setBold(True)
        title.setFont(title_font)
        layout.addWidget(title)

        subtitle = QLabel("GBFF → CDS/AA → OrthoFinder → single-copy orthogroups → CDS")
        subtitle.setStyleSheet("color: #555;")
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
        self.source_combo.addItem("GenBank (GBFF/GB/GBK)", "gbff")
        self.source_combo.addItem("Proteome FASTA (protein-only)", "proteome")
        self.source_combo.currentIndexChanged.connect(self._on_source_changed)
        source_form.addRow("Source type:", self.source_combo)

        self.input_edit = QLineEdit()
        self.input_edit.setPlaceholderText("Folder with .gb/.gbff/.gbk files")
        self.input_browse_button = QPushButton("Browse")
        self.input_browse_button.clicked.connect(self._browse_input)
        source_form.addRow("Input folder:", _path_row(self.input_edit, self.input_browse_button))
        source_layout.addLayout(source_form)

        self.scan_label = QLabel("No folder selected.")
        self.scan_label.setWordWrap(True)
        self.scan_label.setStyleSheet("color: #666;")
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
        mode_form.addRow("Existing results:", self.existing_row)
        mode_layout.addLayout(mode_form)

        self.mode_hint = QLabel(
            "Parse only writes proteomes + CDS + cds_map and does not require OrthoFinder."
        )
        self.mode_hint.setWordWrap(True)
        self.mode_hint.setStyleSheet("color: #666;")
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
        note.setStyleSheet("color: #555;")
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
                label.setStyleSheet("color: #777;")
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
        is_gbff = self.source_combo.currentData() == "gbff"
        self.genetic_code_combo.setEnabled(is_gbff)
        self.skip_pseudo_check.setEnabled(is_gbff)
        self.internal_stop_combo.setEnabled(is_gbff)
        self.keep_stop_check.setEnabled(is_gbff)
        self.input_edit.setPlaceholderText(
            "Folder with .gb/.gbff/.gbk files" if is_gbff else "Folder with .fa/.faa/.fasta files"
        )
        self.scan_label.setText("No folder selected." if not self.input_edit.text() else self.scan_label.text())

    def _on_run_mode_changed(self) -> None:
        mode = self.run_mode_combo.currentData()
        self.existing_row.setVisible(mode == "existing")
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
        path = QFileDialog.getExistingDirectory(self, "Select input folder", self.input_edit.text() or os.path.expanduser("~"))
        if path:
            self.input_edit.setText(path)
            if not self.output_edit.text():
                self.output_edit.setText(os.path.join(self.workdir or tempfile.gettempdir(), f"orthofinder_{_timestamp()}"))
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
        folder = Path(self.input_edit.text())
        if not folder.exists():
            self.scan_label.setText("⚠ Folder not found.")
            self.scan_label.setStyleSheet("color: #b00020;")
            return
        suffixes = self._input_suffixes()
        files = [p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in suffixes]
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
        folder = Path(self.input_edit.text())
        if not folder.is_dir():
            if show:
                QMessageBox.warning(self, "Input missing", "Please select an existing input folder.")
            return False
        if not any(p.is_file() and p.suffix.lower() in self._input_suffixes() for p in folder.iterdir()):
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
        source = "GenBank (GBFF)" if self.source_combo.currentData() == "gbff" else "Proteome FASTA"
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
            card.setFrameShape(QFrame.StyledPanel)
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(12, 8, 12, 8)
            number = QLabel("—")
            number.setAlignment(Qt.AlignCenter)
            font = QFont()
            font.setPointSize(18)
            font.setBold(True)
            number.setFont(font)
            caption_label = QLabel(caption)
            caption_label.setAlignment(Qt.AlignCenter)
            caption_label.setStyleSheet("color: #555;")
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
        self.og_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.og_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        layout.addWidget(self.og_table, 1)

        self.output_info = QLabel("Output files will appear here.")
        self.output_info.setWordWrap(True)
        self.output_info.setStyleSheet("color: #555;")
        layout.addWidget(self.output_info)

        actions = QHBoxLayout()
        self.open_output_button = QPushButton("Open output folder")
        self.open_output_button.clicked.connect(self._open_output_folder)
        self.export_button = QPushButton("Export summary CSV")
        self.export_button.clicked.connect(self._export_summary)
        actions.addWidget(self.open_output_button)
        actions.addWidget(self.export_button)
        actions.addStretch()
        layout.addLayout(actions)

    def setup_control_panel(self):
        super().setup_control_panel()
        self.run_button.setText("Run pipeline")

    # ------------------------------------------------------------------
    # run lifecycle
    # ------------------------------------------------------------------
    def run_analysis(self):
        if not self._validate_inputs(show=True):
            return
        mode = self.run_mode_combo.currentData()
        if mode != "parse_only" and not self._validate_runtime_inputs(show=True):
            return

        output_dir = Path(self.output_edit.text())
        output_dir.mkdir(parents=True, exist_ok=True)
        source_mode = self.source_combo.currentData()
        default_table = self.genetic_code_combo.currentData() or 1

        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running pipeline...")
        self.tab_widget.setCurrentIndex(2)  # Console
        self.add_console_message("Starting OrthoFinder pipeline...", "info")

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
            workdir=str(output_dir),
        )
        self.analysis_thread.progress.connect(self._on_progress)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.result_summary.connect(self._on_summary)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.start()

    def _on_progress(self, message: str) -> None:
        self.progress_bar.setFormat(message)
        self.add_console_message(message, "info")

    def _on_summary(self, summary: dict) -> None:
        self._last_summary = summary

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        if self._last_summary:
            self.update_results(self._last_summary)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("Pipeline completed.", "info")

    def analysis_error(self, error_message: str):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(f"Pipeline failed: {error_message}", "error")
        QMessageBox.critical(self, "Pipeline failed", error_message)

    def stop_analysis(self):
        if self.analysis_thread and self.analysis_thread.isRunning():
            self.analysis_thread.stop()
            self.analysis_thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Stopped by user.", "warning")

    def display_results(self, output_files):
        if self._last_summary:
            self.update_results(self._last_summary)

    def update_results(self, summary: dict):
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
