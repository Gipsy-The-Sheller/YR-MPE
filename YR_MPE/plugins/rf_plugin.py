# rf_plugin.py
#
# Copyright (c) 2026 Zhi-Jie Xu
#
# This file is part of YR-MPE.
#
# YR-MPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Pairwise Robinson-Foulds distances between a set of trees.

IQ-TREE 3 computes the all-to-all split distance with
``-t trees --tree-dist-all`` and writes ``<prefix>.rfdist``.  This plugin
shows the distribution as a box plot (splitting within-tree-set from
between-tree-set pairs when the input carries tree-set metadata) and the
pairwise distances as a lower triangular matrix.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .components.stdout_summary import SUMMARY_PREFIX, emit_summary
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread


# ----------------------------------------------------------------------
# pure helpers (tested without Qt)
# ----------------------------------------------------------------------
def parse_rfdist(path: str) -> Tuple[List[str], List[List[int]]]:
    """Parse an IQ-TREE ``.rfdist`` file into labels and a square matrix."""
    lines = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if text and not text.startswith("#"):
                lines.append(text)
    if not lines:
        raise ValueError("The RF distance file is empty.")

    header = lines[0].split()
    if len(header) < 2:
        raise ValueError("Malformed RF distance header: %r" % lines[0])
    count = int(header[0])
    if len(lines) < count + 1:
        raise ValueError(
            "Expected %d distance rows, found %d." % (count, len(lines) - 1)
        )

    labels: List[str] = []
    matrix: List[List[int]] = []
    for line in lines[1:count + 1]:
        parts = line.split()
        if len(parts) < count + 1:
            raise ValueError("Malformed RF distance row: %r" % line)
        labels.append(parts[0])
        matrix.append([int(round(float(value))) for value in parts[1:count + 1]])

    if len(matrix) != count:
        raise ValueError("RF distance matrix is not square.")
    return labels, matrix


def pairwise_rows(labels: Sequence[str], matrix: Sequence[Sequence[int]]):
    """Return ``(label_i, label_j, rf)`` for the lower triangle."""
    rows = []
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            rows.append((labels[i], labels[j], int(matrix[i][j])))
    return rows


def classify_pairs(groups: Sequence[Optional[dict]]) -> List[str]:
    """Label each pair as within-tree-set, between-tree-sets, or other.

    ``groups`` has one entry per tree: ``None`` for a plain tree file, or a
    dict with ``group_id``, ``group_name`` and ``kind`` (``treeset`` or
    ``tree``).  Tree sets selected from the workspace provide the metadata;
    a file opened from disk has no grouping and falls back to ``All pairs``.
    """
    has_treeset = any(
        group and group.get("kind") == "treeset" for group in groups
    )
    categories = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            first, second = groups[i], groups[j]
            if not has_treeset:
                categories.append("All pairs")
                continue
            first_set = first and first.get("kind") == "treeset"
            second_set = second and second.get("kind") == "treeset"
            if first_set and second_set:
                if first.get("group_id") == second.get("group_id"):
                    categories.append("Within tree set")
                else:
                    categories.append("Between tree sets")
            else:
                categories.append("Between objects")
    return categories


def normalized_rf(rf: int, ntaxa: int) -> Optional[float]:
    """Normalise RF by the maximum split distance for binary trees.

    ``2 * (n - 3)`` is the maximum RF between two fully resolved unrooted
    trees with ``n >= 4`` taxa.  Three-taxon trees are trifurcating and
    differ by at most two splits; two-taxon trees cannot differ.
    """
    if ntaxa >= 4:
        maximum = 2 * (ntaxa - 3)
    elif ntaxa == 3:
        maximum = 2
    else:
        maximum = 0
    if not maximum:
        return None
    return rf / float(maximum)


def count_newick_trees(path: str) -> int:
    """Count non-comment, non-empty lines (one Newick tree per line)."""
    count = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                text = line.strip()
                if text and not text.startswith("#"):
                    count += 1
    except OSError:
        return 0
    return count


def _unroot_tree(tree) -> bool:
    """Convert a rooted binary Newick tree to its unrooted representation.

    IQ-TREE treats a tree whose root has two children as rooted and a tree
    whose root has three or more children as unrooted.  Mixing the two in
    ``--tree-dist-all`` makes IQ-TREE report a sentinel distance (12 for a
    four-taxon test), so every tree is normalised before the comparison.
    """
    root = tree.root
    if len(root.clades) != 2:
        return False
    first, second = root.clades
    if len(second.clades) >= 1:
        keep, move = second, first
    elif len(first.clades) >= 1:
        keep, move = first, second
    else:
        return False
    keep.clades.append(move)
    tree.root = keep
    keep.branch_length = None
    return True


def normalize_tree_file(source: str, target: str):
    """Write an all-unrooted copy of a Newick file.

    Returns a small summary dict.  Lines that Biopython cannot parse are
    passed through unchanged and counted in ``failed``.
    """
    from io import StringIO

    from Bio import Phylo

    summary = {"parsed": 0, "unrooted": 0, "failed": 0}
    with open(source, "r", encoding="utf-8", errors="replace") as source_handle:
        with open(target, "w", encoding="utf-8", newline="") as target_handle:
            for line in source_handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                try:
                    tree = Phylo.read(StringIO(text), "newick")
                except Exception:
                    summary["failed"] += 1
                    target_handle.write(text + "\n")
                    continue
                summary["parsed"] += 1
                if _unroot_tree(tree):
                    summary["unrooted"] += 1
                buffer = StringIO()
                Phylo.write(tree, buffer, "newick")
                target_handle.write(buffer.getvalue().strip() + "\n")
    return summary


def count_taxa(path: str) -> int:
    """Read the first parseable Newick line and return its terminal count."""
    try:
        from io import StringIO

        from Bio import Phylo

        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                try:
                    tree = Phylo.read(StringIO(text), "newick")
                except Exception:
                    continue
                return len(tree.get_terminals())
    except Exception:
        pass
    return 0


def _short_labels(labels: Sequence[str]) -> List[str]:
    return [
        label if len(label) <= 24 else label[:21] + "..."
        for label in labels
    ]


# ----------------------------------------------------------------------
# thread
# ----------------------------------------------------------------------
class RFComparisonThread(BaseProcessThread):
    """Run IQ-TREE ``--tree-dist-all`` and parse the resulting matrix."""

    result_ready = pyqtSignal(object)

    def __init__(self, iqtree_path: str, tree_file: str, labels: Sequence[str],
                 output_dir: str, threads: int = 2, workdir: Optional[str] = None,
                 groups: Optional[Sequence[Optional[dict]]] = None):
        super().__init__(iqtree_path, [tree_file], [], workdir=workdir)
        self.iqtree_path = Path(iqtree_path)
        self.tree_file = Path(tree_file)
        self.labels = list(labels)
        self.groups = list(groups or [])
        self.output_dir = Path(output_dir)
        self.threads = max(1, int(threads))

    def get_tool_name(self) -> str:
        return "IQ-TREE 3"

    def execute_commands(self):
        if not self.iqtree_path.is_file():
            raise FileNotFoundError(
                "IQ-TREE executable is not configured (config.json → IQ-TREE 3)."
            )
        self.output_dir.mkdir(parents=True, exist_ok=True)
        normalized_path = self.output_dir / "rf_input_unrooted.nwk"
        summary = normalize_tree_file(str(self.tree_file), str(normalized_path))
        if summary["unrooted"]:
            self.console_output.emit(
                "Converted %d rooted tree(s) to unrooted topologies."
                % summary["unrooted"],
                "info",
            )
        if summary["failed"]:
            self.console_output.emit(
                "Warning: %d tree line(s) could not be parsed by Biopython; "
                "kept unchanged." % summary["failed"],
                "warning",
            )
        tree_input = normalized_path if summary["parsed"] else self.tree_file
        prefix = self.output_dir / "rf"
        cmd = [
            str(self.iqtree_path),
            "-t", str(tree_input),
            "--tree-dist-all",
            "-T", str(self.threads),
            "--prefix", str(prefix),
            "--redo",
        ]
        self.progress.emit("Computing pairwise Robinson-Foulds distances...")
        result = self.execute_command(cmd, cwd=str(self.output_dir))

        rfdist = Path(str(prefix) + ".rfdist")
        if not rfdist.is_file():
            tail = ((result.stdout or "") + "\n" + (result.stderr or "")).strip()
            if len(tail) > 1500:
                tail = tail[-1500:]
            raise RuntimeError(
                "IQ-TREE did not write %s.\n%s" % (rfdist, tail)
            )

        file_labels, matrix = parse_rfdist(str(rfdist))
        labels = list(self.labels)
        if len(labels) != len(matrix):
            labels = file_labels

        ntaxa = count_taxa(str(tree_input))
        groups = list(self.groups)
        if len(groups) != len(matrix):
            groups = []
        payload = {
            "labels": labels,
            "groups": groups,
            "matrix": matrix,
            "ntaxa": ntaxa,
            "rfdist": str(rfdist),
            "output_dir": str(self.output_dir),
        }
        self.result_ready.emit(payload)
        self.finished.emit([str(rfdist)], [])


# ----------------------------------------------------------------------
# plugin
# ----------------------------------------------------------------------
class RFComparisonPlugin(BasePlugin):
    """Pairwise Robinson-Foulds comparison of a tree set."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self._import_payload = import_data
        self._import_path = ""
        self._labels: List[str] = []
        self._groups: List[dict] = []
        self._payload = None
        self._out_dir = ""
        self.analysis_thread = None

        super().__init__(import_from, import_data, **kwargs)

        self.run_button.setText("Compute RF distances")

    def init_plugin_info(self):
        self.plugin_name = "Pairwise Robinson-Foulds Comparison"
        self.tool_name = "IQ-TREE 3"
        self.citation = []
        self.citation_keys = ["wong2026iqtree", "robinson1981comparison"]
        self.input_types = {
            "Newick tree file": ["nwk", "newick", "tree", "trees", "tre"]
        }
        self.output_types = {"RF matrix": ["tsv", "png"]}
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    def config(self):
        # The inherited check locates IQ-TREE 3 in config.json.  It is not
        # required until Run; keeping the dialog open lets the user inspect
        # help even when the runtime is missing.
        return True

    def control_bindings(self):
        """GUI controls mapped to the IQ-TREE tree-distance command."""
        return [
            {"id": "input", "widget": "tree_edit",
             "label": "Tree file / TreeSet", "kind": "path",
             "cli": ["-t", "{value}"],
             "notes": "One Newick tree per line; a workspace TreeSet also "
                      "populates this field."},
            {"id": "threads", "widget": "threads_spin",
             "label": "Threads", "kind": "int",
             "cli": ["-T", "{value}"],
             "notes": "IQ-TREE --tree-dist-all thread count."},
        ]

    # ------------------------------------------------------------------
    # input tab
    # ------------------------------------------------------------------
    def setup_input_tab(self):
        outer = QVBoxLayout(self.input_tab)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        outer.addWidget(scroll)
        content = QWidget()
        scroll.setWidget(content)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        title = QLabel("Tree comparison")
        title.setStyleSheet("font-size: 18px; font-weight: 600;")
        layout.addWidget(title)
        subtitle = QLabel(
            "Compute all pairwise Robinson-Foulds distances with IQ-TREE 3 "
            "(--tree-dist-all). The result is a distance distribution and a "
            "lower-triangular RF matrix."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        layout.addWidget(subtitle)

        trees_group = QGroupBox("Trees")
        form = QFormLayout(trees_group)
        self.tree_edit = QLineEdit()
        self.tree_edit.setPlaceholderText("Newick file, one tree per line")
        browse = QPushButton("Browse")
        browse.clicked.connect(self._browse_tree_file)
        row = QHBoxLayout()
        row.addWidget(self.tree_edit, 1)
        row.addWidget(browse)
        form.addRow("Tree file:", row)
        self.tree_info = QLabel("No tree file selected.")
        self.tree_info.setWordWrap(True)
        self.tree_info.setStyleSheet("color: palette(mid);")
        form.addRow("", self.tree_info)
        layout.addWidget(trees_group)

        engine_group = QGroupBox("IQ-TREE")
        eform = QFormLayout(engine_group)
        iqtree = self._iqtree_path()
        self.iqtree_label = QLabel(iqtree or "Not configured — set the IQ-TREE 3 path in config.json.")
        self.iqtree_label.setWordWrap(True)
        self.iqtree_label.setStyleSheet(
            "color: palette(mid);" if iqtree else "color: palette(mid);"
        )
        eform.addRow("Executable:", self.iqtree_label)
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, max(1, os.cpu_count() or 1))
        self.threads_spin.setValue(min(4, max(1, (os.cpu_count() or 2) // 2)))
        self.threads_spin.setToolTip("Threads passed to IQ-TREE with -T")
        eform.addRow("Threads:", self.threads_spin)
        layout.addWidget(engine_group)

        note = QLabel(
            "RF counts split differences and is unnormalised. Normalised RF "
            "is RF / 2(n−3) for n ≥ 4 taxa. If TreeSets are selected, the "
            "overview separates within-tree-set from between-tree-set pairs. "
            "Rooted trees are converted to unrooted topologies before the "
            "comparison. Trees should share the same taxon set; IQ-TREE "
            "handles missing taxa but the distance is then less directly "
            "comparable."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color: palette(mid);")
        layout.addWidget(note)
        layout.addStretch()

        self.tree_edit.editingFinished.connect(self._on_tree_path_edited)
        self._apply_import_payload()

    def _apply_import_payload(self):
        payload = self._import_payload
        path = ""
        labels = []
        groups = []
        if isinstance(payload, dict):
            path = str(payload.get("tree_file") or payload.get("path") or "")
            labels = list(payload.get("labels") or [])
            groups = list(payload.get("groups") or [])
        elif isinstance(payload, str):
            path = payload
        if not path:
            return
        self._import_path = path
        self._labels = [str(label) for label in labels]
        self._groups = [dict(group) for group in groups if isinstance(group, dict)]
        self.tree_edit.setText(path)
        self.tree_info.setText("Imported from Selected.")
        self._update_tree_info()

    def _browse_tree_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select a Newick tree file",
            "",
            "Tree files (*.nwk *.newick *.tree *.trees *.tre);;All files (*)",
        )
        if not path:
            return
        self._import_path = ""
        self._labels = []
        self._groups = []
        self.tree_edit.setText(path)
        self._update_tree_info()

    def _on_tree_path_edited(self):
        path = self.tree_edit.text().strip()
        if path != self._import_path:
            self._labels = []
            self._groups = []
        self._update_tree_info()

    def _tree_labels(self, count: int) -> List[str]:
        if len(self._labels) == count:
            return list(self._labels)
        return ["Tree %d" % (index + 1) for index in range(count)]

    def _tree_groups(self, count: int) -> List[Optional[dict]]:
        if len(self._groups) == count:
            return list(self._groups)
        return []

    def _update_tree_info(self):
        path = self.tree_edit.text().strip()
        if not path:
            self.tree_info.setText("No tree file selected.")
            return
        if not Path(path).is_file():
            self.tree_info.setText("File not found.")
            return
        count = count_newick_trees(path)
        if count < 2:
            self.tree_info.setText(
                "%d tree found. RF comparison needs at least two trees." % count
            )
            return
        labels = self._tree_labels(count)
        preview = ", ".join(labels[:4])
        if count > 4:
            preview += ", ..."
        self.tree_info.setText("%d trees: %s" % (count, preview))

    def _iqtree_path(self):
        config = self._find_plugin_config()
        if not config:
            return ""
        raw = Path(str(config.get("path", "")))
        if not raw:
            return ""
        if not raw.is_absolute():
            raw = Path(self.plugin_path) / str(config.get("path", "")).lstrip("/\\")
        return str(raw) if raw.is_file() else ""

    # ------------------------------------------------------------------
    # output tab
    # ------------------------------------------------------------------
    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.summary_label = QLabel("No RF comparison yet.")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.summary_label)

        self.view_tabs = QTabWidget()

        overview = QWidget()
        overview_layout = QVBoxLayout(overview)
        overview_layout.setContentsMargins(0, 0, 0, 0)
        self.overview_figure = Figure(figsize=(5, 3.2), dpi=100)
        self.overview_canvas = FigureCanvas(self.overview_figure)
        self.overview_canvas.setMinimumHeight(240)
        overview_layout.addWidget(self.overview_canvas)
        self.view_tabs.addTab(overview, "Overview")

        matrix = QWidget()
        matrix_layout = QVBoxLayout(matrix)
        matrix_layout.setContentsMargins(0, 0, 0, 0)
        self.matrix_figure = Figure(figsize=(5, 3.8), dpi=100)
        self.matrix_canvas = FigureCanvas(self.matrix_figure)
        self.matrix_canvas.setMinimumHeight(260)
        matrix_layout.addWidget(self.matrix_canvas)
        self.view_tabs.addTab(matrix, "RF matrix")

        pairs = QWidget()
        pairs_layout = QVBoxLayout(pairs)
        pairs_layout.setContentsMargins(0, 0, 0, 0)
        self.pairs_table = QTableWidget(0, 5)
        self.pairs_table.setHorizontalHeaderLabels(
            ["Tree 1", "Tree 2", "Comparison", "RF", "Normalised RF"]
        )
        self.pairs_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.pairs_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.pairs_table.setSortingEnabled(True)
        self.pairs_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        pairs_layout.addWidget(self.pairs_table)
        self.view_tabs.addTab(pairs, "Pairs")

        layout.addWidget(self.view_tabs, 1)

        actions = QHBoxLayout()
        self.export_button = QPushButton("Export tables...")
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self._export_tables)
        actions.addWidget(self.export_button)
        self.save_button = QPushButton("Save figures...")
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self._save_figures)
        actions.addWidget(self.save_button)
        actions.addStretch(1)
        layout.addLayout(actions)

        self._draw_placeholder()

    def _draw_placeholder(self):
        self.overview_figure.clear()
        axis = self.overview_figure.add_subplot(111)
        axis.text(0.5, 0.5, "Run the comparison to see the RF distribution",
                  ha="center", va="center", transform=axis.transAxes,
                  color="#6b7280")
        axis.set_axis_off()
        self.overview_canvas.draw_idle()

        self.matrix_figure.clear()
        axis = self.matrix_figure.add_subplot(111)
        axis.text(0.5, 0.5, "Run the comparison to see the RF matrix",
                  ha="center", va="center", transform=axis.transAxes,
                  color="#6b7280")
        axis.set_axis_off()
        self.matrix_canvas.draw_idle()

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    def run_analysis(self):
        path = self.tree_edit.text().strip()
        if not Path(path).is_file():
            QMessageBox.warning(self, "Input missing", "Please select a tree file.")
            return
        count = count_newick_trees(path)
        if count < 2:
            QMessageBox.warning(
                self, "Need more trees",
                "RF comparison needs at least two trees (one per line).",
            )
            return
        iqtree = self._iqtree_path()
        if not iqtree:
            QMessageBox.warning(
                self, "IQ-TREE missing",
                "IQ-TREE 3 is not configured. Add its path to config.json.",
            )
            return

        if self.workdir and Path(self.workdir).is_dir():
            base = Path(self.workdir)
        else:
            base = Path(tempfile.mkdtemp(prefix="yr_mpe_rf_"))
        self._out_dir = str(base / "rf_comparison")
        labels = self._tree_labels(count)

        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Computing pairwise RF distances...")
        self.tab_widget.setCurrentIndex(2)
        self.add_console_message(
            "Starting RF comparison for %d trees with IQ-TREE 3." % count,
            "info",
        )

        self.analysis_thread = RFComparisonThread(
            iqtree, path, labels, self._out_dir,
            threads=self.threads_spin.value(),
            workdir=self.workdir,
            groups=self._tree_groups(count),
        )
        self.analysis_thread.progress.connect(self._on_progress)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.start()

    def _on_progress(self, message):
        self.progress_bar.setFormat(message)
        self.add_console_message(message, "info")

    def _on_result(self, payload):
        self._payload = payload

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self._display_results()
        if self._payload:
            labels = self._payload["labels"]
            matrix = self._payload["matrix"]
            groups = self._payload.get("groups") or [None] * len(labels)
            if len(groups) != len(labels):
                groups = [None] * len(labels)
            rows = pairwise_rows(labels, matrix)
            categories = classify_pairs(groups)
            values = [row[2] for row in rows]
            import statistics
            by_category = {}
            for category in sorted(set(categories)):
                subset = [value for value, kind in zip(values, categories)
                          if kind == category]
                if subset:
                    by_category[category] = {
                        "pairs": len(subset),
                        "median_rf": statistics.median(subset),
                        "min_rf": min(subset),
                        "max_rf": max(subset),
                    }
            line = emit_summary(
                "rf_comparison",
                trees=len(labels),
                pairs=len(rows),
                median_rf=statistics.median(values) if values else None,
                min_rf=min(values) if values else None,
                max_rf=max(values) if values else None,
                categories=by_category,
            )
            self.add_console_message(SUMMARY_PREFIX + line, "output")
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message("RF comparison completed.", "info")

    def analysis_error(self, message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("RF comparison failed: %s" % message, "error")
        QMessageBox.critical(self, "RF comparison failed", message)

    def stop_analysis(self):
        if self.analysis_thread and self.analysis_thread.isRunning():
            self.analysis_thread.stop()
            self.analysis_thread.wait()
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Stopped by user.", "warning")

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------
    def _display_results(self):
        if not self._payload:
            return
        labels = self._payload["labels"]
        matrix = self._payload["matrix"]
        ntaxa = int(self._payload.get("ntaxa") or 0)
        groups = self._payload.get("groups") or [None] * len(labels)
        if len(groups) != len(labels):
            groups = [None] * len(labels)
        rows = pairwise_rows(labels, matrix)
        categories = classify_pairs(groups)
        values = [row[2] for row in rows]
        if values:
            ordered = sorted(values)
            median = ordered[len(ordered) // 2] if len(ordered) % 2 else (
                ordered[len(ordered) // 2 - 1] + ordered[len(ordered) // 2]
            ) / 2.0
            text = (
                "%d trees · %d pairs · median RF %.1f · range %d–%d%s"
                % (
                    len(labels), len(values), median, min(values), max(values),
                    " · %d taxa" % ntaxa if ntaxa else "",
                )
            )
            counts = {}
            for category in categories:
                counts[category] = counts.get(category, 0) + 1
            parts = []
            if counts.get("Within tree set"):
                parts.append("within %d" % counts["Within tree set"])
            if counts.get("Between tree sets"):
                parts.append("between %d" % counts["Between tree sets"])
            if parts:
                text += " · " + ", ".join(parts)
            self.summary_label.setText(text)
        else:
            self.summary_label.setText("%d trees · no pairs" % len(labels))

        self._fill_pairs_table(rows, categories, ntaxa)
        self._draw_overview(values, categories)
        self._draw_matrix(labels, matrix)
        self.export_button.setEnabled(True)
        self.save_button.setEnabled(True)

    def _fill_pairs_table(self, rows, categories, ntaxa):
        self.pairs_table.setSortingEnabled(False)
        self.pairs_table.setRowCount(len(rows))
        for index, ((first, second, rf), category) in enumerate(
            zip(rows, categories)
        ):
            self.pairs_table.setItem(index, 0, QTableWidgetItem(first))
            self.pairs_table.setItem(index, 1, QTableWidgetItem(second))
            self.pairs_table.setItem(index, 2, QTableWidgetItem(category))
            rf_item = QTableWidgetItem()
            rf_item.setData(Qt.DisplayRole, int(rf))
            self.pairs_table.setItem(index, 3, rf_item)
            normalised = normalized_rf(rf, ntaxa)
            norm_item = QTableWidgetItem()
            norm_item.setData(
                Qt.DisplayRole,
                "" if normalised is None else "%.4f" % normalised,
            )
            self.pairs_table.setItem(index, 4, norm_item)
        self.pairs_table.setSortingEnabled(True)
        self.pairs_table.sortItems(3, Qt.DescendingOrder)

    def _style_axis(self, axis):
        axis.tick_params(colors="#374151")
        axis.xaxis.label.set_color("#374151")
        axis.yaxis.label.set_color("#374151")
        axis.title.set_color("#111827")
        for spine in axis.spines.values():
            spine.set_color("#d1d5db")

    def _draw_overview(self, values, categories):
        self.overview_figure.clear()
        axis = self.overview_figure.add_subplot(111)
        if values:
            preferred = [
                "Within tree set", "Between tree sets",
                "Between objects", "All pairs",
            ]
            present = [name for name in preferred if name in set(categories)]
            buckets = {name: [] for name in present}
            for value, category in zip(values, categories):
                buckets[category].append(value)
            colors = {
                "Within tree set": ("#bfdbfe", "#1d4ed8"),
                "Between tree sets": ("#fde68a", "#b45309"),
                "Between objects": ("#ddd6fe", "#6d28d9"),
                "All pairs": ("#bfdbfe", "#1d4ed8"),
            }
            positions = list(range(1, len(present) + 1))
            box = axis.boxplot(
                [buckets[name] for name in present],
                positions=positions, widths=0.45, patch_artist=True,
                showfliers=False, zorder=2,
            )
            for patch, name in zip(box["boxes"], present):
                fill, edge = colors[name]
                patch.set_facecolor(fill)
                patch.set_edgecolor(edge)
            for index, name in enumerate(present):
                _, edge = colors[name]
                for line in box["whiskers"][2 * index:2 * index + 2]:
                    line.set_color(edge)
                for line in box["caps"][2 * index:2 * index + 2]:
                    line.set_color(edge)
                box["medians"][index].set_color(edge)
                points = buckets[name]
                jitter = (np.linspace(-0.08, 0.08, len(points))
                          if len(points) > 1 else np.array([0.0]))
                axis.scatter(
                    np.full(len(points), positions[index]) + jitter, points,
                    s=18, color=edge, alpha=0.55, zorder=3,
                )
            axis.set_xticks(positions)
            axis.set_xticklabels(
                [
                    name if name == "All pairs"
                    else "%s\n(n=%d)" % (name, len(buckets[name]))
                    for name in present
                ]
            )
            axis.set_ylabel("Robinson-Foulds distance")
            if present == ["All pairs"]:
                axis.set_title("Pairwise RF distances (%d pairs)" % len(values))
            else:
                axis.set_title("RF distances by comparison type (%d pairs)"
                               % len(values))
            axis.set_xlim(0.5, len(present) + 0.5)
        self._style_axis(axis)
        self.overview_figure.tight_layout()
        self.overview_canvas.draw_idle()

    def _draw_matrix(self, labels, matrix):
        self.matrix_figure.clear()
        axis = self.matrix_figure.add_subplot(111)
        array = np.array(matrix, dtype=float)
        masked = np.ma.masked_where(
            np.triu(np.ones_like(array, dtype=bool), 0), array
        )
        image = axis.imshow(masked, cmap="Blues", interpolation="nearest")
        short = _short_labels(labels)
        axis.set_xticks(range(len(labels)))
        axis.set_yticks(range(len(labels)))
        axis.set_xticklabels(short, rotation=45, ha="right", fontsize=8)
        axis.set_yticklabels(short, fontsize=8)
        axis.set_title("RF matrix (lower triangle)")
        threshold = array.max() / 2.0 if array.size else 0
        if len(labels) <= 18:
            for i in range(len(labels)):
                for j in range(i):
                    value = int(array[i][j])
                    axis.text(
                        j, i, str(value), ha="center", va="center",
                        fontsize=8,
                        color="white" if value > threshold else "#1f2937",
                    )
        self._style_axis(axis)
        self.matrix_figure.colorbar(image, ax=axis, shrink=0.8,
                                    label="RF distance")
        self.matrix_figure.tight_layout()
        self.matrix_canvas.draw_idle()

    # ------------------------------------------------------------------
    # export
    # ------------------------------------------------------------------
    def _export_tables(self):
        if not self._payload:
            return
        start = self._out_dir or str(Path.home())
        directory = QFileDialog.getExistingDirectory(
            self, "Export RF tables", start
        )
        if not directory:
            return
        labels = self._payload["labels"]
        matrix = self._payload["matrix"]
        ntaxa = int(self._payload.get("ntaxa") or 0)
        clean = [label.replace("\t", " ").replace("\n", " ") for label in labels]

        matrix_path = Path(directory) / "rf_matrix.tsv"
        with open(matrix_path, "w", encoding="utf-8", newline="") as handle:
            handle.write("Tree\t" + "\t".join(clean) + "\n")
            for index, label in enumerate(clean):
                values = "\t".join(str(int(value)) for value in matrix[index])
                handle.write("%s\t%s\n" % (label, values))

        groups = self._payload.get("groups") or [None] * len(clean)
        if len(groups) != len(clean):
            groups = [None] * len(clean)
        categories = classify_pairs(groups)
        pairs_path = Path(directory) / "rf_pairs.tsv"
        with open(pairs_path, "w", encoding="utf-8", newline="") as handle:
            handle.write("Tree 1\tTree 2\tComparison\tRF\tNormalised RF\n")
            for (first, second, rf), category in zip(
                pairwise_rows(clean, matrix), categories
            ):
                normalised = normalized_rf(rf, ntaxa)
                handle.write(
                    "%s\t%s\t%s\t%d\t%s\n"
                    % (first, second, category, rf,
                       "" if normalised is None else "%.6f" % normalised)
                )
        self.add_console_message(
            "RF tables exported to %s and %s" % (matrix_path, pairs_path),
            "info",
        )

    def _save_figures(self):
        if not self._payload:
            return
        start = self._out_dir or str(Path.home())
        directory = QFileDialog.getExistingDirectory(
            self, "Save RF figures", start
        )
        if not directory:
            return
        overview_path = Path(directory) / "rf_overview.png"
        matrix_path = Path(directory) / "rf_matrix.png"
        self.overview_figure.savefig(overview_path, dpi=200, bbox_inches="tight")
        self.matrix_figure.savefig(matrix_path, dpi=200, bbox_inches="tight")
        self.add_console_message(
            "Figures saved to %s and %s" % (overview_path, matrix_path),
            "info",
        )


class RFComparisonPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return RFComparisonPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
