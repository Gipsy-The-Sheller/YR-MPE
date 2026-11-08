"""Qt shell and workers for the gCF, sCF and ASTER diagnostics plugins."""

from __future__ import annotations

import math
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt, QTimer, pyqtSignal
from PyQt5.QtGui import QColor
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
    QSpinBox,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from . import concordance_common as common


class ConcordanceThread(BaseProcessThread):
    """Run IQ-TREE and, for diagnostics, ASTER on a fixed species tree."""

    result_ready = pyqtSignal(object)

    def __init__(self, mode, iqtree, astral, caster, species_tree, gene_trees,
                 alignment, outdir, prefix, threads, quartets=1000, workdir=None):
        super().__init__(iqtree, [species_tree], [], workdir=workdir)
        self.mode = mode
        self.iqtree = iqtree
        self.astral = astral
        self.caster = caster
        self.species_tree = species_tree
        self.gene_trees = gene_trees
        self.alignment = alignment
        self.outdir = outdir
        self.prefix_path = str(Path(outdir) / prefix)
        self.threads = int(threads)
        self.quartets = int(quartets)
        self.payload = {}

    def get_tool_name(self):
        return "Concordance factor analysis"

    def _run(self, cmd):
        self.console_output.emit(" ".join(str(part) for part in cmd), "command")
        result = self.execute_command(
            [str(part) for part in cmd],
            cwd=self.workdir if self.workdir else None,
        )
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or "").strip()
            raise RuntimeError(
                "Command failed (%s): %s" % (
                    " ".join(str(part) for part in cmd),
                    detail[-1200:] if detail else "no stderr",
                )
            )
        return result

    def execute_commands(self):
        try:
            Path(self.outdir).mkdir(parents=True, exist_ok=True)
            tree_option = (
                "-te" if self.mode in ("scf", "site_diagnostics") else "-t"
            )
            iq_cmd = [
                self.iqtree, tree_option, self.species_tree,
                "--prefix", self.prefix_path,
                "-T", str(self.threads), "-redo",
            ]
            if self.mode in ("gcf", "diagnostics"):
                iq_cmd.extend(["--gcf", self.gene_trees])
            if self.mode in ("scf", "site_diagnostics"):
                iq_cmd.extend(["--scfl", str(self.quartets), "-s", self.alignment])
            self._run(iq_cmd)

            cf_mode = "scf" if self.mode in ("scf", "site_diagnostics") else "gcf"
            columns, rows, summary = common.cf_rows(self.prefix_path, cf_mode)
            wastral_status = ""
            if self.mode == "site_diagnostics":
                caster_out = str(Path(self.outdir) / "caster_site.tre")
                self._run([
                    self.caster, "-C", "-c", self.species_tree,
                    "-u", "2", "--length", "CULength",
                    "-t", str(self.threads), "-o", caster_out,
                    self.alignment,
                ])
                caster_rows = common.parse_astral_tree(caster_out)
                columns, rows, summary = common.merge_site_diagnostics_rows(
                    columns, rows, caster_rows
                )
                summary["caster"] = "CASTER-site"
            if self.mode == "diagnostics":
                astral_out = str(Path(self.outdir) / "astral4.tre")
                self._run([
                    self.astral, "astral4", "-C", "-c", self.species_tree,
                    "-u", "2", "-t", str(self.threads), "-o", astral_out,
                    self.gene_trees,
                ])
                astral_rows = common.parse_astral_tree(astral_out)
                wastral_rows = None
                ok, reason, detail = common.check_wastral_labels(self.gene_trees)
                if ok:
                    wastral_out = str(Path(self.outdir) / "wastral.tre")
                    try:
                        self._run([
                            self.astral, "wastral", "-C", "-c", self.species_tree,
                            "-u", "2", "-t", str(self.threads), "-o", wastral_out,
                            self.gene_trees,
                        ])
                        wastral_rows = common.parse_astral_tree(wastral_out)
                        wastral_status = "available"
                    except Exception as exc:
                        wastral_status = "rejected: %s" % exc
                        self.console_output.emit(
                            "wASTRAL skipped: %s" % exc, "warning"
                        )
                else:
                    wastral_status = "rejected: %s" % reason
                    self.console_output.emit(
                        "wASTRAL rejected: %s (trees=%s, missing=%s, composite=%s, "
                        "nonnumeric=%s)" % (
                            reason, detail.get("trees"), detail.get("missing"),
                            detail.get("composite"), detail.get("nonnumeric"),
                        ),
                        "warning",
                    )
                columns, rows, summary = common.merge_diagnostics_rows(
                    columns, rows, astral_rows, wastral_rows, wastral_status
                )
                summary["astral"] = "ASTRAL-IV"
                summary["wastral"] = wastral_status or "available"

            table_path = str(Path(self.outdir) / ("%s.tsv" % self.mode))
            common.write_table(table_path, columns, rows)
            try:
                taxa = common.tree_taxa(self.species_tree)
            except Exception:
                taxa = []
            artifacts = [
                table_path,
                self.prefix_path + ".cf.stat",
                self.prefix_path + ".cf.branch",
                self.prefix_path + ".cf.tree",
            ]
            if self.mode == "site_diagnostics":
                artifacts.append(str(Path(self.outdir) / "caster_site.tre"))
            self.payload = {
                "table_path": table_path,
                "columns": list(columns),
                "rows": [dict(row) for row in rows],
                "summary": summary,
                "artifacts": artifacts,
                "taxa": taxa,
                "mode": self.mode,
            }
            self.progress.emit("Concordance factors ready")
            self.result_ready.emit(self.payload)
            self.finished.emit(artifacts, [])
        except Exception as exc:
            self.error.emit("Concordance analysis failed: %s" % exc)


class ConcordancePluginBase(BasePlugin):
    """Common UI and result handling for the three concordance plugins."""

    MODE = "gcf"
    TITLE = "Gene concordance factors (gCF)"
    NEEDS_GENE_TREES = True
    NEEDS_ALIGNMENT = False
    RUNS_ASTRAL = False
    RUNS_CASTER = False

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self._pending_payload = import_data
        super().__init__(import_from=import_from, import_data=import_data, **kwargs)
        self._apply_pending_payload()

    def init_plugin_info(self):
        self.plugin_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
        self.plugin_name = self.TITLE
        self.tool_name = "IQ-TREE 3"
        self.citation_keys = ["minh2020gcf"]
        self.input_types = {"Newick": ["nwk", "newick", "tre", "tree"]}
        self.output_types = {"TSV": ".tsv"}

    def control_bindings(self):
        bindings = [
            {"id": "species_tree", "widget": "species_tree_edit",
             "label": "Species tree", "kind": "path", "cli": []},
            {"id": "output_dir", "widget": "output_dir_edit",
             "label": "Output directory", "kind": "path", "cli": []},
            {"id": "threads", "widget": "threads_spin",
             "label": "Threads", "kind": "int", "cli": ["-T", "{value}"]},
        ]
        if self.NEEDS_GENE_TREES:
            bindings.append({
                "id": "gene_trees", "widget": "gene_trees_edit",
                "label": "Gene trees", "kind": "path", "cli": ["--gcf", "{value}"],
            })
        if self.NEEDS_ALIGNMENT:
            bindings.append({
                "id": "alignment", "widget": "alignment_edit",
                "label": "Alignment", "kind": "path",
                "cli": ["--scfl", "{quartets}", "-s", "{value}"],
            })
        return bindings

    # ------------------------------------------------------------------
    # input
    # ------------------------------------------------------------------
    def setup_input_tab(self):
        outer = QVBoxLayout(self.input_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        title = QLabel(self.TITLE)
        title.setStyleSheet("font-size: 17px; font-weight: 600;")
        outer.addWidget(title)
        subtitle = QLabel(
            "Fixed species tree concordance analysis. Results are written as a "
            "TSV table; each split is named by its smaller taxon side."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        outer.addWidget(subtitle)

        group = QGroupBox("Inputs")
        form = QFormLayout(group)
        self.species_tree_edit = self._path_row(
            form, "Species tree:", self._browse_species,
            "Newick species tree with all taxa",
        )
        if self.NEEDS_GENE_TREES:
            self.gene_trees_edit = self._path_row(
                form, "Gene trees:", self._browse_gene_trees,
                "Newick gene trees, one per line (support labels required for wASTRAL)",
            )
        if self.NEEDS_ALIGNMENT:
            self.alignment_edit = self._path_row(
                form, "Alignment:", self._browse_alignment,
                "Aligned FASTA/PHYLIP/NEXUS",
            )
        self.output_dir_edit = self._path_row(
            form, "Output directory:", self._browse_output,
            "Directory for the TSV table and IQ-TREE artifacts",
        )
        outer.addWidget(group)

        options = QGroupBox("Options")
        options_form = QFormLayout(options)
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, max(1, os.cpu_count() or 1))
        self.threads_spin.setValue(min(4, max(1, (os.cpu_count() or 2) // 2)))
        options_form.addRow("Threads:", self.threads_spin)
        if self.MODE in ("scf", "site_diagnostics"):
            self.quartets_spin = QSpinBox()
            self.quartets_spin.setRange(10, 100000)
            self.quartets_spin.setValue(1000)
            self.quartets_spin.setToolTip(
                "Number of quartets resampled for likelihood sCF (--scfl)"
            )
            options_form.addRow("sCF quartets:", self.quartets_spin)
        outer.addWidget(options)

        if self.MODE == "site_diagnostics":
            interp_group = QGroupBox("Interpretation")
            interp_layout = QVBoxLayout(interp_group)
            interp = QLabel(
                "CASTER q1-q3 are site-based quartet scores and may be "
                "negative because of negative quartet weights. sCF is "
                "sensitive to multiple hits."
            )
            interp.setWordWrap(True)
            interp.setStyleSheet("color: palette(mid);")
            interp_layout.addWidget(interp)
            outer.addWidget(interp_group)

        if self.MODE == "diagnostics":
            check_group = QGroupBox("wASTRAL input check")
            check_layout = QVBoxLayout(check_group)
            self.wastral_check_label = QLabel(
                "Select a gene-tree file to check branch labels."
            )
            self.wastral_check_label.setWordWrap(True)
            self.wastral_check_label.setStyleSheet("color: palette(mid);")
            check_layout.addWidget(self.wastral_check_label)
            check_note = QLabel(
                "wASTRAL needs one finite numeric label on every internal "
                "gene-tree branch. Missing, non-numeric or composite labels "
                "such as /1/100 cause wASTRAL to be skipped; the gCF and "
                "ASTRAL-IV columns are still produced."
            )
            check_note.setWordWrap(True)
            check_note.setStyleSheet("color: palette(mid);")
            check_layout.addWidget(check_note)
            outer.addWidget(check_group)
            self._wastral_check_timer = QTimer(self)
            self._wastral_check_timer.setSingleShot(True)
            self._wastral_check_timer.setInterval(300)
            self._wastral_check_timer.timeout.connect(self._update_wastral_check)
            self.gene_trees_edit.textChanged.connect(self._schedule_wastral_check)
            self._update_wastral_check()

        outer.addStretch(1)

        if self.workdir:
            self.output_dir_edit.setText(
                str(Path(self.workdir) / ("concordance_%s" % self.MODE))
            )
        self._apply_pending_payload()

    def _path_row(self, form, label, slot, placeholder):
        edit = QLineEdit()
        edit.setPlaceholderText(placeholder)
        button = QPushButton("Browse...")
        button.clicked.connect(slot)
        row = QHBoxLayout()
        row.addWidget(edit, 1)
        row.addWidget(button)
        form.addRow(label, row)
        return edit

    def _browse_species(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select species tree", "", "Newick (*.nwk *.newick *.tre *.tree);;All files (*)"
        )
        if path:
            self.species_tree_edit.setText(path)

    def _browse_gene_trees(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select gene trees", "", "Newick (*.nwk *.newick *.tre *.tree);;All files (*)"
        )
        if path:
            self.gene_trees_edit.setText(path)

    def _browse_alignment(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select alignment", "", "Alignments (*.fas *.fasta *.phy *.nex *.nexus);;All files (*)"
        )
        if path:
            self.alignment_edit.setText(path)

    def _browse_output(self):
        path = QFileDialog.getExistingDirectory(self, "Select output directory")
        if path:
            self.output_dir_edit.setText(path)

    def _schedule_wastral_check(self):
        timer = getattr(self, "_wastral_check_timer", None)
        if timer is not None:
            timer.start()

    def _update_wastral_check(self):
        label = getattr(self, "wastral_check_label", None)
        if label is None:
            return
        path = self.gene_trees_edit.text().strip()
        if not path:
            label.setText("Select a gene-tree file to check branch labels.")
            label.setStyleSheet("color: palette(mid);")
            return
        if not Path(path).is_file():
            label.setText("Gene-tree file not found.")
            label.setStyleSheet("color: #b45309;")
            return
        ok, reason, detail = common.check_wastral_labels(path)
        if ok:
            trees = int(detail.get("trees", 0))
            label.setText(
                "OK — %d gene %s; all internal branches carry finite "
                "numeric labels. wASTRAL is enabled."
                % (trees, "tree" if trees == 1 else "trees")
            )
            label.setStyleSheet("color: #2e7d32;")
        else:
            label.setText(
                "wASTRAL will be rejected — %s. The gCF and ASTRAL-IV "
                "columns are still produced." % reason
            )
            label.setStyleSheet("color: #b45309;")

    def _apply_pending_payload(self):
        payload = getattr(self, "_pending_payload", None)
        if not isinstance(payload, dict):
            return
        mapping = {
            "species_tree": "species_tree_edit",
            "gene_trees": "gene_trees_edit",
            "alignment": "alignment_edit",
        }
        for key, attr in mapping.items():
            value = payload.get(key)
            if value and hasattr(self, attr):
                getattr(self, attr).setText(str(value))
        output = payload.get("output_dir")
        if output and hasattr(self, "output_dir_edit"):
            self.output_dir_edit.setText(str(output))

    # ------------------------------------------------------------------
    # output
    # ------------------------------------------------------------------
    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.summary_label = QLabel("No concordance result yet.")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.summary_label)

        self.view_tabs = QTabWidget()
        self.table_view = QTableWidget(0, 0)
        self.table_view.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table_view.horizontalHeader().setStretchLastSection(True)
        self.view_tabs.addTab(self.table_view, "Table")

        overview = QWidget()
        overview_layout = QVBoxLayout(overview)
        overview_layout.setContentsMargins(0, 0, 0, 0)
        self.overview_figure = Figure(figsize=(6.0, 3.4), dpi=100)
        self.overview_canvas = FigureCanvas(self.overview_figure)
        self.overview_canvas.setMinimumHeight(260)
        overview_layout.addWidget(self.overview_canvas)
        self.view_tabs.addTab(overview, "Overview")

        simplex = QWidget()
        simplex_layout = QVBoxLayout(simplex)
        simplex_layout.setContentsMargins(0, 0, 0, 0)
        self.simplex_note = QLabel("")
        self.simplex_note.setWordWrap(True)
        self.simplex_note.setStyleSheet("color: palette(mid);")
        simplex_layout.addWidget(self.simplex_note)
        self.simplex_figure = Figure(figsize=(6.0, 4.0), dpi=100)
        self.simplex_canvas = FigureCanvas(self.simplex_figure)
        self.simplex_canvas.setMinimumHeight(260)
        simplex_layout.addWidget(self.simplex_canvas, 1)
        self.view_tabs.addTab(simplex, "Simplex")
        self._draw_simplex_placeholder()

        bipartition = QWidget()
        bipartition_layout = QVBoxLayout(bipartition)
        bipartition_layout.setContentsMargins(0, 0, 0, 0)
        self.bipartition_note = QLabel(
            "Blue = taxon in the smaller side of that split. Columns S1...Sn "
            "follow the Table tab; hover a column header for the full split."
        )
        self.bipartition_note.setWordWrap(True)
        self.bipartition_note.setStyleSheet("color: palette(mid);")
        bipartition_layout.addWidget(self.bipartition_note)
        self.bipartition_key = QLabel("")
        self.bipartition_key.setWordWrap(True)
        self.bipartition_key.setStyleSheet("color: palette(mid);")
        bipartition_layout.addWidget(self.bipartition_key)
        self.bipartition_table = QTableWidget(0, 0)
        self.bipartition_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.bipartition_table.setSelectionMode(QAbstractItemView.NoSelection)
        self.bipartition_table.verticalHeader().setVisible(False)
        self.bipartition_table.setWordWrap(False)
        bipartition_layout.addWidget(self.bipartition_table, 1)
        self.view_tabs.addTab(bipartition, "Bipartition info")
        layout.addWidget(self.view_tabs, 1)

        if self.MODE == "diagnostics":
            note = QLabel(
                "qCF/qDF1/qDF2 are ASTRAL-IV q1/q2/q3 quartet concordance "
                "factors (not support values). Blank CU or x_gcf cells are "
                "not estimable (gCF ≤ 1/3 or = 100%); ASTER nan values are "
                "shown blank. Splits are named by the smaller taxon side."
            )
            note.setWordWrap(True)
            note.setStyleSheet("color: palette(mid);")
            layout.addWidget(note)
        elif self.MODE == "site_diagnostics":
            note = QLabel(
                "q1/q2/q3 are CASTER-site quartet scores and may be negative "
                "because of negative quartet weights. sCF is sensitive to "
                "multiple hits. Blank cells indicate a split missing from "
                "that program."
            )
            note.setWordWrap(True)
            note.setStyleSheet("color: palette(mid);")
            layout.addWidget(note)

        controls = QHBoxLayout()
        self.export_button = QPushButton("Send table to workspace")
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self._export_table)
        controls.addWidget(self.export_button)
        self.save_graph_button = QPushButton("Save graph")
        self.save_graph_button.setEnabled(False)
        self.save_graph_button.setToolTip(
            "Save the currently visible Overview or Simplex figure."
        )
        self.save_graph_button.clicked.connect(self._save_graph)
        controls.addWidget(self.save_graph_button)
        controls.addStretch(1)
        layout.addLayout(controls)
        self.view_tabs.currentChanged.connect(self._update_graph_button)
        self._draw_placeholder()

    def _draw_placeholder(self):
        self.overview_figure.clear()
        axis = self.overview_figure.add_subplot(111)
        axis.text(0.5, 0.5, "Run the analysis to see the overview",
                  ha="center", va="center", transform=axis.transAxes, color="#6b7280")
        axis.set_axis_off()
        self.overview_canvas.draw_idle()

    def _display_payload(self, payload):
        columns = payload.get("columns") or []
        rows = payload.get("rows") or []
        self.table_view.clear()
        self.table_view.setColumnCount(len(columns))
        self.table_view.setRowCount(len(rows))
        self.table_view.setHorizontalHeaderLabels(columns)
        for r, row in enumerate(rows):
            for c, column in enumerate(columns):
                self.table_view.setItem(r, c, QTableWidgetItem(str(row.get(column, ""))))
        self.table_view.resizeColumnsToContents()
        summary = payload.get("summary") or {}
        self.summary_label.setText(
            "%d splits · %s%s" % (
                len(rows),
                summary.get("source", "IQ-TREE 3"),
                (" · wASTRAL: " + str(summary.get("wastral_status")))
                if summary.get("wastral_status") else "",
            )
        )
        self._display_bipartition_info(payload)
        self._draw_overview(rows)
        self._draw_simplex(rows)
        self.export_button.setEnabled(bool(payload.get("table_path")))
        self.view_tabs.setCurrentIndex(0)

    def _display_bipartition_info(self, payload):
        rows = payload.get("rows") or []
        taxa = list(payload.get("taxa") or [])
        if not taxa:
            names = set()
            for row in rows:
                names.update(
                    part for part in str(row.get("split", "")).split(",")
                    if part
                )
            taxa = sorted(names)

        self.bipartition_table.clear()
        self.bipartition_table.setColumnCount(len(rows) + 1)
        self.bipartition_table.setRowCount(len(taxa))
        headers = ["Taxon"] + ["S%d" % (index + 1)
                               for index in range(len(rows))]
        self.bipartition_table.setHorizontalHeaderLabels(headers)
        self.bipartition_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Interactive
        )

        for row_index, taxon in enumerate(taxa):
            item = QTableWidgetItem(taxon)
            item.setFlags(Qt.ItemIsEnabled)
            item.setToolTip(taxon)
            self.bipartition_table.setItem(row_index, 0, item)

        sizes = []
        for column_index, row in enumerate(rows):
            split_text = str(row.get("split", ""))
            members = {
                part for part in split_text.split(",") if part
            }
            sizes.append(len(members))
            header_item = self.bipartition_table.horizontalHeaderItem(
                column_index + 1
            )
            if header_item is not None:
                header_item.setToolTip(split_text)
            for row_index, taxon in enumerate(taxa):
                item = QTableWidgetItem("")
                item.setFlags(Qt.ItemIsEnabled)
                if taxon in members:
                    item.setBackground(QColor("#4e699a"))
                self.bipartition_table.setItem(
                    row_index, column_index + 1, item
                )

        self.bipartition_table.resizeColumnToContents(0)
        for column_index in range(len(rows)):
            self.bipartition_table.setColumnWidth(column_index + 1, 34)
        key = " · ".join(
            "S%d (%d taxa)" % (index + 1, size)
            for index, size in enumerate(sizes)
        )
        self.bipartition_key.setText(key or "No splits.")

    def _draw_simplex_placeholder(self):
        self.simplex_figure.clear()
        axis = self.simplex_figure.add_subplot(111)
        axis.text(0.5, 0.5, "Run the analysis to see the simplex",
                  ha="center", va="center", transform=axis.transAxes,
                  color="#6b7280")
        axis.set_axis_off()
        self.simplex_canvas.draw_idle()

    def _simplex_note_text(self):
        if self.MODE == "gcf":
            return ("Ternary = closed composition of gCF, gDF1 and gDF2; "
                    "colour encodes gDFP (%).")
        if self.MODE == "scf":
            return "Ternary = closed composition of sCF, sDF1 and sDF2."
        if self.MODE == "diagnostics":
            return ("gCF panel is normalized over (gCF, gDF1, gDF2) with "
                    "gDFP as colour; qCF panel uses ASTRAL-IV q1-q3. Panels "
                    "are separate because DF1/DF2 order may differ.")
        return ("sCF and CASTER q1-q3 are closed compositions; panels are "
                "separate because DF1/DF2 order may differ. Negative CASTER "
                "q values cannot be plotted and are omitted.")

    def _draw_simplex(self, rows):
        self.simplex_note.setText(self._simplex_note_text())
        self.simplex_figure.clear()
        if not rows:
            axis = self.simplex_figure.add_subplot(111)
            axis.text(0.5, 0.5, "No rows", ha="center", va="center",
                      transform=axis.transAxes, color="#6b7280")
            axis.set_axis_off()
            self.simplex_canvas.draw_idle()
            return

        if self.MODE == "gcf":
            panels = [(
                rows, ("gCF", "gDF1", "gDF2"),
                ("gCF", "gDF1", "gDF2"), "gCF simplex", "gDFP",
            )]
        elif self.MODE == "scf":
            panels = [(
                rows, ("sCF", "sDF1", "sDF2"),
                ("sCF", "sDF1", "sDF2"), "sCF simplex", None,
            )]
        elif self.MODE == "site_diagnostics":
            panels = [
                (rows, ("sCF", "sDF1", "sDF2"),
                 ("sCF", "sDF1", "sDF2"), "sCF simplex", None),
                (rows, ("q1_caster", "q2_caster", "q3_caster"),
                 ("q1", "q2", "q3"), "CASTER q1-q3 simplex", None),
            ]
        else:
            panels = [
                (rows, ("gCF", "gDF1", "gDF2"),
                 ("gCF", "gDF1", "gDF2"), "gCF simplex", "gDFP"),
                (rows, ("qCF", "qDF1", "qDF2"),
                 ("qCF", "qDF1", "qDF2"), "ASTRAL-IV qCF simplex", None),
            ]

        if len(panels) == 1:
            axes = [self.simplex_figure.add_subplot(111)]
        else:
            grid = self.simplex_figure.add_gridspec(1, 2)
            axes = [
                self.simplex_figure.add_subplot(grid[0, 0]),
                self.simplex_figure.add_subplot(grid[0, 1]),
            ]
        for axis, panel in zip(axes, panels):
            self._draw_simplex_panel(axis, *panel)
        self.simplex_figure.tight_layout()
        self.simplex_canvas.draw_idle()

    def _draw_simplex_panel(self, axis, rows, keys, labels, title, color_key):
        root3 = math.sqrt(3.0)
        axis.plot([0, 1, 0.5, 0], [0, 0, root3 / 2, 0],
                  color="#374151", linewidth=1.0)
        axis.plot([0, 0.75], [0, root3 / 4], color="#9ca3af",
                  linestyle="--", linewidth=0.8)
        axis.plot([1, 0.25], [0, root3 / 4], color="#9ca3af",
                  linestyle="--", linewidth=0.8)
        axis.plot([0.5, 0.5], [root3 / 2, 0], color="#9ca3af",
                  linestyle="--", linewidth=0.8)
        axis.text(0.5, root3 / 2 + 0.035, labels[0],
                  ha="center", va="bottom", fontsize=8)
        axis.text(-0.02, -0.04, labels[1],
                  ha="right", va="top", fontsize=8)
        axis.text(1.02, -0.04, labels[2],
                  ha="left", va="top", fontsize=8)
        axis.set_xlim(-0.18, 1.18)
        axis.set_ylim(-0.12, 0.98)
        axis.set_aspect("equal")
        axis.axis("off")

        points = []
        color_values = []
        omitted = 0
        for row in rows:
            values = [_float(row.get(key)) for key in keys]
            if any(value is None for value in values):
                continue
            if any(value < 0 for value in values):
                omitted += 1
                continue
            total = sum(values)
            if total <= 0:
                continue
            a, b, c = [value / total for value in values]
            points.append((c + 0.5 * a, a * root3 / 2))
            if color_key:
                color_values.append(_float(row.get(color_key)) or 0.0)

        if points:
            xs = [point[0] for point in points]
            ys = [point[1] for point in points]
            if color_key and max(color_values) > 0:
                scatter = axis.scatter(
                    xs, ys, c=color_values, cmap="Blues", vmin=0, vmax=100,
                    s=28, edgecolors="#4e699a", linewidths=0.5, alpha=0.9,
                    zorder=3,
                )
                colorbar = self.simplex_figure.colorbar(
                    scatter, ax=axis, fraction=0.046, pad=0.03
                )
                colorbar.set_label("gDFP (%)", fontsize=7)
                colorbar.ax.tick_params(labelsize=6)
            else:
                axis.scatter(
                    xs, ys, color="#4e699a", s=28, edgecolors="white",
                    linewidths=0.5, alpha=0.9, zorder=3,
                )
        else:
            axis.text(0.5, 0.45, "No non-negative complete vectors",
                      ha="center", va="center", color="#6b7280",
                      fontsize=8)

        if omitted:
            title = "%s (%d negative omitted)" % (title, omitted)
        axis.set_title(title, fontsize=8)

    def _draw_overview(self, rows):
        self.overview_figure.clear()
        if not rows:
            axis = self.overview_figure.add_subplot(111)
            axis.text(0.5, 0.5, "No rows", ha="center", va="center",
                      transform=axis.transAxes, color="#6b7280")
            axis.set_axis_off()
            self.overview_canvas.draw_idle()
            return
        if self.MODE in ("gcf", "scf"):
            self._draw_single_cf_overview(rows)
        elif self.MODE == "site_diagnostics":
            grid = self.overview_figure.add_gridspec(1, 2)
            scf_bars = self.overview_figure.add_subplot(grid[0, 0])
            qcf_bars = self.overview_figure.add_subplot(grid[0, 1])
            self._draw_max_factor_bars(
                scf_bars, rows, ("sCF", "sDF1", "sDF2"),
                concordant_label="sCF max", discordant_label="sDF max",
                title="sCF vs sDF per split",
            )
            self._draw_max_factor_bars(
                qcf_bars, rows, ("q1_caster", "q2_caster", "q3_caster"),
                concordant_label="q1 max", discordant_label="q2/q3 max",
                title="CASTER q1 vs q2/q3",
            )
        else:
            grid = self.overview_figure.add_gridspec(1, 2)
            gcf_bars = self.overview_figure.add_subplot(grid[0, 0])
            qcf_bars = self.overview_figure.add_subplot(grid[0, 1])
            self._draw_max_factor_bars(
                gcf_bars, rows, ("gCF", "gDF1", "gDF2"),
                concordant_label="gCF max", discordant_label="gDF max",
                title="gCF vs gDF per split",
            )
            self._draw_max_factor_bars(
                qcf_bars, rows, ("qCF", "qDF1", "qDF2"),
                concordant_label="qCF max", discordant_label="qDF max",
                title="qCF vs qDF per split",
            )
        self.overview_figure.tight_layout()
        self.overview_canvas.draw_idle()

    def _draw_single_cf_overview(self, rows):
        key = "gCF" if self.MODE == "gcf" else "sCF"
        values = [_float(row.get(key)) or 0.0 for row in rows]
        positions = list(range(1, len(rows) + 1))
        axis = self.overview_figure.add_subplot(111)
        axis.bar(positions, values, color="#4c78a8", width=0.8)
        axis.axhline(100.0 / 3.0, color="#d1495b", linestyle="--",
                     linewidth=1, label="1/3")
        axis.set_ylim(0, 105)
        axis.set_xlabel("Split (smaller side, sorted)")
        axis.set_ylabel("%s (%%)" % key)
        axis.legend(loc="lower right")
        axis.set_title("%s per split" % key)
        axis.tick_params(labelsize=8)

    def _draw_max_factor_bars(self, axis, rows, factors,
                              concordant_label, discordant_label, title):
        counts = [0, 0]
        usable = 0
        for row in rows:
            values = [_float(row.get(key)) for key in factors]
            if not any(value is not None for value in values):
                continue
            usable += 1
            first = values[0]
            if first is not None and first >= max(
                value for value in values if value is not None
            ):
                counts[0] += 1
            else:
                counts[1] += 1
        if not usable:
            axis.text(0.5, 0.5, "No values", ha="center", va="center",
                      transform=axis.transAxes, color="#6b7280")
            axis.set_axis_off()
            return
        positions = [0, 1]
        axis.bar(positions, counts, color=["#4c78a8", "#f58518"], width=0.62)
        ceiling = max(counts) or 1
        axis.set_ylim(0, ceiling * 1.28)
        axis.set_xticks(positions)
        axis.set_xticklabels(
            [concordant_label, discordant_label], fontsize=7
        )
        axis.set_ylabel("Number of splits", fontsize=8)
        axis.set_title("%s (n=%d)" % (title, usable), fontsize=9)
        for position, count in zip(positions, counts):
            axis.text(position, count + ceiling * 0.04, str(count),
                      ha="center", va="bottom", fontsize=8)
        axis.tick_params(axis="y", labelsize=7)

    # ------------------------------------------------------------------
    # run/export
    # ------------------------------------------------------------------
    def _validate_paths(self):
        species = self.species_tree_edit.text().strip()
        if not Path(species).is_file():
            QMessageBox.warning(self, "Input missing", "Please select a species tree.")
            return None
        gene = getattr(self, "gene_trees_edit", None)
        alignment = getattr(self, "alignment_edit", None)
        payload = {"species_tree": species}
        if gene is not None:
            if not Path(gene.text().strip()).is_file():
                QMessageBox.warning(self, "Input missing", "Please select gene trees.")
                return None
            payload["gene_trees"] = gene.text().strip()
        if alignment is not None:
            if not Path(alignment.text().strip()).is_file():
                QMessageBox.warning(self, "Input missing", "Please select an alignment.")
                return None
            payload["alignment"] = alignment.text().strip()
        return payload

    def run_analysis(self):
        payload = self._validate_paths()
        if payload is None:
            return
        outdir = self.output_dir_edit.text().strip()
        if not outdir:
            outdir = tempfile.mkdtemp(prefix="yr_mpe_concordance_")
        prefix = "%s_%s" % (self.MODE, Path(payload["species_tree"]).stem)
        astral = common.resolve_tool(self.plugin_path, "astral") if self.RUNS_ASTRAL else ""
        if self.RUNS_ASTRAL and not astral:
            QMessageBox.warning(
                self, "ASTER missing",
                "ASTRAL-IV/wASTRAL needs the configured ASTER executable.",
            )
            return
        caster = common.resolve_tool(self.plugin_path, "caster-site") if self.RUNS_CASTER else ""
        if self.RUNS_CASTER and not caster:
            QMessageBox.warning(
                self, "CASTER-site missing",
                "This diagnostic needs the configured CASTER-site executable.",
            )
            return
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setFormat("Running concordance analysis...")
        self.tab_widget.setCurrentIndex(2)
        self.analysis_thread = ConcordanceThread(
            self.MODE, self.tool_path, astral, caster,
            payload["species_tree"], payload.get("gene_trees", ""),
            payload.get("alignment", ""), outdir, prefix,
            self.threads_spin.value(),
            quartets=getattr(self, "quartets_spin", None).value()
            if hasattr(self, "quartets_spin") else 1000,
            workdir=self.workdir,
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.start()

    def _on_result(self, payload):
        self._payload = payload

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        payload = getattr(self, "_payload", None)
        if not payload:
            self.add_console_message("No concordance payload returned.", "error")
            return
        self._display_payload(payload)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message(
            "Concordance table written to %s" % payload.get("table_path"), "info"
        )

    def analysis_error(self, message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message(message, "error")
        QMessageBox.critical(self, "Concordance analysis failed", message)

    def stop_analysis(self):
        thread = getattr(self, "analysis_thread", None)
        if thread is not None and thread.isRunning():
            for signal_name in ("progress", "result_ready", "finished", "error",
                                "console_output"):
                signal = getattr(thread, signal_name, None)
                if signal is not None:
                    try:
                        signal.disconnect()
                    except (TypeError, RuntimeError):
                        pass
            thread.terminate()
            thread.wait(2000)
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

    def _build_result_spec(self):
        payload = getattr(self, "_payload", None)
        if not payload:
            raise ValueError("No concordance table is available.")
        table_path = payload["table_path"]
        item = DatasetItemSpec(
            loci_name="%s_table" % self.MODE,
            item_type="table",
            file_path=table_path,
            format="tsv",
            data={
                "table_type": self.MODE,
                "columns": payload.get("columns", []),
                "rows": payload.get("rows", []),
                "summary": payload.get("summary", {}),
                "taxa": payload.get("taxa", []),
                "species_tree": self.species_tree_edit.text().strip(),
                "gene_trees": getattr(self, "gene_trees_edit", None).text().strip()
                if hasattr(self, "gene_trees_edit") else "",
                "alignment": getattr(self, "alignment_edit", None).text().strip()
                if hasattr(self, "alignment_edit") else "",
                "wastral_status": payload.get("summary", {}).get("wastral_status", ""),
            },
            metadata={"plugin": self.plugin_name, "output_dir": os.path.dirname(table_path)},
        )
        dataset = DatasetSpec(
            name="%s_%s" % (self.MODE.upper(), Path(table_path).stem),
            items=[item],
            description="%s concordance table (%d splits)" % (
                self.TITLE, len(payload.get("rows", []))
            ),
            source="concordance_%s" % self.MODE,
            summary=dict(payload.get("summary", {})),
            metadata={
                "plugin": self.plugin_name,
                "table": table_path,
                "species_tree": self.species_tree_edit.text().strip(),
            },
        )
        return ResultSpec(
            source="concordance_%s" % self.MODE,
            datasets=[dataset],
            artifacts=list(payload.get("artifacts", [])),
            summary=dict(payload.get("summary", {})),
        )

    def _update_graph_button(self, index=None):
        if not hasattr(self, "save_graph_button"):
            return
        if index is None:
            index = self.view_tabs.currentIndex()
        self.save_graph_button.setEnabled(index in (1, 2))

    def _save_graph(self):
        index = self.view_tabs.currentIndex()
        if index not in (1, 2):
            return
        figure = self.overview_figure if index == 1 else self.simplex_figure
        stem = "overview" if index == 1 else "simplex"
        path, selected = QFileDialog.getSaveFileName(
            self,
            "Save graph",
            str(Path.home() / ("yr_mpe_%s.png" % stem)),
            "PNG image (*.png);;PDF (*.pdf);;SVG (*.svg)",
        )
        if not path:
            return
        if not Path(path).suffix:
            suffix = ".pdf" if "PDF" in selected else (
                ".svg" if "SVG" in selected else ".png"
            )
            path += suffix
        try:
            figure.savefig(
                path, dpi=300, bbox_inches="tight", facecolor="white"
            )
        except Exception as exc:
            QMessageBox.warning(self, "Save failed", str(exc))
            return
        self.add_console_message("Graph saved to %s" % path, "info")

    def _export_table(self):
        try:
            result = self._build_result_spec()
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
            return
        if not self.result_has_consumer():
            QMessageBox.information(
                self, "Export unavailable",
                "Open this plugin from the YR-MPE workspace to send the table.",
            )
            return
        if self.export_result(result):
            self.add_console_message("Concordance table sent to the workspace.", "info")


def _float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
