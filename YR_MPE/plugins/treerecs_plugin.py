"""TreeRecs 1.2 gene/species tree reconciliation plugin.

TreeRecs takes one or more gene trees and a species tree, then returns gene
trees corrected by a duplication/loss reconciliation. The plugin follows the
standard result contract: everything it produces is emitted as a ResultSpec.
"""

from __future__ import annotations

import glob
import os
import tempfile
from pathlib import Path

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)
from PyQt5.QtSvg import QSvgWidget

from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from .components.nexus_tree import read_tree_file


class TreeRecsThread(BaseProcessThread):
    """Run treerecs.exe and collect the reconciliation output."""

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path, gene_tree, species_tree, outdir, parameters,
                 workdir=None):
        super().__init__(
            tool_path, [gene_tree, species_tree], parameters, workdir=workdir
        )
        self.gene_tree = gene_tree
        self.species_tree = species_tree
        self.outdir = outdir

    def get_tool_name(self):
        return "TreeRecs Reconciliation"

    def execute_commands(self):
        try:
            os.makedirs(self.outdir, exist_ok=True)
            cmd = [
                self.tool_path,
                "-g", self.gene_tree,
                "-s", self.species_tree,
                "-o", self.outdir,
                "-O", "newick:svg",
                "-f", "-q",
                *self.parameters,
            ]
            self.console_output.emit(" ".join(cmd), "command")
            result = self.execute_command(
                cmd, cwd=self.workdir if self.workdir else None
            )
            if result.returncode != 0:
                detail = (result.stderr or "").strip()
                raise RuntimeError(
                    "TreeRecs failed%s" % (": " + detail if detail else "")
                )

            files = sorted(
                path for path in glob.glob(os.path.join(self.outdir, "*"))
                if os.path.isfile(path)
            )
            reconciled = [
                path for path in files
                if path.endswith("_recs.nwk")
                or path.endswith("_recs.tree")
                or path.endswith(".recphyloxml")
            ]
            if not reconciled:
                raise RuntimeError(
                    "TreeRecs produced no reconciled tree; see the console output."
                )

            self.progress.emit("Reconciliation completed")
            self.result_ready.emit({"output_files": files})
            self.finished.emit(files, [])
        except Exception as exc:
            self.error.emit("TreeRecs exception: %s" % exc)


class TreeRecsPlugin(BasePlugin):
    """TreeRecs gene/species tree reconciliation plugin."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        super().__init__(import_from, import_data, **kwargs)
        if import_data is not None:
            self.handle_import_data(import_data)

    def init_plugin_info(self):
        self.plugin_name = "TreeRecs 1.2"
        self.tool_name = "TreeRecs"
        self.citation = []  # citation text lives in YR_MPE/citations/*.bib
        self.citation_keys = ["penel2020treerecs"]
        self.input_types = {"Newick": ["nwk", "newick", "tre", "tree"]}
        self.output_types = {"Newick": [".nwk"], "RecPhyloXML": [".recphyloxml"]}
        self.plugin_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '..'
        )

    def control_bindings(self):
        """GUI controls mapped to the TreeRecs command line."""
        return [
            {"id": "gene_tree", "widget": "gene_tree_edit",
             "label": "Gene tree", "kind": "path",
             "cli": ["-g", "{value}"]},
            {"id": "species_tree", "widget": "species_tree_edit",
             "label": "Species tree", "kind": "path",
             "cli": ["-s", "{value}"]},
            {"id": "duplication_cost", "widget": "dup_cost_spin",
             "label": "Duplication cost", "kind": "int",
             "cli": ["-d", "{value}"]},
            {"id": "loss_cost", "widget": "loss_cost_spin",
             "label": "Loss cost", "kind": "int",
             "cli": ["-l", "{value}"]},
            {"id": "reroot", "widget": "reroot_check",
             "label": "Best root by reconciliation cost", "kind": "bool",
             "cli": ["-r"]},
            {"id": "threshold", "widget": "threshold_edit",
             "label": "Contraction threshold", "kind": "text",
             "cli": ["-t", "{value}"], "omit_values": ["none", "None", ""],
             "notes": "Leave as none to use TreeRecs' default."},
            {"id": "prefix", "widget": "prefix_combo",
             "label": "Species prefix mode", "kind": "text",
             "cli": ["-p", "{value}"], "omit_values": ["N", "n", ""],
             "notes": "Only Y enables separator-based species mapping."},
            {"id": "separator", "widget": "separator_edit",
             "label": "Species separator", "kind": "text",
             "cli": ["-c", "{value}"], "omit_values": ["_", ""]},
            {"id": "events", "widget": "events_check",
             "label": "Write orthology/paralogy events", "kind": "bool",
             "cli": ["--fevent"]},
            {"id": "mapping", "widget": "mapping_edit",
             "label": "Gene-to-species mapping", "kind": "path",
             "cli": ["-S", "{value}"]},
        ]

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def setup_input_tab(self):
        layout = QVBoxLayout(self.input_tab)
        layout.setContentsMargins(12, 12, 12, 12)

        trees_group = QGroupBox("Trees")
        trees_layout = QFormLayout(trees_group)

        self.gene_tree_edit = QLineEdit()
        self.gene_tree_edit.setPlaceholderText(
            "One gene tree, or a file with one Newick tree per line"
        )
        gene_row = QHBoxLayout()
        gene_row.addWidget(self.gene_tree_edit, 1)
        gene_browse = QPushButton("Browse...")
        gene_browse.clicked.connect(
            lambda: self.browse_file(self.gene_tree_edit, "Gene tree")
        )
        gene_row.addWidget(gene_browse)
        trees_layout.addRow("Gene trees:", gene_row)

        self.species_tree_edit = QLineEdit()
        self.species_tree_edit.setPlaceholderText("Species tree (Newick)")
        species_row = QHBoxLayout()
        species_row.addWidget(self.species_tree_edit, 1)
        species_browse = QPushButton("Browse...")
        species_browse.clicked.connect(
            lambda: self.browse_file(self.species_tree_edit, "Species tree")
        )
        species_row.addWidget(species_browse)
        trees_layout.addRow("Species tree:", species_row)

        self.swap_button = QPushButton("Swap gene / species")
        self.swap_button.clicked.connect(self.swap_trees)
        trees_layout.addRow("", self.swap_button)

        self.mapping_edit = QLineEdit()
        self.mapping_edit.setPlaceholderText(
            "Optional gene-to-species mapping file"
        )
        mapping_row = QHBoxLayout()
        mapping_row.addWidget(self.mapping_edit, 1)
        mapping_browse = QPushButton("Browse...")
        mapping_browse.clicked.connect(
            lambda: self.browse_file(self.mapping_edit, "Mapping file")
        )
        mapping_row.addWidget(mapping_browse)
        trees_layout.addRow("Mapping:", mapping_row)

        layout.addWidget(trees_group)

        options_group = QGroupBox("Reconciliation")
        options_form = QFormLayout(options_group)

        self.reroot_check = QCheckBox("Find the best root by reconciliation cost")
        options_form.addRow("", self.reroot_check)

        self.dup_cost_spin = QSpinBox()
        self.dup_cost_spin.setRange(0, 100)
        self.dup_cost_spin.setValue(2)
        options_form.addRow("Duplication cost:", self.dup_cost_spin)

        self.loss_cost_spin = QSpinBox()
        self.loss_cost_spin.setRange(0, 100)
        self.loss_cost_spin.setValue(1)
        options_form.addRow("Loss cost:", self.loss_cost_spin)

        self.threshold_edit = QLineEdit("none")
        self.threshold_edit.setToolTip(
            "none, all, a support value, value+e, or quantiles[N]"
        )
        options_form.addRow("Support threshold:", self.threshold_edit)

        self.prefix_combo = QLineEdit()
        self.prefix_combo.setPlaceholderText("Y or N")
        self.prefix_combo.setText("N")
        options_form.addRow("Species name is prefix:", self.prefix_combo)

        self.separator_edit = QLineEdit("_")
        options_form.addRow("Name separator:", self.separator_edit)

        self.events_check = QCheckBox("Write orthology/paralogy events")
        options_form.addRow("", self.events_check)

        layout.addWidget(options_group)

        self.status_label = QLabel(
            "Load a gene tree and a species tree, then run TreeRecs."
        )
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)
        layout.addStretch(1)

    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.output_info = QLabel("No reconciliation yet.")
        self.output_info.setWordWrap(True)
        layout.addWidget(self.output_info)

        self.output_tabs = QTabWidget()
        svg_page = QWidget()
        svg_layout = QVBoxLayout(svg_page)
        svg_layout.setContentsMargins(0, 0, 0, 0)
        self.svg_preview = QSvgWidget()
        self.svg_preview.setMinimumHeight(320)
        svg_layout.addWidget(self.svg_preview, 1)
        self.output_tabs.addTab(svg_page, "Reconciliation (SVG)")

        text_page = QWidget()
        text_layout = QVBoxLayout(text_page)
        text_layout.setContentsMargins(0, 0, 0, 0)
        self.output_preview = QPlainTextEdit()
        self.output_preview.setReadOnly(True)
        self.output_preview.setPlaceholderText("Reconciled gene tree (Newick)")
        text_layout.addWidget(self.output_preview, 1)
        self.output_tabs.addTab(text_page, "Newick")
        layout.addWidget(self.output_tabs, 1)

        self.import_to_platform_btn = QPushButton("Send to workspace")
        self.import_to_platform_btn.setVisible(False)
        self.import_to_platform_btn.clicked.connect(self.import_to_platform)
        layout.addWidget(self.import_to_platform_btn)

    def browse_file(self, line_edit, label):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select %s" % label, "",
            "Tree files (*.nwk *.newick *.tre *.tree *.nex *.nexus);;All files (*)"
        )
        if path:
            line_edit.setText(path)

    def swap_trees(self):
        gene = self.gene_tree_edit.text().strip()
        species = self.species_tree_edit.text().strip()
        self.gene_tree_edit.setText(species)
        self.species_tree_edit.setText(gene)

    def set_input_trees(self, paths):
        """Populate the two tree fields from workspace item file paths.

        With exactly two paths the first is the gene tree; with more, the
        last path is the species tree and the others are combined into one
        multi-tree gene file.
        """
        paths = [str(path) for path in paths if path]
        if len(paths) < 2:
            self.status_label.setText(
                "Select at least two trees: gene tree(s) and one species tree."
            )
            return
        species = paths[-1]
        genes = paths[:-1]
        if len(genes) == 1:
            gene_file = genes[0]
        else:
            gene_file = self.create_temp_file(suffix='.nwk')
            with open(gene_file, "w", encoding="utf-8") as handle:
                for path in genes:
                    try:
                        content = Path(path).read_text(
                            encoding="utf-8", errors="replace"
                        )
                    except OSError:
                        continue
                    handle.write(content.rstrip() + "\n")
        self.gene_tree_edit.setText(gene_file)
        self.species_tree_edit.setText(species)
        self.status_label.setText(
            "Gene trees: %d file(s); species tree: %s"
            % (len(genes), os.path.basename(species))
        )

    def handle_import_data(self, import_data):
        if isinstance(import_data, dict):
            gene = import_data.get("gene_tree_file")
            species = import_data.get("species_tree_file")
            if gene and species:
                self.gene_tree_edit.setText(str(gene))
                self.species_tree_edit.setText(str(species))
                return
            datasets = import_data.get("dataset_items") or []
            paths = [
                getattr(item, "file_path", "")
                for item in datasets
                if getattr(item, "item_type", "") in ("phylogeny", "treeset")
            ]
            self.set_input_trees(paths)
        elif isinstance(import_data, (list, tuple)):
            self.set_input_trees(import_data)

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    def get_parameters(self):
        params = [
            "-d", str(self.dup_cost_spin.value()),
            "-l", str(self.loss_cost_spin.value()),
        ]
        if self.reroot_check.isChecked():
            params.append("-r")

        # Treerecs 1.2 throws std::invalid_argument ("stod") when given
        # "-t none", which is also its no-contraction default. Only pass a
        # threshold the user actually set.
        threshold = self.threshold_edit.text().strip()
        if threshold and threshold.lower() != "none":
            params.extend(["-t", threshold])

        # Passing "-p N"/"-c _" explicitly switches Treerecs into separator
        # mapping mode; leaves such as A/B/C/D then fail to map. Both options
        # are therefore only emitted when the user asks for them.
        prefix = self.prefix_combo.text().strip().upper()
        if prefix == "Y":
            params.extend(["-p", "Y"])

        separator = self.separator_edit.text().strip()
        if separator and separator != "_":
            params.extend(["-c", separator])

        if self.events_check.isChecked():
            params.append("--fevent")

        mapping = self.mapping_edit.text().strip()
        if mapping:
            params.extend(["-S", mapping])
        return params

    def run_analysis(self):
        gene_tree = self.gene_tree_edit.text().strip()
        species_tree = self.species_tree_edit.text().strip()
        if not gene_tree or not os.path.exists(gene_tree):
            QMessageBox.warning(self, "Warning", "Please provide a gene tree file.")
            return
        if not species_tree or not os.path.exists(species_tree):
            QMessageBox.warning(self, "Warning", "Please provide a species tree file.")
            return
        if not self.tool_path or not os.path.exists(self.tool_path):
            QMessageBox.critical(
                self, "Error",
                "TreeRecs executable file not found! Please check config.json."
            )
            return

        outdir = tempfile.mkdtemp(
            prefix='treerecs_', dir=self.workdir if self.workdir else None
        )
        self.current_outdir = outdir

        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)

        self.analysis_thread = TreeRecsThread(
            self.tool_path, gene_tree, species_tree, outdir,
            self.get_parameters(), workdir=self.workdir,
        )
        self.analysis_thread.progress.connect(self.progress_bar.setFormat)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self.analysis_finished)
        self.analysis_thread.error.connect(self.analysis_error)
        self.analysis_thread.console_output.connect(self.add_console_message)
        self.analysis_thread.start()

    def _on_result(self, payload):
        self._result_payload = payload

    def analysis_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

        self.current_output_files = output_files
        self.display_results(output_files)
        self.tab_widget.setCurrentIndex(1)
        self.add_console_message(
            "Reconciliation completed! Found %d file(s)" % len(output_files),
            "info",
        )
        if (self.import_from in ["YR_MPEA", "DATASET_MANAGER"]
                or self.result_has_consumer()):
            self.import_to_platform_btn.setVisible(True)
        else:
            self.import_to_platform_btn.setVisible(False)

        QMessageBox.information(self, "Completed", "TreeRecs reconciliation completed!")

    def analysis_error(self, error_message):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.add_console_message("Reconciliation failed: %s" % error_message, "error")
        QMessageBox.critical(self, "Error", "Reconciliation failed: %s" % error_message)

    def stop_analysis(self):
        thread = getattr(self, "analysis_thread", None)
        if thread is not None and thread.isRunning():
            for signal_name in ("progress", "result_ready", "finished",
                                "error", "console_output"):
                signal = getattr(thread, signal_name, None)
                if signal is None:
                    continue
                try:
                    signal.disconnect()
                except (TypeError, RuntimeError):
                    pass
            thread.stop()
            thread.wait()
            self.add_console_message("Reconciliation stopped by user.", "warning")
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)

    # ------------------------------------------------------------------
    # results
    # ------------------------------------------------------------------
    def _reconciled_files(self, output_files):
        return [
            path for path in output_files
            if os.path.basename(path).endswith("_recs.nwk")
            or os.path.basename(path).endswith("_recs.tree")
        ]

    def _reconciled_svg(self, output_files):
        return next(
            (path for path in output_files
             if os.path.basename(path).endswith("_recs.svg")),
            None,
        )

    def display_results(self, output_files):
        reconciled = self._reconciled_files(output_files)
        svg = self._reconciled_svg(output_files)
        if not reconciled:
            self.output_info.setText(
                "No reconciled tree found. Generated %d file(s)." % len(output_files)
            )
            self.output_preview.setPlainText(
                "\n".join(os.path.basename(path) for path in output_files)
            )
            self.output_tabs.setCurrentIndex(1)
            return
        first = reconciled[0]
        try:
            trees = read_tree_file(first)
        except OSError as exc:
            self.output_info.setText("Cannot read %s: %s" % (os.path.basename(first), exc))
            self.output_preview.clear()
            return
        if not trees:
            self.output_info.setText(
                "No tree could be parsed from %s." % os.path.basename(first)
            )
            self.output_preview.clear()
            return
        self.output_info.setText(
            "Reconciled tree: %s · %d output file(s)"
            % (os.path.basename(first), len(output_files))
        )
        self.output_preview.setPlainText("\n".join(trees))
        if svg:
            self.svg_preview.load(svg)
            self.output_tabs.setCurrentIndex(0)
        else:
            self.output_tabs.setCurrentIndex(1)

    def _build_result_spec(self):
        files = list(getattr(self, "current_output_files", []) or [])
        if not files:
            raise ValueError("No TreeRecs output is available.")

        items = []
        artifacts = []
        for path in files:
            name = os.path.basename(path)
            if not (name.endswith("_recs.nwk") or name.endswith("_recs.tree")):
                artifacts.append(str(path))
                continue
            trees = read_tree_file(path)
            if not trees:
                raise ValueError("No tree could be parsed from %s." % name)
            if len(trees) > 1:
                items.append(DatasetItemSpec(
                    loci_name=os.path.splitext(name)[0],
                    item_type="treeset",
                    file_path=str(path),
                    format="newick",
                    data={"data": [{"content": tree} for tree in trees]},
                    metadata={"source": self.plugin_name, "tool": "TreeRecs"},
                ))
            else:
                items.append(DatasetItemSpec(
                    loci_name=os.path.splitext(name)[0],
                    item_type="phylogeny",
                    file_path=str(path),
                    format="newick",
                    data={"newick": trees[0]},
                    metadata={"source": self.plugin_name, "tool": "TreeRecs"},
                ))

        if not items:
            raise ValueError("TreeRecs produced no reconciled tree item.")

        dataset = DatasetSpec(
            name="TreeRecs reconciliation",
            items=items,
            description="%d reconciled tree item(s)" % len(items),
            source="treerecs",
            summary={"items": len(items)},
            metadata={
                "gene_tree": self.gene_tree_edit.text().strip(),
                "species_tree": self.species_tree_edit.text().strip(),
                "parameters": self.get_parameters(),
            },
        )
        return ResultSpec(source="treerecs", datasets=[dataset], artifacts=artifacts)

    def import_to_platform(self):
        if not getattr(self, "current_output_files", None):
            QMessageBox.warning(self, "No result", "Run TreeRecs first.")
            return
        if not self.result_has_consumer():
            QMessageBox.warning(
                self, "No workspace",
                "This window is not connected to a YR-MPE workspace."
            )
            return
        try:
            result = self._build_result_spec()
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        if self.export_result(result):
            self.import_to_platform_btn.setEnabled(False)
            self.import_to_platform_btn.setText("Sent to workspace")
            QMessageBox.information(
                self, "Success", "TreeRecs result sent to the current workspace."
            )


class TreeRecsPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None):
        return TreeRecsPlugin(
            import_from=import_from, import_data=import_data, workdir=workdir
        )
