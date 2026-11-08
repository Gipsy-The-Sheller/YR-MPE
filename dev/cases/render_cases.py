"""Render UI assets for the three YR-MPE case studies.

The screenshots are taken from the real Qt widgets, not mocked up.  Data are
the open-task fixtures (Sanger amplicons and linked organelle GenBank files)
and a MrBayes two-run trace for the diagnostics panels.
"""

from __future__ import annotations

import os
import sys
import time
import uuid
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "windows" if sys.platform == "win32" else "offscreen")

from PyQt5.QtCore import Qt  # noqa: E402
from PyQt5.QtGui import QFont, QFontDatabase  # noqa: E402
from PyQt5.QtWidgets import QApplication, QTableWidgetItem  # noqa: E402
from PyQt5.QtTest import QTest  # noqa: E402

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from Bio import SeqIO  # noqa: E402

from YR_MPE.plugins.icytree import IcyTreePlugin  # noqa: E402
from YR_MPE.plugins.mafft_plugin import MAFFTAdvancedDialog  # noqa: E402
from YR_MPE.plugins.minitracer_plugin import MiniTracerPlugin  # noqa: E402
from YR_MPE.plugins.ml_distance_plugin import MLDistancePlugin  # noqa: E402
from YR_MPE.plugins.model_finder_plugin import ModelFinderPlugin  # noqa: E402
from YR_MPE.plugins.mrbayes_plugin import MrBayesPlugin  # noqa: E402
from YR_MPE.plugins.phylobayes_plugin import PhyloBayesPlugin  # noqa: E402
from YR_MPE.plugins.seqdbg_plugin import SeqDBGWindow  # noqa: E402
from YR_MPE.plugins.seqmatrix import Lg_SeqMatrix  # noqa: E402
from YR_MPE.sequence_editor import SequenceAlignmentViewer  # noqa: E402

TASK02 = Path("E:/YR-MPE_tutorial/benchmark/open_tasks/task02")
ORGANELLE = Path("E:/YR-MPE_tutorial/YR-MPE_dev/dev/gbff_to_orthofinder/testdata/gb")
BAYES = HERE / "bayes"
ASSETS = HERE / "assets"
BAYES.mkdir(parents=True, exist_ok=True)
ASSETS.mkdir(parents=True, exist_ok=True)


FONT_PATH = Path("E:/Codes/YRTools/fonts/MiSans-Medium.ttf")


def apply_ui_font(app=None) -> None:
    """Use the host MiSans font for every screenshot."""
    app = app or QApplication.instance()
    if app is None or not FONT_PATH.exists():
        return
    font_id = QFontDatabase.addApplicationFont(str(FONT_PATH))
    families = QFontDatabase.applicationFontFamilies(font_id) if font_id != -1 else []
    if families:
        app_font = QFont(families[0], 10)
        app_font.setWeight(QFont.Medium)
        app.setFont(app_font)


def wait(ms: int = 300) -> None:
    QTest.qWait(ms)


def save(widget, name: str) -> None:
    widget.show()
    QApplication.processEvents()
    wait(500)
    path = ASSETS / name
    widget.grab().save(str(path))
    print("saved", path)


def load_fasta(path: Path):
    return list(SeqIO.parse(str(path), "fasta"))


def render_distance_case() -> None:
    """ITS2 ML distance matrix in the current distance plugin."""
    matrix = HERE / "single_gene" / "its2.mldist"
    widget = MLDistancePlugin()
    widget.resize(1120, 760)
    widget.display_results([str(matrix)])
    widget.tab_widget.setCurrentIndex(1)
    save(widget, "F_its2_distance_matrix.png")
    widget.close()


def render_seqmatrix() -> None:
    """A 12-taxon x 3-locus Sanger matrix with the missing ITS cell visible."""
    loci = ("18S", "ITS", "28S")
    records = {locus: {rec.id: str(rec.seq) for rec in load_fasta(TASK02 / "input" / f"{locus}.fasta")}
               for locus in loci}
    taxa = sorted(set().union(*(set(values) for values in records.values())))

    widget = Lg_SeqMatrix()
    widget.resize(1180, 660)

    widget.sequence_data = {}
    widget.seqmatrix_data = {}
    for locus in loci:
        entries = []
        for taxon, sequence in records[locus].items():
            sequence_id = str(uuid.uuid4())
            entries.append({"uuid": sequence_id, "name": taxon, "sequence": sequence})
            widget.seqmatrix_data[sequence_id] = {"name": taxon, "sequence": sequence}
        widget.sequence_data.setdefault("Sanger amplicons", {})[f"{locus}.fasta"] = entries

    widget.seq_matrix.setColumnCount(1 + len(loci))
    widget.seq_matrix.setRowCount(len(taxa))
    widget.seq_matrix.setHorizontalHeaderItem(0, QTableWidgetItem("Sequence Name"))
    for column, locus in enumerate(loci, start=1):
        widget.seq_matrix.setHorizontalHeaderItem(column, QTableWidgetItem(locus))
    for row, taxon in enumerate(taxa):
        widget.seq_matrix.setItem(row, 0, QTableWidgetItem(taxon))
        for column, locus in enumerate(loci, start=1):
            if taxon not in records[locus]:
                continue
            sequence_id = next(
                entry["uuid"] for entry in widget.sequence_data["Sanger amplicons"][f"{locus}.fasta"]
                if entry["name"] == taxon
            )
            item = QTableWidgetItem(taxon)
            item.setData(Qt.UserRole, sequence_id)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            widget.seq_matrix.setItem(row, column, item)

    widget.seq_matrix.resizeColumnsToContents()
    widget.seq_matrix.setColumnWidth(0, 175)
    for column in range(1, 4):
        widget.seq_matrix.setColumnWidth(column, 155)
    widget.update_tree_view()
    widget.tree_explorer.collapseAll()
    if widget.tree_explorer.topLevelItemCount():
        widget.tree_explorer.topLevelItem(0).setExpanded(True)
    widget.update_highlighting()
    save(widget, "A_seqmatrix_sanger.png")
    widget.close()


def _seqmatrix_state(state: str, name: str) -> None:
    taxa = ["E. alderi", "E. fulvus", "O. cellarius", "P. liuae"]
    accessions = ["MK299689", "MK299693", "AY014116", "MN056414"]
    loci = ["18S", "ITS"]
    widget = Lg_SeqMatrix()
    widget.resize(980, 520)
    widget.seq_matrix.setColumnCount(1 + len(loci))
    widget.seq_matrix.setRowCount(len(taxa))
    widget.seq_matrix.setHorizontalHeaderItem(0, QTableWidgetItem("Sequence Name"))
    for column, locus in enumerate(loci, start=1):
        widget.seq_matrix.setHorizontalHeaderItem(column, QTableWidgetItem(locus))
    for row, taxon in enumerate(taxa):
        widget.seq_matrix.setItem(row, 0, QTableWidgetItem(taxon))

    if state == "ncbi":
        widget.sequence_data["Downloaded"] = {"NCBI Downloads": []}
        for row, accession in enumerate(accessions):
            for column in (1, 2):
                sequence_id = str(uuid.uuid4())
                widget.sequence_data["Downloaded"]["NCBI Downloads"].append(
                    {"uuid": sequence_id, "name": accession, "sequence": "ACGT"}
                )
                item = QTableWidgetItem(accession)
                item.setData(Qt.UserRole, sequence_id)
                widget.seq_matrix.setItem(row, column, item)
    else:
        duplicate_id = str(uuid.uuid4())
        widget.sequence_data["Local"] = {
            "amplicons.fasta": [
                {"uuid": duplicate_id, "name": "MK299689", "sequence": "ACGT"}
            ]
        }
        for row, accession in enumerate(accessions):
            for column in (1, 2):
                item = QTableWidgetItem(accession)
                item.setData(Qt.UserRole, str(uuid.uuid4()))
                widget.seq_matrix.setItem(row, column, item)
        for row, column in ((0, 1), (2, 2)):
            item = QTableWidgetItem("MK299689")
            item.setData(Qt.UserRole, duplicate_id)
            widget.seq_matrix.setItem(row, column, item)

    widget.seq_matrix.setColumnWidth(0, 175)
    widget.seq_matrix.setColumnWidth(1, 150)
    widget.seq_matrix.setColumnWidth(2, 150)
    widget.update_tree_view()
    widget.update_highlighting()
    save(widget, name)
    widget.close()


def render_seqmatrix_states() -> None:
    _seqmatrix_state("ncbi", "D_seqmatrix_ncbi_blue.png")
    _seqmatrix_state("duplicate", "E_seqmatrix_duplicate_red.png")


def render_mafft_dialog() -> None:
    dialog = MAFFTAdvancedDialog()
    dialog.resize(640, 480)
    dialog.adjust_direction_combo.setCurrentIndex(1)  # --adjustdirection
    dialog.adjust_direction_combo.setToolTip(
        "Handle potentially reversed sequencing reads (--adjustdirection)"
    )
    save(dialog, "B_mafft_adjustdirection.png")
    dialog.close()


def render_seqdbg() -> None:
    files = sorted(ORGANELLE.glob("*.gb"))
    widget = SeqDBGWindow()
    widget.resize(1420, 900)
    widget.show()
    QApplication.processEvents()
    widget.on_files_added([str(path) for path in files])

    deadline = time.time() + 180
    while widget.parser_thread is not None and widget.parser_thread.isRunning() and time.time() < deadline:
        QApplication.processEvents()
        time.sleep(0.05)
    QApplication.processEvents()
    wait(800)

    # Pick a gene with a well-populated adjacency table so the Sankey panel is
    # informative in the screenshot.
    if widget.table_viewer.stats_data:
        target = None
        for row, stat in enumerate(widget.table_viewer.stats_data):
            if stat.gene_name.lower() in {"cox1", "cob", "nad1", "nad5"}:
                target = row
                break
        if target is None:
            target = min(3, len(widget.table_viewer.stats_data) - 1)
        widget.table_viewer.on_table_cell_clicked(target, 0)
        QApplication.processEvents()
        wait(500)

    save(widget, "A_seqdbg_organelle.png")
    widget.close()


def render_phylobayes() -> None:
    alignment = BAYES / "prmA.phy"
    widget = PhyloBayesPlugin(workdir=str(BAYES))
    widget.resize(1080, 720)
    widget.file_path_edit.setText(str(alignment))
    widget.is_alignment_file.setChecked(True)
    widget.use_CAT_model.setChecked(True)
    widget.use_inf_CAT.setChecked(True)
    widget.use_Gamma_model.setChecked(True)
    widget.gamma_categories_spinbox.setValue(4)
    widget.subst_model_combo.setCurrentText("GTR")
    widget.n_parallel_chains.setValue(2)
    widget.n_mcmc_generations.setValue(10000)
    widget.sampling_frequency.setValue(10)
    widget.input_tabs.setCurrentIndex(0)
    widget.tab_widget.setCurrentIndex(0)
    save(widget, "A_phylobayes_infcat.png")
    widget.close()


def render_mrbayes() -> None:
    """MrBayes setup for the two-run chain used by the diagnostics panels."""
    alignment = Path("E:/YR-MPE_tutorial/YR-MPE_dev/dev/fig3/mrbayes/locus001.nex")
    widget = MrBayesPlugin(workdir=str(alignment.parent))
    widget.resize(1120, 820)
    widget.file_path_edit.setText(str(alignment))
    widget.datatype_combo.setCurrentText("DNA")
    widget.dna_rate_num_spinbox.setValue(6)
    widget.state_freq_pr_combo_dna.setCurrentText("estimated(dirichlet)")
    widget.rate_hetero_combo.setCurrentText("Gamma (+G)")
    widget.gamma_categories_spinbox.setValue(4)
    widget.generation_spinbox.setText("30000")
    widget.sampling_frequency_spinbox.setText("100")
    widget.run_num_spinbox.setValue(2)
    widget.chain_num_spinbox.setValue(4)
    widget.tab_widget.setCurrentIndex(0)
    save(widget, "A_mrbayes_setup.png")
    widget.close()


def _load_mcmc(widget: MiniTracerPlugin, paths) -> None:
    for path in paths:
        path = str(path)
        widget.mcmc_files.append(path)
        widget.mcmc_data[path] = widget.parse_trace_file(path)
    widget.update_file_table()
    widget.update_stats_table()


def _select_param(widget: MiniTracerPlugin, name: str) -> int:
    for row in range(widget.stats_table.rowCount()):
        item = widget.stats_table.item(row, 0)
        if item and item.text() == name:
            widget.stats_table.selectRow(row)
            return row
    return 0


def render_minitracer_single() -> None:
    trace = Path("E:/YR-MPE_tutorial/YR-MPE_dev/dev/fig3/mrbayes/locus001.nex.run1.p")
    widget = MiniTracerPlugin()
    widget.resize(1220, 820)
    _load_mcmc(widget, [trace])
    widget.file_table.selectRow(0)
    widget.update_stats_table()
    row = _select_param(widget, "lnLike")
    widget.show_trace_plot_for_selection(0, row)
    widget.tab_widget.setCurrentIndex(0)
    save(widget, "B_minitracer_trace_density.png")
    widget.close()


def render_minitracer_compare() -> None:
    run1 = Path("E:/YR-MPE_tutorial/YR-MPE_dev/dev/fig3/mrbayes/locus001.nex.run1.p")
    run2 = Path("E:/YR-MPE_tutorial/YR-MPE_dev/dev/fig3/mrbayes/locus001.nex.run2.p")
    widget = MiniTracerPlugin()
    widget.resize(1220, 820)
    _load_mcmc(widget, [run1, run2])
    widget.file_table.clearSelection()
    widget.file_table.selectRow(0)
    widget.file_table.selectRow(1)
    widget.update_stats_table()
    row = _select_param(widget, "lnLike")
    widget.update_compare_chains_analysis(0, 1)
    widget.tab_widget.setCurrentIndex(1)
    save(widget, "C_minitracer_compare_psrf.png")
    widget.close()


def render_icytree(newick_path: Path, name: str) -> None:
    """Load a case tree into the embedded IcyTree viewer."""
    widget = IcyTreePlugin(plugin_path=str(REPO / "YR_MPE"))
    widget.resize(1120, 760)
    widget.show()
    QApplication.processEvents()
    wait(4000)
    widget.load_tree_from_file(newick_path.read_text(encoding="utf-8"))
    QApplication.processEvents()
    wait(4000)
    save(widget, name)
    widget.close()


def render_single_gene() -> None:
    """ITS2 single-gene case: alignment viewer, ModelFinder table and tree."""
    root = HERE / "single_gene"
    alignment = root / "its2_trimmed.fasta"
    report = root / "its2.iqtree"
    tree = root / "its2.treefile"

    viewer = SequenceAlignmentViewer()
    viewer.resize(1220, 820)
    viewer.load_alignment(str(alignment))
    save(viewer, "A_its2_alignment.png")
    viewer.close()

    finder = ModelFinderPlugin()
    finder.resize(1120, 760)
    finder.display_results([str(report)])
    finder.tab_widget.setCurrentIndex(1)
    save(finder, "B_its2_modelfinder.png")
    finder.close()

    render_icytree(tree, "C_its2_tree.png")


def main() -> int:
    app = QApplication(sys.argv)
    apply_ui_font(app)

    render_single_gene()
    render_distance_case()
    render_seqmatrix()
    render_seqmatrix_states()
    render_mafft_dialog()
    render_icytree(TASK02 / "submission" / "tree.newick", "C_sanger_tree.png")
    render_seqdbg()
    render_icytree(Path("E:/YR-MPE_tutorial/benchmark/open_tasks/task03/submission/tree.nwk"),
                   "B_organelle_tree.png")
    render_phylobayes()
    render_mrbayes()
    render_minitracer_single()
    render_minitracer_compare()
    print("assets:", ASSETS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
