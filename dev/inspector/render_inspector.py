"""Standalone preview of the right-side InspectorPanel.

The panel is a plain QWidget; this script builds a temporary workspace store
with sample provenance and renders the panel without starting YR-MPEA.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

if os.environ.get("QT_QPA_PLATFORM") is None:
    os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "offscreen"

from PyQt5.QtCore import Qt  # noqa: E402
from PyQt5.QtGui import QFont, QFontDatabase  # noqa: E402
from PyQt5.QtWidgets import QApplication  # noqa: E402

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

HERE = Path(__file__).resolve()
DEV_ROOT = HERE.parents[2]
sys.path.insert(0, str(DEV_ROOT))

from YR_MPE.platforms.methods.inspector_panel import InspectorPanel  # noqa: E402
from YR_MPE.runtime.workspace_store import WorkspaceStore  # noqa: E402


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def build_sample_store(root: Path) -> WorkspaceStore:
    store = WorkspaceStore(str(root), name="YPhylo demo")
    dataset = store.add_dataset(
        "OrthoFinder_plugin_full",
        description="4 species, 335 CDS, 80 orthogroups",
        source="one_step_orthofinder",
        settings={"summary": {"species": 4, "cds": 335, "single_copy_ogs": 62}},
    )

    gb = _write(root / "input" / "all.gb", "LOCUS       demo\n")
    imported = store.new_import(
        [str(gb)],
        importer="One-Step OrthoFinder",
        parameters={"genetic_code": 5},
        copy=True,
        dataset_id=dataset.id,
    )

    fasta = _write(root / "input" / "OG0000008.cds.fna", ">a\nATGCAA\n>b\nATGCAA\n")
    og1 = store.add_item(
        dataset.id, "sequence", "OG0000008",
        source_path=str(fasta), copy=True, file_format="fasta",
        import_id=imported.id, length=6, sequence_count=2,
    )
    og2 = store.add_item(
        dataset.id, "sequence", "OG0000009",
        source_path=str(fasta), copy=True, file_format="fasta",
        import_id=imported.id, length=6, sequence_count=2,
    )
    og3 = store.add_item(
        dataset.id, "sequence", "OG0000010",
        source_path=str(fasta), copy=True, file_format="fasta",
        import_id=imported.id, length=6, sequence_count=2,
    )

    aln_path = _write(root / "runs" / "mafft" / "OG0000008.aln.fas", ">a\nATGCAA\n>b\nATGCAA\n")
    mafft = store.new_run(
        "mafft", "MAFFT", [og1.id], name="MAFFT run 1",
        tool_version="7.520", parameters={"auto": True},
    )
    aln = store.add_item(
        dataset.id, "alignment", "OG0000008 alignment",
        source_path=str(aln_path), copy=True, file_format="fasta",
        run_id=mafft.id, is_aligned=True, length=6, sequence_count=2,
    )
    aln.parent_id = og1.id
    store.finish_run(mafft.id, output_ids=[aln.id])

    tree1_path = _write(root / "runs" / "iqtree1" / "tree.nwk", "(a:0.1,b:0.1);\n")
    iqtree1 = store.new_run(
        "iqtree", "IQ-TREE", [aln.id], name="IQ-TREE run 1",
        tool_version="3.0.1", parameters={"model": "GTR+F+I+G4", "boot": 1000},
    )
    tree1 = store.add_item(
        dataset.id, "phylogeny", "Tree (GTR+F+I+G4)",
        source_path=str(tree1_path), copy=True, file_format="newick",
        run_id=iqtree1.id, length=0, sequence_count=2,
    )
    tree1.parent_id = aln.id
    store.finish_run(iqtree1.id, output_ids=[tree1.id])

    tree2_path = _write(root / "runs" / "iqtree2" / "tree.nwk", "(a:0.2,b:0.2);\n")
    iqtree2 = store.new_run(
        "iqtree", "IQ-TREE", [aln.id], name="IQ-TREE run 2",
        tool_version="3.0.1", parameters={"model": "HKY+F+G4", "boot": 1000},
    )
    tree2 = store.add_item(
        dataset.id, "phylogeny", "Tree (HKY+F+G4)",
        source_path=str(tree2_path), copy=True, file_format="newick",
        run_id=iqtree2.id, length=0, sequence_count=2,
    )
    tree2.parent_id = aln.id
    store.finish_run(iqtree2.id, output_ids=[tree2.id])

    store.add_stream(
        "OG0000008 pipeline", dataset.id,
        [mafft.id, iqtree1.id, iqtree2.id],
    )
    store.save()
    return store, tree1


def main() -> int:
    app = QApplication(sys.argv)
    font_path = Path("E:/Codes/YRTools/fonts/MiSans-Medium.ttf")
    if font_path.exists():
        font_id = QFontDatabase.addApplicationFont(str(font_path))
        families = QFontDatabase.applicationFontFamilies(font_id) if font_id != -1 else []
        if families:
            app_font = QFont(families[0], 10)
            app_font.setWeight(QFont.Medium)
            app.setFont(app_font)

    root = Path(tempfile.mkdtemp())
    store, tree_item = build_sample_store(root)

    out_dir = HERE.parent / "screenshots"
    out_dir.mkdir(parents=True, exist_ok=True)

    panel = InspectorPanel(store=store)
    panel.resize(400, 900)
    panel.show()
    panel.show_item(tree_item.id)
    for _ in range(6):
        app.processEvents()
    panel.grab().save(str(out_dir / "01_inspector_item.png"))

    panel._metadata_expanded = True
    panel._apply_metadata_mode(animated=False)
    for _ in range(4):
        app.processEvents()
    panel.grab().save(str(out_dir / "04_inspector_metadata_expanded.png"))
    panel._metadata_expanded = False
    panel._apply_metadata_mode(animated=False)

    panel.show_run(tree_item.run_id)
    for _ in range(4):
        app.processEvents()
    panel.grab().save(str(out_dir / "03_inspector_edge.png"))

    panel.show_dataset(tree_item.dataset_id)
    for _ in range(4):
        app.processEvents()
    panel.grab().save(str(out_dir / "02_inspector_dataset.png"))

    panel.close()
    print("wrote", out_dir / "01_inspector_item.png")
    print("wrote", out_dir / "02_inspector_dataset.png")
    print("wrote", out_dir / "03_inspector_edge.png")
    print("wrote", out_dir / "04_inspector_metadata_expanded.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
