"""Persistent, interactive workspace demo for the chain feedback.

The demo state is written to dev/workspace/demo_workspace/. The window stays
open; hover an icon to see the grey chain, click to keep it blue and slide the
inspector in. Pass ``--reset`` to recreate the demo state.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

if os.environ.get("QT_QPA_PLATFORM") is None:
    os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "offscreen"

from PyQt5.QtCore import Qt  # noqa: E402
from PyQt5.QtGui import QFont, QFontDatabase  # noqa: E402
from PyQt5.QtWidgets import QApplication  # noqa: E402

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_SEQUENCE, ITEM_TYPE_MODEL, ITEM_TYPE_DISTANCE,
    ITEM_TYPE_PHYLOGENY,
)
from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import YR_MPEA_Widget  # noqa: E402

DEMO_DIR = HERE.parent / "demo_workspace"


def add_item(manager, dataset_id, item_type, name, chain, parent_id=None, **kwargs):
    item = DatasetItem(item_type=item_type)
    item.dataset_id = dataset_id
    item.loci_name = name
    item.name = name
    item.ownership_uuid = chain
    item.parent_id = parent_id
    for key, value in kwargs.items():
        setattr(item, key, value)
    if not manager.add_item(item, dataset_id):
        raise RuntimeError("could not add %s" % name)
    return manager.get_item(item.id)


def create_demo(manager: DatasetSelectionManager) -> str:
    dataset_id = manager.create_dataset("Persisted chain demo")
    fasta_path = DEMO_DIR / "demo_sequences.fna"
    records = [
        SeqRecord(Seq("ATG" + "AAA" * 12 + "CTG" * 8), id="taxon_A"),
        SeqRecord(Seq("ATG" + "AAG" * 12 + "CTA" * 8), id="taxon_B"),
        SeqRecord(Seq("ATG" + "AAA" * 12 + "CTT" * 8), id="taxon_C"),
    ]
    from Bio import SeqIO
    SeqIO.write(records, str(fasta_path), "fasta")

    root = add_item(
        manager, dataset_id, ITEM_TYPE_SEQUENCE, "Sequence import",
        "demo-chain-a", sequences=list(records), sequence_count=len(records),
        length=len(records[0].seq), file_path=str(fasta_path), format="fasta",
    )
    model_a = add_item(
        manager, dataset_id, ITEM_TYPE_MODEL, "Model A", "demo-chain-a",
        parent_id=root.id, data={"name": "Model A", "model": "GTR+F+I+G4"},
    )
    distance_a = add_item(
        manager, dataset_id, ITEM_TYPE_DISTANCE, "Distance A", "demo-chain-a",
        parent_id=model_a.id, data={"method": "ML"},
    )
    add_item(
        manager, dataset_id, ITEM_TYPE_PHYLOGENY, "Tree A", "demo-chain-a",
        parent_id=distance_a.id, data={"newick": "(taxon_A,taxon_B,taxon_C);"},
    )

    model_b = add_item(
        manager, dataset_id, ITEM_TYPE_MODEL, "Model B", "demo-chain-b",
        data={"name": "Model B", "model": "HKY+F+G4"},
    )
    distance_b = add_item(
        manager, dataset_id, ITEM_TYPE_DISTANCE, "Distance B", "demo-chain-b",
        parent_id=model_b.id, data={"method": "p-distance"},
    )
    add_item(
        manager, dataset_id, ITEM_TYPE_PHYLOGENY, "Tree B", "demo-chain-b",
        parent_id=distance_b.id, data={"newick": "(taxon_A,taxon_C,taxon_B);"},
    )
    return dataset_id


def load_demo(reset: bool = False) -> DatasetSelectionManager:
    DEMO_DIR.mkdir(parents=True, exist_ok=True)
    state_file = DEMO_DIR / "dataset_selection_state.json"
    if reset and state_file.exists():
        state_file.unlink()
    manager = DatasetSelectionManager()
    manager.workspace_path = str(DEMO_DIR)
    if state_file.exists():
        manager._load_state()
    else:
        create_demo(manager)
        manager._save_state()
    return manager


def main() -> int:
    app = QApplication(sys.argv)
    font_path = Path("E:/Codes/YRTools/fonts/MiSans-Medium.ttf")
    if font_path.exists():
        font_id = QFontDatabase.addApplicationFont(str(font_path))
        families = QFontDatabase.applicationFontFamilies(font_id) if font_id != -1 else []
        if families:
            font = QFont(families[0], 10)
            font.setWeight(QFont.Medium)
            app.setFont(font)

    manager = load_demo(reset="--reset" in sys.argv)
    window = YR_MPEA_Widget()
    window.setWindowTitle("YPhylo Workbench — persisted chain demo")
    window.dataset_selection_manager = manager
    window.workspace.dataset_selection_manager = manager
    window.workspace.rebuild_from_manager()
    window.inspector.set_sources(manager, manager.store)
    window.inspector.show_empty()
    window.resize(1250, 720)
    window.show()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
