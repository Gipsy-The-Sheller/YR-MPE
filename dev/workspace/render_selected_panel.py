"""Render the right-side Selected tab with a gCF-style basket."""

from __future__ import annotations

import os
import sys
from pathlib import Path

if os.environ.get("QT_QPA_PLATFORM") is None:
    os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "offscreen"

from PyQt5.QtCore import Qt  # noqa: E402
from PyQt5.QtGui import QFont, QFontDatabase  # noqa: E402
from PyQt5.QtTest import QTest  # noqa: E402
from PyQt5.QtWidgets import QApplication  # noqa: E402

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_TREESET, ITEM_TYPE_PHYLOGENY,
)
from YR_MPE.platforms.mpea import YR_MPEA_Widget  # noqa: E402


def make_item(manager, dataset_id, item_type, name, data):
    item = DatasetItem(item_type=item_type)
    item.dataset_id = dataset_id
    item.loci_name = name
    item.name = name
    item.data = data
    manager.add_item(item, dataset_id)
    return manager.get_item(item.id)


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

    from dev.workspace.preview_live import load_demo
    manager = load_demo(reset=True)
    dataset_id = next(iter(manager.datasets))
    gene_trees = make_item(
        manager, dataset_id, ITEM_TYPE_TREESET, "Gene trees",
        {"data": [{"content": "(A,B);"}, {"content": "(A,C);"}]},
    )
    species_tree = make_item(
        manager, dataset_id, ITEM_TYPE_PHYLOGENY, "Species tree",
        {"newick": "(A,B,C);"},
    )
    manager.add_to_basket(gene_trees.id)
    manager.add_to_basket(species_tree.id)

    window = YR_MPEA_Widget()
    window.setWindowTitle("YPhylo Workbench — Selected basket")
    window.dataset_selection_manager = manager
    window.workspace.dataset_selection_manager = manager
    window.workspace.rebuild_from_manager()
    window.inspector.set_sources(manager, manager.store)
    window.selected_panel.set_sources(manager, manager.store)
    window.resize(1250, 720)
    window.show()
    QTest.qWait(120)
    window._open_selected_panel()
    QTest.qWait(320)
    window.grab().save(str(HERE.parent / "screenshots" / "09_selected_basket.png"))
    window.close()
    print("screenshot saved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
