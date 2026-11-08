"""Capture the click-to-open inspector animation in the full workbench."""

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
    DatasetItem, ITEM_TYPE_MODEL, ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY,
)
from YR_MPE.platforms.mpea import YR_MPEA_Widget  # noqa: E402


def add_item(manager, dataset_id, chain, name, item_type, data):
    item = DatasetItem(item_type=item_type)
    item.dataset_id = dataset_id
    item.loci_name = name
    item.name = name
    item.ownership_uuid = chain
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

    window = YR_MPEA_Widget()
    window.resize(1250, 720)
    window.workspace.workspace_hint.setText("Click a data icon to open the inspector")
    window.show()
    QTest.qWait(120)

    manager = window.dataset_selection_manager
    dataset_id = manager.create_dataset("Inspector preview")
    window.workspace.current_dataset_id = dataset_id
    rows = [
        ("Model A", ITEM_TYPE_MODEL, "chain-a", "file/model.svg", 1, 0,
         {"name": "Model A", "model": "GTR+F+I+G4"}),
        ("Distance A", ITEM_TYPE_DISTANCE, "chain-a", "file/distance.svg", 2, 0,
         {"method": "ML"}),
        ("Tree A", ITEM_TYPE_PHYLOGENY, "chain-a", "file/phylogeny.svg", 3, 0,
         {"newick": "(A,B);"}),
        ("Model B", ITEM_TYPE_MODEL, "chain-b", "file/model.svg", 1, 1,
         {"name": "Model B", "model": "HKY+F+G4"}),
        ("Distance B", ITEM_TYPE_DISTANCE, "chain-b", "file/distance.svg", 2, 1,
         {"method": "p-distance"}),
    ]
    first = None
    for label, item_type, chain, icon, row, col, data in rows:
        item = add_item(manager, dataset_id, chain, label, item_type, data)
        window.workspace._create_item_button(item, "models", row, col, icon)
        first = first or item

    out_dir = HERE.parent / "screenshots"
    out_dir.mkdir(parents=True, exist_ok=True)
    QTest.qWait(60)
    window.grab().save(str(out_dir / "06_inspector_closed.png"))

    window.workspace._on_item_single_click(first.id)
    QTest.qWait(90)
    window.grab().save(str(out_dir / "07_inspector_mid.png"))
    QTest.qWait(220)
    window.grab().save(str(out_dir / "08_inspector_open.png"))

    window.close()
    print("screenshots:", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
