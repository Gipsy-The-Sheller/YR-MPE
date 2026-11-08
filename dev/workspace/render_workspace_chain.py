"""Render hover/active chain feedback in the single-gene workspace."""

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

from YR_MPE.factories import ResourceFactory  # noqa: E402
from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_MODEL, ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY,
)
from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import SingleGeneWorkspace  # noqa: E402


def make_item(manager, dataset_id, chain, name, item_type, data):
    item = DatasetItem(item_type=item_type)
    item.dataset_id = dataset_id
    item.loci_name = name
    item.name = name
    item.ownership_uuid = chain
    item.data = data
    assert manager.add_item(item, dataset_id)
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

    manager = DatasetSelectionManager(disable_auto_save=True)
    dataset_id = manager.create_dataset("Chain preview")
    workspace = SingleGeneWorkspace(
        resource_factory=ResourceFactory(),
        dataset_selection_manager=manager,
    )
    workspace.current_dataset_id = dataset_id
    workspace.workspace_hint.setText("Hover or click a data icon")
    workspace.resize(720, 360)

    data = {
        "model-a": (ITEM_TYPE_MODEL, "Model A", "file/model.svg", 1, 0),
        "dist-a": (ITEM_TYPE_DISTANCE, "Distance A", "file/distance.svg", 2, 0),
        "tree-a": (ITEM_TYPE_PHYLOGENY, "Tree A", "file/phylogeny.svg", 3, 0),
        "model-b": (ITEM_TYPE_MODEL, "Model B", "file/model.svg", 1, 1),
        "dist-b": (ITEM_TYPE_DISTANCE, "Distance B", "file/distance.svg", 2, 1),
        "tree-b": (ITEM_TYPE_PHYLOGENY, "Tree B", "file/phylogeny.svg", 3, 1),
    }
    items = {}
    for name, (item_type, label, icon, row, col) in data.items():
        chain = "chain-a" if name.endswith("-a") else "chain-b"
        item = make_item(manager, dataset_id, chain, label, item_type,
                         {"name": label})
        items[name] = item
        workspace._create_item_button(item, "models", row, col, icon)

    out_dir = HERE.parent / "screenshots"
    out_dir.mkdir(parents=True, exist_ok=True)
    workspace.show()
    app.processEvents()
    workspace.grab().save(str(out_dir / "01_chain_idle.png"))

    workspace._on_item_hovered(items["model-a"].id)
    app.processEvents()
    workspace.grab().save(str(out_dir / "02_chain_hover.png"))

    workspace._clear_hover_chain()
    workspace._on_item_single_click(items["model-a"].id)
    app.processEvents()
    workspace.grab().save(str(out_dir / "03_chain_active.png"))

    workspace._on_item_hovered(items["model-b"].id)
    app.processEvents()
    workspace.grab().save(str(out_dir / "04_other_chain_hover.png"))

    workspace._clear_hover_chain()
    workspace._on_item_single_click(items["model-b"].id)
    app.processEvents()
    workspace.grab().save(str(out_dir / "05_chain_switched.png"))

    # Selecting a child keeps only the child green; its ancestor model is
    # green internally but must render as blue chain context.
    workspace._on_item_single_click(items["dist-a"].id)
    app.processEvents()
    workspace.grab().save(str(out_dir / "06_click_distance.png"))

    workspace.close()
    print("screenshots:", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
