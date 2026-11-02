"""Render the workspace dataset control and Dataset Manager for review."""

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
DEV_ROOT = HERE.parents[2]
sys.path.insert(0, str(DEV_ROOT))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.factories import ResourceFactory  # noqa: E402
from YR_MPE.platforms.methods.dataset_manager import DatasetManager  # noqa: E402
from YR_MPE.platforms.methods.dataset_selection_manager import DatasetSelectionManager  # noqa: E402
from YR_MPE.platforms.mpea import SingleGeneWorkspace  # noqa: E402
from YR_MPE.runtime.results import DatasetItemSpec, DatasetSpec  # noqa: E402


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

    out_dir = HERE.parent / "screenshots"
    out_dir.mkdir(parents=True, exist_ok=True)

    manager = DatasetSelectionManager(disable_auto_save=True)
    workspace = SingleGeneWorkspace(
        resource_factory=ResourceFactory(),
        dataset_selection_manager=manager,
    )

    items = []
    for index in range(1, 6):
        og = "OG%07d" % index
        records = [
            SeqRecord(Seq("ATGCAA" * 20), id="%s_%s" % (species, og))
            for species in ("HM347959", "JN637766", "KF150644", "KJ806277")
        ]
        items.append(
            DatasetItemSpec(
                loci_name=og,
                item_type="sequence",
                sequences=records,
                file_path=str(HERE.parent / "testdata" / ("%s.cds.fna" % og)),
                is_aligned=False,
                data={"orthogroup": og},
            )
        )

    spec = DatasetSpec(
        name="OrthoFinder_plugin_full",
        items=items,
        description="4 species, 335 CDS, 80 orthogroups, 5 shown",
        source="one_step_orthofinder",
        summary={"species": 4, "cds": 335, "orthogroups": 80, "single_copy_ogs": 5},
        metadata={"input": "dev/gbff_to_orthofinder/testdata/gb", "output": "plugin_full"},
    )
    dataset_id = workspace.add_result_spec(spec)[0]
    workspace.current_dataset_id = dataset_id

    workspace.resize(760, 240)
    workspace.show()
    for _ in range(5):
        app.processEvents()
    workspace.grab().save(str(out_dir / "01_workspace_dataset.png"))

    dialog = DatasetManager(dataset_name=spec.name, workspace=workspace)
    dialog.resize(1000, 620)
    dialog.show()
    app.processEvents()
    dialog.grab().save(str(out_dir / "02_dataset_manager.png"))

    dialog.close()
    workspace.close()
    print("wrote", out_dir / "01_workspace_dataset.png")
    print("wrote", out_dir / "02_dataset_manager.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
