"""Render the OrthoFinder plugin panel to PNG files.

The renderer uses the host platform style and font; it does not set a QStyle.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Use the host platform/style; never force Fusion or another QStyle.
# YRTools owns the application stylesheet and font.
if os.environ.get("QT_QPA_PLATFORM") is None:
    os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "offscreen"

from PyQt5.QtCore import Qt  # noqa: E402
from PyQt5.QtGui import QFont, QFontDatabase  # noqa: E402
from PyQt5.QtWidgets import QApplication  # noqa: E402

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

HERE = Path(__file__).resolve()
DEV_ROOT = HERE.parents[2]  # YR-MPE_dev
sys.path.insert(0, str(DEV_ROOT))

from YR_MPE.plugins.orthofinder_plugin import OrthoFinderPlugin  # noqa: E402


def main() -> int:
    app = QApplication(sys.argv)
    # Mirror YRTools: load MiSans Medium and use it as the application font.
    font_path = Path("E:/Codes/YRTools/fonts/MiSans-Medium.ttf")
    if font_path.exists():
        font_id = QFontDatabase.addApplicationFont(str(font_path))
        families = QFontDatabase.applicationFontFamilies(font_id) if font_id != -1 else []
        if families:
            app_font = QFont(families[0], 10)
            app_font.setWeight(QFont.Medium)
            app.setFont(app_font)
            print("font:", app_font.toString())

    out_dir = HERE.parent / "screenshots"
    out_dir.mkdir(parents=True, exist_ok=True)

    plugin = OrthoFinderPlugin(workdir=str(HERE.parent / "testdata"))
    plugin.resize(1280, 900)

    input_dir = HERE.parent / "testdata" / "gb"
    output_dir = HERE.parent / "testdata" / "plugin_full"
    plugin.input_edit.setText(str(input_dir))
    plugin.output_edit.setText(str(output_dir))
    plugin._scan_input()

    # Step 1: Input
    plugin._set_step(0)
    plugin.show()
    app.processEvents()
    plugin.grab().save(str(out_dir / "01_input.png"))

    # Multi-record single-file mode
    multi_dir = HERE.parent / "testdata" / "gb_multirecord"
    multi_dir.mkdir(parents=True, exist_ok=True)
    multi_file = multi_dir / "all.gb"
    if not multi_file.exists():
        with multi_file.open("wb") as handle:
            for name in ("HM347959.gb", "JN637766.gb"):
                handle.write((input_dir / name).read_bytes())
    plugin.source_combo.setCurrentIndex(1)
    plugin.input_edit.setText(str(multi_file))
    plugin.output_edit.setText(str(HERE.parent / "testdata" / "plugin_multirecord"))
    plugin._scan_input()
    plugin._set_step(0)
    app.processEvents()
    plugin.grab().save(str(out_dir / "05_input_multirecord.png"))

    # Restore folder mode for the remaining screenshots
    plugin.source_combo.setCurrentIndex(0)
    plugin.input_edit.setText(str(input_dir))
    plugin.output_edit.setText(str(output_dir))
    plugin._scan_input()

    # Step 2: OrthoFinder
    plugin._set_step(1)
    plugin._check_runtime(show=False)
    app.processEvents()
    plugin.grab().save(str(out_dir / "02_orthofinder.png"))

    # Step 3: Review
    plugin._set_step(2)
    app.processEvents()
    plugin.grab().save(str(out_dir / "03_review.png"))

    # Output tab with sample results (simulate a YR-MPEA host connection so
    # the export control is visible in the screenshot).
    plugin.tab_widget.setCurrentIndex(1)
    plugin._opened_from_yrmpea = True
    plugin.result_ready.connect(lambda result: None)
    plugin.update_results({
        "species": 4,
        "cds": 335,
        "orthogroups": 80,
        "single_copy_ogs": 62,
        "output_dir": str(output_dir),
        "elapsed": 19.2,
        "run_mode": "full",
        "cds_mapping": True,
    })
    app.processEvents()
    plugin.grab().save(str(out_dir / "04_output.png"))

    # Citation tab: parsed from citations/software.bib via config.json
    plugin.tab_widget.setCurrentIndex(3)
    app.processEvents()
    plugin.grab().save(str(out_dir / "06_citation.png"))

    plugin.close()
    print("wrote:")
    for path in sorted(out_dir.glob("*.png")):
        print(" ", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
