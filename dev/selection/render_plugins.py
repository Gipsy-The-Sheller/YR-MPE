"""Render the two SELECTION plugin panels to PNG files.

Uses the host platform style and font; never forces a QStyle.
"""

from __future__ import annotations

import os
import sys
import time
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

from YR_MPE.plugins.batch_dnds_plugin import BatchDnDsPlugin  # noqa: E402
from YR_MPE.plugins.codon_ztest_plugin import CodonZTestPlugin  # noqa: E402

SEQUENCES = {
    "A": "ATG" + "AAA" * 10 + "CTG" * 10 + "TAT" * 5,
    "B": "ATG" + "AAG" * 10 + "CTG" * 10 + "TAC" * 5,
    "C": "ATG" + "AAA" * 10 + "CTT" * 10 + "TAT" * 5,
    "D": "ATG" + "AAG" * 10 + "CTC" * 10 + "TAC" * 5,
}


def write_alignment(path: Path, names=("A", "B", "C", "D")) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for name in names:
            handle.write(">%s\n%s\n" % (name, SEQUENCES[name]))


def wait_for(plugin, timeout=120.0) -> None:
    deadline = time.time() + timeout
    while plugin.is_running and time.time() < deadline:
        QApplication.processEvents()
        time.sleep(0.02)
    for _ in range(50):
        QApplication.processEvents()
        time.sleep(0.01)


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
            print("font:", app_font.toString())

    testdata = HERE.parent / "testdata"
    out_dir = HERE.parent / "screenshots"
    out_dir.mkdir(parents=True, exist_ok=True)
    alignment = testdata / "codon.fas"
    write_alignment(alignment)

    folder = testdata / "alignments"
    for index, names in enumerate((("A", "B", "C"), ("A", "B", "D"), ("B", "C", "D")), 1):
        write_alignment(folder / ("gene%d.fas" % index), names)

    ztest = CodonZTestPlugin(workdir=str(testdata))
    ztest.resize(1180, 780)
    ztest.alignment_edit.setText(str(alignment))
    ztest._scan_alignment()
    ztest.tab_widget.setCurrentIndex(0)
    ztest.show()
    app.processEvents()
    ztest.grab().save(str(out_dir / "01_ztest_input.png"))

    ztest.run_analysis()
    wait_for(ztest)
    ztest.tab_widget.setCurrentIndex(1)
    app.processEvents()
    ztest.grab().save(str(out_dir / "02_ztest_output.png"))

    ztest.view_combo.setCurrentIndex(1)  # MEGA-style matrix
    app.processEvents()
    ztest.grab().save(str(out_dir / "02b_ztest_matrix.png"))

    ztest.view_combo.setCurrentIndex(0)
    ztest.method_combo.setCurrentIndex(ztest.method_combo.findData("yn00"))
    ztest.run_analysis()
    wait_for(ztest)
    ztest.tab_widget.setCurrentIndex(1)
    app.processEvents()
    ztest.grab().save(str(out_dir / "02c_ztest_yn00.png"))
    ztest.is_running = False
    ztest.close()

    batch = BatchDnDsPlugin(workdir=str(testdata))
    batch.resize(1180, 780)
    batch.folder_edit.setText(str(folder))
    batch.output_edit.setText(str(testdata / "mle_out"))
    batch._scan_folder()
    batch.tab_widget.setCurrentIndex(0)
    batch.show()
    app.processEvents()
    batch.grab().save(str(out_dir / "03_batch_input.png"))

    batch.run_analysis()
    wait_for(batch)
    batch.tab_widget.setCurrentIndex(1)
    app.processEvents()
    batch.grab().save(str(out_dir / "04_batch_counting.png"))

    batch.tab_widget.setCurrentIndex(0)
    batch.method_combo.setCurrentIndex(batch.method_combo.findData("yn00"))
    batch.run_analysis()
    wait_for(batch)
    batch.tab_widget.setCurrentIndex(1)
    app.processEvents()
    batch.grab().save(str(out_dir / "05_batch_yn00.png"))

    batch.tab_widget.setCurrentIndex(0)
    batch.method_combo.setCurrentIndex(batch.method_combo.findData("gy"))
    app.processEvents()
    batch.grab().save(str(out_dir / "06_batch_mle_input.png"))

    batch.is_running = False
    batch.close()
    print("screenshots:", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
