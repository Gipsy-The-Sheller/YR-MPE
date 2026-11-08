"""Headless tests for the codon Z-test and batch dN/dS plugins."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.plugins.batch_dnds_plugin import BatchDnDsPlugin  # noqa: E402
from YR_MPE.plugins.codon_ztest_plugin import CodonZTestPlugin  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


def write_alignment(path: Path, names=("A", "B", "C", "D")):
    sequences = {
        "A": "ATG" + "AAA" * 10 + "CTG" * 10,
        "B": "ATG" + "AAG" * 10 + "CTG" * 10,
        "C": "ATG" + "AAA" * 10 + "CTT" * 10,
        "D": "ATG" + "AAG" * 10 + "CTC" * 10,
    }
    with path.open("w", encoding="utf-8") as handle:
        for name in names:
            handle.write(">%s\n%s\n" % (name, sequences[name]))


class SelectionPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.alignment = self.root / "codon.fas"
        write_alignment(self.alignment)

    def _wait(self, plugin, timeout=60000):
        if plugin.analysis_thread:
            plugin.analysis_thread.wait(timeout)
            for _ in range(30):
                APP.processEvents()

    def test_codon_ztest(self):
        plugin = CodonZTestPlugin(workdir=str(self.root))
        try:
            plugin.alignment_edit.setText(str(self.alignment))
            plugin.run_analysis()
            self._wait(plugin)
            self.assertEqual(plugin.table.rowCount(), 6)
            self.assertIn("sequences", plugin.summary_label.text())
            conclusions = [plugin.table.item(row, 11).text() for row in range(6)]
            self.assertEqual(conclusions[0], "purifying selection")
            self.assertEqual(conclusions[1:3], ["not significant"] * 2)
            self.assertTrue(all(c == "saturated/not estimable"
                                for c in conclusions[3:]))
            self.assertEqual(plugin.matrix.rowCount(), 4)
            self.assertEqual(plugin.matrix.columnCount(), 4)
            # MEGA layout: p-values above the diagonal, omega below it.
            self.assertEqual(plugin.matrix.item(0, 0).text(), "—")
            self.assertTrue(plugin.matrix.item(0, 1).text().startswith("0."))
            self.assertEqual(plugin.matrix.item(1, 0).text(), "0")
            self.assertIn("0.05", plugin.view_hint.text())
            self.assertEqual(len(plugin._matrix_tsv().splitlines()), 5)
        finally:
            plugin.is_running = False
            plugin.close()

    def test_codon_yn00(self):
        plugin = CodonZTestPlugin(workdir=str(self.root))
        try:
            plugin.alignment_edit.setText(str(self.alignment))
            plugin.method_combo.setCurrentIndex(plugin.method_combo.findData("yn00"))
            plugin.run_analysis()
            self._wait(plugin)
            self.assertEqual(plugin.table.rowCount(), 6)
            self.assertIn("YN00", plugin.summary_label.text())
            dS_values = [plugin.table.item(row, 7).text() for row in range(6)]
            self.assertTrue(any("±" in value for value in dS_values))
        finally:
            plugin.is_running = False
            plugin.close()

    def test_batch_counting(self):
        folder = self.root / "alignments"
        folder.mkdir()
        write_alignment(folder / "gene1.fas")
        write_alignment(folder / "gene2.fas")
        plugin = BatchDnDsPlugin(workdir=str(self.root))
        try:
            plugin.folder_edit.setText(str(folder))
            plugin.run_analysis()
            self._wait(plugin)
            self.assertGreater(plugin.table.rowCount(), 0)
            self.assertIn("pairwise", plugin.summary_label.text())
        finally:
            plugin.is_running = False
            plugin.close()

    def test_batch_yn00(self):
        folder = self.root / "alignments_yn00"
        folder.mkdir()
        write_alignment(folder / "gene1.fas")
        write_alignment(folder / "gene2.fas")
        plugin = BatchDnDsPlugin(workdir=str(self.root))
        try:
            plugin.folder_edit.setText(str(folder))
            plugin.method_combo.setCurrentIndex(plugin.method_combo.findData("yn00"))
            plugin.run_analysis()
            self._wait(plugin)
            self.assertGreater(plugin.table.rowCount(), 0)
            self.assertEqual(plugin.table.item(0, 3).text(), "YN00")
            self.assertIn("pairwise", plugin.summary_label.text())
        finally:
            plugin.is_running = False
            plugin.close()

    def test_batch_mle_if_iqtree_available(self):
        plugin = BatchDnDsPlugin(workdir=str(self.root))
        try:
            if not plugin._iqtree_path():
                self.skipTest("IQ-TREE 3 is not configured")
            folder = self.root / "mle"
            folder.mkdir()
            write_alignment(folder / "gene.fas", names=("A", "B", "C", "D"))
            plugin.folder_edit.setText(str(folder))
            plugin.output_edit.setText(str(self.root / "mle_out"))
            plugin.method_combo.setCurrentIndex(plugin.method_combo.findData("gy"))
            plugin.run_analysis()
            self._wait(plugin, timeout=120000)
            self.assertEqual(plugin.table.rowCount(), 1)
            self.assertEqual(plugin.table.item(0, 9).text(), "ok")
            self.assertNotEqual(plugin.table.item(0, 6).text(), "NA")
        finally:
            plugin.is_running = False
            plugin.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
