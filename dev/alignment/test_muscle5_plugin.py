"""Tests for the MUSCLE5 plugin migration."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.plugins.muscle5_plugin import Muscle5Plugin, Muscle5Thread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
MUSCLE = REPO / "YR_MPE" / "bin" / "muscle5" / "muscle.exe"


def records(*ids):
    return [SeqRecord(Seq("ATGCAAATGCAA"), id=name) for name in ids]


class Muscle5PluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = Muscle5Plugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "Muscle5")
        self.assertIn("edgar2022muscle5", self.plugin.citation_keys)
        self.assertTrue(self.plugin.tool_path)
        self.assertEqual(self.plugin.get_parameters(), ["-align"])
        self.plugin.mode_combo.setCurrentIndex(1)
        self.assertEqual(self.plugin.get_parameters(), ["-super5"])

    def test_import_list_once(self):
        self.plugin.close()
        self.plugin = Muscle5Plugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)
        self.assertFalse(self.plugin.batch_mode)

    def test_import_batch(self):
        self.plugin.close()
        self.plugin = Muscle5Plugin(
            import_from="YR_MPEA",
            import_data=[records("a", "b"), records("c", "d")],
            workdir=str(self.root),
        )
        self.assertTrue(self.plugin.batch_mode)
        self.assertEqual(len(self.plugin.imported_files), 2)

    def test_result_spec(self):
        output = self.root / "locus_m5.fas"
        output.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        self.plugin.alignment_output_files = [str(output)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "muscle5")
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "alignment")
        self.assertEqual(item.sequence_count, 2)

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import muscle5_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(MUSCLE.is_file(), "MUSCLE5 runtime not available")
class Muscle5ThreadTests(unittest.TestCase):
    def test_muscle_writes_aligned_fasta(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGCAAATGCAA\n>b\nATGCAAATGCAA\n>c\nATGCAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = Muscle5Thread(
            str(MUSCLE), [str(source)], ["-align"],
            workdir=str(root), suffix="_m5",
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        output = Path(payloads[0]["output_files"][0])
        self.assertTrue(output.is_file())
        self.assertIn(">a", output.read_text(encoding="utf-8").lower())


if __name__ == "__main__":
    unittest.main(verbosity=2)
