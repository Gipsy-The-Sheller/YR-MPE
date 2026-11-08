"""Tests for the TrimAl plugin migration."""

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

from YR_MPE.plugins.trimal_plugin import TrimAlPlugin, TrimAlThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
TRIMAL = REPO / "YR_MPE" / "bin" / "trimAl" / "bin" / "trimal.exe"


def records(*ids):
    return [SeqRecord(Seq("ATGCAA"), id=name) for name in ids]


class TrimAlPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = TrimAlPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "TrimAl")
        self.assertIn("capellagutierrez2009trimal", self.plugin.citation_keys)
        self.assertIn("-automated1", self.plugin.get_parameters())
        self.assertEqual(self.plugin.run_button.text(), "Run TrimAl")

    def test_import_once(self):
        self.plugin.close()
        self.plugin = TrimAlPlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)
        self.assertFalse(self.plugin.batch_mode)

    def test_result_requires_fasta_output(self):
        output = self.root / "x_trimal.fas"
        output.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        self.plugin.alignment_output_files = [str(output)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "trimal")
        self.plugin.output_format_combo.setCurrentIndex(2)
        with self.assertRaises(ValueError):
            self.plugin._build_result_spec()

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import trimal_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(TRIMAL.is_file(), "TrimAl runtime not available")
class TrimAlThreadTests(unittest.TestCase):
    def test_trimal_writes_trimmed_fasta(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGAAAATGCAA\n>b\nATGAAAATGCAA\n>c\nATGAAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = TrimAlThread(
            str(TRIMAL), [str(source)], ["-automated1"],
            workdir=str(root), suffix="_trimal",
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
