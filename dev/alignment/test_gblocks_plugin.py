"""Tests for the GBlocks plugin migration."""

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

from YR_MPE.plugins.gblocks_plugin import GBlocksPlugin, GBlocksThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
GBLOCKS = REPO / "YR_MPE" / "bin" / "Gblocks" / "Gblocks.exe"


def records(*ids):
    return [SeqRecord(Seq("ATGAAAATGCAA"), id=name) for name in ids]


class GBlocksPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = GBlocksPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata_and_parameters(self):
        self.assertEqual(self.plugin.tool_name, "GBlocks")
        self.assertIn("talavera2007improvement", self.plugin.citation_keys)
        params = self.plugin.get_parameters()
        self.assertTrue(any(p.startswith("-b6=") for p in params))
        self.assertEqual(self.plugin.run_button.text(), "Run GBlocks")

    def test_import_once(self):
        self.plugin.close()
        self.plugin = GBlocksPlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)

    def test_result_spec(self):
        output = self.root / "x_gblocks.fas"
        output.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        self.plugin.alignment_output_files = [str(output)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "gblocks")
        self.assertEqual(result.datasets[0].items[0].sequence_count, 2)

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import gblocks_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(GBLOCKS.is_file(), "GBlocks runtime not available")
class GBlocksThreadTests(unittest.TestCase):
    def test_gblocks_writes_trimmed_fasta(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGAAAATGCAA\n>b\nATGAAAATGCAA\n>c\nATGAAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = GBlocksThread(
            str(GBLOCKS), [str(source)],
            ["-t=d", "-b3=8", "-b4=10", "-b6=n"],
            workdir=str(root), suffix="_gblocks",
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
