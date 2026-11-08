"""Tests for the ML distance plugin migration."""

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

from YR_MPE.plugins.ml_distance_plugin import (  # noqa: E402
    MLDistancePlugin, MLDistanceThread,
)

APP = QApplication.instance() or QApplication(sys.argv)
IQTREE = REPO / "YR_MPE" / "bin" / "iqtree3" / "iqtree3.exe"


def records(*ids):
    return [SeqRecord(Seq("ATGCAAATGCAA"), id=name) for name in ids]


class MLDistancePluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = MLDistancePlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "IQ-TREE 3")
        self.assertIn("wong2026iqtree", self.plugin.citation_keys)
        self.assertEqual(self.plugin.run_button.text(), "Run ML distance")

    def test_import_once(self):
        self.plugin.close()
        self.plugin = MLDistancePlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)

    def test_result_spec(self):
        matrix = self.root / "input.mldist"
        matrix.write_text(
            "2\nA 0.0 0.1\nB 0.1 0.0\n", encoding="utf-8"
        )
        self.plugin.current_output_files = [str(matrix)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "ml_distance")
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "distance")
        self.assertEqual(item.data["labels"], ["A", "B"])
        self.assertEqual(item.data["matrix"][0][1], 0.1)

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import ml_distance_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(IQTREE.is_file(), "IQ-TREE runtime not available")
class MLDistanceThreadTests(unittest.TestCase):
    def test_mldist_written(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGCAAATGCAA\n>b\nATGCAAATGCAA\n>c\nATGCAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = MLDistanceThread(
            str(IQTREE), [str(source)],
            ["-m", "JC", "-nt", "1", "-redo"],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        matrix = Path(payloads[0]["output_files"][0])
        self.assertTrue(matrix.is_file())
        self.assertTrue(matrix.read_text(encoding="utf-8").strip())


if __name__ == "__main__":
    unittest.main(verbosity=2)
