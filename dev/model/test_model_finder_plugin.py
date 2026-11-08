"""Tests for the ModelFinder plugin migration."""

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

from YR_MPE.plugins.model_finder_plugin import (  # noqa: E402
    ModelFinderPlugin, ModelFinderThread,
)

APP = QApplication.instance() or QApplication(sys.argv)
IQTREE = REPO / "YR_MPE" / "bin" / "iqtree3" / "iqtree3.exe"


def records(*ids):
    return [SeqRecord(Seq("ATGCAAATGCAA"), id=name) for name in ids]


class ModelFinderPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = ModelFinderPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "IQ-Tree 3")
        self.assertIn("kalyaanamoorthy2017modelfinder", self.plugin.citation_keys)
        self.assertIn("-mset", self.plugin.get_parameters())
        self.assertEqual(self.plugin.run_button.text(), "Run ModelFinder")

    def test_import_list_once(self):
        self.plugin.close()
        self.plugin = ModelFinderPlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)

    def test_result_spec_from_iqtree(self):
        report = self.root / "input.fas.iqtree"
        report.write_text(
            "List of models sorted by BIC scores:\n\n"
            "Model  LogL  AIC  w-AIC  AICc  w-AICc  BIC  w-BIC\n"
            "JC  -10.0  20.0 +  0.5  22.0 +  0.5  24.0 +  0.6\n"
            "K2P  -9.0  19.0 +  0.5  21.0 +  0.5  23.0 +  0.4\n\n",
            encoding="utf-8",
        )
        self.plugin.analysis_output_files = [str(report)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "modelfinder")
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "model")
        self.assertEqual(item.file_path, str(report))
        self.assertEqual(item.data["best_model"], "JC")

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import model_finder_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(IQTREE.is_file(), "IQ-TREE runtime not available")
class ModelFinderThreadTests(unittest.TestCase):
    def test_modelfinder_writes_iqtree(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGCAAATGCAA\n>b\nATGCAAATGCAA\n>c\nATGCAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = ModelFinderThread(
            str(IQTREE), [str(source)],
            ["-st", "DNA", "-mset", "ALL", "-nt", "1"],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        report = Path(payloads[0]["output_files"][0])
        self.assertTrue(report.is_file())
        self.assertIn("List of models", report.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
