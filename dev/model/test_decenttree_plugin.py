"""Tests for the DecentTree plugin migration."""

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

from YR_MPE.plugins.decenttree_plugin import DecentTreePlugin, DecentTreeThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
DECENT = REPO / "YR_MPE" / "bin" / "decenttree" / "decenttree.exe"


class DecentTreePluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = DecentTreePlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "DecentTree")
        self.assertIn("wang2023decenttree", self.plugin.citation_keys)
        self.assertEqual(self.plugin.run_button.text(), "Build tree")

    def test_result_spec(self):
        tree = self.root / "tree.nwk"
        tree.write_text("(A:0.1,B:0.1,C:0.1);\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "decenttree")
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "phylogeny")
        self.assertEqual(item.format, "newick")
        self.assertIn("newick", item.data)

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import decenttree_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(DECENT.is_file(), "DecentTree runtime not available")
class DecentTreeThreadTests(unittest.TestCase):
    def test_tree_written(self):
        root = Path(tempfile.mkdtemp())
        matrix = root / "input.mldist"
        matrix.write_text(
            "4\n"
            "a 0 0.1 0.2 0.3\n"
            "b 0.1 0 0.4 0.5\n"
            "c 0.2 0.4 0 0.6\n"
            "d 0.3 0.5 0.6 0\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = DecentTreeThread(
            str(DECENT), [str(matrix)], ["-t", "BIONJ-R", "-no-banner"],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        tree = Path(payloads[0]["output_files"][0])
        self.assertTrue(tree.is_file())
        self.assertIn("(", tree.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
