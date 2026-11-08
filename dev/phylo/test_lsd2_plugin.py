"""Tests for the LSD2 plugin migration."""

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

from YR_MPE.plugins.lsd2_plugin import LSD2Plugin  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class LSD2PluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = LSD2Plugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "LSD2")
        self.assertIn("to2016fast", self.plugin.citation_keys)

    def test_result_spec_and_finish_preview(self):
        tree = self.root / "tree.nwk_result.nwk"
        tree.write_text("(A:1.0,B:1.0,(C:1.0,D:1.0):1.0);\n", encoding="utf-8")
        result = self.root / "tree.nwk.result"
        result.write_text("rate\n1.0\n", encoding="utf-8")

        self.plugin.on_lsd2_finished([str(tree), str(result)], [])
        self.assertEqual(len(self.plugin.current_output_files), 2)
        self.assertIn("A:1.0", self.plugin.output_preview.toPlainText())

        spec = self.plugin._build_result_spec()
        self.assertEqual(spec.validate(), [])
        self.assertEqual(spec.source, "lsd2")
        item = spec.datasets[0].items[0]
        self.assertEqual(item.item_type, "phylogeny")
        self.assertTrue(item.data.get("dated"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
