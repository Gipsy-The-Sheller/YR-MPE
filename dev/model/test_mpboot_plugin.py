"""Tests for the MPBoot plugin migration."""

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

from YR_MPE.plugins.mpboot_plugin import MPBootPlugin, MPBootThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
MPBOOT = REPO / "YR_MPE" / "bin" / "mpboot" / "mpboot.exe"


class MPBootPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = MPBootPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "MPBoot")
        self.assertIn("hoang2018mpboot", self.plugin.citation_keys)
        self.assertEqual(self.plugin.run_button.text(), "Run MPBoot")

    def test_result_spec(self):
        tree = self.root / "input.fas.treefile"
        tree.write_text("(a,b,c);\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "mpboot")
        self.assertEqual(result.datasets[0].items[0].item_type, "phylogeny")

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import mpboot_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(MPBOOT.is_file(), "MPBoot runtime not available")
class MPBootThreadTests(unittest.TestCase):
    def test_treefile_written(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n>b\nATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTTCT\n>c\nATGGCTGCTGCTGTTGCTGCTGCTGCTGCTGCTGCT\n>d\nATGGCTGCTGCTGCTGCTCCTGCTGCTGATGCTGCT\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = MPBootThread(
            str(MPBOOT), [str(source)], [],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        files = [Path(path).name for path in payloads[0]["output_files"]]
        self.assertTrue(any(name.endswith(".treefile") for name in files))


if __name__ == "__main__":
    unittest.main(verbosity=2)
