"""Smoke tests for the rewritten initializer (no network access)."""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt  # noqa: E402
from PyQt5.QtWidgets import QApplication  # noqa: E402

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import Initializer as initializer_module  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class InitializerTests(unittest.TestCase):
    def setUp(self):
        self._info = initializer_module.QMessageBox.information
        self._warning = initializer_module.QMessageBox.warning
        initializer_module.QMessageBox.information = lambda *args, **kwargs: None
        initializer_module.QMessageBox.warning = lambda *args, **kwargs: None
        self.widget = initializer_module.Initializer()

    def tearDown(self):
        initializer_module.QMessageBox.information = self._info
        initializer_module.QMessageBox.warning = self._warning
        self.widget.close()

    def test_template_and_sources(self):
        self.assertEqual(len(self.widget.tools_data), 20)
        names = [tool["name"] for tool in self.widget.tools_data]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(self.widget.source_combo.count(), 3)
        self.assertEqual(self.widget.source_combo.itemText(0), "GitHub")
        self.assertEqual(self.widget.source_combo.itemText(1), "Gitee")
        self.assertEqual(self.widget.source_combo.itemText(2), "Custom URL")

    def test_bundled_paths_resolve(self):
        self.widget.auto_config()
        configured = {tool["name"]: tool["path"] for tool in self.widget.configured_tools}
        self.assertEqual(configured["MAFFT"], "/mafft/mafft.bat")
        self.assertEqual(configured["IQ-TREE 3"], "/iqtree3/iqtree3.exe")
        self.assertEqual(configured["MrBayes-MPI-BEAGLE"], "/mb/mb.exe")
        self.assertEqual(configured["caster-site"], "/aster/caster-site.exe")

    def test_archive_extraction(self):
        root = Path(tempfile.mkdtemp())
        (root / "YR_MPE").mkdir()
        archive = root / "bins.zip"
        with zipfile.ZipFile(archive, "w") as bundle:
            bundle.writestr("bundle/bin/mafft/mafft.bat", "@echo off\n")
            bundle.writestr("bundle/bin/iqtree3/iqtree3.exe", "x")
        self.widget.plugin_path = root
        self.widget.install_bundle(str(archive))
        self.assertTrue((root / "YR_MPE/bin/mafft/mafft.bat").is_file())
        self.assertTrue((root / "YR_MPE/bin/iqtree3/iqtree3.exe").is_file())

    def test_manifest_url_uses_version(self):
        url = self.widget._resolved_url(
            "https://example.org/v{version}/YR-MPE-bins-win64.zip"
        )
        self.assertNotIn("{version}", url)
        self.assertIn(self.widget._version(), url)


if __name__ == "__main__":
    unittest.main(verbosity=2)
