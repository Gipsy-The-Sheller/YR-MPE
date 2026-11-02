"""Tests for OrthoFinder drag-and-drop input handling."""

from __future__ import annotations

import os
import shutil
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

from YR_MPE.plugins.orthofinder_plugin import OrthoFinderPlugin  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
GB_DIR = REPO / "dev" / "gbff_to_orthofinder" / "testdata" / "gb"


class TestPlugin(OrthoFinderPlugin):
    def config(self):
        return True


class DragDropTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = TestPlugin(import_from="YR_MPEA", workdir=str(self.root))

    def tearDown(self):
        self.plugin.is_running = False
        self.plugin.close()
        shutil.rmtree(self.root, ignore_errors=True)
        for path in getattr(self.plugin, "_drop_temp_dirs", []):
            shutil.rmtree(path, ignore_errors=True)

    def test_genbank_folder(self):
        ok = self.plugin._apply_dropped_paths([GB_DIR], show_dialogs=False)
        self.assertTrue(ok)
        self.assertEqual(self.plugin.source_combo.currentData(), "gbff")
        self.assertEqual(Path(self.plugin.input_edit.text()), GB_DIR)
        self.assertTrue(self.plugin.scan_label.text())

    def test_multirecord_file(self):
        file_path = GB_DIR / "HM347959.gb"
        ok = self.plugin._apply_dropped_paths([file_path], show_dialogs=False)
        self.assertTrue(ok)
        self.assertEqual(self.plugin.source_combo.currentData(), "gbff_records")
        self.assertEqual(Path(self.plugin.input_edit.text()), file_path)

    def test_multiple_genbank_files(self):
        files = [GB_DIR / "HM347959.gb", GB_DIR / "JN637766.gb"]
        ok = self.plugin._apply_dropped_paths(files, show_dialogs=False)
        self.assertTrue(ok)
        self.assertEqual(self.plugin.source_combo.currentData(), "gbff")
        dropped_dir = Path(self.plugin.input_edit.text())
        self.assertTrue(dropped_dir.is_dir())
        self.assertEqual({path.name for path in dropped_dir.iterdir()}, {path.name for path in files})

    def test_proteome_folder(self):
        folder = self.root / "proteomes"
        folder.mkdir()
        (folder / "sp1.faa").write_text(">a\nMQ\n", encoding="utf-8")
        ok = self.plugin._apply_dropped_paths([folder], show_dialogs=False)
        self.assertTrue(ok)
        self.assertEqual(self.plugin.source_combo.currentData(), "proteome")

    def test_single_proteome_file_becomes_temp_folder(self):
        file_path = self.root / "sp1.faa"
        file_path.write_text(">a\nMQ\n", encoding="utf-8")
        ok = self.plugin._apply_dropped_paths([file_path], show_dialogs=False)
        self.assertTrue(ok)
        self.assertEqual(self.plugin.source_combo.currentData(), "proteome")
        self.assertTrue(Path(self.plugin.input_edit.text()).is_dir())

    def test_empty_folder_rejected(self):
        folder = self.root / "empty"
        folder.mkdir()
        ok = self.plugin._apply_dropped_paths([folder], show_dialogs=False)
        self.assertFalse(ok)

    def test_mixed_files_rejected(self):
        gb = GB_DIR / "HM347959.gb"
        fa = self.root / "sp1.faa"
        fa.write_text(">a\nMQ\n", encoding="utf-8")
        ok = self.plugin._apply_dropped_paths([gb, fa], show_dialogs=False)
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main(verbosity=2)
