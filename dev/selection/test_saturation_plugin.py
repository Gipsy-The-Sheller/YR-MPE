"""Tests for the substitution-saturation worker thread."""

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

from YR_MPE.plugins.substitution_saturation_plugin import (  # noqa: E402
    SaturationThread, SubstitutionSaturationPlugin,
)

APP = QApplication.instance() or QApplication(sys.argv)

ALIGNMENT = (
    ">A\nATGCAAATGCAAATGCAAATGCAA\n"
    ">B\nATGCAAATGCAAATGCAAATGTAA\n"
    ">C\nATGCAAATGCAAATGCAAATGCAT\n"
    ">D\nATGCAAATGCAAATGCAAATGCTT\n"
)


class SaturationTests(unittest.TestCase):
    def test_metadata(self):
        plugin = SubstitutionSaturationPlugin(workdir=tempfile.mkdtemp())
        try:
            self.assertTrue(plugin.run_button.text())
            self.assertIn("kuck2014bacoca", plugin.citation_keys)
        finally:
            plugin.close()

    def test_thread_returns_metrics(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(ALIGNMENT, encoding="utf-8")
        payloads = []
        errors = []
        thread = SaturationThread(
            [str(source)], {"gap_treatment": "pairwise"}, workdir=str(root)
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        result = payloads[0][0]
        self.assertEqual(result["filename"], "input.fas")
        self.assertIn("c_value", result)
        self.assertIn("p_values", result)


if __name__ == "__main__":
    unittest.main(verbosity=2)
