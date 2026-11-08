"""Tests for the two-chain harmonic-mean comparison panel."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

import numpy as np  # noqa: E402

from YR_MPE.plugins.components.harmonic_compare import (  # noqa: E402
    HarmonicMeanPanel,
)

APP = QApplication.instance() or QApplication(sys.argv)


class HarmonicCompareTests(unittest.TestCase):
    def test_clearly_different_chains_get_stars(self):
        rng = np.random.default_rng(5)
        panel = HarmonicMeanPanel()
        panel.set_data(
            ["chain-1", "chain-2"],
            [rng.normal(-1200, 1.0, 400), rng.normal(-1210, 1.0, 400)],
            "lnL",
        )
        self.assertLess(panel.p_value, 0.001)
        self.assertEqual(panel.significance, "***")
        self.assertTrue(any(text.get_text() == "***" for text in panel.axes.texts))

    def test_similar_chains_are_not_significant(self):
        rng = np.random.default_rng(6)
        base = rng.normal(-1200, 4.0, 800)
        panel = HarmonicMeanPanel()
        panel.set_data(
            ["chain-1", "chain-2"],
            [base, base + rng.normal(0, 0.2, 800)],
            "lnL",
        )
        self.assertGreater(panel.p_value, 0.05)
        self.assertEqual(panel.significance, "ns")


if __name__ == "__main__":
    unittest.main(verbosity=2)
