"""Tests for MiniTracer's log-space harmonic-mean statistic."""

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

from YR_MPE.plugins.components.methods.mcmc_utils import (  # noqa: E402
    calculate_harmonic_mean, calculate_psrf,
)
from YR_MPE.plugins.components.summary_table import (  # noqa: E402
    SummaryTable, format_number,
)
from YR_MPE.plugins.components.tracecomp import TraceComp  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


def expected_log_hm(values):
    scaled = -np.asarray(values, dtype=float)
    maximum = np.max(scaled)
    return -(maximum + np.log(np.mean(np.exp(scaled - maximum))))


class HarmonicMeanTests(unittest.TestCase):
    def test_is_exp_harmonic_log(self):
        values = np.array([1.0, 2.0, 4.0])
        harmonic, stderr = calculate_harmonic_mean(values)
        self.assertAlmostEqual(harmonic, expected_log_hm(values), places=9)
        self.assertGreater(stderr, 0.0)

    def test_negative_log_likelihood_values(self):
        values = np.array([-1200.0, -1201.0, -1199.0])
        harmonic, stderr = calculate_harmonic_mean(values)
        self.assertAlmostEqual(harmonic, expected_log_hm(values), places=9)
        self.assertGreater(stderr, 0.0)

    def test_empty_is_undefined(self):
        self.assertEqual(calculate_harmonic_mean([]), (None, None))
        self.assertEqual(calculate_harmonic_mean([np.nan]), (None, None))

    def test_number_format(self):
        self.assertEqual(format_number(1.23456), "1.2346")
        self.assertEqual(format_number(123456.789), "1.23e+05")
        self.assertEqual(format_number(None), "N/A")

    def test_summary_form_values(self):
        table = SummaryTable()
        table.set_data(np.array([1.0, 2.0, 4.0]))
        self.assertEqual(len(table.values), 7)
        self.assertEqual(table.values["Mean"], "2.3333")
        self.assertNotEqual(table.values["Harmonic mean (log)"], "N/A")
        self.assertNotEqual(table.values["Harmonic mean SE"], "N/A")

        table.set_data(np.array([-1200.0, -1201.0, -1199.0]))
        self.assertNotEqual(table.values["Harmonic mean (log)"], "N/A")


class PSRFTests(unittest.TestCase):
    def test_identical_chains_are_near_one(self):
        rng = np.random.default_rng(7)
        first = rng.normal(size=500)
        # Identical chains have no between-chain variance; the finite-sample
        # correction leaves R-hat just below one.
        psrf = calculate_psrf(first, first.copy())
        self.assertIsNotNone(psrf)
        self.assertAlmostEqual(psrf, 1.0, delta=0.005)

    def test_separated_chains_are_larger(self):
        rng = np.random.default_rng(11)
        first = rng.normal(0.0, 1.0, 500)
        second = rng.normal(2.0, 1.0, 500)
        psrf = calculate_psrf(first, second)
        self.assertIsNotNone(psrf)
        self.assertGreater(psrf, 1.5)

    def test_constant_equal_and_different_chains(self):
        self.assertEqual(calculate_psrf([1.0] * 50, [1.0] * 50), 1.0)
        self.assertTrue(np.isinf(calculate_psrf([1.0] * 50, [2.0] * 50)))

    def test_short_and_unequal_chains(self):
        self.assertIsNone(calculate_psrf([1.0], [1.0, 2.0]))
        psrf = calculate_psrf(np.arange(100.0), np.arange(90.0))
        self.assertIsNotNone(psrf)
        self.assertAlmostEqual(psrf, 1.0, delta=0.02)

    def test_tracecomp_table_shows_psrf(self):
        widget = TraceComp()
        first = np.linspace(0.0, 1.0, 50)
        second = np.linspace(0.0, 1.0, 50)
        widget.set_data(first, np.arange(50), second, np.arange(50))
        self.assertEqual(widget.table.rowCount(), 6)
        self.assertEqual(widget.table.item(5, 0).text(), "PSRF (R-hat)")
        self.assertNotEqual(widget.table.item(5, 1).text(), "")
        self.assertIn("PSRF", widget.ax.get_title())
        widget.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
