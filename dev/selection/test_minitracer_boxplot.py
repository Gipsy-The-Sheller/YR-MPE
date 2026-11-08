"""Tests for the MiniTracer multi-parameter box plot."""

from __future__ import annotations

import os
import sys
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from YR_MPE.plugins.components.boxplot import BoxPlot  # noqa: E402
from YR_MPE.plugins.components.skyline import SkylinePlot  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class BoxPlotTests(unittest.TestCase):
    def test_one_box_per_parameter_and_zero_line(self):
        plot = BoxPlot()
        plot.set_data(
            ["rate_1", "rate_2", "rate_3"],
            [np.array([0.1, 0.2, 0.3]), np.array([0.0, 0.1, 0.2]),
             np.array([-0.1, -0.2, -0.3])],
            burnin_fraction=0.0,
        )
        self.assertEqual(len(plot.axes.patches), 3)
        zero_lines = [
            line for line in plot.axes.lines
            if len(line.get_ydata()) >= 2
            and all(value == 0 for value in line.get_ydata())
        ]
        self.assertTrue(zero_lines)

    def test_burnin_is_applied(self):
        plot = BoxPlot()
        values = np.arange(10, dtype=float)
        plot.set_data(["x"], [values], burnin_fraction=0.5)
        self.assertEqual(len(plot.axes.patches), 1)


class SkylinePlotTests(unittest.TestCase):
    def test_sample_lines_and_hpd(self):
        plot = SkylinePlot()
        rng = np.random.default_rng(1)
        plot.set_data(
            ["a", "b", "c"],
            [rng.normal(5, 1, 200), rng.normal(4, 1, 200),
             rng.normal(3, 1, 200)],
            burnin_fraction=0.0,
        )
        thin = [line for line in plot.axes.lines
                if line.get_alpha() == 0.3]
        self.assertTrue(thin)
        dashed = [line for line in plot.axes.lines
                  if line.get_linestyle() == "--"]
        self.assertGreaterEqual(len(dashed), 2)
        self.assertEqual(len(plot.axes.get_xticklabels()), 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
