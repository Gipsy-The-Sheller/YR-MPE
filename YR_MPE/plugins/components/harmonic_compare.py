"""Compare two chains by log harmonic mean with an approximate test."""

from __future__ import annotations

import math

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QLabel, QSizePolicy, QVBoxLayout, QWidget
from scipy import stats

from .methods.mcmc_utils import calculate_ESS, calculate_harmonic_mean


class HarmonicMeanPanel(QWidget):
    """Interval plot for two log harmonic means and their significance bar."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.labels = []
        self.values = []
        self.stderrs = []
        self.delta = None
        self.p_value = None
        self.significance = ""
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.summary = QLabel("Load two chains and select a common parameter.")
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.figure = Figure(figsize=(3.2, 2.4))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setMinimumHeight(170)
        self.canvas.setMaximumHeight(260)
        layout.addWidget(self.canvas)
        self.axes = self.figure.add_subplot(111)

    def clear(self):
        self.labels = []
        self.values = []
        self.stderrs = []
        self.delta = None
        self.p_value = None
        self.significance = ""
        self.summary.setText("Load two chains and select a common parameter.")
        self.axes.clear()
        self.canvas.draw_idle()

    @staticmethod
    def _stars(p_value):
        if p_value is None:
            return "ns"
        if p_value < 0.001:
            return "***"
        if p_value < 0.01:
            return "**"
        if p_value < 0.05:
            return "*"
        return "ns"

    def set_data(self, labels, samples_list, parameter_name):
        if len(samples_list) < 2:
            self.clear()
            return
        summaries = []
        for samples in samples_list:
            values = np.asarray(samples, dtype=float)
            values = values[np.isfinite(values)]
            if len(values) < 2:
                self.clear()
                return
            effective_n = calculate_ESS(values)
            harmonic, stderr = calculate_harmonic_mean(
                values, effective_n=effective_n
            )
            summaries.append((harmonic, stderr))
        if any(harmonic is None for harmonic, _ in summaries):
            self.clear()
            return

        self.labels = list(labels)
        self.values = [item[0] for item in summaries]
        self.stderrs = [item[1] for item in summaries]
        self.delta = self.values[0] - self.values[1]
        se_delta = math.hypot(self.stderrs[0], self.stderrs[1])
        if se_delta > 0:
            z_value = abs(self.delta) / se_delta
            self.p_value = 2.0 * stats.norm.sf(z_value)
        else:
            self.p_value = 1.0
        self.significance = self._stars(self.p_value)
        p_text = ("p < 1e-4" if self.p_value < 1e-4
                  else "p = %.4g" % self.p_value)
        self.summary.setText(
            "%s · Δ log HM = %.4f · %s (%s) · approximate delta-method test"
            % (parameter_name, self.delta, p_text, self.significance)
        )
        self.redraw()

    def redraw(self):
        self.axes.clear()
        if len(self.values) < 2:
            self.canvas.draw_idle()
            return

        x = np.arange(len(self.values))
        errors = [1.96 * value for value in self.stderrs]
        self.axes.errorbar(
            x, self.values, yerr=errors, fmt="o", color="black",
            capsize=4, markersize=4, linewidth=1.0,
        )
        self.axes.set_xticks(x)
        self.axes.set_xticklabels(self.labels)
        self.axes.set_ylabel("log harmonic mean")

        low = min(self.values[i] - errors[i] for i in range(len(self.values)))
        high = max(self.values[i] + errors[i] for i in range(len(self.values)))
        span = max(high - low, 1e-9)
        bar_y = high + 0.18 * span
        tick = 0.04 * span
        self.axes.plot([x[0], x[0], x[1], x[1]],
                       [bar_y - tick, bar_y, bar_y, bar_y - tick],
                       color="black", linewidth=1.0)
        self.axes.text(
            (x[0] + x[1]) / 2, bar_y + 0.03 * span, self.significance,
            ha="center", va="bottom", fontsize=11,
        )
        self.axes.set_ylim(low - 0.15 * span, bar_y + 0.25 * span)
        self.figure.tight_layout()
        self.canvas.draw_idle()
