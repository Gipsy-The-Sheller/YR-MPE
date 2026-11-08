"""Posterior skyline view: one thin line per MCMC sample across parameters."""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QSizePolicy, QVBoxLayout, QWidget

from .methods.mcmc_utils import calculate_95HPD


class SkylinePlot(QWidget):
    """Connect the selected parameters across each posterior sample."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.labels = []
        self.samples = []
        self.burnin_fraction = 0.1
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        self.axes = self.figure.add_subplot(111)

    def set_data(self, labels, samples, burnin_fraction=0.1):
        self.labels = list(labels)
        self.samples = [np.asarray(values, dtype=float) for values in samples]
        self.burnin_fraction = max(0.0, min(0.999, float(burnin_fraction)))
        self.redraw()

    def clear(self):
        self.labels = []
        self.samples = []
        self.axes.clear()
        self.canvas.draw_idle()

    def redraw(self):
        self.axes.clear()
        if len(self.samples) < 2:
            self.canvas.draw_idle()
            return

        arrays = []
        for values in self.samples:
            start = int(len(values) * self.burnin_fraction)
            arrays.append(values[start:])
        length = min(len(values) for values in arrays)
        if length < 1:
            self.canvas.draw_idle()
            return
        arrays = [values[:length] for values in arrays]

        x = np.arange(len(self.labels))
        # Keep the plot responsive for very long chains while retaining the
        # shape of the posterior cloud.
        stride = max(1, length // 1200)
        for index in range(0, length, stride):
            y = [values[index] for values in arrays]
            if not np.all(np.isfinite(y)):
                continue
            self.axes.plot(
                x, y, color="0.25", alpha=0.3, linewidth=0.4, zorder=1
            )

        lower = []
        upper = []
        medians = []
        for values in arrays:
            finite = values[np.isfinite(values)]
            if len(finite) < 2:
                lower.append(np.nan)
                upper.append(np.nan)
                medians.append(np.nan)
                continue
            hpd = calculate_95HPD(finite)
            lower.append(hpd[0])
            upper.append(hpd[1])
            medians.append(float(np.median(finite)))

        self.axes.plot(x, lower, color="black", linestyle="--", linewidth=1.1,
                       label="95% HPD", zorder=3)
        self.axes.plot(x, upper, color="black", linestyle="--", linewidth=1.1,
                       zorder=3)
        self.axes.plot(x, medians, color="black", linewidth=1.6,
                       label="Median", zorder=4)
        self.axes.fill_between(x, lower, upper, color="0.5", alpha=0.08,
                               zorder=2)

        all_values = np.concatenate([v for v in arrays if len(v)])
        finite_all = all_values[np.isfinite(all_values)]
        if len(finite_all) and finite_all.min() <= 0.0 <= finite_all.max():
            self.axes.axhline(0.0, color="black", linewidth=0.8)

        self.axes.set_xticks(x)
        self.axes.set_xticklabels(self.labels)
        if any(len(label) > 10 for label in self.labels):
            self.axes.tick_params(axis="x", labelrotation=30)
            for tick in self.axes.get_xticklabels():
                tick.set_horizontalalignment("right")
        self.axes.set_ylabel("Value")
        self.axes.set_xlabel("Parameter")
        self.axes.legend(loc="best", frameon=False, fontsize=8)
        self.figure.tight_layout()
        self.canvas.draw_idle()
