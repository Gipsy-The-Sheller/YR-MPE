"""Box-plot view for several selected MCMC parameters (Tracer-style)."""

from __future__ import annotations

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QSizePolicy, QVBoxLayout, QWidget


class BoxPlot(QWidget):
    """One box per selected parameter, after burn-in."""

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
        if not self.samples:
            self.canvas.draw_idle()
            return

        plotted = []
        tick_labels = []
        for label, values in zip(self.labels, self.samples):
            start = int(len(values) * self.burnin_fraction)
            values = values[start:]
            if len(values):
                plotted.append(values)
                tick_labels.append(label)
        if not plotted:
            self.canvas.draw_idle()
            return

        box = self.axes.boxplot(
            plotted, labels=tick_labels, widths=0.5,
            showfliers=True, patch_artist=True,
            boxprops={"facecolor": "white", "edgecolor": "black"},
            whiskerprops={"color": "black"},
            capprops={"color": "black"},
            flierprops={"marker": "o", "markerfacecolor": "black",
                        "markeredgecolor": "black", "markersize": 3},
        )
        for median in box["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)

        values = np.concatenate(plotted)
        if values.min() <= 0.0 <= values.max():
            self.axes.axhline(0.0, color="black", linewidth=0.8)

        self.axes.set_ylabel("Value")
        self.axes.set_xlabel("Parameter")
        if any(len(label) > 10 for label in tick_labels):
            self.axes.tick_params(axis="x", labelrotation=30)
            for tick in self.axes.get_xticklabels():
                tick.set_horizontalalignment("right")
        self.figure.tight_layout()
        self.canvas.draw_idle()
