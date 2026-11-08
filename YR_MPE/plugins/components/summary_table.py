#!/usr/bin/env python3
"""
Compact summary-statistics form for a single MCMC chain.

A four-column grid (label/value/label/value) replaces QTableWidget: it uses
less vertical space and has no header or scroll frame.
"""

import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QGridLayout, QLabel, QVBoxLayout, QWidget

from .methods.mcmc_utils import (
    calculate_ESS, calculate_95HPD, calculate_harmonic_mean,
)


def format_number(value):
    """Four decimals, or two-digit scientific notation for large values."""
    if value is None:
        return "N/A"
    if abs(value) > 100000:
        return f"{value:.2e}"
    return f"{value:.4f}"


class SummaryTable(QWidget):
    """Compact statistics form for one MCMC parameter."""

    METRICS = [
        "Number of samples",
        "Mean",
        "Harmonic mean (log)",
        "Harmonic mean SE",
        "Range",
        "95% HPD interval",
        "ESS",
    ]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.values = {}
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)
        self.setLayout(layout)

        grid = QGridLayout()
        grid.setContentsMargins(2, 2, 2, 2)
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(2)
        layout.addLayout(grid)

        font = QFont()
        font.setPointSize(9)
        self.setFont(font)

        self.value_labels = {}
        for index, metric in enumerate(self.METRICS):
            row, pair = divmod(index, 2)
            label = QLabel(metric)
            label.setStyleSheet("color: palette(mid);")
            label.setFont(font)
            value = QLabel("")
            value.setFont(font)
            value.setTextInteractionFlags(Qt.TextSelectableByMouse)
            if metric.startswith("Harmonic"):
                label.setToolTip(
                    "Computed as log(HM(exp(x))): likelihood traces are "
                    "reported in log space."
                )
            grid.addWidget(label, row, pair * 2)
            grid.addWidget(value, row, pair * 2 + 1)
            self.value_labels[metric] = value
            self.values[metric] = ""

        layout.addStretch(1)

    def set_data(self, samples):
        """Set one parameter's samples and refresh all values."""
        if samples is None or len(samples) == 0:
            for metric in self.METRICS:
                self.values[metric] = ""
                self.value_labels[metric].setText("")
            return

        samples = np.asarray(samples, dtype=float)
        n_samples = len(samples)
        mean_val = np.mean(samples)
        min_val = np.min(samples)
        max_val = np.max(samples)
        hpd_interval = calculate_95HPD(samples)
        ess = calculate_ESS(samples)
        harmonic, harmonic_se = calculate_harmonic_mean(samples)

        values = {
            "Number of samples": str(n_samples),
            "Mean": format_number(mean_val),
            "Harmonic mean (log)": format_number(harmonic),
            "Harmonic mean SE": format_number(harmonic_se),
            "Range": "[%s, %s]" % (format_number(min_val), format_number(max_val)),
            "95% HPD interval": "[%s, %s]" % (
                format_number(hpd_interval[0]), format_number(hpd_interval[1])
            ),
            "ESS": format_number(ess),
        }
        for metric in self.METRICS:
            self.values[metric] = values[metric]
            self.value_labels[metric].setText(values[metric])
