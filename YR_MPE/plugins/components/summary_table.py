#!/usr/bin/env python3
"""
Summary Table component for displaying MCMC chain statistics.
"""

import numpy as np
from PyQt5.QtWidgets import QWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QSizePolicy
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from .methods.mcmc_utils import calculate_ESS, calculate_95HPD


class SummaryTable(QWidget):
    """
    QWidget for displaying summary statistics of a single MCMC chain.
    Shows number of samples, mean, range, 95% HPD interval, and ESS.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        """Initialize the user interface with compact table layout."""
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)
        
        # Create table widget for statistics
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(['Statistic', 'Value'])
        self.table.setRowCount(5)
        # Disable vertical header (index) to avoid display issues
        self.table.verticalHeader().setVisible(False)
        
        # Set the statistic labels in the first column
        metrics = [
            'Number of samples',
            'Mean',
            'Range',
            '95% HPD interval',
            'ESS'
        ]
        for i, metric in enumerate(metrics):
            item = QTableWidgetItem(metric)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)  # Make non-editable
            self.table.setItem(i, 0, item)
        
        # Apply compact styling similar to TraceComp
        font = QFont()
        font.setPointSize(9)  # Smaller font size
        self.table.setFont(font)
        self.table.horizontalHeader().setFont(font)
        self.table.verticalHeader().setFont(font)
        
        # Adjust row heights and column widths
        for i in range(5):
            self.table.setRowHeight(i, 25)  # Reduced row height
        self.table.setColumnWidth(0, 120)  # Adjust column widths
        self.table.setColumnWidth(1, 120)
        
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)  # Make table read-only
        
        # Set padding to 0px as requested
        self.table.setStyleSheet("QTableWidget { padding: 0px; } "
                                "QTableWidgetItem { padding: 0px; }")
        
        layout.addWidget(self.table)
        
    def set_data(self, samples):
        """
        Set the data and calculate summary statistics.
        
        Parameters:
        samples: array-like, MCMC samples
        """
        if samples is None or len(samples) == 0:
            # Clear the table
            for i in range(5):
                self.table.setItem(i, 1, QTableWidgetItem(""))
            return
            
        samples = np.asarray(samples, dtype=float)
        n_samples = len(samples)
        mean_val = np.mean(samples)
        min_val = np.min(samples)
        max_val = np.max(samples)
        hpd_interval = calculate_95HPD(samples)
        ess = calculate_ESS(samples)
        
        # Format values for display
        values = [
            str(n_samples),
            f"{mean_val:.6f}",
            f"[{min_val:.6f}, {max_val:.6f}]",
            f"[{hpd_interval[0]:.6f}, {hpd_interval[1]:.6f}]",
            f"{ess:.2f}"
        ]
        
        for i, value in enumerate(values):
            item = QTableWidgetItem(value)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(i, 1, item)