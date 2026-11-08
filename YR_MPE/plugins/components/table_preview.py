"""Read and display a TSV result table in a simple Qt dialog."""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)


def read_tsv(path: str) -> tuple[List[str], List[List[str]]]:
    """Read a tab-separated table, skipping comment lines."""
    columns: List[str] = []
    rows: List[List[str]] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.rstrip("\n\r")
            if not text or text.startswith("#"):
                continue
            parts = text.split("\t")
            if not columns:
                columns = parts
            else:
                rows.append(parts)
    if not columns:
        raise ValueError("The table is empty.")
    return columns, rows


class TablePreviewDialog(QDialog):
    """A read-only table viewer used by the workspace and Inspector."""

    def __init__(self, path: str, title: str = "Result table", parent=None):
        super().__init__(parent)
        self.path = path
        self.setWindowTitle(title)
        self.setAttribute(Qt.WA_DeleteOnClose, True)
        self.resize(900, 560)
        layout = QVBoxLayout(self)

        heading = QLabel(Path(path).name)
        heading.setStyleSheet("font-weight: 600;")
        layout.addWidget(heading)

        try:
            columns, rows = read_tsv(path)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Table unavailable", str(exc))
            columns, rows = ["Error"], [[str(exc)]]

        self.table = QTableWidget(len(rows), len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSortingEnabled(False)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        for row_index, row in enumerate(rows):
            for column_index in range(len(columns)):
                value = row[column_index] if column_index < len(row) else ""
                self.table.setItem(row_index, column_index, QTableWidgetItem(value))
        layout.addWidget(self.table, 1)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        buttons.addWidget(close)
        layout.addLayout(buttons)


def show_table(parent, path: str, title: str = "Result table") -> TablePreviewDialog:
    dialog = TablePreviewDialog(path, title=title, parent=parent)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    return dialog
