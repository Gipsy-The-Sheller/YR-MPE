"""Tests for the toolbar feature search."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import QEvent, Qt
from PyQt5.QtGui import QColor, QIcon, QKeyEvent, QPixmap
from PyQt5.QtWidgets import (
    QAction, QApplication, QMenu, QToolBar, QToolButton, QVBoxLayout, QWidget,
)

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.feature_search import (  # noqa: E402
    ACTION_INDEX_ROLE, MENU_ICON_ROLE, MENU_NAME_ROLE, SEARCH_ROLE,
    FeatureSearch,
)
from YR_MPE.platforms.methods.hover_toolbar import (  # noqa: E402
    install_hover_expand,
)

APP = QApplication.instance() or QApplication(sys.argv)


def icon(color: str) -> QIcon:
    pixmap = QPixmap(16, 16)
    pixmap.fill(QColor(color))
    return QIcon(pixmap)


class FeatureSearchTests(unittest.TestCase):
    def setUp(self):
        self.host = QWidget()
        layout = QVBoxLayout(self.host)
        self.toolbar = QToolBar()
        layout.addWidget(self.toolbar)

        self.align_button = QToolButton()
        self.align_button.setText("ALIGN")
        self.align_button.setIcon(icon("red"))
        self.align_button.setPopupMode(QToolButton.InstantPopup)
        self.toolbar.addWidget(self.align_button)
        self.expander = install_hover_expand(self.align_button)

        self.menu = QMenu(self.align_button)
        self.align_button.setMenu(self.menu)

        self.mafft = QAction(icon("blue"), "MAFFT", self.menu)
        self.mafft.setToolTip("align nucleotide sequences")
        self.menu.addAction(self.mafft)

        trim_menu = QMenu("Trim by...", self.menu)
        self.trimal = QAction("TrimAl", trim_menu)
        self.trimal.setToolTip("trim gaps")
        trim_menu.addAction(self.trimal)
        self.menu.addMenu(trim_menu)

        disabled = QAction("Planned feature", self.menu)
        disabled.setEnabled(False)
        self.menu.addAction(disabled)

        history = QAction("C:/Users/example/workspace", self.menu)
        history.setData("C:/Users/example/workspace")
        self.menu.addAction(history)

        self.triggered = []
        self.mafft.triggered.connect(lambda: self.triggered.append("mafft"))
        self.trimal.triggered.connect(lambda: self.triggered.append("trimal"))

        self.search = FeatureSearch(self.host)
        self.toolbar.addWidget(self.search)
        self.search.bind_buttons([self.align_button])

    def tearDown(self):
        self.host.close()

    def test_entries_skip_disabled_and_paths(self):
        labels = [entry.label for entry in self.search.entries()]
        self.assertEqual(labels, ["MAFFT", "TrimAl"])
        self.assertEqual(self.search.entries()[0].menu_label, "ALIGN")
        self.assertEqual(self.search.entries()[1].path, "ALIGN → Trim by... → TrimAl")

    def test_model_roles_and_search_text(self):
        index = self.search.model.index(0, 0)
        self.assertEqual(index.data(Qt.DisplayRole), "MAFFT")
        self.assertEqual(index.data(MENU_NAME_ROLE), "ALIGN")
        self.assertIsInstance(index.data(MENU_ICON_ROLE), QIcon)
        self.assertIn("nucleotide", index.data(SEARCH_ROLE))
        self.assertEqual(index.data(ACTION_INDEX_ROLE), 0)

    def test_filter_and_activate(self):
        self.search.line_edit.setText("mafft")
        self.search._on_text_edited("mafft")
        self.assertEqual(self.search.completer.completionCount(), 1)
        index = self.search.completer.completionModel().index(0, 0)
        self.search._activate_index(index)
        self.assertEqual(self.triggered, ["mafft"])
        self.assertEqual(self.search.line_edit.text(), "")

    def test_enter_triggers_first_match_when_popup_hidden(self):
        self.search.line_edit.setText("trimal")
        self.search._on_text_edited("trimal")
        self.search.completer.popup().hide()
        self.search._on_return_pressed()
        self.assertEqual(self.triggered, ["trimal"])

    def test_escape_clears_when_popup_hidden(self):
        self.search.line_edit.setText("mafft")
        self.search.completer.popup().hide()
        event = QKeyEvent(QEvent.KeyPress, Qt.Key_Escape, Qt.NoModifier)
        handled = self.search.eventFilter(self.search.line_edit, event)
        self.assertTrue(handled)
        self.assertEqual(self.search.line_edit.text(), "")

    def test_search_field_expands_in_toolbar(self):
        self.host.resize(1000, 100)
        self.host.show()
        APP.processEvents()
        self.assertGreater(self.search.width(), 300)
        self.assertGreater(self.search.line_edit.width(), 250)

    def test_search_hides_while_button_or_menu_is_expanded(self):
        from PyQt5.QtTest import QTest

        self.host.resize(1000, 100)
        self.host.show()
        APP.processEvents()
        self.assertGreater(self.search.width(), 0)

        self.expander.expansionChanged.emit(True)
        APP.processEvents()
        self.assertEqual(self.search.width(), 0)

        self.expander.expansionChanged.emit(False)
        QTest.qWait(120)
        self.assertGreater(self.search.width(), 0)

        self.menu.aboutToShow.emit()
        APP.processEvents()
        self.assertEqual(self.search.width(), 0)

        self.menu.aboutToHide.emit()
        QTest.qWait(120)
        self.assertGreater(self.search.width(), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
