"""Toolbar-wide feature search with an icon-rich completion popup."""

from __future__ import annotations

import re
import time
from dataclasses import dataclass
from typing import List, Optional

from PyQt5.QtCore import QEvent, QModelIndex, QRect, QSize, Qt, QTimer
from PyQt5.QtGui import (
    QColor, QFontMetrics, QIcon, QPainter, QPen, QPixmap,
    QStandardItem, QStandardItemModel,
)
from PyQt5.QtWidgets import (
    QAction, QApplication, QCompleter, QHBoxLayout, QLineEdit, QListView,
    QSizePolicy, QStyle, QStyledItemDelegate, QStyleOptionViewItem,
    QToolBar, QToolButton, QWidget,
)

from .hover_toolbar import HoverExpander

SEARCH_ROLE = Qt.UserRole + 1
MENU_NAME_ROLE = Qt.UserRole + 2
MENU_ICON_ROLE = Qt.UserRole + 3
ACTION_INDEX_ROLE = Qt.UserRole + 4

_FEATURE_COLOR = QColor("#1f2328")
_MENU_COLOR = QColor("#8a8f98")


def _clean_label(text: str) -> str:
    value = re.sub(r"(?<!&)&(?!&)", "", str(text or ""))
    value = value.replace("&&", "&").strip()
    return re.sub(r"\s+", " ", value)


def _looks_like_path_action(action: QAction) -> bool:
    """Skip workspace-history entries, which are paths rather than features."""
    data = action.data()
    if not isinstance(data, str) or len(data) < 4:
        return False
    return bool(re.match(r"^(?:[A-Za-z]:[\\/]|/|\\\\|~[\\/])", data.strip()))


def _search_icon(color: str = "#6b7280") -> QIcon:
    pixmap = QPixmap(16, 16)
    pixmap.fill(Qt.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing, True)
    painter.setPen(QPen(QColor(color), 1.8))
    painter.drawEllipse(2, 2, 8, 8)
    painter.drawLine(9, 9, 14, 14)
    painter.end()
    return QIcon(pixmap)


@dataclass
class FeatureEntry:
    label: str
    menu_label: str
    menu_icon: QIcon
    feature_icon: QIcon
    path: str
    description: str
    action: QAction


class FeatureCompleter(QCompleter):
    """Filter on the composite search role but insert only the feature name."""

    def splitPath(self, path):
        return [str(path)]

    def pathFromIndex(self, index: QModelIndex) -> str:
        return str(index.data(Qt.DisplayRole) or "")


class FeatureItemDelegate(QStyledItemDelegate):
    """Paint feature icon + name, then a quieter 'in <menu>' suffix."""

    def sizeHint(self, option, index):
        return QSize(360, 30)

    def paint(self, painter, option, index):
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        opt.text = ""
        style = opt.widget.style() if opt.widget is not None else QApplication.style()
        style.drawPrimitive(QStyle.PE_PanelItemViewItem, opt, painter, opt.widget)

        painter.save()
        selected = bool(opt.state & QStyle.State_Selected)
        rect = opt.rect.adjusted(2, 1, -2, -1)
        feature_color = (
            opt.palette.highlightedText().color()
            if selected else opt.palette.text().color()
        )
        menu_color = QColor(feature_color if selected else _MENU_COLOR)
        if selected:
            menu_color.setAlpha(190)

        x = rect.left() + 6
        feature_icon = index.data(Qt.DecorationRole)
        if isinstance(feature_icon, QIcon) and not feature_icon.isNull():
            feature_icon.paint(
                painter, QRect(x, rect.center().y() - 9, 18, 18)
            )
        x += 26

        metrics = QFontMetrics(opt.font)
        name = str(index.data(Qt.DisplayRole) or "")
        menu_name = str(index.data(MENU_NAME_ROLE) or "")
        menu_icon = index.data(MENU_ICON_ROLE)
        in_text = "in"
        gap = 6
        menu_icon_size = 14
        suffix_width = (
            metrics.horizontalAdvance(in_text) + gap
            + menu_icon_size + 4 + metrics.horizontalAdvance(menu_name)
        )
        suffix_x = max(x + 24, rect.right() - suffix_width - 6)
        name_width = max(24, suffix_x - x - 8)
        name = metrics.elidedText(name, Qt.ElideRight, name_width)

        painter.setPen(feature_color)
        painter.drawText(
            QRect(x, rect.top(), name_width, rect.height()),
            Qt.AlignVCenter | Qt.AlignLeft | Qt.TextSingleLine,
            name,
        )

        x = suffix_x
        painter.setPen(menu_color)
        painter.drawText(
            QRect(x, rect.top(), metrics.horizontalAdvance(in_text), rect.height()),
            Qt.AlignVCenter | Qt.AlignLeft | Qt.TextSingleLine,
            in_text,
        )
        x += metrics.horizontalAdvance(in_text) + gap

        if isinstance(menu_icon, QIcon) and not menu_icon.isNull():
            menu_icon.paint(
                painter, QRect(x, rect.center().y() - menu_icon_size // 2,
                               menu_icon_size, menu_icon_size)
            )
            x += menu_icon_size + 4

        menu_width = max(0, rect.right() - x - 6)
        menu_name = metrics.elidedText(menu_name, Qt.ElideRight, menu_width)
        painter.drawText(
            QRect(x, rect.top(), menu_width, rect.height()),
            Qt.AlignVCenter | Qt.AlignLeft | Qt.TextSingleLine,
            menu_name,
        )
        painter.restore()


class FeaturePopup(QListView):
    """Keep our delegate when QCompleter installs its own on show."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._feature_delegate = None

    def set_feature_delegate(self, delegate):
        self._feature_delegate = delegate
        super().setItemDelegate(delegate)

    def setItemDelegate(self, delegate):
        if self._feature_delegate is not None and not isinstance(
                delegate, FeatureItemDelegate):
            delegate = self._feature_delegate
        super().setItemDelegate(delegate)


class FeatureSearch(QWidget):
    """A search field that fills the free toolbar width.

    The index is built from the live QToolBar menus, so newly added actions are
    picked up when the field receives focus.
    """

    def __init__(self, host: QWidget, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.host = host
        self._entries: List[FeatureEntry] = []
        self._bound_buttons = set()
        self._expanded_buttons = set()
        self._open_menus = 0
        self._popup_watch = False
        self._popup_seen = False
        self._popup_deadline = 0.0
        self._visibility_timer = QTimer(self)
        self._visibility_timer.setSingleShot(True)
        self._visibility_timer.timeout.connect(self._apply_visibility)

        self.setObjectName("featureSearchHost")
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setMinimumWidth(240)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 0, 8, 0)
        layout.setSpacing(0)

        self.line_edit = QLineEdit(self)
        self.line_edit.setObjectName("featureSearch")
        self.line_edit.setPlaceholderText("Search features…")
        self.line_edit.setClearButtonEnabled(True)
        self.line_edit.setMinimumHeight(30)
        self.line_edit.setMaximumHeight(30)
        self.line_edit.setAccessibleName("Search features")
        self.line_edit.addAction(_search_icon(), QLineEdit.LeadingPosition)
        self.line_edit.setStyleSheet(
            "QLineEdit#featureSearch {"
            " border: 1px solid palette(mid);"
            " border-radius: 6px;"
            " padding: 3px 8px;"
            " background: palette(base);"
            "}"
            "QLineEdit#featureSearch:focus {"
            " border: 1px solid palette(highlight);"
            "}"
        )
        layout.addWidget(self.line_edit, 1)

        self.model = QStandardItemModel(self)
        self.completer = FeatureCompleter(self)
        self.completer.setModel(self.model)
        self.completer.setCompletionRole(SEARCH_ROLE)
        self.completer.setFilterMode(Qt.MatchContains)
        self.completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.completer.setCompletionMode(QCompleter.PopupCompletion)
        self.completer.setMaxVisibleItems(9)
        self.completer.setWrapAround(False)
        self.popup_view = FeaturePopup()
        self.completer.setPopup(self.popup_view)
        self.popup_view.set_feature_delegate(FeatureItemDelegate(self.popup_view))
        self.line_edit.setCompleter(self.completer)
        self._style_popup()

        self.line_edit.textEdited.connect(self._on_text_edited)
        self.line_edit.returnPressed.connect(self._on_return_pressed)
        self.line_edit.installEventFilter(self)
        self.completer.activated[QModelIndex].connect(self._activate_index)

        self.refresh()

    # ------------------------------------------------------------------
    # toolbar-button visibility
    # ------------------------------------------------------------------
    def bind_buttons(self, buttons):
        """Hide the search while any of *buttons* is expanded or open."""
        for button in buttons:
            if id(button) in self._bound_buttons:
                continue
            self._bound_buttons.add(id(button))
            expander = button.findChild(HoverExpander)
            if expander is not None:
                expander.expansionChanged.connect(
                    lambda expanded, btn=button:
                        self._set_button_expanded(btn, expanded)
                )
            menu = button.menu()
            if menu is not None:
                menu.aboutToShow.connect(
                    lambda btn=button: self._set_menu_open(btn, True)
                )
                menu.aboutToHide.connect(
                    lambda btn=button: self._set_menu_open(btn, False)
                )
            else:
                # Some QToolButtons build their popup from actions() and do
                # not expose menu(); poll the application popup while open.
                button.pressed.connect(
                    lambda btn=button: self._watch_button_popup(btn)
                )

    def _set_button_expanded(self, button, expanded: bool):
        if expanded:
            self._expanded_buttons.add(button)
        else:
            self._expanded_buttons.discard(button)
        self._schedule_visibility()

    def _set_menu_open(self, button, opened: bool):
        self._open_menus = max(0, self._open_menus + (1 if opened else -1))
        self._schedule_visibility()

    def _watch_button_popup(self, button):
        self._popup_watch = True
        self._popup_seen = False
        self._popup_deadline = time.monotonic() + 0.8
        self._set_collapsed(True)
        QTimer.singleShot(80, self._poll_popup)

    def _poll_popup(self):
        popup = QApplication.activePopupWidget()
        if popup is not None:
            self._popup_seen = True
            self._popup_watch = True
            self._set_collapsed(True)
            QTimer.singleShot(120, self._poll_popup)
            return
        if not self._popup_seen and time.monotonic() < self._popup_deadline:
            # The button's action menu can open after the first poll.
            self._set_collapsed(True)
            QTimer.singleShot(80, self._poll_popup)
            return
        self._popup_watch = False
        self._schedule_visibility()

    def _schedule_visibility(self):
        hide = (bool(self._expanded_buttons) or self._open_menus > 0
                or self._popup_watch)
        if hide:
            self._visibility_timer.stop()
            self.completer.popup().hide()
            self._set_collapsed(True)
        else:
            # Avoid flicker when moving between adjacent menu buttons.
            self._visibility_timer.start(80)

    def _apply_visibility(self):
        if (not self._expanded_buttons and self._open_menus == 0
                and not self._popup_watch):
            self._set_collapsed(False)

    def _set_collapsed(self, collapsed: bool):
        # QToolBar re-shows child widgets during layout; zero width is stable
        # and removes the search from the row while a menu is expanded.
        if collapsed:
            self.setMinimumWidth(0)
            self.setMaximumWidth(0)
            self.line_edit.setVisible(False)
        else:
            self.setMaximumWidth(16777215)
            self.setMinimumWidth(240)
            self.line_edit.setVisible(True)

    # ------------------------------------------------------------------
    # index
    # ------------------------------------------------------------------
    def refresh(self):
        """Rebuild from the current toolbar menus."""
        text = self.line_edit.text()
        self.model.clear()
        self._entries = []
        seen_actions = set()
        seen_buttons = set()

        for toolbar in self.host.findChildren(QToolBar):
            for button in toolbar.findChildren(QToolButton):
                if id(button) in seen_buttons:
                    continue
                if button.parent() is not toolbar:
                    continue
                label = _clean_label(button.text())
                if not label:
                    continue
                seen_buttons.add(id(button))
                menu_icon = button.icon()
                menu = button.menu()
                actions = menu.actions() if menu is not None else button.actions()
                self._walk_actions(
                    actions, label, menu_icon, [label], seen_actions
                )

        if text.strip():
            self.completer.setCompletionPrefix(text.strip())
            self.completer.complete()
            self._style_popup()
        else:
            self.completer.popup().hide()
        return len(self._entries)

    def entries(self) -> List[FeatureEntry]:
        return list(self._entries)

    def _walk_actions(self, actions, menu_label, menu_icon, path,
                      seen_actions):
        for action in actions:
            if action.isSeparator():
                continue
            if id(action) in seen_actions:
                continue
            label = _clean_label(action.text())
            if not label:
                continue
            submenu = action.menu()
            if submenu is not None:
                if not action.isEnabled():
                    continue
                self._walk_actions(
                    submenu.actions(), menu_label, menu_icon,
                    path + [label], seen_actions,
                )
                continue
            if not action.isEnabled():
                continue
            if _looks_like_path_action(action):
                continue
            seen_actions.add(id(action))

            feature_icon = action.icon()
            if feature_icon.isNull():
                feature_icon = menu_icon
            description = " ".join(filter(None, [
                action.toolTip(), action.statusTip(),
            ]))
            full_path = " → ".join(path + [label])
            search_text = " ".join([
                label, menu_label, full_path, description,
            ]).lower()
            # Also match squashed spellings such as "iqtree" for IQ-TREE.
            search_text += " " + re.sub(r"[^a-z0-9]+", "", search_text)
            entry = FeatureEntry(
                label=label,
                menu_label=menu_label,
                menu_icon=menu_icon,
                feature_icon=feature_icon,
                path=full_path,
                description=description,
                action=action,
            )
            index = len(self._entries)
            self._entries.append(entry)

            item = QStandardItem(feature_icon, label)
            item.setEditable(False)
            item.setData(search_text, SEARCH_ROLE)
            item.setData(menu_label, MENU_NAME_ROLE)
            item.setData(menu_icon, MENU_ICON_ROLE)
            item.setData(index, ACTION_INDEX_ROLE)
            item.setData(full_path, Qt.ToolTipRole)
            self.model.appendRow(item)

    # ------------------------------------------------------------------
    # interaction
    # ------------------------------------------------------------------
    def _style_popup(self):
        popup = self.popup_view
        popup.setObjectName("featureSearchPopup")
        popup.setUniformItemSizes(True)
        popup.setMinimumWidth(max(380, min(720, self.line_edit.width())))
        popup.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        popup.setAccessibleName("Feature search results")

    def _on_text_edited(self, text: str):
        text = text.strip()
        if text:
            self.completer.setCompletionPrefix(text)
            self.completer.complete()
            self._style_popup()
        else:
            self.completer.setCompletionPrefix("")
            self.completer.popup().hide()

    def _on_return_pressed(self):
        index = self.completer.popup().currentIndex()
        if not index.isValid():
            completion = self.completer.completionModel()
            if completion.rowCount() == 0:
                return
            index = completion.index(0, 0)
        self._activate_index(index)

    def _activate_index(self, index: QModelIndex):
        if not index.isValid():
            return
        entry_index = index.data(ACTION_INDEX_ROLE)
        if entry_index is None:
            return
        try:
            entry = self._entries[int(entry_index)]
        except (IndexError, TypeError, ValueError):
            return
        self.line_edit.clear()
        self.completer.popup().hide()
        entry.action.trigger()
        # QCompleter writes the selected label after activation; clear once
        # more after its click/key handler has finished.
        QTimer.singleShot(0, self.line_edit.clear)

    def eventFilter(self, obj, event):
        if obj is self.line_edit:
            if event.type() == QEvent.FocusIn:
                self.refresh()
            elif event.type() == QEvent.KeyPress:
                if event.key() == Qt.Key_Escape and not self.completer.popup().isVisible():
                    self.line_edit.clear()
                    return True
        return super().eventFilter(obj, event)
