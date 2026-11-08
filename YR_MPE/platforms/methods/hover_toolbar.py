"""Icon-first toolbar buttons that expand their label on hover."""

from __future__ import annotations

from PyQt5.QtCore import (
    QEasingCurve, QEvent, QObject, QPropertyAnimation, Qt, pyqtSignal,
)
from PyQt5.QtWidgets import QToolButton

COLLAPSED_WIDTH = 56


class HoverExpander(QObject):
    """Animate a QToolButton's width so its text appears on hover."""

    expansionChanged = pyqtSignal(bool)

    def __init__(self, button: QToolButton, collapsed: int = COLLAPSED_WIDTH):
        super().__init__(button)
        self.button = button
        self.collapsed = collapsed
        self._animation = QPropertyAnimation(button, b"maximumWidth", self)
        self._animation.setDuration(180)

        button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        button.setMinimumWidth(collapsed)
        button.setMaximumWidth(collapsed)
        # InstantPopup buttons do not need a menu arrow; hiding it keeps the
        # collapsed state icon-centred.
        style = button.styleSheet() or ""
        button.setStyleSheet(
            style + "\nQToolButton::menu-indicator { image: none; width: 0px; }"
        )
        button.installEventFilter(self)

    def expanded_width(self) -> int:
        text = self.button.text().strip()
        if not text:
            return self.collapsed
        icon_width = self.button.iconSize().width() or 45
        text_width = self.button.fontMetrics().horizontalAdvance(text)
        # QToolButton's sizeHint includes the style's icon/text margins; use it
        # so the label is not elided mid-word at the end of the animation.
        return max(
            self.collapsed + 16,
            icon_width + text_width + 32,
            self.button.sizeHint().width(),
        )

    def _animate(self, target: int, duration: int, curve: QEasingCurve):
        self._animation.stop()
        self._animation.setStartValue(self.button.maximumWidth())
        self._animation.setEndValue(target)
        self._animation.setDuration(duration)
        self._animation.setEasingCurve(curve)
        self._animation.start()

    def eventFilter(self, obj, event):
        if obj is self.button:
            if event.type() == QEvent.Enter:
                self.expansionChanged.emit(True)
                self._animate(
                    self.expanded_width(), 180, QEasingCurve.OutCubic
                )
            elif event.type() == QEvent.Leave:
                self.expansionChanged.emit(False)
                self._animate(
                    self.collapsed, 150, QEasingCurve.InCubic
                )
        return super().eventFilter(obj, event)


def install_hover_expand(button: QToolButton, collapsed: int = COLLAPSED_WIDTH):
    """Make *button* icon-only until the pointer enters it."""
    if not button.text().strip():
        return None
    return HoverExpander(button, collapsed)
