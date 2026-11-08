"""Qt workflow previewer: data circles, tool squares and orthogonal edges."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

from PyQt5.QtCore import QPointF, QRectF, Qt, QTimer, pyqtSignal
from PyQt5.QtGui import (
    QBrush, QColor, QFont, QIcon, QPainter, QPainterPath, QPen, QPolygonF,
)
from PyQt5.QtWidgets import (
    QDialog, QGraphicsItem, QGraphicsScene, QGraphicsView, QHBoxLayout,
    QLabel, QPushButton, QVBoxLayout,
)

from ...runtime.workflow_graph import WorkflowGraph, WorkflowNode, layered_layout

ICON_ROOT = Path(__file__).resolve().parents[2] / "icons"
DATA_RADIUS = 36.0
TOOL_HALF = 34.0
LABEL_HEIGHT = 26.0
NODE_FONT = QFont("Sans Serif", 9)


def _icon(icon_name: str) -> QIcon:
    path = ICON_ROOT / str(icon_name or "file/sequence.svg")
    if not path.exists():
        path = ICON_ROOT / "file" / "sequence.svg"
    return QIcon(str(path))


class WorkflowNodeItem(QGraphicsItem):
    """A data circle or tool rounded square with an icon and caption."""

    def __init__(self, node: WorkflowNode):
        super().__init__()
        self.node = node
        self._hover = False
        self.setAcceptHoverEvents(True)
        self.setToolTip(
            "%s\nKind: %s\n%s"
            % (node.label, node.kind,
               "\n".join("%s: %s" % (key, value)
                         for key, value in node.metadata.items()
                         if value not in ("", None)))
        )
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        if node.kind == "tool":
            self.setCursor(Qt.PointingHandCursor)
        self.setZValue(2)

    def _half_width(self) -> float:
        return DATA_RADIUS if self.node.kind == "data" else TOOL_HALF

    def boundingRect(self) -> QRectF:
        half = self._half_width()
        return QRectF(-half - 5, -half - 5, 2 * half + 10,
                      half + LABEL_HEIGHT + 12)

    def shape(self) -> QPainterPath:
        half = self._half_width()
        path = QPainterPath()
        if self.node.kind == "data":
            path.addEllipse(QPointF(0, 0), half, half)
        else:
            path.addRoundedRect(
                QRectF(-half, -half, 2 * half, 2 * half), 6, 6
            )
        return path

    def edge_point(self, side: str) -> QPointF:
        half = self._half_width()
        if side == "left":
            return QPointF(-half, 0)
        if side == "right":
            return QPointF(half, 0)
        if side == "top":
            return QPointF(0, -half)
        return QPointF(0, half)

    def hoverEnterEvent(self, event):
        self._hover = True
        self.update()
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event):
        self._hover = False
        self.update()
        super().hoverLeaveEvent(event)

    def paint(self, painter: QPainter, option, widget=None):
        painter.setRenderHint(QPainter.Antialiasing, True)
        half = self._half_width()
        border = QColor("#111111") if not self._hover else QColor("#2563eb")
        pen = QPen(border, 5.0)
        brush = QBrush(QColor("#ffffff"))
        painter.setPen(pen)
        painter.setBrush(brush)
        if self.node.kind == "data":
            painter.drawEllipse(QPointF(0, 0), half, half)
            icon_size = 32
        else:
            painter.drawRoundedRect(
                QRectF(-half, -half, 2 * half, 2 * half), 6, 6
            )
            icon_size = 42
        pixmap = _icon(self.node.icon).pixmap(icon_size, icon_size)
        if not pixmap.isNull():
            painter.drawPixmap(
                int(-pixmap.width() / 2),
                int(-pixmap.height() / 2),
                pixmap,
            )
        painter.setPen(QPen(QColor("#333333")))
        painter.setFont(NODE_FONT)
        metrics = painter.fontMetrics()
        label = metrics.elidedText(
            self.node.label or self.node.id, Qt.ElideMiddle, 150
        )
        painter.drawText(
            QRectF(-78, half + 5, 156, LABEL_HEIGHT),
            Qt.AlignHCenter | Qt.AlignTop,
            label,
        )


class WorkflowEdgeItem(QGraphicsItem):
    """Orthogonal edge with a filled arrowhead."""

    def __init__(self, source: WorkflowNodeItem, target: WorkflowNodeItem):
        super().__init__()
        self.source = source
        self.target = target
        self._path = self._build_path()
        self.setZValue(0)

    def _build_path(self) -> QPainterPath:
        start = self.source.scenePos() + self.source.edge_point("right")
        end = self.target.scenePos() + self.target.edge_point("left")
        path = QPainterPath(start)
        if end.x() >= start.x():
            middle = (start.x() + end.x()) / 2.0
            path.lineTo(middle, start.y())
            path.lineTo(middle, end.y())
            path.lineTo(end)
        else:
            # Backward edge: leave right, route below/above, enter from left.
            offset = 54.0
            middle = max(start.x(), end.x()) + offset
            path.lineTo(start.x() + offset, start.y())
            path.lineTo(start.x() + offset, end.y() - offset)
            path.lineTo(end.x() - offset, end.y() - offset)
            path.lineTo(end.x() - offset, end.y())
            path.lineTo(end)
        return path

    def boundingRect(self) -> QRectF:
        return self._path.boundingRect().adjusted(-8, -8, 8, 8)

    def paint(self, painter: QPainter, option, widget=None):
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.setPen(QPen(QColor("#111111"), 3.2,
                            Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin))
        painter.setBrush(Qt.NoBrush)
        painter.drawPath(self._path)
        end = self._path.currentPosition()
        before = self._path.pointAtPercent(0.98)
        direction = end - before
        length = max(direction.manhattanLength(), 1.0)
        ux, uy = direction.x() / length, direction.y() / length
        size = 12.0
        base = end - QPointF(ux * size, uy * size)
        normal = QPointF(-uy, ux) * (size * 0.52)
        polygon = QPolygonF([
            end,
            base + normal,
            base - normal,
        ])
        painter.setBrush(QBrush(QColor("#111111")))
        painter.setPen(Qt.NoPen)
        painter.drawPolygon(polygon)


class WorkflowPreviewWidget(QGraphicsView):
    """Scrollable/zoomable workflow canvas."""

    node_clicked = pyqtSignal(str)

    def __init__(self, graph: Optional[WorkflowGraph] = None, parent=None):
        super().__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.setRenderHints(
            QPainter.Antialiasing | QPainter.TextAntialiasing
            | QPainter.SmoothPixmapTransform
        )
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorViewCenter)
        self.setBackgroundBrush(QBrush(QColor("#ffffff")))
        self.setFrameShape(QGraphicsView.NoFrame)
        self._node_items: Dict[str, WorkflowNodeItem] = {}
        self._graph: Optional[WorkflowGraph] = None
        if graph is not None:
            self.set_graph(graph)

    def set_graph(self, graph: WorkflowGraph):
        scene = self.scene()
        scene.clear()
        self._node_items.clear()
        self._graph = graph
        if not graph.nodes:
            item = scene.addText("No workflow data yet.")
            item.setDefaultTextColor(QColor("#777777"))
            item.setPos(0, 0)
            self.setSceneRect(item.boundingRect().adjusted(-30, -30, 30, 30))
            QTimer.singleShot(0, self.fit_view)
            return

        positions = layered_layout(graph)
        for node in graph.nodes:
            item = WorkflowNodeItem(node)
            item.setPos(*positions.get(node.id, (0, 0)))
            scene.addItem(item)
            self._node_items[node.id] = item
        for edge in graph.edges:
            source = self._node_items.get(edge.source)
            target = self._node_items.get(edge.target)
            if source is None or target is None:
                continue
            scene.addItem(WorkflowEdgeItem(source, target))
        rect = scene.itemsBoundingRect().adjusted(-40, -40, 40, 40)
        scene.setSceneRect(rect)
        QTimer.singleShot(0, self.fit_view)

    def fit_view(self):
        if self.scene() is None:
            return
        rect = self.scene().itemsBoundingRect()
        if rect.isEmpty():
            return
        self.fitInView(rect.adjusted(-20, -20, 20, 20),
                       Qt.KeepAspectRatio)

    def zoom_by(self, factor: float):
        current = self.transform().m11()
        target = current * factor
        if target < 0.25 or target > 3.5:
            return
        self.scale(factor, factor)

    def mousePressEvent(self, event):
        item = self.itemAt(event.pos())
        while item is not None and not isinstance(item, WorkflowNodeItem):
            item = item.parentItem()
        if isinstance(item, WorkflowNodeItem) and item.node.kind == "tool":
            self.node_clicked.emit(item.node.id)
            event.accept()
            return
        super().mousePressEvent(event)

    def wheelEvent(self, event):
        self.zoom_by(1.12 if event.angleDelta().y() > 0 else 1 / 1.12)
        event.accept()


class WorkflowPreviewDialog(QDialog):
    """Non-modal dialog wrapper with fit/zoom controls."""

    def __init__(self, graph: WorkflowGraph, parent=None, step_opener=None):
        super().__init__(parent)
        self.step_opener = step_opener
        self.setWindowTitle("Workflow Preview - YPhylo Workbench")
        self.setMinimumSize(900, 520)
        layout = QVBoxLayout(self)
        toolbar = QHBoxLayout()
        title = QLabel(graph.title or "Workflow")
        title.setStyleSheet("font-size: 16px; font-weight: 600;")
        toolbar.addWidget(title)
        toolbar.addStretch(1)
        fit_button = QPushButton("Fit")
        fit_button.clicked.connect(lambda: self.preview.fit_view())
        zoom_in = QPushButton("+")
        zoom_in.setFixedWidth(34)
        zoom_in.clicked.connect(lambda: self.preview.zoom_by(1.2))
        zoom_out = QPushButton("-")
        zoom_out.setFixedWidth(34)
        zoom_out.clicked.connect(lambda: self.preview.zoom_by(1 / 1.2))
        toolbar.addWidget(zoom_out)
        toolbar.addWidget(zoom_in)
        toolbar.addWidget(fit_button)
        layout.addLayout(toolbar)
        self.preview = WorkflowPreviewWidget(graph)
        self.preview.node_clicked.connect(self._on_node_clicked)
        layout.addWidget(self.preview, 1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        row = QHBoxLayout()
        row.addStretch(1)
        row.addWidget(close_button)
        layout.addLayout(row)

    def _on_node_clicked(self, node_id: str):
        if not callable(self.step_opener):
            return
        graph = self.preview._graph
        if graph is None:
            return
        node = graph.node_map().get(node_id)
        if node is None or node.kind != "tool":
            return
        context = dict(node.metadata)
        context.setdefault("step_id", context.get("step", node.id))
        try:
            config = self.step_opener(context)
        except Exception:
            return
        if isinstance(config, dict):
            node.metadata["config"] = config
            item = self.preview._node_items.get(node_id)
            if item is not None:
                item.setToolTip(
                    "%s\nConfigured step; click to adjust parameters."
                    % node.label
                )

