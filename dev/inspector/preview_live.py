"""Live, interactive inspector preview.

Run from the repository root:

    python dev/inspector/preview_live.py

Unlike render_inspector.py, this keeps the window open and wires the panel to
a small mock host so item selection, preview, run switching/rename, metadata
expansion, and operation intents can be exercised interactively.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

if os.environ.get("QT_QPA_PLATFORM") is None:
    os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "xcb"

from PyQt5.QtCore import QTimer, Qt  # noqa: E402
from PyQt5.QtGui import QFont, QFontDatabase  # noqa: E402
from PyQt5.QtWidgets import (  # noqa: E402
    QApplication,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

HERE = Path(__file__).resolve()
DEV_ROOT = HERE.parents[2]
sys.path.insert(0, str(DEV_ROOT))
sys.path.insert(0, str(HERE.parent))

from render_inspector import build_sample_store  # noqa: E402
from YR_MPE.platforms.methods.inspector_panel import InspectorPanel  # noqa: E402


class PreviewDialog(QDialog):
    def __init__(self, title: str, text: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(720, 480)
        layout = QVBoxLayout(self)
        view = QPlainTextEdit()
        view.setReadOnly(True)
        view.setPlainText(text)
        layout.addWidget(view)
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)


class LivePreview(QMainWindow):
    def __init__(self, store, tree_item):
        super().__init__()
        self.store = store
        self.tree_item = tree_item
        self.setWindowTitle("YPhylo Workbench — Inspector preview")
        self.resize(1080, 820)

        self.panel = InspectorPanel(store=store)
        self.panel.item_activated.connect(self._select_item_in_list)
        self.panel.dataset_activated.connect(self._select_dataset_in_list)
        self.panel.preview_requested.connect(self._on_preview)
        self.panel.action_requested.connect(self._on_action)
        self.panel.run_renamed.connect(self._on_run_renamed)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(10, 10, 10, 10)
        title = QLabel("Workspace (mock)")
        title_font = QFont(QApplication.font())
        title_font.setBold(True)
        title_font.setPointSize(title_font.pointSize() + 1)
        title.setFont(title_font)
        left_layout.addWidget(title)
        hint = QLabel("Click a node or a run. The inspector on the right is the real widget.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: palette(mid);")
        left_layout.addWidget(hint)

        self.list = QListWidget()
        self.list.itemClicked.connect(self._on_list_clicked)
        left_layout.addWidget(self.list, 1)
        self._populate_list()

        splitter = QSplitter(Qt.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(left)
        splitter.addWidget(self.panel)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([420, 420])
        self.setCentralWidget(splitter)

        self.panel.show_item(tree_item.id)
        self._select_item_in_list(tree_item.id)

    def _populate_list(self):
        dataset = self.store.get_dataset(self.tree_item.dataset_id)
        root = QListWidgetItem("Dataset: %s" % dataset.name)
        root.setData(Qt.UserRole, ("dataset", dataset.id))
        self.list.addItem(root)

        for item in self.store.items_in_dataset(dataset.id):
            if item.run_id:
                continue
            node = QListWidgetItem("  %s  [%s]" % (item.name, item.type))
            node.setData(Qt.UserRole, ("item", item.id))
            self.list.addItem(node)

        for run in self.store.index.runs.values():
            if not (set(run.input_ids) | set(run.output_ids)) & set(dataset.item_ids):
                continue
            node = QListWidgetItem("  %s  [run]" % run.name)
            node.setData(Qt.UserRole, ("run", run.id))
            self.list.addItem(node)
            for output_id in run.output_ids:
                output = self.store.get_item(output_id)
                if output is None:
                    continue
                child = QListWidgetItem("      %s  [%s]" % (output.name, output.type))
                child.setData(Qt.UserRole, ("item", output.id))
                self.list.addItem(child)

    def _on_list_clicked(self, node):
        payload = node.data(Qt.UserRole)
        if not payload:
            return
        kind, record_id = payload
        if kind == "dataset":
            self.panel.show_dataset(record_id)
        elif kind == "run":
            self.panel.show_run(record_id)
        else:
            self.panel.show_item(record_id)

    def _select_item_in_list(self, item_id):
        for index in range(self.list.count()):
            node = self.list.item(index)
            payload = node.data(Qt.UserRole)
            if payload == ("item", item_id):
                self.list.setCurrentItem(node)
                return

    def _select_dataset_in_list(self, dataset_id):
        for index in range(self.list.count()):
            node = self.list.item(index)
            if node.data(Qt.UserRole) == ("dataset", dataset_id):
                self.list.setCurrentItem(node)
                return

    def _on_preview(self, item_id, kind):
        item = self.store.get_item(item_id)
        if item is None:
            return
        text = ""
        if item.file and self.store.path(item.file).exists():
            text = self.store.path(item.file).read_text(encoding="utf-8", errors="replace")
        if not text:
            text = "No file content; metadata: %s" % item.metadata
        if item.type in ("phylogeny", "treeset"):
            title = "IcyTree preview (mock)"
        elif item.type in ("sequence", "alignment"):
            title = "Sequence preview (mock)"
        else:
            title = "%s preview (mock)" % item.type
        PreviewDialog(title, text, self).exec_()

    def _on_action(self, action_id, target_id, params):
        QMessageBox.information(
            self,
            "Inspector operation",
            "action: %s\ntarget: %s\nparameters: %s"
            % (action_id, target_id, params),
        )

    def _on_run_renamed(self, run_id, name):
        for index in range(self.list.count()):
            node = self.list.item(index)
            if node.data(Qt.UserRole) == ("run", run_id):
                node.setText("  %s  [run]" % name)
                return


def main() -> int:
    app = QApplication(sys.argv)
    font_path = Path("E:/Codes/YRTools/fonts/MiSans-Medium.ttf")
    if font_path.exists():
        font_id = QFontDatabase.addApplicationFont(str(font_path))
        families = QFontDatabase.applicationFontFamilies(font_id) if font_id != -1 else []
        if families:
            app_font = QFont(families[0], 10)
            app_font.setWeight(QFont.Medium)
            app.setFont(app_font)

    root = Path(tempfile.mkdtemp())
    store, tree_item = build_sample_store(root)
    window = LivePreview(store, tree_item)
    window.show()
    if "--smoke" in sys.argv:
        QTimer.singleShot(800, app.quit)
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
