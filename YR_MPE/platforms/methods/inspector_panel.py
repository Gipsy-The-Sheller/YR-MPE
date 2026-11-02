"""Right-side inspector for YR-MPEA / YPhylo Workbench.

The panel is deliberately host-agnostic:

- it reads items and provenance from a ``DatasetSelectionManager`` and/or
  ``WorkspaceStore``;
- it emits intent signals instead of opening host dialogs;
- it does not set a QStyle or font, so YRTools owns the appearance.

V1 is read-only for data items: preview/open and jump-to-item are wired.
Edges expose their command and a rerun intent; the host decides whether the
rerun tool is available yet.
"""

from __future__ import annotations

import json
from datetime import datetime
from typing import Any, Dict, List, Optional

from PyQt5.QtCore import QEasingCurve, QPropertyAnimation, QSize, Qt, pyqtSignal
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QStyle,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

TYPE_LABELS = {
    "sequence": "Sequence",
    "alignment": "Alignment",
    "model": "Model",
    "distance": "Distance",
    "phylogeny": "Tree",
    "treeset": "Tree set",
    "chain": "MCMC chain",
    "variant": "Variant",
    "coalescent": "Coalescent",
    "clock": "Clock",
}

TYPE_SHORT = {
    "sequence": "Seq",
    "alignment": "Aln",
    "model": "Model",
    "distance": "Dist",
    "phylogeny": "Tree",
    "treeset": "Trees",
    "chain": "Chain",
    "variant": "Var",
    "coalescent": "Coal",
    "clock": "Clock",
}

PREVIEW_LABELS = {
    "tree": "Preview (IcyTree)",
    "alignment": "Preview sequences",
    "model": "Preview model",
    "distance": "Preview distance matrix",
    "chain": "Preview chain",
}


def _attr(obj: Any, *names: str, default=None):
    for name in names:
        if hasattr(obj, name):
            value = getattr(obj, name)
            if value not in (None, ""):
                return value
    return default


def _text(value: Any, limit: int = 240) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        try:
            value = json.dumps(value, ensure_ascii=False)
        except Exception:
            value = str(value)
    value = str(value)
    if len(value) > limit:
        return value[: limit - 1] + "…"
    return value


def _display_time(value: Any) -> str:
    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d %H:%M")
    if isinstance(value, str):
        try:
            return datetime.fromisoformat(value).strftime("%Y-%m-%d %H:%M")
        except Exception:
            pass
    return _text(value)


class InspectorPanel(QWidget):
    """Selection-driven provenance/preview panel."""

    item_activated = pyqtSignal(str)
    dataset_activated = pyqtSignal(str)
    preview_requested = pyqtSignal(str, str)          # item_id, kind
    run_renamed = pyqtSignal(str, str)                # run_id, new_name
    action_requested = pyqtSignal(str, str, object)   # action_id, target_id, params

    def __init__(self, manager=None, store=None, parent=None):
        super().__init__(parent)
        self.manager = manager
        self.store = store
        self.current_item_id: Optional[str] = None
        self.current_run_id: Optional[str] = None
        self.current_dataset_id: Optional[str] = None
        self._item_nodes: Dict[str, QTreeWidgetItem] = {}
        self._run_nodes: Dict[str, QTreeWidgetItem] = {}
        self._metadata_expanded = False
        self._has_operations = False
        self._has_run = False
        self._metadata_collapsed_height = 132
        self._metadata_animation: Optional[QPropertyAnimation] = None

        self.setObjectName("inspector")
        self._icons = {}
        self._build_ui()
        self.show_empty()

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def set_sources(self, manager=None, store=None):
        self.manager = manager if manager is not None else self.manager
        self.store = store if store is not None else self.store
        if self.current_item_id:
            self.show_item(self.current_item_id)
        elif self.current_run_id:
            self.show_run(self.current_run_id)
        elif self.current_dataset_id:
            self.show_dataset(self.current_dataset_id)
        else:
            self.show_empty()

    def show_empty(self):
        self.current_item_id = None
        self.current_run_id = None
        self.current_dataset_id = None
        self.title_label.setText("Inspector")
        self.subtitle_label.setText("Select a data node or dataset")
        self._has_operations = False
        self._has_run = False
        self._metadata_expanded = False
        self._clear_operations()
        self._set_form_rows([])
        self.workflow_tree.clear()
        self.metadata_toggle.setEnabled(False)
        self._apply_metadata_mode(animated=False)
        self._set_secondary_visible(False)

    def show_item(self, item_id: str):
        item = self._get_item(item_id)
        if item is None:
            self.show_empty()
            return

        self.current_item_id = item_id
        self.current_run_id = None
        self.current_dataset_id = _attr(item, "dataset_id")
        item_type = _attr(item, "item_type", "type", default="sequence")
        name = _attr(item, "loci_name", "name", default=item_id)
        dataset_name = self._dataset_name(self.current_dataset_id)
        run = self._run_for_item(item)

        self.title_label.setText(name)
        self.subtitle_label.setText(
            "%s · %s" % (TYPE_LABELS.get(item_type, item_type), dataset_name or "no dataset")
        )

        self._fill_item_operations(item, item_type, run)
        self._fill_item_metadata(item, item_type, dataset_name, run)
        self._fill_workflow(self.current_dataset_id, item_id=item_id, run_id=None)
        self._fill_runs(item)
        self.metadata_toggle.setEnabled(True)
        self._apply_metadata_mode(animated=False)

    def show_run(self, run_id: str):
        store = self._get_store()
        run = store.get_run(run_id) if store is not None else None
        if run is None:
            self.show_empty()
            return

        self.current_item_id = None
        self.current_run_id = run_id
        self.current_dataset_id = self._dataset_for_run(run)
        dataset_name = self._dataset_name(self.current_dataset_id)

        self.title_label.setText(run.name)
        self.subtitle_label.setText(
            "Run · %s · %s" % (run.method, dataset_name or "no dataset")
        )
        self._has_run = False
        self._fill_run_operations(run)
        self._fill_run_metadata(run, dataset_name)
        self._fill_workflow(self.current_dataset_id, item_id=None, run_id=run_id)
        self.run_section.setVisible(False)
        self.metadata_toggle.setEnabled(True)
        self._apply_metadata_mode(animated=False)

    def show_dataset(self, dataset_id: str):
        dataset = self._get_dataset(dataset_id)
        if dataset is None:
            self.show_empty()
            return

        self.current_item_id = None
        self.current_run_id = None
        self.current_dataset_id = dataset_id
        name = _attr(dataset, "name", default=dataset_id)
        self.title_label.setText(name)
        self.subtitle_label.setText("Dataset")
        self._has_operations = False
        self._has_run = False
        self._clear_operations()
        self._fill_dataset_metadata(dataset)
        self._fill_workflow(dataset_id, item_id=None, run_id=None)
        self.run_section.setVisible(False)
        self.metadata_toggle.setEnabled(True)
        self._apply_metadata_mode(animated=False)

    def refresh(self):
        if self.current_item_id:
            self.show_item(self.current_item_id)
        elif self.current_run_id:
            self.show_run(self.current_run_id)
        elif self.current_dataset_id:
            self.show_dataset(self.current_dataset_id)
        else:
            self.show_empty()

    # ------------------------------------------------------------------
    # semantic contract for a future YRTools agent
    # ------------------------------------------------------------------
    def agent_contract(self) -> Dict[str, Any]:
        """Describe state and callable operations without touching Qt events."""
        return {
            "schema": "yrmpe.ui.inspector.v1",
            "id": self.objectName(),
            "state": {
                "item_id": self.current_item_id,
                "run_id": self.current_run_id,
                "dataset_id": self.current_dataset_id,
                "metadata_expanded": self._metadata_expanded,
            },
            "actions": [
                {
                    "id": "inspector.show_item",
                    "description": "Select and display a data item",
                    "parameters": {"item_id": {"type": "string", "required": True}},
                },
                {
                    "id": "inspector.show_run",
                    "description": "Select and display an analysis edge",
                    "parameters": {"run_id": {"type": "string", "required": True}},
                },
                {
                    "id": "inspector.show_dataset",
                    "description": "Display a dataset",
                    "parameters": {"dataset_id": {"type": "string", "required": True}},
                },
                {
                    "id": "inspector.preview",
                    "description": "Open the registered preview for an item",
                    "parameters": {
                        "item_id": {"type": "string", "required": True},
                        "kind": {"type": "string", "required": False},
                    },
                },
                {
                    "id": "inspector.rename_run",
                    "description": "Rename an analysis run",
                    "parameters": {
                        "run_id": {"type": "string", "required": True},
                        "name": {"type": "string", "required": True},
                    },
                },
                {
                    "id": "inspector.set_metadata_expanded",
                    "description": "Expand or collapse the metadata section",
                    "parameters": {"expanded": {"type": "boolean", "required": True}},
                },
                {
                    "id": "inspector.rerun_run",
                    "description": "Rerun an analysis edge with edited parameters",
                    "parameters": {
                        "run_id": {"type": "string", "required": True},
                        "parameters": {"type": "object", "required": False},
                    },
                },
                {
                    "id": "inspector.invoke_operation",
                    "description": "Invoke a visible operation by id",
                    "parameters": {
                        "operation": {"type": "string", "required": True},
                        "target_id": {"type": "string", "required": True},
                    },
                },
            ],
        }

    def agent_invoke(self, action_id: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Dispatch a semantic action. GUI buttons use the same handlers."""
        params = params or {}
        try:
            if action_id == "inspector.show_item":
                self.show_item(params["item_id"])
            elif action_id == "inspector.show_run":
                self.show_run(params["run_id"])
            elif action_id == "inspector.show_dataset":
                self.show_dataset(params["dataset_id"])
            elif action_id == "inspector.preview":
                item_id = params["item_id"]
                kind = params.get("kind") or self._kind_for_item(item_id)
                self.preview_requested.emit(item_id, kind)
            elif action_id == "inspector.rename_run":
                store = self._get_store()
                if store is None:
                    raise RuntimeError("Workspace store is not available")
                store.rename_run(params["run_id"], params["name"])
                self.run_renamed.emit(params["run_id"], params["name"])
                self.refresh()
            elif action_id == "inspector.set_metadata_expanded":
                self._metadata_expanded = bool(params["expanded"])
                self._apply_metadata_mode(animated=False)
            elif action_id == "inspector.rerun_run":
                self.action_requested.emit(
                    "rerun_run",
                    params["run_id"],
                    {"parameters": params.get("parameters", {})},
                )
            elif action_id == "inspector.invoke_operation":
                self._invoke_operation(params["operation"], params["target_id"])
            else:
                raise KeyError("unknown action %r" % action_id)
            return {"ok": True}
        except Exception as exc:
            return {"ok": False, "error": str(exc)}

    def _invoke_operation(self, operation: str, target_id: str):
        if operation == "show":
            self.item_activated.emit(target_id)
        elif operation == "preview":
            self.preview_requested.emit(target_id, self._kind_for_item(target_id))
        elif operation in ("import_lsd2", "import_mcmctree"):
            self.action_requested.emit(operation, target_id, {})
        elif operation == "view_command":
            self._view_command(target_id)
        elif operation == "rerun_run":
            self.action_requested.emit("rerun_run", target_id, {})
        else:
            raise KeyError("unknown operation %r" % operation)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        self.setMinimumWidth(300)
        self.setMaximumWidth(560)
        self.setStyleSheet(
            'QPushButton[class="inspectorLink"] { text-align: left; border: none; '
            "color: palette(link); padding: 1px 0; }"
            'QPushButton[class="inspectorLink"]:hover { text-decoration: underline; }'
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        self.title_label = QLabel("Inspector")
        font = QFont(QApplication.font())
        font.setPointSize(font.pointSize() + 2)
        font.setBold(True)
        self.title_label.setFont(font)
        self.title_label.setWordWrap(True)
        layout.addWidget(self.title_label)

        self.subtitle_label = QLabel("")
        self.subtitle_label.setWordWrap(True)
        self.subtitle_label.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.subtitle_label)

        # Operations --------------------------------------------------
        self.operations_label = self._section_label("Operations")
        layout.addWidget(self.operations_label)
        self.operations_widget = QWidget()
        self.operations_layout = QVBoxLayout(self.operations_widget)
        self.operations_layout.setContentsMargins(0, 0, 0, 0)
        self.operations_layout.setSpacing(1)
        layout.addWidget(self.operations_widget)

        # Metadata ----------------------------------------------------
        self.metadata_label = self._section_label("Metadata")
        layout.addWidget(self.metadata_label)
        self.metadata_widget = QWidget()
        self.metadata_form = QFormLayout(self.metadata_widget)
        self.metadata_form.setContentsMargins(0, 0, 0, 0)
        self.metadata_form.setHorizontalSpacing(6)
        self.metadata_form.setVerticalSpacing(3)
        self.metadata_form.setLabelAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.metadata_form.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.metadata_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)

        self.metadata_scroll = QScrollArea()
        self.metadata_scroll.setObjectName("inspector.metadata.scroll")
        self.metadata_scroll.setWidgetResizable(True)
        self.metadata_scroll.setFrameShape(QFrame.NoFrame)
        self.metadata_scroll.setWidget(self.metadata_widget)
        self.metadata_scroll.setMaximumHeight(self._metadata_collapsed_height)
        self.metadata_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.metadata_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        layout.addWidget(self.metadata_scroll, 1)

        self.metadata_toggle = QToolButton()
        self.metadata_toggle.setObjectName("inspector.metadata.toggle")
        self.metadata_toggle.setArrowType(Qt.DownArrow)
        self.metadata_toggle.setAutoRaise(True)
        self.metadata_toggle.setFixedHeight(16)
        self.metadata_toggle.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.metadata_toggle.setToolTip("Expand metadata")
        self.metadata_toggle.clicked.connect(self._toggle_metadata)
        layout.addWidget(self.metadata_toggle)

        self._metadata_animation = QPropertyAnimation(self.metadata_scroll, b"maximumHeight", self)
        self._metadata_animation.setDuration(180)
        self._metadata_animation.setEasingCurve(QEasingCurve.OutCubic)
        self._metadata_animation.finished.connect(self._on_metadata_animation_finished)

        # Workflow ----------------------------------------------------
        self.workflow_label = self._section_label("Workflow")
        layout.addWidget(self.workflow_label)
        self.workflow_tree = QTreeWidget()
        self.workflow_tree.setObjectName("inspector.workflow")
        self.workflow_tree.setColumnCount(2)
        self.workflow_tree.setHeaderLabels(["Data / method", "Type"])
        self.workflow_tree.setAlternatingRowColors(True)
        self.workflow_tree.setUniformRowHeights(True)
        self.workflow_tree.setIndentation(10)
        self.workflow_tree.setIconSize(QSize(16, 16))
        self.workflow_tree.header().setStretchLastSection(False)
        self.workflow_tree.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.workflow_tree.header().setSectionResizeMode(1, QHeaderView.Fixed)
        self.workflow_tree.setColumnWidth(1, 64)
        self.workflow_tree.itemClicked.connect(self._on_tree_clicked)
        self.workflow_tree.itemDoubleClicked.connect(self._on_tree_double_clicked)
        layout.addWidget(self.workflow_tree, 1)

        # Runs --------------------------------------------------------
        self.run_section = QWidget()
        self.run_section.setObjectName("inspector.runs")
        run_layout = QVBoxLayout(self.run_section)
        run_layout.setContentsMargins(0, 0, 0, 0)
        run_layout.setSpacing(3)
        run_layout.addWidget(self._section_label("Runs"))
        self.run_combo = QComboBox()
        self.run_combo.setObjectName("inspector.runs.combo")
        self.run_combo.currentIndexChanged.connect(self._on_run_combo)
        run_layout.addWidget(self.run_combo)

        rename_row = QHBoxLayout()
        self.run_name_edit = QLineEdit()
        self.run_name_edit.setObjectName("inspector.runs.name")
        self.run_name_edit.setPlaceholderText("Run name")
        self.rename_button = QPushButton("Rename")
        self.rename_button.setObjectName("inspector.runs.rename")
        self.rename_button.clicked.connect(self._on_rename_run)
        rename_row.addWidget(self.run_name_edit, 1)
        rename_row.addWidget(self.rename_button)
        run_layout.addLayout(rename_row)
        layout.addWidget(self.run_section)

    def _section_label(self, text):
        label = QLabel(text)
        font = QFont(QApplication.font())
        font.setBold(True)
        label.setFont(font)
        return label

    def _link_button(self, text, object_name=None):
        button = QPushButton(text)
        button.setProperty("class", "inspectorLink")
        button.setFlat(True)
        button.setCursor(Qt.PointingHandCursor)
        if object_name:
            button.setObjectName(object_name)
        return button

    # ------------------------------------------------------------------
    # metadata collapse / expand
    # ------------------------------------------------------------------
    def _toggle_metadata(self):
        if self.metadata_form.rowCount() == 0:
            return
        self._metadata_expanded = not self._metadata_expanded
        self._apply_metadata_mode(animated=True)

    def _apply_metadata_mode(self, animated=False):
        expanded = self._metadata_expanded
        self.metadata_toggle.setArrowType(Qt.UpArrow if expanded else Qt.DownArrow)
        self.metadata_toggle.setToolTip("Collapse metadata" if expanded else "Expand metadata")
        self.metadata_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarAlwaysOn if expanded else Qt.ScrollBarAlwaysOff
        )
        self._set_secondary_visible(not expanded)

        if self._metadata_animation is not None:
            self._metadata_animation.stop()

        if not animated:
            self.metadata_scroll.setMaximumHeight(
                16777215 if expanded else self._metadata_collapsed_height
            )
            return

        start = max(0, self.metadata_scroll.height())
        if expanded:
            end = max(240, self.height() - 80)
        else:
            end = self._metadata_collapsed_height
        self._metadata_animation.setStartValue(start)
        self._metadata_animation.setEndValue(end)
        self._metadata_animation.start()

    def _on_metadata_animation_finished(self):
        if self._metadata_expanded:
            self.metadata_scroll.setMaximumHeight(16777215)

    def _set_secondary_visible(self, visible):
        self._set_operations_visible(visible and self._has_operations)
        self.workflow_label.setVisible(visible)
        self.workflow_tree.setVisible(visible)
        self.run_section.setVisible(visible and self._has_run)

    # ------------------------------------------------------------------
    # data access
    # ------------------------------------------------------------------
    def _get_store(self):
        if self.store is not None:
            return self.store
        if self.manager is not None and getattr(self.manager, "store", None):
            return self.manager.store
        return None

    def _get_item(self, item_id: str):
        if self.manager is not None:
            item = self.manager.get_item(item_id)
            if item is not None:
                return item
        store = self._get_store()
        if store is not None:
            return store.get_item(item_id)
        return None

    def _get_dataset(self, dataset_id: str):
        if self.manager is not None:
            dataset = self.manager.get_dataset(dataset_id)
            if dataset is not None:
                return dataset
        store = self._get_store()
        if store is not None:
            return store.get_dataset(dataset_id)
        return None

    def _dataset_name(self, dataset_id: str) -> str:
        dataset = self._get_dataset(dataset_id) if dataset_id else None
        return _attr(dataset, "name", default="") if dataset else ""

    def _dataset_items(self, dataset_id: str) -> List[Any]:
        if self.manager is not None:
            items = self.manager.get_items_by_dataset(dataset_id)
            if items:
                return items
        store = self._get_store()
        if store is not None:
            return store.items_in_dataset(dataset_id)
        return []

    def _run_for_item(self, item):
        run_id = _attr(item, "run_id", default=None)
        store = self._get_store()
        if store is None:
            return None
        if run_id:
            return store.get_run(run_id)
        item_id = _attr(item, "id", default="")
        for run in store.index.runs.values():
            if item_id in run.output_ids:
                return run
        return None

    def _dataset_for_run(self, run):
        store = self._get_store()
        if store is None:
            return None
        for dataset in store.index.datasets.values():
            if set(run.input_ids) & set(dataset.item_ids):
                return dataset.id
            if set(run.output_ids) & set(dataset.item_ids):
                return dataset.id
        return None

    def _stream_for_run(self, run):
        store = self._get_store()
        if run is None or store is None:
            return None
        for stream in store.index.streams.values():
            if run.id in stream.run_ids:
                return stream
        return None

    # ------------------------------------------------------------------
    # rendering
    # ------------------------------------------------------------------
    def _preview_kind(self, item_type: str):
        if item_type in ("phylogeny", "treeset"):
            return "tree"
        if item_type in ("sequence", "alignment"):
            return "alignment"
        if item_type == "model":
            return "model"
        if item_type == "distance":
            return "distance"
        if item_type == "chain":
            return "chain"
        return None

    def _fill_item_operations(self, item, item_type, run):
        self._has_operations = True
        self._clear_operations()
        self._add_operation("Show in dataset", self._on_show, "inspector.operation.show")
        kind = self._preview_kind(item_type)
        if kind is not None:
            self._add_operation(
                PREVIEW_LABELS.get(kind, "Preview"),
                self._on_preview,
                "inspector.operation.preview",
            )
        if kind == "tree":
            self._add_operation(
                "Import to LSD2", lambda: self._emit_operation("import_lsd2"),
                "inspector.operation.import_lsd2",
            )
            self._add_operation(
                "Import to MCMCTree", lambda: self._emit_operation("import_mcmctree"),
                "inspector.operation.import_mcmctree",
            )
        self._set_operations_visible(True)

    def _fill_run_operations(self, run):
        self._has_operations = True
        self._clear_operations()
        self._add_operation(
            "View command", lambda: self._view_command(run.id),
            "inspector.operation.view_command",
        )
        self._add_operation(
            "Modify && rerun", lambda: self._emit_operation("rerun_run", run.id),
            "inspector.operation.rerun",
        )
        self._set_operations_visible(True)

    def _set_operations_visible(self, visible):
        self.operations_label.setVisible(bool(visible))
        self.operations_widget.setVisible(bool(visible))

    def _clear_operations(self):
        while self.operations_layout.count():
            item = self.operations_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.hide()
                widget.setParent(None)
                widget.deleteLater()

    def _add_operation(self, text, slot, object_name=None):
        button = self._link_button(text, object_name)
        button.clicked.connect(slot)
        self.operations_layout.addWidget(button)

    def _emit_operation(self, operation, target_id=None):
        target_id = target_id or self.current_item_id or self.current_run_id
        if target_id:
            self.action_requested.emit(operation, target_id, {})

    def _fill_item_metadata(self, item, item_type: str, dataset_name: str, run):
        rows = [
            ("Type", TYPE_LABELS.get(item_type, item_type)),
            ("Dataset", dataset_name),
            ("File", _attr(item, "file_path", "file", default="")),
            ("Format", _attr(item, "format", default="")),
            ("Length", _attr(item, "length", default="")),
            ("Sequences", _attr(item, "sequence_count", default="")),
            ("Aligned", "yes" if _attr(item, "is_aligned", default=False) else "no"),
            ("Import", _attr(item, "import_id", default="")),
            ("Parent", _attr(item, "parent_id", default="")),
            ("Created", _display_time(_attr(item, "created_at", default=""))),
        ]
        if run is not None:
            link = self._link_button("Open edge: %s" % run.name, "inspector.metadata.edge")
            link.clicked.connect(lambda: self.show_run(run.id))
            rows.append(("Edge", link))
        self._set_form_rows(rows)

    def _fill_run_metadata(self, run, dataset_name: str):
        stream = self._stream_for_run(run)
        rows = [
            ("Method", run.method),
            ("Tool", run.tool + (" " + run.tool_version if run.tool_version else "")),
            ("Status", run.status),
            ("Dataset", dataset_name),
            ("Created", _display_time(run.created_at)),
            ("Completed", _display_time(run.completed_at)),
            ("Inputs", ", ".join(run.input_ids)),
            ("Outputs", ", ".join(run.output_ids)),
            ("Parameters", run.parameters),
            ("Command", " ".join(run.command) if run.command else ""),
        ]
        if stream is not None:
            rows.append(("Stream", stream.name))
        self._set_form_rows(rows)

    def _fill_dataset_metadata(self, dataset):
        rows = [
            ("Name", _attr(dataset, "name", default="")),
            ("Description", _attr(dataset, "description", default="")),
            ("Source", _attr(dataset, "source", default="")),
            ("Loci", _attr(dataset, "loci_count", default="")),
            ("Max taxa", _attr(dataset, "taxa_count", default="")),
        ]
        settings = _attr(dataset, "settings", default={}) or {}
        for key, value in settings.get("summary", {}).items():
            rows.append((key, _text(value)))
        self._set_form_rows(rows)

    def _set_form_rows(self, rows):
        while self.metadata_form.rowCount():
            self.metadata_form.removeRow(0)
        for field, value in rows:
            if value in (None, ""):
                continue
            label = QLabel(str(field))
            label.setStyleSheet("color: palette(mid);")
            if isinstance(value, QWidget):
                self.metadata_form.addRow(label, value)
            else:
                full = _text(value)
                display = full
                if field in ("Parent", "Inputs", "Outputs") and len(display) > 30:
                    display = display[:14] + "…" + display[-12:]
                elif len(display) > 42 and not any(ch.isspace() for ch in display):
                    display = display[:20] + "…" + display[-18:]
                field_label = QLabel(display)
                field_label.setWordWrap(True)
                field_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                field_label.setToolTip(full)
                self.metadata_form.addRow(label, field_label)

    # ------------------------------------------------------------------
    # workflow tree
    # ------------------------------------------------------------------
    def _icon(self, name):
        if name not in self._icons:
            style = self.style()
            mapping = {
                "dataset": QStyle.SP_DirIcon,
                "item": QStyle.SP_FileIcon,
                "import": QStyle.SP_ArrowDown,
                "run": QStyle.SP_MediaPlay,
            }
            self._icons[name] = style.standardIcon(mapping[name])
        return self._icons[name]

    def _fill_workflow(self, dataset_id: str, item_id: Optional[str] = None,
                       run_id: Optional[str] = None):
        self.workflow_tree.clear()
        self._item_nodes = {}
        self._run_nodes = {}
        if not dataset_id:
            return
        dataset = self._get_dataset(dataset_id)
        items = self._dataset_items(dataset_id)
        if dataset is None and not items:
            return

        item_ids = {_attr(item, "id", default="") for item in items}
        root = QTreeWidgetItem([
            _attr(dataset, "name", default="Dataset") if dataset else "Dataset",
            "DS",
        ])
        root.setData(0, Qt.UserRole, ("dataset", dataset_id))
        root.setIcon(0, self._icon("dataset"))
        root.setExpanded(True)
        self.workflow_tree.addTopLevelItem(root)

        for item in items:
            record_id = _attr(item, "id", default="")
            node = QTreeWidgetItem([
                _attr(item, "loci_name", "name", default=record_id),
                TYPE_SHORT.get(_attr(item, "item_type", "type", default=""), ""),
            ])
            node.setData(0, Qt.UserRole, ("item", record_id))
            node.setIcon(0, self._icon("item"))
            node.setToolTip(0, _text(_attr(item, "file_path", "file", default=""), 500))
            self._item_nodes[record_id] = node

        store = self._get_store()
        import_nodes = {}
        if store is not None:
            for record in store.index.imports.values():
                if record.dataset_id != dataset_id:
                    continue
                node = QTreeWidgetItem([
                    "Import: %s" % record.importer,
                    "IMP",
                ])
                node.setData(0, Qt.UserRole, ("import", record.id))
                node.setIcon(0, self._icon("import"))
                node.setToolTip(0, "\n".join(record.source_paths))
                import_nodes[record.id] = node
                root.addChild(node)

            for record in store.index.runs.values():
                if not (set(record.input_ids) | set(record.output_ids)) & item_ids:
                    continue
                node = QTreeWidgetItem([
                    "Run: %s" % record.name,
                    "RUN",
                ])
                node.setData(0, Qt.UserRole, ("run", record.id))
                node.setIcon(0, self._icon("run"))
                node.setToolTip(0, "Status: %s" % record.status)
                self._run_nodes[record.id] = node

        for item in items:
            record_id = _attr(item, "id", default="")
            node = self._item_nodes[record_id]
            producing_run = _attr(item, "run_id", default=None)
            import_id = _attr(item, "import_id", default=None)
            parent_id = _attr(item, "parent_id", default=None)
            if producing_run and producing_run in self._run_nodes:
                self._run_nodes[producing_run].addChild(node)
            elif import_id and import_id in import_nodes:
                import_nodes[import_id].addChild(node)
            elif parent_id and parent_id in self._item_nodes:
                self._item_nodes[parent_id].addChild(node)
            else:
                root.addChild(node)

        for record in (store.index.runs.values() if store is not None else []):
            if record.id not in self._run_nodes:
                continue
            node = self._run_nodes[record.id]
            anchor = None
            for input_id in record.input_ids:
                if input_id in self._item_nodes:
                    anchor = self._item_nodes[input_id]
                    break
            if anchor is not None:
                anchor.addChild(node)
            else:
                root.addChild(node)

        for node in import_nodes.values():
            node.setExpanded(node.childCount() <= 30)
        for node in self._run_nodes.values():
            node.setExpanded(node.childCount() <= 30)

        selected = None
        if run_id and run_id in self._run_nodes:
            selected = self._run_nodes[run_id]
        elif item_id and item_id in self._item_nodes:
            selected = self._item_nodes[item_id]
        if selected is not None:
            self.workflow_tree.setCurrentItem(selected)
            self.workflow_tree.scrollToItem(selected)
            self._bold_ancestors(selected)
        else:
            self.workflow_tree.setCurrentItem(root)
        self.workflow_tree.expandItem(root)

    def _bold_ancestors(self, node):
        bold = QFont(self.workflow_tree.font())
        bold.setBold(True)
        current = node.parent()
        while current is not None:
            for column in range(self.workflow_tree.columnCount()):
                current.setFont(column, bold)
            current = current.parent()

    def _on_tree_clicked(self, node, column):
        payload = node.data(0, Qt.UserRole)
        if not payload:
            return
        kind, record_id = payload
        if kind == "item":
            self.show_item(record_id)
            self.item_activated.emit(record_id)
        elif kind == "dataset":
            self.show_dataset(record_id)
            self.dataset_activated.emit(record_id)
        elif kind == "run":
            self.show_run(record_id)

    def _on_tree_double_clicked(self, node, column):
        payload = node.data(0, Qt.UserRole)
        if payload and payload[0] == "item":
            self.preview_requested.emit(payload[1], self._kind_for_item(payload[1]))

    def _kind_for_item(self, item_id: str) -> str:
        item = self._get_item(item_id)
        if item is None:
            return ""
        return self._preview_kind(_attr(item, "item_type", "type", default="")) or ""

    # ------------------------------------------------------------------
    # operations
    # ------------------------------------------------------------------
    def _on_show(self):
        if self.current_item_id:
            self.item_activated.emit(self.current_item_id)

    def _on_preview(self):
        if not self.current_item_id:
            return
        self.preview_requested.emit(
            self.current_item_id, self._kind_for_item(self.current_item_id)
        )

    def _view_command(self, run_id):
        store = self._get_store()
        run = store.get_run(run_id) if store is not None else None
        if run is None:
            return
        dialog = QDialog(self)
        dialog.setWindowTitle("Run command")
        dialog.resize(720, 420)
        layout = QVBoxLayout(dialog)
        text = QPlainTextEdit()
        text.setReadOnly(True)
        lines = []
        if run.command:
            lines.append("$ " + " ".join(run.command))
        else:
            lines.append("(no command recorded)")
        if run.parameters:
            lines.append("")
            lines.append("parameters:")
            lines.append(json.dumps(run.parameters, indent=2, ensure_ascii=False))
        text.setPlainText("\n".join(lines))
        layout.addWidget(text)
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec_()

    # ------------------------------------------------------------------
    # runs
    # ------------------------------------------------------------------
    def _siblings(self, item) -> List[Any]:
        item_type = _attr(item, "item_type", "type", default="")
        dataset_id = _attr(item, "dataset_id", default="")
        parent_id = _attr(item, "parent_id", default=None)
        candidates = []
        if self.manager is not None:
            candidates = list(self.manager.items.values())
        else:
            store = self._get_store()
            if store is not None:
                candidates = list(store.index.items.values())
        siblings = [
            other for other in candidates
            if _attr(other, "dataset_id", default="") == dataset_id
            and _attr(other, "item_type", "type", default="") == item_type
            and _attr(other, "parent_id", default=None) == parent_id
        ]
        siblings.sort(key=lambda other: _display_time(_attr(other, "created_at", default="")))
        return siblings

    def _fill_runs(self, item):
        run = self._run_for_item(item)
        siblings = self._siblings(item)
        if run is None:
            siblings = []
        else:
            siblings = [other for other in siblings if self._run_for_item(other) is not None]
            if not any(_attr(other, "id", default="") == _attr(item, "id", default="") for other in siblings):
                siblings.insert(0, item)
        self._has_run = run is not None
        self.run_section.setVisible(run is not None)

        self.run_combo.blockSignals(True)
        self.run_combo.clear()
        current_index = -1
        for sibling in siblings:
            sibling_run = self._run_for_item(sibling)
            label = _attr(sibling, "loci_name", "name", default="")
            if sibling_run is not None:
                label = "%s — %s" % (sibling_run.name, label)
            self.run_combo.addItem(label, _attr(sibling, "id", default=""))
            if _attr(sibling, "id", default="") == _attr(item, "id", default=""):
                current_index = self.run_combo.count() - 1
        if current_index >= 0:
            self.run_combo.setCurrentIndex(current_index)
        self.run_combo.blockSignals(False)

        self.run_name_edit.setText(run.name if run is not None else "")
        self.run_name_edit.setEnabled(run is not None)
        self.rename_button.setEnabled(run is not None)

    def _on_run_combo(self, index):
        if index < 0:
            return
        item_id = self.run_combo.itemData(index)
        if item_id:
            self.item_activated.emit(item_id)

    def _on_rename_run(self):
        if self.current_item_id:
            item = self._get_item(self.current_item_id)
            run = self._run_for_item(item) if item is not None else None
        else:
            store = self._get_store()
            run = store.get_run(self.current_run_id) if store is not None else None
        store = self._get_store()
        if run is None or store is None:
            return
        name = self.run_name_edit.text().strip()
        if not name or name == run.name:
            return
        store.rename_run(run.id, name)
        self.run_renamed.emit(run.id, name)
        self.refresh()
