"""Right-side inspector for YR-MPEA / YPhylo Workbench.

The panel is deliberately host-agnostic:

- it reads items and provenance from a ``DatasetSelectionManager`` and/or
  ``WorkspaceStore``;
- it emits intent signals (``item_activated``, ``preview_requested``,
  ``run_renamed``) instead of opening dialogs itself;
- it does not set a QStyle or a font, so YRTools owns the appearance.

V1 is read-only: preview/open and jump-to-item are wired; inline parameter
editing and import actions are intentionally left to a later version.
"""

from __future__ import annotations

import json
from datetime import datetime
from typing import Any, Dict, List, Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
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
    "tree": "Open in IcyTree",
    "alignment": "Open sequences",
    "model": "Open model",
    "distance": "Open distance matrix",
    "chain": "Open chain",
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
    return _text(value)


class InspectorPanel(QWidget):
    """Selection-driven provenance/preview panel."""

    item_activated = pyqtSignal(str)
    dataset_activated = pyqtSignal(str)
    preview_requested = pyqtSignal(str, str)  # item_id, kind
    run_renamed = pyqtSignal(str, str)        # run_id, new_name

    def __init__(self, manager=None, store=None, parent=None):
        super().__init__(parent)
        self.manager = manager
        self.store = store
        self.current_item_id: Optional[str] = None
        self.current_dataset_id: Optional[str] = None

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
        elif self.current_dataset_id:
            self.show_dataset(self.current_dataset_id)
        else:
            self.show_empty()

    def show_empty(self):
        self.current_item_id = None
        self.current_dataset_id = None
        self.title_label.setText("Inspector")
        self.subtitle_label.setText("Select a data node or dataset")
        self.preview_group.setEnabled(False)
        self.preview_hint.setText("No selection.")
        self.preview_button.setText("Open preview")
        self.show_button.setEnabled(False)
        self.meta_table.setRowCount(0)
        self.workflow_tree.clear()
        self.run_group.setVisible(False)

    def show_item(self, item_id: str):
        item = self._get_item(item_id)
        if item is None:
            self.show_empty()
            return

        self.current_item_id = item_id
        self.current_dataset_id = _attr(item, "dataset_id")
        item_type = _attr(item, "item_type", "type", default="sequence")
        name = _attr(item, "loci_name", "name", default=item_id)
        dataset_name = self._dataset_name(self.current_dataset_id)

        self.title_label.setText(name)
        self.subtitle_label.setText(
            "%s · %s" % (TYPE_LABELS.get(item_type, item_type), dataset_name or "no dataset")
        )

        kind = self._preview_kind(item_type)
        self.preview_group.setEnabled(True)
        self.preview_button.setText(PREVIEW_LABELS.get(kind, "Open preview"))
        self.preview_button.setEnabled(kind is not None)
        self.show_button.setEnabled(True)
        self.preview_hint.setText(self._preview_hint(item, item_type))

        self._fill_item_metadata(item, item_type, dataset_name)
        self._fill_workflow(self.current_dataset_id, item_id)
        self._fill_runs(item)

    def show_dataset(self, dataset_id: str):
        dataset = self._get_dataset(dataset_id)
        if dataset is None:
            self.show_empty()
            return

        self.current_item_id = None
        self.current_dataset_id = dataset_id
        name = _attr(dataset, "name", default=dataset_id)
        self.title_label.setText(name)
        self.subtitle_label.setText("Dataset")

        self.preview_group.setEnabled(False)
        self.preview_hint.setText("Select an item to preview it.")
        self.preview_button.setEnabled(False)
        self.show_button.setEnabled(False)

        self._fill_dataset_metadata(dataset)
        self._fill_workflow(dataset_id, None)
        self.run_group.setVisible(False)

    def refresh(self):
        if self.current_item_id:
            self.show_item(self.current_item_id)
        elif self.current_dataset_id:
            self.show_dataset(self.current_dataset_id)
        else:
            self.show_empty()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        self.setMinimumWidth(300)
        self.setMaximumWidth(560)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

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

        # Preview -----------------------------------------------------
        self.preview_group = QGroupBox("Preview")
        preview_layout = QVBoxLayout(self.preview_group)
        self.preview_hint = QLabel("No selection.")
        self.preview_hint.setWordWrap(True)
        self.preview_hint.setStyleSheet("color: palette(mid);")
        preview_layout.addWidget(self.preview_hint)

        preview_buttons = QHBoxLayout()
        self.preview_button = QPushButton("Open preview")
        self.preview_button.clicked.connect(self._on_preview)
        self.show_button = QPushButton("Show in workspace")
        self.show_button.clicked.connect(self._on_show)
        preview_buttons.addWidget(self.preview_button)
        preview_buttons.addWidget(self.show_button)
        preview_layout.addLayout(preview_buttons)
        layout.addWidget(self.preview_group)

        # Metadata ----------------------------------------------------
        metadata_group = QGroupBox("Metadata")
        metadata_layout = QVBoxLayout(metadata_group)
        self.meta_table = QTableWidget(0, 2)
        self.meta_table.setHorizontalHeaderLabels(["Field", "Value"])
        self.meta_table.verticalHeader().setVisible(False)
        self.meta_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.meta_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.meta_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.meta_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.meta_table.setMinimumHeight(180)
        self.meta_table.setMaximumHeight(260)
        metadata_layout.addWidget(self.meta_table)
        layout.addWidget(metadata_group)

        # Workflow ----------------------------------------------------
        workflow_group = QGroupBox("Workflow")
        workflow_layout = QVBoxLayout(workflow_group)
        self.workflow_tree = QTreeWidget()
        self.workflow_tree.setColumnCount(2)
        self.workflow_tree.setHeaderLabels(["Data / method", "Type"])
        self.workflow_tree.setAlternatingRowColors(True)
        self.workflow_tree.setUniformRowHeights(True)
        self.workflow_tree.header().setStretchLastSection(False)
        self.workflow_tree.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.workflow_tree.header().setSectionResizeMode(1, QHeaderView.Fixed)
        self.workflow_tree.setColumnWidth(1, 64)
        self.workflow_tree.itemClicked.connect(self._on_tree_clicked)
        self.workflow_tree.itemDoubleClicked.connect(self._on_tree_double_clicked)
        workflow_layout.addWidget(self.workflow_tree)
        layout.addWidget(workflow_group, 1)

        # Runs --------------------------------------------------------
        self.run_group = QGroupBox("Runs")
        run_layout = QVBoxLayout(self.run_group)
        self.run_combo = QComboBox()
        self.run_combo.currentIndexChanged.connect(self._on_run_combo)
        run_layout.addWidget(self.run_combo)

        rename_row = QHBoxLayout()
        self.run_name_edit = QLineEdit()
        self.run_name_edit.setPlaceholderText("Run name")
        self.rename_button = QPushButton("Rename")
        self.rename_button.clicked.connect(self._on_rename_run)
        rename_row.addWidget(self.run_name_edit, 1)
        rename_row.addWidget(self.rename_button)
        run_layout.addLayout(rename_row)
        layout.addWidget(self.run_group)

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

    def _stream_for_item(self, item):
        run = self._run_for_item(item)
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

    def _preview_hint(self, item, item_type: str) -> str:
        if item_type in ("sequence", "alignment"):
            count = _attr(item, "sequence_count", default=0)
            length = _attr(item, "length", default=0)
            return "%s sequences · %s sites" % (count, length)
        if item_type in ("phylogeny", "treeset"):
            return "Newick tree"
        if item_type == "model":
            return "Model selection result"
        if item_type == "distance":
            return "Pairwise distance matrix"
        if item_type == "chain":
            return "MCMC chain output"
        return TYPE_LABELS.get(item_type, item_type)

    def _fill_item_metadata(self, item, item_type: str, dataset_name: str):
        rows = [
            ("Type", TYPE_LABELS.get(item_type, item_type)),
            ("Dataset", dataset_name),
            ("File", _attr(item, "file_path", "file", default="")),
            ("Format", _attr(item, "format", default="")),
            ("Length", _attr(item, "length", default="")),
            ("Sequences", _attr(item, "sequence_count", default="")),
            ("Aligned", "yes" if _attr(item, "is_aligned", default=False) else "no"),
            ("Import", _attr(item, "import_id", default="")),
            ("Run", _attr(item, "run_id", default="")),
            ("Parent", _attr(item, "parent_id", default="")),
            ("Created", _display_time(_attr(item, "created_at", default=""))),
        ]
        run = self._run_for_item(item)
        if run is not None:
            rows.extend([
                ("Run name", run.name),
                ("Method", run.method),
                ("Tool", run.tool + (" " + run.tool_version if run.tool_version else "")),
                ("Status", run.status),
            ])
            stream = self._stream_for_item(item)
            if stream is not None:
                rows.append(("Stream", stream.name))
        metadata = _attr(item, "metadata", default={}) or {}
        for key, value in metadata.items():
            if key in ("source_path", "sha256"):
                rows.append((key, _text(value)))
        self._set_table_rows(rows)

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
        self._set_table_rows(rows)

    def _set_table_rows(self, rows):
        self.meta_table.setRowCount(0)
        for field, value in rows:
            if value in (None, ""):
                continue
            row = self.meta_table.rowCount()
            self.meta_table.insertRow(row)
            field_item = QTableWidgetItem(str(field))
            value_item = QTableWidgetItem(_text(value))
            value_item.setToolTip(_text(value, limit=2000))
            self.meta_table.setItem(row, 0, field_item)
            self.meta_table.setItem(row, 1, value_item)

    # ------------------------------------------------------------------
    # workflow tree
    # ------------------------------------------------------------------
    def _fill_workflow(self, dataset_id: str, current_item_id: Optional[str]):
        self.workflow_tree.clear()
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
        root.setExpanded(True)
        self.workflow_tree.addTopLevelItem(root)

        item_nodes = {}
        for item in items:
            item_id = _attr(item, "id", default="")
            node = QTreeWidgetItem([
                _attr(item, "loci_name", "name", default=item_id),
                TYPE_SHORT.get(_attr(item, "item_type", "type", default=""), ""),
            ])
            node.setData(0, Qt.UserRole, ("item", item_id))
            node.setToolTip(0, _text(_attr(item, "file_path", "file", default=""), 500))
            item_nodes[item_id] = node

        store = self._get_store()
        import_nodes = {}
        run_nodes = {}
        if store is not None:
            for record in store.index.imports.values():
                if record.dataset_id != dataset_id:
                    continue
                node = QTreeWidgetItem([
                    "Import: %s" % record.importer,
                    "IMP",
                ])
                node.setData(0, Qt.UserRole, ("import", record.id))
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
                node.setToolTip(0, "Status: %s" % record.status)
                run_nodes[record.id] = node

        # Attach items under their producing run/import, then parent.
        for item in items:
            item_id = _attr(item, "id", default="")
            node = item_nodes[item_id]
            run_id = _attr(item, "run_id", default=None)
            import_id = _attr(item, "import_id", default=None)
            parent_id = _attr(item, "parent_id", default=None)
            if run_id and run_id in run_nodes:
                run_nodes[run_id].addChild(node)
            elif import_id and import_id in import_nodes:
                import_nodes[import_id].addChild(node)
            elif parent_id and parent_id in item_nodes:
                item_nodes[parent_id].addChild(node)
            else:
                root.addChild(node)

        # Attach each run under its first input, otherwise under the dataset.
        for record in (store.index.runs.values() if store is not None else []):
            if record.id not in run_nodes:
                continue
            node = run_nodes[record.id]
            anchor = None
            for input_id in record.input_ids:
                if input_id in item_nodes:
                    anchor = item_nodes[input_id]
                    break
            if anchor is not None:
                anchor.addChild(node)
            else:
                root.addChild(node)

        # Expand the useful top levels, but not hundreds of leaves.
        for record in import_nodes.values():
            record.setExpanded(record.childCount() <= 30)
        for record in run_nodes.values():
            record.setExpanded(record.childCount() <= 30)

        if current_item_id and current_item_id in item_nodes:
            current = item_nodes[current_item_id]
            self.workflow_tree.setCurrentItem(current)
            self.workflow_tree.scrollToItem(current)
            self._bold_ancestors(current)
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
            self.item_activated.emit(record_id)
        elif kind == "dataset":
            self.dataset_activated.emit(record_id)

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
    # preview / jump
    # ------------------------------------------------------------------
    def _on_preview(self):
        if not self.current_item_id:
            return
        self.preview_requested.emit(self.current_item_id, self._kind_for_item(self.current_item_id))

    def _on_show(self):
        if self.current_item_id:
            self.item_activated.emit(self.current_item_id)

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
            # Root/imported data has no analysis run; do not treat every
            # same-type root item as an alternative run.
            siblings = []
        else:
            siblings = [other for other in siblings if self._run_for_item(other) is not None]
            if not any(_attr(other, "id", default="") == _attr(item, "id", default="") for other in siblings):
                siblings.insert(0, item)
        self.run_group.setVisible(run is not None)

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
        if not self.current_item_id:
            return
        item = self._get_item(self.current_item_id)
        run = self._run_for_item(item) if item is not None else None
        store = self._get_store()
        if run is None or store is None:
            return
        name = self.run_name_edit.text().strip()
        if not name or name == run.name:
            return
        store.rename_run(run.id, name)
        self.run_renamed.emit(run.id, name)
        self.refresh()
