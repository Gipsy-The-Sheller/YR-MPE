"""Right-side Selected basket.

The Details panel answers "what is this node?".  The Selected panel answers
"what am I about to send to a tool?" and can hold objects from several chains
and types (for example a gene-tree set plus a species tree for gCF/sCF).
"""

from __future__ import annotations

from typing import Iterable, List, Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QBrush, QColor
from PyQt5.QtWidgets import (
    QAbstractItemView, QGroupBox, QHBoxLayout, QLabel, QPushButton,
    QToolButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from YR_MPE.factories import ResourceFactory

from .dataset_item_button import ITEM_MIME_TYPE
from .dataset_models import (
    ITEM_TYPE_ALIGNMENT, ITEM_TYPE_CHAIN, ITEM_TYPE_DISTANCE,
    ITEM_TYPE_MODEL, ITEM_TYPE_PHYLOGENY, ITEM_TYPE_SEQUENCE,
    ITEM_TYPE_TABLE, ITEM_TYPE_TREESET,
)

TYPE_LABELS = {
    ITEM_TYPE_SEQUENCE: "Sequence",
    ITEM_TYPE_ALIGNMENT: "Alignment",
    ITEM_TYPE_MODEL: "Model",
    ITEM_TYPE_DISTANCE: "Distance",
    ITEM_TYPE_PHYLOGENY: "Tree",
    ITEM_TYPE_TREESET: "Tree set",
    ITEM_TYPE_CHAIN: "Chain",
    ITEM_TYPE_TABLE: "Table",
}


def _attr(obj, *names, default=""):
    for name in names:
        if hasattr(obj, name):
            value = getattr(obj, name)
            if value:
                return value
    return default


class SelectedDropButton(QToolButton):
    """Toolbar button that accepts workspace objects dragged onto it."""

    item_dropped = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)

    @staticmethod
    def _item_id(mime):
        if not mime.hasFormat(ITEM_MIME_TYPE):
            return ""
        return bytes(mime.data(ITEM_MIME_TYPE)).decode("utf-8")

    def dragEnterEvent(self, event):
        if self._item_id(event.mimeData()):
            event.setDropAction(Qt.CopyAction)
            event.accept()
            self.setDown(True)
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if self._item_id(event.mimeData()):
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dragLeaveEvent(self, event):
        self.setDown(False)
        super().dragLeaveEvent(event)

    def dropEvent(self, event):
        item_id = self._item_id(event.mimeData())
        self.setDown(False)
        if not item_id:
            event.ignore()
            return
        event.setDropAction(Qt.CopyAction)
        event.accept()
        self.item_dropped.emit(item_id)


class SelectedPanel(QWidget):
    """Multi-object basket and type-aware 'send to' recommendations."""

    send_requested = pyqtSignal(str, object)  # tool_id, item_ids
    item_activated = pyqtSignal(str)
    changed = pyqtSignal()

    def __init__(self, manager=None, store=None, parent=None, task_queue=None):
        super().__init__(parent)
        self.manager = manager
        self.store = store
        self.task_queue = task_queue
        self._job_inspectors = {}
        self._items: List[object] = []
        self.resource_factory = ResourceFactory()
        self.setObjectName("selectedPanel")
        self._build_ui()
        self.set_task_queue(task_queue)
        self.refresh()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)

        self.title_label = QLabel("Selected & Queue")
        self.title_label.setStyleSheet("font-size: 17px; font-weight: 600;")
        layout.addWidget(self.title_label)
        self.subtitle_label = QLabel(
            "Selected objects form the next input set. Running and finished "
            "analyses stay here so long jobs can continue after their window "
            "is closed."
        )
        self.subtitle_label.setWordWrap(True)
        self.subtitle_label.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.subtitle_label)

        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setRootIsDecorated(True)
        self.tree.setIndentation(14)
        self.tree.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tree.setSelectionMode(QAbstractItemView.SingleSelection)
        self.tree.itemDoubleClicked.connect(self._on_double_click)
        layout.addWidget(self.tree, 1)

        row = QHBoxLayout()
        self.remove_button = QPushButton("Remove")
        self.remove_button.clicked.connect(self._on_remove)
        self.clear_button = QPushButton("Clear all")
        self.clear_button.clicked.connect(self._on_clear)
        row.addWidget(self.remove_button)
        row.addWidget(self.clear_button)
        row.addStretch(1)
        layout.addLayout(row)

        self.send_group = QGroupBox("Send to")
        self.send_layout = QVBoxLayout(self.send_group)
        self.send_layout.setContentsMargins(6, 6, 6, 6)
        self.send_layout.setSpacing(4)
        self.send_tree = QTreeWidget()
        self.send_tree.setHeaderHidden(True)
        self.send_tree.setRootIsDecorated(True)
        self.send_tree.setIndentation(14)
        self.send_tree.setMinimumHeight(150)
        self.send_tree.setMaximumHeight(320)
        self.send_tree.itemClicked.connect(self._on_send_tree_clicked)
        self.send_layout.addWidget(self.send_tree)
        layout.addWidget(self.send_group)

        self.jobs_group = QWidget()
        jobs_layout = QVBoxLayout(self.jobs_group)
        jobs_layout.setContentsMargins(0, 0, 0, 0)
        jobs_layout.setSpacing(4)

        self.queue_title_label = QLabel("Queue")
        title_font = self.queue_title_label.font()
        title_font.setBold(True)
        self.queue_title_label.setFont(title_font)
        jobs_layout.addWidget(self.queue_title_label)

        self.jobs_tree = QTreeWidget()
        self.jobs_tree.setColumnCount(2)
        self.jobs_tree.setHeaderLabels(["Task", "Status"])
        self.jobs_tree.setRootIsDecorated(False)
        self.jobs_tree.setIndentation(0)
        self.jobs_tree.setAlternatingRowColors(True)
        self.jobs_tree.setSelectionMode(QAbstractItemView.SingleSelection)
        self.jobs_tree.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.jobs_tree.setMinimumHeight(96)
        self.jobs_tree.setMaximumHeight(180)
        self.jobs_tree.header().setStretchLastSection(False)
        self.jobs_tree.header().resizeSection(0, 240)
        self.jobs_tree.header().resizeSection(1, 82)
        self.jobs_tree.itemDoubleClicked.connect(self._on_job_double_click)
        self.jobs_tree.itemSelectionChanged.connect(self._update_job_buttons)
        jobs_layout.addWidget(self.jobs_tree)

        self.jobs_empty_label = QLabel("No queued analyses.")
        self.jobs_empty_label.setStyleSheet("color: palette(mid);")
        jobs_layout.addWidget(self.jobs_empty_label)

        self.queue_notice_label = QLabel("")
        self.queue_notice_label.setWordWrap(True)
        self.queue_notice_label.setVisible(False)
        jobs_layout.addWidget(self.queue_notice_label)

        job_buttons = QHBoxLayout()
        self.inspect_job_button = QPushButton("Inspect")
        self.inspect_job_button.setToolTip("Inspect this task and preview its result")
        self.inspect_job_button.setEnabled(False)
        self.inspect_job_button.clicked.connect(self._inspect_current_job)
        self.cancel_job_button = QPushButton("Cancel")
        self.cancel_job_button.setEnabled(False)
        self.cancel_job_button.clicked.connect(self._cancel_current_job)
        job_buttons.addWidget(self.inspect_job_button)
        job_buttons.addWidget(self.cancel_job_button)
        job_buttons.addStretch(1)
        jobs_layout.addLayout(job_buttons)
        layout.addWidget(self.jobs_group)

        self.empty_label = QLabel(
            "No objects selected. Right-click an icon → Add to Selected, "
            "or drag it onto the Selected button."
        )
        self.empty_label.setWordWrap(True)
        self.empty_label.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.empty_label)

        self._body_widgets = [
            self.tree, self.remove_button, self.clear_button,
            self.send_group, self.jobs_group, self.empty_label,
        ]
        self._collapsed = False

    def set_collapsed(self, collapsed: bool):
        collapsed = bool(collapsed)
        if collapsed == self._collapsed:
            return
        if collapsed:
            for widget in self._body_widgets:
                widget.setVisible(False)
        else:
            self.tree.setVisible(True)
            self.remove_button.setVisible(True)
            self.clear_button.setVisible(True)
            jobs = self._jobs()
            self.empty_label.setVisible(not self._items and not jobs)
            self.send_group.setVisible(self.send_tree.topLevelItemCount() > 0)
            self.jobs_group.setVisible(True)
            self.queue_notice_label.setVisible(
                bool(self.queue_notice_label.text())
            )
        self._collapsed = collapsed

    def collapsed_height(self) -> int:
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        return max(
            48,
            margins.top() + self.title_label.sizeHint().height() + spacing
            + self.subtitle_label.sizeHint().height() + margins.bottom(),
        )

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def set_sources(self, manager=None, store=None):
        self.manager = manager if manager is not None else self.manager
        self.store = store if store is not None else self.store
        self.refresh()

    def set_task_queue(self, queue):
        self.task_queue = queue
        if queue is not None:
            for signal_name in ("jobs_changed", "job_updated", "job_finished"):
                signal = getattr(queue, signal_name, None)
                if signal is None:
                    continue
                try:
                    signal.disconnect(self.refresh_jobs)
                except (TypeError, RuntimeError):
                    pass
                signal.connect(self.refresh_jobs)
            finished = getattr(queue, "job_finished", None)
            if finished is not None:
                try:
                    finished.disconnect(self._on_job_finished_notice)
                except (TypeError, RuntimeError):
                    pass
                finished.connect(self._on_job_finished_notice)
            added = getattr(queue, "job_added", None)
            if added is not None:
                try:
                    added.disconnect(self._clear_job_notice)
                except (TypeError, RuntimeError):
                    pass
                added.connect(self._clear_job_notice)
        self.refresh_jobs()

    def _clear_job_notice(self, *args):
        self.queue_notice_label.clear()
        self.queue_notice_label.setVisible(False)

    def _on_job_finished_notice(self, job_id):
        job = self._job_by_id(job_id)
        if job is None:
            return
        status = getattr(job, "status_label", "Finished")
        title = getattr(job, "title", "Analysis")
        elapsed = ""
        if callable(getattr(job, "elapsed_text", None)):
            try:
                elapsed = " · %s" % job.elapsed_text()
            except Exception:
                elapsed = ""
        state = getattr(job, "state", "")
        color = {
            "finished": "#15803d",
            "error": "#b91c1c",
            "cancelled": "#6b7280",
        }.get(state, "#b91c1c")
        self.queue_notice_label.setText("%s: %s%s" % (status, title, elapsed))
        self.queue_notice_label.setStyleSheet(
            "color: %s; font-weight: 600;" % color
        )
        self.queue_notice_label.setVisible(not self._collapsed)
        self.queue_notice_label.setToolTip(
            "Select the task and use Inspect / Preview to open its result."
        )

    def _job_by_id(self, job_id):
        queue = self.task_queue
        getter = getattr(queue, "get_job", None) if queue is not None else None
        return getter(job_id) if callable(getter) else None

    def _jobs(self):
        queue = self.task_queue
        if queue is None or not callable(getattr(queue, "jobs", None)):
            return []
        try:
            return list(queue.jobs())
        except Exception:
            return []

    TYPE_GROUPS = [
        (ITEM_TYPE_ALIGNMENT, "Alignments", "file/sequence.svg"),
        (ITEM_TYPE_SEQUENCE, "Sequences", "file/sequence.svg"),
        (ITEM_TYPE_PHYLOGENY, "Trees", "file/phylogeny.svg"),
        (ITEM_TYPE_TREESET, "Tree sets", "file/treeset.svg"),
        (ITEM_TYPE_MODEL, "Models", "file/model.svg"),
        (ITEM_TYPE_DISTANCE, "Distances", "file/distance.svg"),
        (ITEM_TYPE_CHAIN, "Chains", "file/chain.svg"),
        (ITEM_TYPE_TABLE, "Tables", "file/table.svg"),
    ]

    def refresh(self):
        self._items = list(self.manager.get_basket_items()) if self.manager else []
        self.tree.clear()

        grouped = {}
        for item in self._items:
            grouped.setdefault(
                _attr(item, "item_type", "type", default=""), []
            ).append(item)

        ordered_types = [entry[0] for entry in self.TYPE_GROUPS]
        ordered_types += [key for key in grouped if key not in ordered_types]
        specs = {key: (label, icon) for key, label, icon in self.TYPE_GROUPS}

        for item_type in ordered_types:
            items = grouped.get(item_type)
            if not items:
                continue
            label, icon_path = specs.get(
                item_type, (TYPE_LABELS.get(item_type, item_type) + "s", "file/dataset.svg")
            )
            parent = QTreeWidgetItem(["%s (%d)" % (label, len(items))])
            parent.setFlags(parent.flags() & ~Qt.ItemIsSelectable)
            font = parent.font(0)
            font.setBold(True)
            parent.setFont(0, font)
            parent.setIcon(0, self.resource_factory.get_icon(icon_path))
            self.tree.addTopLevelItem(parent)

            for item in sorted(
                items,
                key=lambda record: _attr(record, "loci_name", "name", default=record.id).lower(),
            ):
                name = _attr(item, "loci_name", "name", default=item.id)
                dataset_id = getattr(item, "dataset_id", "") or ""
                dataset_record = (self.manager.datasets.get(dataset_id)
                                  if self.manager else None)
                dataset = (getattr(dataset_record, "name", "")
                           if dataset_record else dataset_id)
                child = QTreeWidgetItem([name])
                child.setData(0, Qt.UserRole, item.id)
                child.setIcon(0, self.resource_factory.get_icon(icon_path))
                child.setToolTip(0, "%s · %s" % (name, dataset))
                parent.addChild(child)
        self.tree.expandAll()
        self._update_recommendations([_attr(item, "item_type", "type")
                                      for item in self._items])
        self.refresh_jobs()
        jobs = self._jobs()
        self.empty_label.setVisible(
            not self._items and not jobs and not self._collapsed
        )
        self.remove_button.setEnabled(bool(self._items))
        self.clear_button.setEnabled(bool(self._items))

    def refresh_jobs(self):
        """Render the queue without changing the Selected basket."""
        if not hasattr(self, "jobs_tree"):
            return
        current = self._current_job_id()
        jobs = self._jobs()
        self.jobs_tree.clear()
        self.queue_title_label.setText(
            "Queue (%d)" % len(jobs) if jobs else "Queue"
        )
        colors = {
            "queued": "#6b7280",
            "running": "#2563eb",
            "finished": "#15803d",
            "error": "#b91c1c",
            "cancelled": "#9ca3af",
        }
        for job in jobs:
            status = getattr(job, "status_label", "") or str(
                getattr(job, "state", "")
            ).title()
            title = getattr(job, "title", "") or getattr(
                job, "tool_name", "Analysis"
            )
            row = QTreeWidgetItem([str(title), str(status)])
            row.setData(0, Qt.UserRole, getattr(job, "job_id", ""))
            color = QColor(colors.get(getattr(job, "state", ""), "#6b7280"))
            row.setForeground(1, QBrush(color))
            font = row.font(1)
            font.setBold(True)
            row.setFont(1, font)
            tooltip = []
            command = getattr(job, "command", None)
            if command:
                tooltip.append(" ".join(str(part) for part in command))
            progress = getattr(job, "progress", "")
            if progress:
                tooltip.append(progress)
            error = getattr(job, "error_message", "")
            if error:
                tooltip.append(error)
            if tooltip:
                row.setToolTip(0, "\n".join(tooltip))
                row.setToolTip(1, "\n".join(tooltip))
            self.jobs_tree.addTopLevelItem(row)
        if current:
            for index in range(self.jobs_tree.topLevelItemCount()):
                row = self.jobs_tree.topLevelItem(index)
                if row.data(0, Qt.UserRole) == current:
                    self.jobs_tree.setCurrentItem(row)
                    break
        self.jobs_empty_label.setVisible(not jobs)
        if not jobs:
            self.queue_notice_label.clear()
            self.queue_notice_label.setVisible(False)
        else:
            self.queue_notice_label.setVisible(
                bool(self.queue_notice_label.text()) and not self._collapsed
            )
        self.jobs_group.setVisible(not self._collapsed)
        self._update_job_buttons()

    def _current_job_id(self) -> Optional[str]:
        item = self.jobs_tree.currentItem() if hasattr(self, "jobs_tree") else None
        if item is None:
            return None
        job_id = item.data(0, Qt.UserRole)
        return str(job_id) if job_id else None

    def _current_job(self):
        job_id = self._current_job_id()
        if not job_id or self.task_queue is None:
            return None
        getter = getattr(self.task_queue, "get_job", None)
        return getter(job_id) if callable(getter) else None

    def _update_job_buttons(self):
        job = self._current_job()
        self.inspect_job_button.setEnabled(job is not None)
        self.cancel_job_button.setEnabled(
            bool(job is not None and getattr(job, "can_cancel", False))
        )

    def _on_job_double_click(self, item, column):
        if item is not None and item.data(0, Qt.UserRole):
            self._inspect_job(str(item.data(0, Qt.UserRole)))

    def _inspect_current_job(self):
        job_id = self._current_job_id()
        if job_id:
            self._inspect_job(job_id)

    def _inspect_job(self, job_id):
        if self.task_queue is None:
            return
        dialog = self._job_inspectors.get(job_id)
        try:
            if dialog is not None and dialog.isVisible():
                dialog.raise_()
                dialog.activateWindow()
                return
        except RuntimeError:
            self._job_inspectors.pop(job_id, None)
        from .task_job_inspector import TaskJobInspector
        dialog = TaskJobInspector(self.task_queue, job_id, parent=self)
        dialog.setAttribute(Qt.WA_DeleteOnClose, True)
        dialog.destroyed.connect(
            lambda *args, key=job_id: self._job_inspectors.pop(key, None)
        )
        self._job_inspectors[job_id] = dialog
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _cancel_current_job(self):
        job = self._current_job()
        if job is None:
            return
        if self.task_queue is not None and callable(
            getattr(self.task_queue, "cancel", None)
        ):
            self.task_queue.cancel(job.job_id)
        self.refresh_jobs()

    def recommended_operations(self):
        """Return send-to operations for the current basket as data.

        The recommendation rules are driven by the item types in Selected:

        * tree / tree set -> ASTRAL-IV, wASTRAL, RF
        * tree + tree set -> TreeRecs, RF
        * tree -> LSD2
        * sequence (FASTA) / alignment (PHYLIP, NEXUS) -> aligners and
          alignment-based tree, model, distance and selection tools
        * distance matrix -> DecentTree

        Each operation is always present with ``available`` so an agent can
        explain why a tool cannot consume the current selection. The UI uses
        the same list and hides unavailable entries.
        """
        counts = {}
        for item in self._items:
            item_type = _attr(item, "item_type", "type")
            counts[item_type] = counts.get(item_type, 0) + 1

        def count(item_type):
            return counts.get(item_type, 0)

        has_sequence = count(ITEM_TYPE_SEQUENCE) > 0
        has_alignment = count(ITEM_TYPE_ALIGNMENT) > 0
        has_distance = count(ITEM_TYPE_DISTANCE) > 0
        has_phylogeny = count(ITEM_TYPE_PHYLOGENY) > 0
        has_treeset = count(ITEM_TYPE_TREESET) > 0
        has_chain = count(ITEM_TYPE_CHAIN) > 0
        multiple_trees = has_treeset or count(ITEM_TYPE_PHYLOGENY) >= 2
        sequence_data = has_sequence or has_alignment

        groups = [
            ("Species Tree", [
                ("astral4", "ASTRAL-IV", "software/astral.svg",
                 multiple_trees, ["tree set or >=2 trees"]),
                ("wastral", "wASTRAL", "software/astral.svg",
                 multiple_trees, ["tree set or >=2 trees"]),
                ("treerecs", "TreeRecs 1.2", "coalescent.svg",
                 has_phylogeny and has_treeset, ["tree + tree set"]),
                ("gcf", "Gene concordance (gCF)", "coalescent.svg",
                 has_phylogeny and multiple_trees, ["tree + tree set"]),
                ("scf", "Site concordance (sCF)", "coalescent.svg",
                 has_phylogeny and has_alignment, ["tree + alignment"]),
                ("site_diagnostics", "Site diagnostics (sCF + CASTER)",
                 "software/caster.svg", has_phylogeny and has_alignment,
                 ["tree + alignment"]),
                ("msc_diagnostics", "MSC diagnostics (gCF + ASTER)",
                 "software/astral.svg", has_phylogeny and multiple_trees,
                 ["tree + tree set"]),
                ("rf", "Robinson-Foulds comparison",
                 "software/robinsonfoluds.svg", multiple_trees,
                 ["tree set or >=2 trees"]),
            ]),
            ("Molecular Clock", [
                ("lsd2", "LSD2", "software/lsd2.svg",
                 has_phylogeny, ["tree"]),
            ]),
            ("Alignment", [
                ("mafft", "MAFFT", "software/mafft.svg",
                 sequence_data, ["sequence or alignment"]),
                ("muscle5", "Muscle5", "software/muscle.svg",
                 sequence_data, ["sequence or alignment"]),
                ("clustal_omega", "Clustal Omega", "software/clustalo.svg",
                 sequence_data, ["sequence or alignment"]),
                ("macse", "MACSE2", "software/macse.svg",
                 sequence_data, ["sequence or alignment"]),
            ]),
            ("Trimming", [
                ("trimal", "TrimAl", "software/trimal.svg",
                 has_alignment, ["alignment"]),
                ("gblocks", "GBlocks", "software/gblocks.svg",
                 has_alignment, ["alignment"]),
            ]),
            ("Tree Inference", [
                ("iqtree", "IQ-TREE 3", "software/iqtree.svg",
                 has_alignment, ["alignment"]),
                ("mrbayes", "MrBayes", "software/mrbayes.svg",
                 has_alignment, ["alignment"]),
                ("raxmlng", "RAxML-NG", "software/raxml.svg",
                 has_alignment, ["alignment"]),
                ("mpboot", "MPBoot", "mp.svg",
                 has_alignment, ["alignment"]),
                ("phylobayes", "PhyloBayes-MPI", "software/phylobayes.svg",
                 has_alignment, ["alignment"]),
                ("caster_site", "CASTER-site", "software/caster.svg",
                 has_alignment, ["alignment"]),
                ("best", "BEST / MrBayes MSC", "software/mrbayes.svg",
                 has_alignment, ["alignment"]),
            ]),
            ("Model Selection", [
                ("modelfinder", "ModelFinder", "find_model.svg",
                 has_alignment, ["alignment"]),
                ("modeltest_ng", "ModelTest-NG", "find_model.svg",
                 has_alignment, ["alignment"]),
                ("model_parameter", "Model Parameter Estimation",
                 "find_model.svg", has_alignment, ["alignment"]),
            ]),
            ("Distances & Tests", [
                ("ml_distance", "ML pairwise distance", "dist.svg",
                 has_alignment, ["alignment"]),
                ("uncorrected_distance", "Uncorrected p-distance",
                 "dist.svg", has_alignment, ["alignment"]),
                ("codon_ztest", "Codon-based Z-test", "find_model.svg",
                 has_alignment, ["alignment"]),
                ("saturation", "Substitution saturation (BaCoCa)",
                 "saturation.svg", has_alignment, ["alignment"]),
            ]),
            ("Distance Tree", [
                ("decenttree", "DecentTree", "bionj.svg",
                 has_distance, ["distance matrix"]),
            ]),
            ("Diagnostics", [
                ("minitracer", "MiniTracer", "software/minitracer.svg",
                 has_chain, ["chain"]),
            ]),
        ]
        return [
            {
                "group": group_name,
                "operations": [
                    {
                        "id": tool_id,
                        "label": label,
                        "icon": icon_path,
                        "available": bool(available),
                        "requires": requires,
                    }
                    for tool_id, label, icon_path, available, requires in entries
                ],
            }
            for group_name, entries in groups
        ]

    def _update_recommendations(self, types: Iterable[str] = ()):
        self.send_tree.clear()
        if not self._items:
            self.send_group.setVisible(False)
            self._update_body_visibility()
            return
        visible_groups = 0
        for group in self.recommended_operations():
            entries = [entry for entry in group["operations"]
                       if entry["available"]]
            if not entries:
                continue
            visible_groups += 1
            parent = QTreeWidgetItem([group["group"]])
            parent.setFlags(parent.flags() & ~Qt.ItemIsSelectable)
            font = parent.font(0)
            font.setBold(True)
            parent.setFont(0, font)
            self.send_tree.addTopLevelItem(parent)
            for entry in entries:
                child = QTreeWidgetItem([entry["label"]])
                child.setData(0, Qt.UserRole, entry["id"])
                child.setIcon(0, self.resource_factory.get_icon(entry["icon"]))
                child.setToolTip(
                    0, "Send the Selected objects to %s" % entry["label"]
                )
                parent.addChild(child)
        self.send_tree.expandAll()
        self.send_group.setVisible(bool(visible_groups) and not self._collapsed)
        self._update_body_visibility()

    def _update_body_visibility(self):
        if self._collapsed:
            return
        self.send_group.setVisible(self.send_tree.topLevelItemCount() > 0)

    def _on_send_tree_clicked(self, item, column):
        tool_id = item.data(0, Qt.UserRole)
        if tool_id:
            self._on_send(tool_id)

    # ------------------------------------------------------------------
    # actions
    # ------------------------------------------------------------------
    def _current_item_id(self) -> Optional[str]:
        item = self.tree.currentItem()
        if item is None:
            return None
        item_id = item.data(0, Qt.UserRole)
        return str(item_id) if item_id else None

    def _on_double_click(self, item, column):
        item_id = item.data(0, Qt.UserRole) if item is not None else None
        if item_id:
            self.item_activated.emit(str(item_id))

    def _on_remove(self):
        item_id = self._current_item_id()
        if item_id and self.manager:
            self.manager.remove_from_basket(item_id)
            self.refresh()
            self.changed.emit()

    def _on_clear(self):
        if self.manager:
            self.manager.clear_basket()
            self.refresh()
            self.changed.emit()

    def _on_send(self, tool_id: str):
        self.send_requested.emit(tool_id, [item.id for item in self._items])
