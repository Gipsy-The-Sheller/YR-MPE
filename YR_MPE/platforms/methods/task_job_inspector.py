"""Read-only inspector for one task-queue job.

The inspector is deliberately separate from the plugin window.  A running job
may have no visible window at all, and a finished job should still expose its
status, command, console tail and artifacts before the user reopens the tool.
"""

from __future__ import annotations

import os
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QDialog, QFormLayout, QHBoxLayout, QLabel, QPlainTextEdit, QPushButton,
    QTabWidget, QVBoxLayout, QWidget,
)

from YR_MPE.runtime.task_queue import (
    CANCELLED, ERROR, FINISHED, QUEUED, RUNNING,
)


STATE_COLORS = {
    QUEUED: "#6b7280",
    RUNNING: "#2563eb",
    FINISHED: "#15803d",
    ERROR: "#b91c1c",
    CANCELLED: "#9ca3af",
}


def _tail(text: str, limit: int = 20000) -> str:
    text = text or ""
    return text[-limit:]


def _read_preview(paths, limit: int = 200000) -> str:
    for path in paths or []:
        try:
            target = Path(path)
            if not target.is_file() or target.stat().st_size == 0:
                continue
            data = target.read_text(encoding="utf-8", errors="replace")
            if len(data) > limit:
                data = data[:limit] + "\n... (preview truncated)"
            return data
        except (OSError, UnicodeError):
            continue
    return ""


class TaskJobInspector(QDialog):
    """Non-blocking view of a task's state, console and outputs."""

    def __init__(self, queue, job_id: str, parent=None):
        super().__init__(parent)
        self.queue = queue
        self.job_id = job_id
        self.setWindowTitle("Task Inspect")
        self.setMinimumSize(720, 540)
        self.setAttribute(Qt.WA_DeleteOnClose, True)
        self._build_ui()
        for signal_name in ("job_updated", "job_finished", "jobs_changed"):
            signal = getattr(queue, signal_name, None)
            if signal is not None:
                try:
                    signal.connect(self.refresh)
                except (TypeError, RuntimeError):
                    pass
        self.refresh()

    def _job(self):
        getter = getattr(self.queue, "get_job", None)
        return getter(self.job_id) if callable(getter) else None

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)

        self.title_label = QLabel("Task")
        self.title_label.setStyleSheet("font-size: 17px; font-weight: 600;")
        layout.addWidget(self.title_label)

        self.status_label = QLabel("")
        self.status_label.setStyleSheet("font-weight: 600;")
        layout.addWidget(self.status_label)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_overview(), "Overview")
        self.tabs.addTab(self._build_console(), "Console")
        self.tabs.addTab(self._build_output(), "Output")
        layout.addWidget(self.tabs, 1)

        buttons = QHBoxLayout()
        self.show_window_button = QPushButton("Show tool window")
        self.show_window_button.clicked.connect(self._show_window)
        self.import_button = QPushButton("Import result")
        self.import_button.clicked.connect(self._import_result)
        self.cancel_button = QPushButton("Cancel task")
        self.cancel_button.clicked.connect(self._cancel_task)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        buttons.addWidget(self.show_window_button)
        buttons.addWidget(self.import_button)
        buttons.addWidget(self.cancel_button)
        buttons.addStretch(1)
        buttons.addWidget(close_button)
        layout.addLayout(buttons)

    def _build_overview(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        form = QFormLayout()
        self.tool_value = QLabel("-")
        self.created_value = QLabel("-")
        self.started_value = QLabel("-")
        self.finished_value = QLabel("-")
        self.progress_value = QLabel("-")
        self.progress_value.setWordWrap(True)
        self.error_value = QLabel("-")
        self.error_value.setWordWrap(True)
        form.addRow("Tool:", self.tool_value)
        form.addRow("Created:", self.created_value)
        form.addRow("Started:", self.started_value)
        form.addRow("Finished:", self.finished_value)
        form.addRow("Progress:", self.progress_value)
        form.addRow("Error:", self.error_value)
        layout.addLayout(form)

        layout.addWidget(QLabel("Command"))
        self.command_text = QPlainTextEdit()
        self.command_text.setReadOnly(True)
        self.command_text.setMaximumHeight(72)
        layout.addWidget(self.command_text)

        layout.addWidget(QLabel("Artifacts"))
        self.artifacts_text = QPlainTextEdit()
        self.artifacts_text.setReadOnly(True)
        self.artifacts_text.setMaximumHeight(120)
        layout.addWidget(self.artifacts_text)
        return page

    def _build_console(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        self.console_text = QPlainTextEdit()
        self.console_text.setReadOnly(True)
        self.console_text.setStyleSheet(
            "font-family: Consolas, monospace; background: #272822; color: #f8f8f2;"
        )
        layout.addWidget(self.console_text)
        return page

    def _build_output(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        self.output_text = QPlainTextEdit()
        self.output_text.setReadOnly(True)
        layout.addWidget(self.output_text)
        return page

    def refresh(self, *args):
        job = self._job()
        if job is None:
            self.title_label.setText("Task not found")
            self.status_label.setText("")
            return
        self.title_label.setText(str(getattr(job, "title", "Task")))
        state = str(getattr(job, "state", ""))
        status = str(getattr(job, "status_label", state.title()))
        elapsed = ""
        if callable(getattr(job, "elapsed_text", None)):
            try:
                elapsed = " · %s" % job.elapsed_text()
            except Exception:
                elapsed = ""
        self.setWindowTitle("Task Inspect - %s" % getattr(job, "title", ""))
        self.status_label.setText("%s%s" % (status, elapsed))
        self.status_label.setStyleSheet(
            "font-weight: 600; color: %s;" % STATE_COLORS.get(state, "#374151")
        )
        self.tool_value.setText(str(getattr(job, "tool_name", "")))
        created = getattr(job, "created_at", 0)
        started = getattr(job, "started_at", 0)
        finished = getattr(job, "finished_at", 0)
        self.created_value.setText(self._time(created))
        self.started_value.setText(self._time(started))
        self.finished_value.setText(self._time(finished))
        progress = str(getattr(job, "progress", "") or "")
        self.progress_value.setText(progress or "-")
        error = str(getattr(job, "error_message", "") or "")
        self.error_value.setText(error or "-")
        command = getattr(job, "command", None) or []
        self.command_text.setPlainText(" ".join(str(part) for part in command))
        artifacts = list(getattr(job, "output_files", []) or [])
        self.artifacts_text.setPlainText("\n".join(artifacts) or "-")
        plugin = getattr(job, "plugin", None)
        console = str(getattr(job, "console_text", "") or "")
        if not console and plugin is not None:
            try:
                from YR_MPE.runtime.task_queue import _safe_console
                console = _safe_console(plugin)
            except (RuntimeError, AttributeError):
                console = ""
        self.console_text.setPlainText(_tail(console or error or "-"))

        preview = _read_preview(artifacts)
        if not preview and plugin is not None:
            try:
                from YR_MPE.runtime.task_queue import _plugin_outputs
                preview = _read_preview(_plugin_outputs(plugin))
            except (RuntimeError, AttributeError):
                preview = ""
        self.output_text.setPlainText(preview or "No result file is available yet.")

        window = getattr(job, "window", None)
        try:
            has_window = window is not None and callable(
                getattr(window, "isVisible", None)
            )
        except RuntimeError:
            has_window = False
        self.show_window_button.setEnabled(bool(has_window))
        try:
            can_import = bool(
                state == FINISHED
                and plugin is not None
                and callable(getattr(plugin, "import_to_platform", None))
            )
        except RuntimeError:
            can_import = False
        self.import_button.setEnabled(can_import)
        self.cancel_button.setEnabled(
            bool(getattr(job, "can_cancel", False))
        )

    @staticmethod
    def _time(timestamp):
        if not timestamp:
            return "-"
        try:
            import datetime
            return datetime.datetime.fromtimestamp(timestamp).strftime(
                "%Y-%m-%d %H:%M:%S"
            )
        except Exception:
            return "-"

    def _show_window(self):
        show = getattr(self.queue, "show_job_window", None)
        if callable(show):
            show(self.job_id)

    def _import_result(self):
        importer = getattr(self.queue, "import_job", None)
        if callable(importer):
            result = importer(self.job_id)
            if not result.get("ok"):
                self.error_value.setText(str(result.get("error", "Import failed")))
            self.refresh()

    def _cancel_task(self):
        cancel = getattr(self.queue, "cancel", None)
        if callable(cancel):
            cancel(self.job_id)
            self.refresh()
