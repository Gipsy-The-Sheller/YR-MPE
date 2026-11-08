"""Host-level task queue for long-running plugin analyses.

The queue owns the plugin object for the lifetime of a job.  A plugin window can
be closed while its worker thread continues; the queue keeps both the window
and the plugin alive so completion handlers can still update the hidden output
page.  Completion and error callbacks run through a message-box guard so a
queued job never blocks the GUI with a modal ``QMessageBox``.

The queue is intentionally small: one FIFO queue, a configurable concurrency
limit, and in-memory job records.  Persistence belongs to the workspace run
store; this layer only coordinates execution and UI state.
"""

from __future__ import annotations

import datetime as _datetime
import os
import time
import traceback
import uuid
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from PyQt5.QtCore import QObject, QThread, Qt, QTimer, pyqtSignal
from PyQt5.QtWidgets import QApplication, QMessageBox


QUEUED = "queued"
RUNNING = "running"
FINISHED = "finished"
ERROR = "error"
CANCELLED = "cancelled"

TERMINAL_STATES = frozenset((FINISHED, ERROR, CANCELLED))

STATE_LABELS = {
    QUEUED: "Queued",
    RUNNING: "Running",
    FINISHED: "Finished",
    ERROR: "Error",
    CANCELLED: "Cancelled",
}


def _clock() -> float:
    return time.time()


def _format_time(timestamp: Optional[float]) -> str:
    if not timestamp:
        return ""
    return _datetime.datetime.fromtimestamp(timestamp).strftime("%Y-%m-%d %H:%M:%S")


def _safe_window(plugin):
    """Return the top-level window for a plugin without assuming its type."""
    try:
        window = plugin.window()
        if window is not None:
            return window
    except Exception:
        pass
    return plugin if getattr(plugin, "isWidgetType", lambda: False)() else None


def _window_visible(window) -> bool:
    if window is None:
        return False
    try:
        return bool(window.isVisible())
    except (RuntimeError, AttributeError):
        return False


def _safe_cli_preview(plugin) -> List[str]:
    preview = getattr(plugin, "cli_preview", None)
    if not callable(preview):
        return []
    try:
        value = preview()
        if isinstance(value, (list, tuple)):
            return [str(item) for item in value]
        if value:
            return [str(value)]
    except Exception:
        pass
    return []


def _safe_console(plugin, limit: int = 20000) -> str:
    try:
        widget = getattr(plugin, "console_text", None)
        if widget is not None and callable(getattr(widget, "toPlainText", None)):
            try:
                return widget.toPlainText()[-limit:]
            except Exception:
                pass
        rows = getattr(plugin, "console_output", None) or []
    except (RuntimeError, AttributeError):
        return ""
    lines = []
    for row in rows[-400:]:
        if isinstance(row, dict):
            lines.append(str(row.get("formatted") or row.get("message") or ""))
        else:
            lines.append(str(row))
    return "\n".join(lines)[-limit:]


def _plugin_outputs(plugin) -> List[str]:
    for attr in (
        "current_output_files", "analysis_output_files", "alignment_output_files",
        "output_files",
    ):
        try:
            value = getattr(plugin, attr, None)
        except (RuntimeError, AttributeError):
            continue
        if isinstance(value, (list, tuple)):
            rows = [str(path) for path in value if path]
            if rows:
                return rows
    return []


class _MessageCapture:
    """Temporarily turn QMessageBox statics into recorders.

    Plugin completion handlers often end with ``QMessageBox.information`` or
    ``critical``.  A queued job is not allowed to block on those dialogs; the
    message is recorded for the queue inspector instead.  The patch is active
    only for the synchronous call that starts or completes a job, so normal GUI
    actions are unaffected.
    """

    _METHODS = (
        ("information", QMessageBox.Ok),
        ("warning", QMessageBox.Ok),
        ("critical", QMessageBox.Ok),
        ("question", QMessageBox.No),
        ("about", None),
    )

    def __init__(self):
        self.messages: List[Dict[str, str]] = []
        self._originals = {}
        self._exec_original = None

    def __enter__(self):
        for name, result in self._METHODS:
            original = getattr(QMessageBox, name, None)
            if original is None:
                continue
            self._originals[name] = original
            setattr(QMessageBox, name, staticmethod(
                self._make_recorder(name, result)
            ))
        self._exec_original = getattr(QMessageBox, "exec_", None)
        if self._exec_original is not None:
            try:
                QMessageBox.exec_ = self._fake_exec
            except Exception:
                self._exec_original = None
        return self

    def __exit__(self, exc_type, exc, tb):
        for name, original in self._originals.items():
            try:
                setattr(QMessageBox, name, staticmethod(original))
            except Exception:
                pass
        if self._exec_original is not None:
            try:
                QMessageBox.exec_ = self._exec_original
            except Exception:
                pass
        self._originals.clear()
        return False

    def _record(self, kind: str, args, kwargs):
        strings = [str(value) for value in args if isinstance(value, str)]
        title = strings[0] if strings else str(kwargs.get("title", ""))
        text = strings[1] if len(strings) > 1 else str(kwargs.get("text", ""))
        self.messages.append({
            "kind": kind,
            "title": title,
            "text": text,
        })

    def _make_recorder(self, kind: str, result):
        def recorder(*args, **kwargs):
            self._record(kind, args, kwargs)
            return result
        return recorder

    def _fake_exec(self, box):
        try:
            text = box.text() if callable(getattr(box, "text", None)) else ""
            title = box.windowTitle() if callable(
                getattr(box, "windowTitle", None)
            ) else ""
            self.messages.append({
                "kind": "dialog",
                "title": str(title),
                "text": str(text),
            })
        except Exception:
            pass
        return QMessageBox.Ok


def _extract_outputs(args: Sequence[Any]) -> Tuple[List[str], List[str]]:
    files: List[str] = []
    html: List[str] = []
    if args:
        value = args[0]
        if isinstance(value, (list, tuple)):
            files = [str(item) for item in value if item]
        elif isinstance(value, str) and value:
            files = [value]
    if len(args) > 1 and isinstance(args[1], (list, tuple)):
        html = [str(item) for item in args[1] if item]
    return files, html


def _first_callable(plugin, names: Iterable[str]):
    for name in names:
        value = getattr(plugin, name, None)
        if callable(value):
            return value
    return None


class TaskJob:
    """A single queued plugin run."""

    def __init__(self, plugin, title: str = ""):
        self.job_id = "job-%s" % uuid.uuid4().hex[:12]
        self.plugin = plugin
        self.title = title or str(
            getattr(plugin, "plugin_name", "")
            or getattr(plugin, "tool_name", "")
            or plugin.__class__.__name__
        )
        inputs = getattr(plugin, "imported_files", None) or []
        if inputs and not title:
            try:
                self.title += " · %s" % os.path.basename(str(inputs[0]))
            except Exception:
                pass
        self.tool_name = str(
            getattr(plugin, "tool_name", "")
            or getattr(plugin, "plugin_name", "")
            or plugin.__class__.__name__
        )
        self.state = QUEUED
        self.progress = ""
        self.error_message = ""
        self.messages: List[Dict[str, str]] = []
        self.output_files: List[str] = []
        self.console_text = ""
        self.result_payload = None
        self.command: List[str] = _safe_cli_preview(plugin)
        self.created_at = _clock()
        self.started_at: Optional[float] = None
        self.finished_at: Optional[float] = None
        self.window = None
        self.threads: List[QThread] = []
        self.multi_thread = False
        self._completed = 0
        self._had_error = False
        self._closed = False
        self._cancel_requested = False

    @property
    def status_label(self) -> str:
        return STATE_LABELS.get(self.state, self.state.title())

    @property
    def is_terminal(self) -> bool:
        return self.state in TERMINAL_STATES

    @property
    def can_cancel(self) -> bool:
        return self.state in (QUEUED, RUNNING)

    def elapsed_text(self) -> str:
        end = self.finished_at or _clock()
        start = self.started_at or self.created_at
        seconds = max(0, int(end - start))
        if seconds < 60:
            return "%ds" % seconds
        minutes, seconds = divmod(seconds, 60)
        if minutes < 60:
            return "%dm %02ds" % (minutes, seconds)
        hours, minutes = divmod(minutes, 60)
        return "%dh %02dm" % (hours, minutes)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "job_id": self.job_id,
            "title": self.title,
            "tool": self.tool_name,
            "state": self.state,
            "status": self.status_label,
            "progress": self.progress,
            "error": self.error_message,
            "messages": list(self.messages),
            "command": list(self.command),
            "artifacts": list(self.output_files),
            "console_tail": self.console_text[-4000:],
            "created_at": _format_time(self.created_at),
            "started_at": _format_time(self.started_at),
            "finished_at": _format_time(self.finished_at),
            "window_visible": _window_visible(self.window),
        }


class TaskQueueManager(QObject):
    """FIFO task queue with one GUI-thread plugin per running slot."""

    jobs_changed = pyqtSignal()
    job_added = pyqtSignal(str)
    job_updated = pyqtSignal(str)
    job_finished = pyqtSignal(str)

    def __init__(self, parent=None, max_concurrent: int = 1):
        super().__init__(parent)
        self.max_concurrent = max(1, int(max_concurrent or 1))
        self._jobs: List[TaskJob] = []
        self._pending: List[TaskJob] = []
        self._active: List[TaskJob] = []

    # ------------------------------------------------------------------
    # queries
    # ------------------------------------------------------------------
    def jobs(self) -> List[TaskJob]:
        return list(self._jobs)

    def get_job(self, job_id: str) -> Optional[TaskJob]:
        for job in self._jobs:
            if job.job_id == job_id:
                return job
        return None

    def active_jobs(self) -> List[TaskJob]:
        return list(self._active)

    def attach_plugin(self, plugin, window=None):
        """Tell a plugin that its Run button should use this queue."""
        setter = getattr(plugin, "set_task_queue", None)
        if callable(setter):
            try:
                setter(self, window=window)
            except TypeError:
                setter(self)
        else:
            plugin.task_queue = self
            plugin.task_window = window

    # ------------------------------------------------------------------
    # submission / scheduling
    # ------------------------------------------------------------------
    def submit(self, plugin, title: str = "", window=None) -> str:
        if plugin is None or not callable(getattr(plugin, "run_analysis", None)):
            raise ValueError("plugin does not provide run_analysis()")
        for job in self._jobs:
            if job.plugin is plugin and not job.is_terminal:
                return job.job_id
        if getattr(plugin, "is_running", False):
            raise RuntimeError("plugin is already running outside the queue")

        job = TaskJob(plugin, title=title)
        job.window = window or _safe_window(plugin)
        if job.window is not None:
            try:
                job.window.setAttribute(Qt.WA_DeleteOnClose, False)
            except Exception:
                pass
        self._jobs.append(job)
        self._pending.append(job)
        self.job_added.emit(job.job_id)
        self.jobs_changed.emit()
        self._pump()
        return job.job_id

    def _pump(self):
        while self._pending and len(self._active) < self.max_concurrent:
            job = self._pending.pop(0)
            if job.state != QUEUED:
                continue
            self._start_job(job)
        self.jobs_changed.emit()

    def _start_job(self, job: TaskJob):
        job.state = RUNNING
        job.started_at = _clock()
        job.finished_at = None
        job.error_message = ""
        job.messages = []
        job.output_files = []
        job.console_text = ""
        job.result_payload = None
        job._closed = False
        job._completed = 0
        job._had_error = False
        job.window = job.window or _safe_window(job.plugin)
        if job.window is not None:
            try:
                job.window.setAttribute(Qt.WA_DeleteOnClose, False)
            except Exception:
                pass
        job.command = _safe_cli_preview(job.plugin)
        self._active.append(job)
        self.job_updated.emit(job.job_id)
        self.jobs_changed.emit()

        capture = _MessageCapture()
        try:
            with capture:
                job.plugin.run_analysis()
            job.messages.extend(capture.messages)
        except Exception as exc:
            job.messages.extend(capture.messages)
            job.error_message = str(exc)
            self._fail_job(job, "Failed to start: %s" % exc)
            return

        if job._cancel_requested:
            self._cancel_running(job)
            return

        threads = self._plugin_threads(job.plugin)
        if not threads:
            if getattr(job.plugin, "is_running", False):
                # A synchronous plugin is unusual, but do not report success
                # before it actually returns.  Poll its state without blocking
                # the event loop.
                QTimer.singleShot(250, lambda: self._watch_synchronous(job))
                return
            details = self._message_text(job.messages)
            self._fail_job(
                job,
                details or "The plugin did not start a worker thread.",
            )
            return

        job.threads = threads
        job.multi_thread = bool(
            callable(getattr(job.plugin, "thread_finished", None))
            and _first_callable(
                job.plugin, ("analysis_finished", "on_analysis_finished",
                             "handle_finished")
            ) is None
        ) or bool(
            len(threads) > 1
            and callable(getattr(job.plugin, "thread_finished", None))
        )
        self._wire_threads(job)
        QTimer.singleShot(1500, lambda: self._check_wired_job(job))
        self.job_updated.emit(job.job_id)
        self.jobs_changed.emit()

    def _check_wired_job(self, job: TaskJob):
        """Catch a worker that finished before its signal was connected."""
        if job._closed or job._cancel_requested or job.state != RUNNING:
            return
        if not job.threads or any(thread.isRunning() for thread in job.threads):
            return
        if job._completed:
            return
        job.output_files = job.output_files or _plugin_outputs(job.plugin)
        if job._had_error or job.error_message:
            self._fail_job(job, job.error_message or "Analysis failed.")
        else:
            self._finish_job(job)

    def _watch_synchronous(self, job: TaskJob):
        if job._closed or job._cancel_requested:
            return
        if getattr(job.plugin, "is_running", False):
            QTimer.singleShot(250, lambda: self._watch_synchronous(job))
            return
        files = _plugin_outputs(job.plugin)
        job.output_files = files
        if job._had_error or job.error_message:
            self._fail_job(job, job.error_message or "Analysis failed.")
        else:
            self._finish_job(job)

    def _plugin_threads(self, plugin) -> List[QThread]:
        found: List[QThread] = []

        def add(value):
            if isinstance(value, QThread) and value not in found:
                found.append(value)

        for attr in ("analysis_thread", "alignment_thread", "thread"):
            add(getattr(plugin, attr, None))
        for attr in ("analysis_threads", "threads", "worker_threads"):
            values = getattr(plugin, attr, None)
            if isinstance(values, (list, tuple)):
                for value in values:
                    add(value)
        try:
            for value in vars(plugin).values():
                if isinstance(value, QThread):
                    add(value)
        except Exception:
            pass
        running = [thread for thread in found if thread.isRunning()]
        return running or found

    # ------------------------------------------------------------------
    # worker callbacks
    # ------------------------------------------------------------------
    def _disconnect_thread_signals(self, thread: QThread):
        for name in ("finished", "error"):
            signal = getattr(thread, name, None)
            if signal is None:
                continue
            try:
                signal.disconnect()
            except (TypeError, RuntimeError):
                pass

    def _wire_threads(self, job: TaskJob):
        plugin = job.plugin
        for index, thread in enumerate(job.threads):
            self._disconnect_thread_signals(thread)

            progress = getattr(thread, "progress", None)
            if progress is not None:
                try:
                    progress.connect(
                        lambda text, _job=job: self._on_progress(_job, text)
                    )
                except (TypeError, RuntimeError):
                    pass

            result = getattr(thread, "result_ready", None)
            if result is not None:
                try:
                    result.connect(
                        lambda payload, _job=job: self._on_result(_job, payload)
                    )
                except (TypeError, RuntimeError):
                    pass

            finished = getattr(thread, "finished", None)
            if finished is not None:
                try:
                    finished.connect(
                        lambda *args, _job=job, _index=index:
                        self._on_thread_finished(_job, _index, args)
                    )
                except (TypeError, RuntimeError):
                    pass

            error = getattr(thread, "error", None)
            if error is not None:
                try:
                    error.connect(
                        lambda message, _job=job, _index=index:
                        self._on_thread_error(_job, _index, message)
                    )
                except (TypeError, RuntimeError):
                    pass

    def _on_progress(self, job: TaskJob, text):
        if job._closed or job._cancel_requested:
            return
        job.progress = str(text or "")
        self.job_updated.emit(job.job_id)
        self.jobs_changed.emit()

    def _on_result(self, job: TaskJob, payload):
        if job._closed or job._cancel_requested:
            return
        job.result_payload = payload

    def _on_thread_finished(self, job: TaskJob, index: int, args):
        if job._closed or job._cancel_requested:
            return
        files, html_files = _extract_outputs(args)
        for path in files:
            if path and path not in job.output_files:
                job.output_files.append(path)
        capture = _MessageCapture()
        handler = None
        if job.multi_thread:
            handler = getattr(job.plugin, "thread_finished", None)
        else:
            handler = _first_callable(
                job.plugin, ("analysis_finished", "on_analysis_finished",
                             "handle_finished")
            )
        if callable(handler):
            try:
                with capture:
                    if job.multi_thread:
                        handler(files, html_files, index + 1)
                    else:
                        handler(files, html_files)
            except Exception:
                job._had_error = True
                job.error_message = traceback.format_exc(limit=4)
        job.messages.extend(capture.messages)

        if job.multi_thread:
            job._completed += 1
            if job._completed >= len(job.threads):
                if job._had_error:
                    self._fail_job(job, job.error_message or "One or more chains failed.")
                else:
                    self._finish_job(job)
        else:
            if job._had_error:
                self._fail_job(job, job.error_message or "Analysis failed.")
            else:
                self._finish_job(job)

    def _on_thread_error(self, job: TaskJob, index: int, message):
        if job._closed or job._cancel_requested:
            return
        job._had_error = True
        job.error_message = str(message or "Analysis failed.")
        capture = _MessageCapture()
        if job.multi_thread:
            handler = getattr(job.plugin, "thread_error", None)
            try:
                if callable(handler):
                    with capture:
                        handler(job.error_message, index + 1)
            except Exception:
                job.error_message += "\n" + traceback.format_exc(limit=4)
            job.messages.extend(capture.messages)
            job._completed += 1
            if job._completed >= len(job.threads):
                self._fail_job(job, job.error_message)
        else:
            handler = _first_callable(
                job.plugin, ("analysis_error", "on_analysis_error",
                             "handle_error")
            )
            try:
                if callable(handler):
                    with capture:
                        handler(job.error_message)
            except Exception:
                job.error_message += "\n" + traceback.format_exc(limit=4)
            job.messages.extend(capture.messages)
            self._fail_job(job, job.error_message)

    # ------------------------------------------------------------------
    # terminal state / cancellation
    # ------------------------------------------------------------------
    def _remove_active(self, job: TaskJob):
        if job in self._active:
            self._active.remove(job)

    def _finish_job(self, job: TaskJob):
        if job._closed:
            return
        job._closed = True
        job.state = FINISHED
        job.finished_at = _clock()
        job.progress = ""
        if not job.output_files:
            job.output_files = _plugin_outputs(job.plugin)
        job.console_text = _safe_console(job.plugin)
        self._remove_active(job)
        self.job_updated.emit(job.job_id)
        self.job_finished.emit(job.job_id)
        self.jobs_changed.emit()
        self._pump()

    def _fail_job(self, job: TaskJob, message: str):
        if job._closed:
            return
        job._closed = True
        job.state = ERROR
        job.error_message = str(message or "Analysis failed.")
        job.finished_at = _clock()
        job.progress = ""
        if not job.output_files:
            job.output_files = _plugin_outputs(job.plugin)
        job.console_text = _safe_console(job.plugin)
        self._remove_active(job)
        self.job_updated.emit(job.job_id)
        self.job_finished.emit(job.job_id)
        self.jobs_changed.emit()
        self._pump()

    def cancel(self, job_id: str) -> bool:
        job = self.get_job(job_id)
        if job is None or job.is_terminal:
            return False
        job._cancel_requested = True
        if job.state == QUEUED:
            if job in self._pending:
                self._pending.remove(job)
            job._closed = True
            job.state = CANCELLED
            job.finished_at = _clock()
            job.progress = ""
            self.job_updated.emit(job.job_id)
            self.job_finished.emit(job.job_id)
            self.jobs_changed.emit()
            self._pump()
            return True
        self._cancel_running(job)
        return True

    def _cancel_running(self, job: TaskJob):
        if job._closed:
            return
        job._cancel_requested = True
        guard = _MessageCapture()
        try:
            with guard:
                stop = getattr(job.plugin, "stop_analysis", None)
                if callable(stop):
                    stop()
        except Exception:
            pass
        finally:
            job.messages.extend(guard.messages)
        for thread in job.threads:
            if thread.isRunning():
                try:
                    thread.stop()
                    thread.wait(3000)
                except Exception:
                    pass
        job._closed = True
        job.state = CANCELLED
        job.finished_at = _clock()
        job.progress = ""
        try:
            job.console_text = _safe_console(job.plugin)
        except RuntimeError:
            pass
        self._remove_active(job)
        self.job_updated.emit(job.job_id)
        self.job_finished.emit(job.job_id)
        self.jobs_changed.emit()
        self._pump()

    # ------------------------------------------------------------------
    # UI helpers
    # ------------------------------------------------------------------
    def show_job_window(self, job_id: str) -> bool:
        job = self.get_job(job_id)
        if job is None:
            return False
        window = job.window or _safe_window(job.plugin)
        job.window = window or job.window
        if window is not None:
            try:
                window.show()
                window.raise_()
                window.activateWindow()
                return True
            except Exception:
                pass
        try:
            job.plugin.show()
            job.plugin.raise_()
            job.plugin.activateWindow()
            return True
        except Exception:
            return False

    def import_job(self, job_id: str) -> Dict[str, Any]:
        job = self.get_job(job_id)
        if job is None:
            return {"ok": False, "error": "unknown job"}
        if job.state != FINISHED:
            return {"ok": False, "error": "job is not finished"}
        importer = getattr(job.plugin, "import_to_platform", None)
        if not callable(importer):
            return {"ok": False, "error": "plugin does not expose import_to_platform()"}
        guard = _MessageCapture()
        try:
            with guard:
                importer()
        except Exception as exc:
            return {"ok": False, "error": str(exc)}
        job.messages.extend(guard.messages)
        return {"ok": True, "messages": list(guard.messages)}

    def job_payload(self, job_id: str) -> Optional[Dict[str, Any]]:
        job = self.get_job(job_id)
        return job.to_dict() if job is not None else None

    def state_payload(self) -> Dict[str, Any]:
        return {
            "schema": "yrmpe.jobs.v1",
            "max_concurrent": self.max_concurrent,
            "running": len(self._active),
            "queued": len(self._pending),
            "jobs": [job.to_dict() for job in self._jobs],
        }

    @staticmethod
    def _message_text(messages: Iterable[Dict[str, str]]) -> str:
        rows = []
        for row in messages:
            text = str(row.get("text") or "").strip()
            if text:
                rows.append(text)
        return "\n".join(rows)


__all__ = [
    "TaskQueueManager", "TaskJob", "QUEUED", "RUNNING", "FINISHED",
    "ERROR", "CANCELLED", "STATE_LABELS",
]
