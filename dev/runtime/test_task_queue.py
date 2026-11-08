"""Tests for the host task queue and the Selected & Queue panel."""

from __future__ import annotations

import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import QThread, Qt, pyqtSignal
from PyQt5.QtWidgets import QApplication, QDialog, QVBoxLayout, QWidget

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.methods.selected_panel import SelectedPanel  # noqa: E402
from YR_MPE.runtime.task_queue import (  # noqa: E402
    CANCELLED, ERROR, FINISHED, QUEUED, RUNNING, TaskQueueManager,
)

APP = QApplication.instance() or QApplication(sys.argv)


class FakeWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, list)
    error = pyqtSignal(str)

    def __init__(self, delay=0.03, fail=False):
        super().__init__()
        self.delay = delay
        self.fail = fail
        self._stop = False

    def run(self):
        if self.delay:
            time.sleep(self.delay)
        if self._stop:
            return
        if self.fail:
            self.error.emit("synthetic failure")
            return
        self.progress.emit("half")
        self.finished.emit(["output.tree"], [])

    def stop(self):
        self._stop = True


class FakePlugin(QWidget):
    def __init__(self, delay=0.03, fail=False):
        super().__init__()
        self.plugin_name = "Fake analysis"
        self.tool_name = "fake"
        self.is_running = False
        self.thread = None
        self.current_output_files = []
        self.console_output = []
        self.delay = delay
        self.fail = fail
        self.completion_dialogs = 0

    def run_analysis(self):
        self.is_running = True
        self.thread = FakeWorker(self.delay, self.fail)
        self.thread.finished.connect(self.analysis_finished)
        self.thread.error.connect(self.analysis_error)
        self.thread.start()

    def analysis_finished(self, files, html):
        from PyQt5.QtWidgets import QMessageBox
        self.is_running = False
        self.current_output_files = files
        self.completion_dialogs += 1
        QMessageBox.information(self, "Finished", "synthetic finish")

    def analysis_error(self, message):
        from PyQt5.QtWidgets import QMessageBox
        self.is_running = False
        QMessageBox.critical(self, "Error", message)

    def stop_analysis(self):
        worker = self.thread
        if worker is not None and worker.isRunning():
            worker.stop()
            worker.wait(3000)
        self.is_running = False

    def cli_preview(self):
        return ["fake", "--run"]


def wait_for(app, predicate, timeout=3.0):
    end = time.time() + timeout
    while time.time() < end:
        app.processEvents()
        if predicate():
            return True
        time.sleep(0.01)
    app.processEvents()
    return predicate()


class HandlePlugin(QWidget):
    """Minimal plugin whose completion methods use the handle_* names."""

    def __init__(self):
        super().__init__()
        self.plugin_name = "Handle analysis"
        self.tool_name = "handle"
        self.is_running = False
        self.thread = None
        self.current_output_files = []
        self.console_output = []
        self.finished_called = 0

    def run_analysis(self):
        self.is_running = True
        self.thread = FakeWorker()
        self.thread.finished.connect(self.handle_finished)
        self.thread.error.connect(self.handle_error)
        self.thread.start()

    def handle_finished(self, files, html):
        self.is_running = False
        self.current_output_files = files
        self.finished_called += 1

    def handle_error(self, message):
        self.is_running = False

    def stop_analysis(self):
        if self.thread is not None:
            self.thread.wait(3000)
        self.is_running = False


class TaskQueueTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="yrmpe-queue-"))
        self.manager = DatasetSelectionManager(disable_auto_save=True)
        self.manager.workspace_path = str(self.root)
        self.queue = TaskQueueManager(max_concurrent=1)
        self.panel = SelectedPanel(self.manager, task_queue=self.queue)
        self.plugin = FakePlugin()
        self.plugin.show()
        self.queue.attach_plugin(self.plugin)

    def tearDown(self):
        job = self.queue.jobs()[-1] if self.queue.jobs() else None
        if job is not None:
            thread = getattr(self.plugin, "thread", None)
            if thread is not None and thread.isRunning():
                thread.stop()
                thread.wait(3000)
        self.panel.close()
        self.plugin.close()

    def test_finished_job_updates_status_and_panel(self):
        job_id = self.queue.submit(self.plugin)
        self.assertEqual(self.queue.get_job(job_id).state, RUNNING)
        self.assertTrue(wait_for(
            APP, lambda: self.queue.get_job(job_id).state == FINISHED
        ))
        job = self.queue.get_job(job_id)
        self.assertEqual(job.output_files, ["output.tree"])
        self.assertEqual(
            [row["text"] for row in job.messages], ["synthetic finish"]
        )
        APP.processEvents()
        self.assertEqual(self.panel.jobs_tree.topLevelItemCount(), 1)
        self.assertEqual(
            self.panel.jobs_tree.topLevelItem(0).text(1), "Finished"
        )
        self.assertIn("Finished", self.panel.queue_notice_label.text())

    def test_handle_finished_is_recognized_by_queue(self):
        plugin = HandlePlugin()
        try:
            job_id = self.queue.submit(plugin)
            self.assertTrue(wait_for(
                APP, lambda: self.queue.get_job(job_id).state == FINISHED
            ))
            self.assertEqual(plugin.finished_called, 1)
            self.assertFalse(plugin.is_running)
            self.assertEqual(plugin.current_output_files, ["output.tree"])
        finally:
            plugin.close()

    def test_error_job_is_reported_without_modal_block(self):
        self.plugin.fail = True
        job_id = self.queue.submit(self.plugin)
        self.assertTrue(wait_for(
            APP, lambda: self.queue.get_job(job_id).state == ERROR
        ))
        job = self.queue.get_job(job_id)
        self.assertIn("synthetic failure", job.error_message)
        APP.processEvents()
        self.assertEqual(
            self.panel.jobs_tree.topLevelItem(0).text(1), "Error"
        )

    def test_second_job_waits_and_can_be_cancelled(self):
        first = FakePlugin(delay=0.4)
        second = FakePlugin(delay=0.01)
        first.show()
        second.show()
        self.queue.attach_plugin(first)
        self.queue.attach_plugin(second)
        first_id = self.queue.submit(first)
        second_id = self.queue.submit(second)
        self.assertEqual(self.queue.get_job(first_id).state, RUNNING)
        self.assertEqual(self.queue.get_job(second_id).state, QUEUED)
        self.assertTrue(self.queue.cancel(second_id))
        self.assertEqual(self.queue.get_job(second_id).state, CANCELLED)
        self.assertTrue(self.queue.cancel(first_id))
        self.assertEqual(self.queue.get_job(first_id).state, CANCELLED)
        first.close()
        second.close()

    def test_closing_dialog_does_not_cancel_job(self):
        dialog = QDialog()
        dialog.setAttribute(Qt.WA_DeleteOnClose, True)
        dialog.setLayout(QVBoxLayout())
        dialog.layout().addWidget(self.plugin)
        dialog.show()
        job_id = self.queue.submit(self.plugin, window=dialog)
        APP.processEvents()
        dialog.close()
        APP.processEvents()
        self.assertEqual(self.queue.get_job(job_id).state, RUNNING)
        self.assertTrue(wait_for(
            APP, lambda: self.queue.get_job(job_id).state == FINISHED
        ))
        self.assertIs(self.queue.get_job(job_id).window, dialog)

    def test_job_payload_for_agent_and_inspector(self):
        job_id = self.queue.submit(self.plugin)
        self.assertTrue(wait_for(
            APP, lambda: self.queue.get_job(job_id).state == FINISHED
        ))
        payload = self.queue.job_payload(job_id)
        self.assertEqual(payload["state"], FINISHED)
        self.assertEqual(payload["command"], ["fake", "--run"])
        self.assertIn("artifacts", payload)


if __name__ == "__main__":
    unittest.main(verbosity=2)
