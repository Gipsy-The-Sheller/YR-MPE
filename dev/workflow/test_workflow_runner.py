"""Tests for the compiled-workflow runner plugin."""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import QApplication, QMessageBox

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import SingleGeneWorkspace  # noqa: E402
from YR_MPE.plugins.workflow_runner_plugin import (  # noqa: E402
    WorkflowRunnerPlugin,
)

APP = QApplication.instance() or QApplication(sys.argv)

MAFFT = REPO / "YR_MPE" / "bin" / "mafft" / "mafft.bat"
TRIMAL = REPO / "YR_MPE" / "bin" / "trimAl" / "bin" / "trimal.exe"
FASTA = ">a\nATGCAAATGCAA\n>b\nATGCAAATGCAA\n>c\nATGAAAATGAAA\n"


class WorkflowRunnerTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="yrmpe-wf-"))
        self.fasta = self.root / "raw.fas"
        self.fasta.write_text(FASTA, encoding="utf-8")
        self.spec = self.root / "workflow.jsonl"
        self.spec.write_text(
            "\n".join([
                json.dumps({
                    "schema": "yrmpe.workflow.v1",
                    "name": "test workflow",
                    "inputs": {"raw": str(self.fasta)},
                }),
                json.dumps({
                    "id": "mafft", "op": "mafft.run",
                    "outputs": {"alignment": "aln.fas"},
                    "params": {"input": "%input.raw.fas%", "threads": 1},
                }),
            ]) + "\n",
            encoding="utf-8",
        )
        self._dialogs = {}
        for name in ("warning", "critical", "information", "question"):
            self._dialogs[name] = getattr(QMessageBox, name)
            setattr(QMessageBox, name,
                    staticmethod(lambda *a, **k: QMessageBox.Ok))

    def tearDown(self):
        for name, value in self._dialogs.items():
            setattr(QMessageBox, name, value)

    def _plugin(self):
        return WorkflowRunnerPlugin(workdir=str(self.root), headless=True)

    def test_compile_and_validate(self):
        plugin = self._plugin()
        try:
            plugin.workflow_path_edit.setText(str(self.spec))
            plan = plugin.load_workflow()
            self.assertIsNotNone(plan)
            self.assertEqual(plugin.validate_workflow(plan), [])
            self.assertEqual(len(plan.steps), 1)
            self.assertIn("raw", plugin.input_edits)
        finally:
            plugin.close()

    def test_unknown_op_is_rejected(self):
        spec = self.root / "bad.jsonl"
        spec.write_text(
            json.dumps({
                "schema": "yrmpe.workflow.v1",
                "inputs": {"raw": str(self.fasta)},
            }) + "\n" + json.dumps({
                "id": "nope", "op": "does_not_exist.run",
                "outputs": {"x": "x.txt"},
            }) + "\n",
            encoding="utf-8",
        )
        plugin = self._plugin()
        try:
            plugin.workflow_path_edit.setText(str(spec))
            plan = plugin.load_workflow()
            problems = plugin.validate_workflow(plan)
            self.assertTrue(any("does_not_exist" in p for p in problems))
        finally:
            plugin.close()

    @unittest.skipUnless(MAFFT.exists(), "bundled MAFFT is not available")
    def test_two_step_run(self):
        with open(self.spec, "a", encoding="utf-8") as handle:
            handle.write(json.dumps({
                "id": "trimal", "op": "trimal.run",
                "outputs": {"trimmed": "trim.fas"},
                "params": {
                    "input": "%mafft.alignment.fas%",
                    "mode": "-gappyout",
                },
            }) + "\n")
        plugin = self._plugin()
        captured = []
        plugin.result_ready.connect(lambda spec: captured.append(spec))
        plugin.workflow_path_edit.setText(str(self.spec))
        plugin.output_dir_edit.setText(str(self.root / "out"))
        plugin.load_workflow()
        plugin.agent_invoke("workflow.run")

        def wait():
            if plugin.is_running:
                QTimer.singleShot(250, wait)
            else:
                APP.quit()

        QTimer.singleShot(400, wait)
        QTimer.singleShot(120000, APP.quit)
        APP.exec_()

        status = plugin.inspect_status()
        self.assertFalse(status["is_running"])
        self.assertEqual(status["steps_done"], 2)
        self.assertEqual([o["step"] for o in plugin._outputs],
                         ["mafft", "trimal"])
        self.assertTrue((self.root / "out" / "workflow_run.json").exists())
        self.assertEqual(len(captured), 1)
        types = [item.item_type
                 for dataset in captured[0].datasets
                 for item in dataset.items]
        self.assertEqual(types.count("alignment"), 2)
        runs = captured[0].runs
        self.assertEqual([run.method for run in runs], ["mafft.run", "trimal.run"])
        self.assertEqual([run.metadata["step_id"] for run in runs],
                         ["mafft", "trimal"])
        self.assertTrue(all(run.output_paths for run in runs))
        self.assertTrue(any(run.command for run in runs))
        self.assertTrue(runs[0].parameters)
        self.assertTrue(runs[1].parameters)
        plugin.close()

    @unittest.skipUnless(MAFFT.exists(), "bundled MAFFT is not available")
    def test_workflow_import_creates_step_runs(self):
        with open(self.spec, "a", encoding="utf-8") as handle:
            handle.write(json.dumps({
                "id": "trimal", "op": "trimal.run",
                "outputs": {"trimmed": "trim.fas"},
                "params": {
                    "input": "%mafft.alignment.fas%",
                    "mode": "-gappyout",
                },
            }) + "\n")
        manager = DatasetSelectionManager()
        manager.set_workspace_path(str(self.root / "workspace"))
        workspace = SingleGeneWorkspace(dataset_selection_manager=manager)
        plugin = self._plugin()
        plugin.result_ready.connect(workspace.add_result_spec)
        plugin.workflow_path_edit.setText(str(self.spec))
        plugin.output_dir_edit.setText(str(self.root / "out"))
        plugin.load_workflow()
        plugin.agent_invoke("workflow.run")

        def wait():
            if plugin.is_running:
                QTimer.singleShot(250, wait)
            else:
                APP.quit()

        QTimer.singleShot(400, wait)
        QTimer.singleShot(120000, APP.quit)
        APP.exec_()

        runs = sorted(
            manager.store.index.runs.values(), key=lambda run: run.method
        )
        self.assertEqual([run.method for run in runs],
                         ["mafft.run", "trimal.run"])
        self.assertEqual([run.metadata["step_id"] for run in runs],
                         ["mafft", "trimal"])
        self.assertEqual(len(runs[0].output_ids), 1)
        self.assertEqual(runs[1].input_ids, runs[0].output_ids)
        for run in runs:
            self.assertTrue(run.artifact_hashes)
            self.assertEqual(run.status, "completed")
        dataset_id = list(manager.store.index.datasets)[0]
        items = manager.get_items_by_dataset(dataset_id)
        self.assertEqual(
            {item.run_id for item in items},
            {run.id for run in runs},
        )
        plugin.close()
        workspace.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
