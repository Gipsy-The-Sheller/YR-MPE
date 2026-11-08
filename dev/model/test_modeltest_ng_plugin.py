"""Tests for the ModelTest-NG plugin."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.factories import PluginFactory  # noqa: E402
from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_MODEL,
)
from YR_MPE.plugins.modeltest_ng_plugin import (  # noqa: E402
    ModelTestNGPlugin, ModelTestNGPluginEntry, ModelTestNGThread,
    detect_data_type, parse_modeltest_output,
)

APP = QApplication.instance() or QApplication(sys.argv)
FIXTURE = REPO / "dev" / "model" / "testdata" / "four_taxa.fas"
TOOL = REPO / "YR_MPE" / "bin" / "modeltest-ng" / "modeltest-ng.exe"


class ModelTestNGParserTests(unittest.TestCase):
    def test_detects_dna(self):
        self.assertEqual(detect_data_type(str(FIXTURE)), "nt")

    def test_aic_and_aicc_best_blocks_stay_distinct(self):
        report = Path(tempfile.mkdtemp()) / "synthetic.out"
        report.write_text(
            "Best model according to BIC\n---\nModel:  TVM\nlnL:  -1.0\n---\n"
            "Best model according to AIC\n---\nModel:  GTR\nlnL:  -2.0\n---\n"
            "Best model according to AICc\n---\nModel:  TrN+I\nlnL:  -3.0\n---\n",
            encoding="utf-8",
        )
        parsed = parse_modeltest_output(str(report))
        self.assertEqual(parsed["best"]["BIC"]["Model"], "TVM")
        self.assertEqual(parsed["best"]["AIC"]["Model"], "GTR")
        self.assertEqual(parsed["best"]["AICc"]["Model"], "TrN+I")

    def test_host_actions_and_model_handoff(self):
        from PyQt5.QtWidgets import QAction
        from YR_MPE.platforms.mpea import YR_MPEA_Widget

        host = YR_MPEA_Widget()
        try:
            labels = {action.text() for action in host.findChildren(QAction)}
            self.assertIn("Find Best Substitution Model (ModelTest-NG)", labels)
            self.assertIn("RAxML-NG", labels)

            manager = host.dataset_selection_manager
            dataset_id = manager.create_dataset("model-handoff")
            item = DatasetItem(item_type=ITEM_TYPE_MODEL)
            item.dataset_id = dataset_id
            item.loci_name = "best"
            item.data = {
                "type": "model_table",
                "best_model": "TrN+I",
                "data": [{"Model": "TrN+I"}],
            }
            self.assertTrue(manager.add_item(item, dataset_id))
            self.assertTrue(manager.add_to_basket(item.id))
            self.assertEqual(host._selected_single_model_name(), "TrN+I")
        finally:
            host.close()

    def test_help_doc_maps_every_cli_control_to_help(self):
        plugin = ModelTestNGPlugin()
        try:
            doc = plugin.help_doc()
            self.assertEqual(doc["schema"], "yrmpe.skill.v1")
            self.assertEqual(doc["id"], "modeltest-ng")
            self.assertTrue(doc["markdown"])
            self.assertEqual(doc["warnings"], [])
            controls = {row["id"]: row for row in doc["controls"]}
            self.assertEqual(controls["data_type"]["choices"]["DNA"], "nt")
            self.assertEqual(controls["topology"]["cli_args"], ["-t", "ml"])
            help_text = subprocess.run(
                [str(TOOL), "--help"], capture_output=True, text=True
            ).stdout
            for row in doc["controls"]:
                for flag in row.get("cli", []):
                    if str(flag).startswith("-"):
                        self.assertIn(flag, help_text)
        finally:
            plugin.close()

    def test_generic_agent_surface(self):
        plugin = ModelTestNGPlugin()
        try:
            plugin.imported_files = [str(FIXTURE)]
            plugin.file_path_edit.setText(str(FIXTURE))
            contract = plugin.agent_contract()
            self.assertEqual(contract["schema"], "yrmpe.plugin.v1")
            self.assertIn("is_running", contract["status"])

            controls = plugin.agent_invoke("plugin.inspect_controls", {})
            self.assertTrue(controls["ok"])
            ids = {row["id"] for row in controls["result"]["controls"]}
            self.assertIn("data_type", ids)

            changed = plugin.agent_invoke(
                "plugin.modify_options",
                {"options": {"data_type": "nt", "threads": 3,
                             "criterion": "AICc"}},
            )
            self.assertTrue(changed["ok"])
            self.assertEqual(
                set(changed["result"]["applied"]),
                {"data_type", "threads", "criterion"},
            )
            self.assertEqual(plugin.data_type_combo.currentText(), "DNA")
            self.assertEqual(plugin.threads_spinbox.value(), 3)
            self.assertEqual(plugin.criterion_combo.currentText(), "AICc")

            status = plugin.agent_invoke("plugin.inspect_status", {})
            self.assertTrue(status["ok"])
            self.assertFalse(status["result"]["is_running"])
        finally:
            plugin.close()

    def test_workflow_mode_locks_inputs_and_applies_config(self):
        plugin = ModelTestNGPlugin()
        try:
            context = {
                "step_id": "s02",
                "op": "modeltest.run",
                "inputs": {"input": "%mafft.alignment.fas%"},
                "outputs": {"model": "%s02.model.out%"},
            }
            plugin.bind_workflow_step(context)
            self.assertEqual(
                plugin.file_path_edit.text(), "%mafft.alignment.fas%"
            )
            self.assertTrue(plugin.file_path_edit.isReadOnly())
            self.assertFalse(plugin.file_browse_btn.isEnabled())
            self.assertEqual(plugin.run_button.text(), "Close && Apply")
            captured = []
            plugin.workflow_configured.connect(captured.append)
            plugin._apply_workflow_step()
            self.assertEqual(len(captured), 1)
            config = captured[0]
            self.assertEqual(config["id"], "s02")
            self.assertEqual(config["op"], "modeltest.run")
            self.assertEqual(
                config["inputs"]["input"], "%mafft.alignment.fas%"
            )
            self.assertEqual(
                config["locked_inputs"][0]["placeholder"],
                "%mafft.alignment.fas%",
            )
        finally:
            plugin.close()

    def test_plugin_metadata_and_factory(self):
        plugin = ModelTestNGPlugin()
        try:
            self.assertEqual(plugin.tool_name, "ModelTest-NG")
            self.assertTrue(plugin.tool_path)
            self.assertEqual(plugin.get_citation_keys(), ["darriba2020modeltest"])
            factory = PluginFactory()
            self.assertIsInstance(
                factory.get_modeltest_plugin(), ModelTestNGPluginEntry
            )
        finally:
            plugin.close()


@unittest.skipUnless(TOOL.is_file(), "ModelTest-NG runtime not available")
class ModelTestNGThreadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(tempfile.mkdtemp())
        cls.payloads = []
        cls.errors = []
        cls.thread = ModelTestNGThread(
            str(TOOL), [str(FIXTURE)], data_type="nt", topology="ml",
            schemes="11", threads=2, workdir=str(cls.root),
        )
        cls.thread.result_ready.connect(cls.payloads.append)
        cls.thread.error.connect(cls.errors.append)
        cls.thread.run()

    def test_run_succeeds(self):
        self.assertEqual(self.errors, [])
        self.assertEqual(len(self.payloads), 1)
        payload = self.payloads[0]
        self.assertEqual(len(payload["output_files"]), 1)
        self.assertTrue(Path(payload["output_files"][0]).is_file())

    def test_parses_all_criteria(self):
        payload = self.payloads[0]
        report = payload["parsed"][payload["output_files"][0]]
        self.assertEqual(set(report["criteria"]), {"BIC", "AIC", "AICc"})
        self.assertTrue(report["criteria"]["BIC"])
        self.assertIn("Model", report["best"]["BIC"])

    def test_result_spec_is_importable(self):
        plugin = ModelTestNGPlugin(workdir=str(self.root))
        try:
            plugin._on_result(self.payloads[0])
            plugin.criterion_combo.setCurrentText("BIC")
            result = plugin._build_result_spec()
            self.assertEqual(result.validate(), [])
            item = result.datasets[0].items[0]
            self.assertEqual(item.item_type, "model")
            self.assertEqual(item.data["type"], "model_table")
            self.assertEqual(
                item.data["best_model"], item.data["data"][0]["Model"]
            )
            self.assertIn("BIC", item.data["criteria"])
        finally:
            plugin.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
