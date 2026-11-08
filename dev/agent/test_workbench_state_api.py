"""Tests for menu/workspace/sidebar/skill discovery APIs."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_ALIGNMENT,
)
from YR_MPE.platforms.mpea import YR_MPEA_Widget  # noqa: E402
from YR_MPE.plugins.modeltest_ng_plugin import ModelTestNGPlugin  # noqa: E402
from YR_MPE.runtime.host_actions import build_workbench_registry  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
DEMO_WORKFLOW = REPO / "dev" / "workflow" / "demo_pipeline.jsonl"


class WorkbenchStateApiTests(unittest.TestCase):
    def setUp(self):
        self.host = YR_MPEA_Widget()
        self.registry = build_workbench_registry(self.host)
        manager = self.host.dataset_selection_manager
        self.dataset_id = manager.create_dataset("state-api")
        self.item = DatasetItem(item_type=ITEM_TYPE_ALIGNMENT)
        self.item.dataset_id = self.dataset_id
        self.item.loci_name = "gene1"
        self.item.format = "fasta"
        self.item.is_aligned = True
        self.item.length = 120
        self.item.sequences = [
            SeqRecord(Seq("ACGT" * 30), id="taxon%d" % index)
            for index in range(4)
        ]
        self.item.sequence_count = 4
        self.item.data = {"source": "test"}
        self.assertTrue(manager.add_item(self.item, self.dataset_id))
        self.assertTrue(manager.add_to_basket(self.item.id))
        self.host.current_dataset_id = self.dataset_id

    def tearDown(self):
        self.host.close()

    def test_menu_tree_exposes_stable_action_ids(self):
        result = self.registry.invoke("workbench.menu", {})
        self.assertTrue(result["ok"])
        actions = result["result"]["actions"]
        labels = {row["label"]: row for row in actions}
        self.assertIn("Find Best Substitution Model (ModelTest-NG)", labels)
        self.assertIn("RAxML-NG", labels)
        raxml = labels["RAxML-NG"]
        self.assertEqual(
            raxml["path"],
            ["PHYLOGENY", "Maximum Likelihood (ML)", "RAxML-NG"],
        )
        self.assertTrue(raxml["id"])
        self.assertIn("maximum-likelihood", raxml["description"])

        # The menu description should be short enough for an agent to use
        # without opening another discovery call.
        self.assertEqual(
            labels["MAFFT"]["description"],
            "batch align nucleotide and amino acid sequences",
        )
        self.assertEqual(
            labels["MACSE"]["description"],
            "codon aligner aware of frameshifts",
        )

        filtered = self.registry.invoke(
            "workbench.menu", {"query": "frameshift"}
        )["result"]["actions"]
        self.assertIn("MACSE", {row["label"] for row in filtered})

    def test_invoke_menu_rejects_unknown_and_disabled(self):
        unknown = self.registry.invoke(
            "workbench.invoke_menu", {"id": "not_a_real_action"}
        )
        self.assertFalse(unknown["ok"])
        disabled = self.registry.invoke(
            "workbench.invoke_menu",
            {"id": "species_tree_gene_concordance_gcf_planned"},
        )
        self.assertFalse(disabled["ok"])

    def test_workspace_state_lists_items_and_selection(self):
        result = self.registry.invoke("workbench.workspace", {})
        self.assertTrue(result["ok"])
        payload = result["result"]
        self.assertEqual(payload["dataset_count"], 1)
        self.assertEqual(payload["item_count"], 1)
        item = payload["items"][0]
        self.assertEqual(item["name"], "gene1")
        self.assertTrue(item["is_aligned"])
        self.assertIn(item["id"], payload["basket_item_ids"])
        self.assertIn(item["id"], payload["selected_item_ids"])
        self.assertEqual(payload["items_by_type"]["alignment"], [item["id"]])

    def test_sidebar_state_shows_recommended_operations(self):
        result = self.registry.invoke("workbench.sidebar", {})
        self.assertTrue(result["ok"])
        payload = result["result"]
        self.assertIn(self.item.id, payload["basket_item_ids"])
        operations = {
            op["id"]: op
            for group in payload["recommended_operations"]
            for op in group["operations"]
        }
        self.assertTrue(operations["iqtree"]["available"])
        self.assertFalse(operations["astral4"]["available"])
        self.host.inspector.show_item(self.item.id)
        payload = self.registry.invoke("workbench.sidebar", {})["result"]
        operation_ids = {
            op["id"] for op in payload["inspector"]["operations"]
        }
        self.assertIn("inspector.operation.show", operation_ids)
        self.assertIn("inspector.operation.preview", operation_ids)

    def test_import_workflow_jsonl_opens_plan(self):
        plan = self.host.import_workflow_file(str(DEMO_WORKFLOW))
        self.assertIsNotNone(plan)
        self.assertEqual(
            [step["id"] for step in plan.steps], ["mafft", "modeltest", "raxml"]
        )
        self.assertEqual(
            self.host.current_workflow_plan.name,
            "Demo: MAFFT -> ModelTest-NG -> RAxML-NG",
        )
        dialog = getattr(self.host, "workflow_preview_dialog", None)
        self.assertIsNotNone(dialog)
        dialog.close()

    def test_show_workflow_preview_opens_dialog(self):
        self.host.show_workflow_preview()
        dialog = getattr(self.host, "workflow_preview_dialog", None)
        self.assertIsNotNone(dialog)
        dialog.close()

    def test_human_file_dialog_is_not_used_by_agent(self):
        result = self.registry.invoke(
            "workbench.invoke_menu", {"id": "align_open_sequence_files"}
        )
        self.assertFalse(result["ok"])
        self.assertIn("workbench.import_files", result["error"])

    def test_workflow_graph_action(self):
        result = self.registry.invoke("workbench.workflow_graph", {})
        self.assertTrue(result["ok"])
        self.assertIn("nodes", result["result"])
        self.assertIn("edges", result["result"])
        self.assertGreaterEqual(result["result"]["node_count"], 1)

    def test_skills_catalog_and_help_doc(self):
        skills = self.registry.invoke("workbench.skills", {})
        self.assertTrue(skills["ok"])
        ids = {row["id"] for row in skills["result"]["skills"]}
        self.assertIn("modeltest-ng", ids)
        self.assertIn("raxml-ng", ids)
        self.assertIn("orthofinder", ids)

        doc = self.registry.invoke(
            "workbench.help_doc", {"id": "raxml-ng"}
        )
        self.assertTrue(doc["ok"])
        self.assertTrue(doc["result"]["markdown"])
        self.assertIn("--model", doc["result"]["markdown"])

    def test_live_plugin_help_doc_has_control_bindings(self):
        plugin = ModelTestNGPlugin()
        try:
            result = self.registry.invoke(
                "workbench.help_doc", {"plugin": "ModelTest-NG"}
            )
            self.assertTrue(result["ok"])
            doc = result["result"]
            self.assertEqual(doc["id"], "modeltest-ng")
            controls = {row["id"]: row for row in doc["controls"]}
            self.assertEqual(controls["threads"]["cli_args"], ["-p", "2"])
            self.assertEqual(controls["data_type"]["resolved"], "auto")
        finally:
            plugin.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
