"""Tests for the standalone InspectorPanel (offscreen Qt)."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QPushButton

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from render_inspector import build_sample_store  # noqa: E402
from YR_MPE.platforms.methods.inspector_panel import InspectorPanel  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class InspectorPanelTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.store, self.tree_item = build_sample_store(self.root)
        self.panel = InspectorPanel(store=self.store)

    def tearDown(self):
        self.panel.close()

    def _operation(self, name):
        return self.panel.findChild(QPushButton, "inspector.operation.%s" % name)

    def test_show_item(self):
        panel = self.panel
        panel.show_item(self.tree_item.id)
        self.assertEqual(panel.title_label.text(), "Tree (GTR+F+I+G4)")
        self.assertIn("Tree", panel.subtitle_label.text())
        self.assertIsNotNone(self._operation("preview"))
        self.assertEqual(self._operation("preview").text(), "Preview (IcyTree)")
        self.assertGreater(panel.metadata_form.rowCount(), 0)
        self.assertFalse(panel.run_section.isHidden())
        self.assertEqual(panel.run_combo.count(), 2)
        self.assertEqual(panel.run_name_edit.text(), "IQ-TREE run 1")
        self.assertGreater(panel.workflow_tree.topLevelItemCount(), 0)
        self.assertTrue(panel.workflow_tree.topLevelItem(0).isExpanded())

    def test_run_icon_uses_tool(self):
        panel = self.panel
        run = self.store.get_run(self.tree_item.run_id)
        self.assertEqual(panel._run_icon_path(run), "software/iqtree.svg")
        mafft = next(
            record for record in self.store.index.runs.values()
            if record.method == "mafft"
        )
        self.assertEqual(panel._run_icon_path(mafft), "software/mafft.svg")

    def test_sequence_has_no_run_switcher(self):
        panel = self.panel
        sequence = next(
            item for item in self.store.index.items.values()
            if item.type == "sequence"
        )
        panel.show_item(sequence.id)
        self.assertEqual(self._operation("preview").text(), "Preview sequences")
        self.assertTrue(panel.run_section.isHidden())

    def test_edge_view_and_rename(self):
        panel = self.panel
        panel.show_run(self.tree_item.run_id)
        self.assertEqual(panel.title_label.text(), "IQ-TREE run 1")
        self.assertIsNotNone(self._operation("view_command"))
        self.assertIsNotNone(self._operation("rerun"))

        renamed = []
        panel.run_renamed.connect(lambda run_id, name: renamed.append((run_id, name)))
        panel.run_name_edit.setText("IQ-TREE final")
        panel._on_rename_run()
        self.assertEqual(self.store.get_run(self.tree_item.run_id).name, "IQ-TREE final")
        self.assertEqual(renamed[0][1], "IQ-TREE final")

    def test_show_dataset(self):
        panel = self.panel
        panel.show_dataset(self.tree_item.dataset_id)
        self.assertEqual(panel.title_label.text(), "OrthoFinder_plugin_full")
        self.assertTrue(panel.run_section.isHidden())
        self.assertTrue(panel.workflow_tree.topLevelItem(0).isExpanded())
        self.assertGreaterEqual(panel.workflow_tree.topLevelItem(0).childCount(), 1)

    def test_tree_click_emits_item(self):
        panel = self.panel
        panel.show_item(self.tree_item.id)
        activated = []
        panel.item_activated.connect(activated.append)
        root = panel.workflow_tree.topLevelItem(0)
        import_node = root.child(0)
        item_node = import_node.child(0)
        panel._on_tree_clicked(item_node, 0)
        self.assertEqual(len(activated), 1)

    def test_metadata_expand_collapse(self):
        panel = self.panel
        panel.show_item(self.tree_item.id)
        self.assertFalse(panel._metadata_expanded)
        self.assertEqual(panel.metadata_scroll.maximumHeight(), panel._metadata_collapsed_height)
        self.assertEqual(panel.metadata_toggle.arrowType(), Qt.DownArrow)

        panel._toggle_metadata()
        self.assertTrue(panel._metadata_expanded)
        self.assertEqual(panel.metadata_toggle.arrowType(), Qt.UpArrow)
        self.assertTrue(panel.workflow_tree.isHidden())
        self.assertEqual(
            panel.metadata_scroll.verticalScrollBarPolicy(), Qt.ScrollBarAlwaysOn
        )

        panel._toggle_metadata()
        self.assertFalse(panel._metadata_expanded)
        self.assertEqual(panel.metadata_toggle.arrowType(), Qt.DownArrow)
        self.assertFalse(panel.workflow_tree.isHidden())

    def test_agent_contract(self):
        panel = self.panel
        contract = panel.agent_contract()
        ids = {action["id"] for action in contract["actions"]}
        self.assertIn("inspector.show_item", ids)
        self.assertIn("inspector.show_run", ids)
        self.assertIn("inspector.preview", ids)
        self.assertIn("inspector.rename_run", ids)
        self.assertIn("inspector.set_metadata_expanded", ids)
        self.assertIn("inspector.rerun_run", ids)
        self.assertIn("inspector.invoke_operation", ids)

        result = panel.agent_invoke(
            "inspector.set_metadata_expanded", {"expanded": True}
        )
        self.assertTrue(result["ok"])
        self.assertTrue(panel._metadata_expanded)

        emitted = []
        panel.action_requested.connect(lambda *args: emitted.append(args))
        result = panel.agent_invoke(
            "inspector.rerun_run",
            {"run_id": self.tree_item.run_id, "parameters": {"boot": 2000}},
        )
        self.assertTrue(result["ok"])
        self.assertEqual(emitted[0][0], "rerun_run")
        self.assertEqual(emitted[0][2]["parameters"]["boot"], 2000)

        result = panel.agent_invoke("inspector.show_item", {"item_id": self.tree_item.id})
        self.assertTrue(result["ok"])
        self.assertEqual(panel.current_item_id, self.tree_item.id)

        result = panel.agent_invoke(
            "inspector.rename_run",
            {"run_id": self.tree_item.run_id, "name": "agent renamed"},
        )
        self.assertTrue(result["ok"])
        self.assertEqual(self.store.get_run(self.tree_item.run_id).name, "agent renamed")


if __name__ == "__main__":
    unittest.main(verbosity=2)
