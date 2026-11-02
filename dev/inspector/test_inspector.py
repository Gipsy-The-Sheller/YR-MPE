"""Tests for the standalone InspectorPanel (offscreen Qt)."""

from __future__ import annotations

import os
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

    def test_show_item(self):
        panel = self.panel
        panel.show_item(self.tree_item.id)
        self.assertEqual(panel.title_label.text(), "Tree (GTR+F+I+G4)")
        self.assertIn("Tree", panel.subtitle_label.text())
        self.assertEqual(panel.preview_button.text(), "Open in IcyTree")
        self.assertTrue(panel.preview_button.isEnabled())
        self.assertGreater(panel.meta_table.rowCount(), 0)
        self.assertTrue(panel.run_group.isVisible() or panel.run_group.isVisibleTo(panel))
        self.assertEqual(panel.run_combo.count(), 2)
        self.assertEqual(panel.run_name_edit.text(), "IQ-TREE run 1")
        self.assertGreater(panel.workflow_tree.topLevelItemCount(), 0)
        self.assertTrue(panel.workflow_tree.topLevelItem(0).isExpanded())

    def test_sequence_has_no_run_switcher(self):
        panel = self.panel
        sequence = next(
            item for item in self.store.index.items.values()
            if item.type == "sequence"
        )
        panel.show_item(sequence.id)
        self.assertEqual(panel.preview_button.text(), "Open sequences")
        self.assertFalse(panel.run_group.isVisible())

    def test_rename_run(self):
        panel = self.panel
        panel.show_item(self.tree_item.id)
        renamed = []
        panel.run_renamed.connect(lambda run_id, name: renamed.append((run_id, name)))
        panel.run_name_edit.setText("IQ-TREE final")
        panel._on_rename_run()
        self.assertEqual(self.store.get_run(self.tree_item.run_id).name, "IQ-TREE final")
        self.assertEqual(len(renamed), 1)
        self.assertEqual(renamed[0][1], "IQ-TREE final")

    def test_show_dataset(self):
        panel = self.panel
        panel.show_dataset(self.tree_item.dataset_id)
        self.assertEqual(panel.title_label.text(), "OrthoFinder_plugin_full")
        self.assertFalse(panel.preview_button.isEnabled())
        self.assertFalse(panel.run_group.isVisible())
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
