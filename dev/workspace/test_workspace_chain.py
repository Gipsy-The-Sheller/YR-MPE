"""Workspace chain feedback and inspector-animation tests (offscreen Qt)."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_MODEL, ITEM_TYPE_DISTANCE, ITEM_TYPE_PHYLOGENY,
    SELECTION_STATE_GREEN, SELECTION_STATE_NONE,
)
from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import SingleGeneWorkspace, YR_MPEA_Widget  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


def make_item(manager, dataset_id, chain_uuid, name, item_type=ITEM_TYPE_MODEL,
              parent_id=None):
    item = DatasetItem(item_type=item_type)
    item.dataset_id = dataset_id
    item.loci_name = name
    item.name = name
    item.ownership_uuid = chain_uuid
    item.parent_id = parent_id
    item.data = {"name": name}
    assert manager.add_item(item, dataset_id)
    return manager.get_item(item.id)


class WorkspaceChainTests(unittest.TestCase):
    def setUp(self):
        self.manager = DatasetSelectionManager(disable_auto_save=True)
        self.dataset_id = self.manager.create_dataset("chain-test")
        self.workspace = SingleGeneWorkspace(dataset_selection_manager=self.manager)
        self.workspace.current_dataset_id = self.dataset_id
        self.a1 = make_item(self.manager, self.dataset_id, "chain-a", "A1")
        self.a2 = make_item(self.manager, self.dataset_id, "chain-a", "A2")
        self.b1 = make_item(self.manager, self.dataset_id, "chain-b", "B1")
        self.buttons = {}
        for index, item in enumerate((self.a1, self.a2, self.b1)):
            button = self.workspace._create_item_button(
                item, "models", 1, index, "file/model.svg"
            )
            self.buttons[item.id] = button

    def tearDown(self):
        self.workspace.close()

    def _mode(self, item):
        return self.buttons[item.id].chain_highlight

    def test_hover_highlights_only_chain(self):
        self.workspace._on_item_hovered(self.a1.id)
        self.assertEqual(self._mode(self.a1), "hover")
        self.assertEqual(self._mode(self.a2), "hover")
        self.assertIsNone(self._mode(self.b1))
        self.workspace._clear_hover_chain()
        self.assertIsNone(self._mode(self.a1))
        self.assertIsNone(self._mode(self.a2))

    def test_click_activates_chain_and_clears_previous(self):
        self.workspace._on_item_single_click(self.a1.id)
        # The directly selected item keeps its green selection; its chain
        # context is blue.
        self.assertIsNone(self._mode(self.a1))
        self.assertEqual(self._mode(self.a2), "active")
        self.assertIsNone(self._mode(self.b1))
        self.assertEqual(self.a1.selection_state, SELECTION_STATE_GREEN)

        self.workspace._on_item_single_click(self.b1.id)
        self.assertIsNone(self._mode(self.a1))
        self.assertIsNone(self._mode(self.a2))
        self.assertIsNone(self._mode(self.b1))
        self.assertEqual(self.b1.selection_state, SELECTION_STATE_GREEN)

    def test_hover_does_not_erase_active_chain(self):
        self.workspace._on_item_single_click(self.a1.id)
        self.workspace._on_item_hovered(self.b1.id)
        self.assertIsNone(self._mode(self.a1))
        self.assertEqual(self._mode(self.a2), "active")
        self.assertEqual(self._mode(self.b1), "hover")
        self.workspace._clear_hover_chain()
        self.assertIsNone(self._mode(self.a1))
        self.assertIsNone(self._mode(self.b1))

    def test_button_click_emits_activation_immediately(self):
        seen = []
        self.buttons[self.a1.id].activated.connect(seen.append)
        self.buttons[self.a1.id].click()
        self.assertEqual(seen, [self.a1.id])

    def test_hover_active_chain_restores_blue_on_leave(self):
        self.workspace._on_item_single_click(self.a1.id)
        self.workspace._on_item_hovered(self.a2.id)
        self.assertIsNone(self._mode(self.a1))
        self.assertEqual(self._mode(self.a2), "hover")
        self.workspace._clear_hover_chain()
        self.assertIsNone(self._mode(self.a1))
        self.assertEqual(self._mode(self.a2), "active")

    def test_only_clicked_item_is_green_even_with_green_ancestor(self):
        parent = make_item(self.manager, self.dataset_id, "chain-c", "Parent model")
        child = make_item(
            self.manager, self.dataset_id, "chain-c", "Child distance",
            item_type=ITEM_TYPE_DISTANCE, parent_id=parent.id,
        )
        parent_button = self.workspace._create_item_button(
            parent, "models", 1, 8, "file/model.svg"
        )
        child_button = self.workspace._create_item_button(
            child, "distances", 2, 8, "file/distance.svg"
        )
        self.workspace._on_item_single_click(child.id)
        # The selection engine keeps the ancestor green for dependency
        # resolution, but only the clicked item may look green.
        self.assertEqual(child.selection_state, SELECTION_STATE_GREEN)
        self.assertEqual(parent.selection_state, SELECTION_STATE_GREEN)
        self.assertIsNone(child_button.chain_highlight)
        self.assertEqual(parent_button.chain_highlight, "active")

    def test_dataset_has_its_own_grid_row(self):
        button = self.workspace._add_dataset_button(
            self.manager.get_dataset(self.dataset_id)
        )
        index = self.workspace.grid_layout.indexOf(button)
        row, _column, _rowspan, _colspan = self.workspace.grid_layout.getItemPosition(index)
        self.assertEqual(row, 6)

    def test_materialize_dataset_items_is_typed_and_idempotent(self):
        tree = make_item(
            self.manager, self.dataset_id, "chain-c", "tree",
            item_type=ITEM_TYPE_PHYLOGENY,
        )
        self.workspace.materialize_dataset_items(self.dataset_id)
        self.assertEqual(len(self.workspace.item_buttons["models"]), 6)
        self.assertEqual(len(self.workspace.item_buttons["phylogenies"]), 1)
        # A repeated import must not duplicate the workspace icons.
        self.workspace.materialize_dataset_items(self.dataset_id)
        self.assertEqual(len(self.workspace.item_buttons["models"]), 6)
        self.assertIn(tree.id, self.workspace._materialized_item_ids)

    def test_deselect_clears_chain(self):
        self.workspace._on_item_single_click(self.a1.id)
        self.workspace._on_item_deselect(self.a1.id)
        self.assertEqual(self.a1.selection_state, SELECTION_STATE_NONE)
        self.assertEqual(self.a2.selection_state, SELECTION_STATE_NONE)
        self.assertIsNone(self._mode(self.a1))
        self.assertIsNone(self._mode(self.a2))
        self.assertIsNone(self.workspace._active_chain_key)


class InspectorAnimationTests(unittest.TestCase):
    def setUp(self):
        self.window = YR_MPEA_Widget()
        self.window.resize(1200, 800)
        self.window.show()
        QTest.qWait(80)

    def tearDown(self):
        self.window.close()

    def test_close_and_open_animation(self):
        self.window.close_inspector_animated()
        QTest.qWait(350)
        self.assertFalse(self.window.inspector.isVisible())
        self.assertFalse(self.window.inspector_button.isChecked())

        self.window._open_details_panel()
        QTest.qWait(350)
        self.assertTrue(self.window.right_splitter.isVisible())
        self.assertTrue(self.window.inspector_button.isChecked())
        self.assertFalse(self.window.selected_button.isChecked())
        self.assertGreaterEqual(self.window.inspectorWidth, 300)

        self.window._open_selected_panel()
        QTest.qWait(120)
        self.assertTrue(self.window.selected_button.isChecked())
        self.assertFalse(self.window.inspector_button.isChecked())
        self.assertEqual(self.window._active_right_panel, "selected")


if __name__ == "__main__":
    unittest.main(verbosity=2)
