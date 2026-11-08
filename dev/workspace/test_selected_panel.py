"""Tests for the multi-object Selected basket and recommendations."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PyQt5.QtCore import QMimeData, Qt
from PyQt5.QtTest import QSignalSpy
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.methods.dataset_item_button import ITEM_MIME_TYPE  # noqa: E402
from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_ALIGNMENT, ITEM_TYPE_DISTANCE,
    ITEM_TYPE_PHYLOGENY, ITEM_TYPE_SEQUENCE, ITEM_TYPE_TREESET,
)
from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.methods.selected_panel import (  # noqa: E402
    SelectedDropButton, SelectedPanel,
)

APP = QApplication.instance() or QApplication(sys.argv)


def make_item(manager, dataset_id, item_type, name, data=None, **attrs):
    item = DatasetItem(item_type=item_type)
    item.dataset_id = dataset_id
    item.loci_name = name
    item.name = name
    item.data = data or {}
    for key, value in attrs.items():
        setattr(item, key, value)
    assert manager.add_item(item, dataset_id)
    return manager.get_item(item.id)


class SelectedPanelTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.manager = DatasetSelectionManager(disable_auto_save=True)
        self.manager.workspace_path = str(self.root)
        self.dataset_id = self.manager.create_dataset("selected-test")
        self.treeset = make_item(
            self.manager, self.dataset_id, ITEM_TYPE_TREESET, "Gene trees",
            {"data": [{"content": "(A,B);"}, {"content": "(A,C);"}]},
        )
        self.tree = make_item(
            self.manager, self.dataset_id, ITEM_TYPE_PHYLOGENY, "Species tree",
            {"newick": "(A,B,C);"},
        )
        self.alignment = make_item(
            self.manager, self.dataset_id, ITEM_TYPE_ALIGNMENT, "Locus",
            sequences=[SeqRecord(Seq("ACGT"), id="A")],
            sequence_count=1, is_aligned=True,
        )
        self.panel = SelectedPanel(self.manager)

    def tearDown(self):
        self.panel.close()

    def _group_labels(self):
        root = self.panel.tree.invisibleRootItem()
        return [root.child(i).text(0) for i in range(root.childCount())]

    def _basket_names(self):
        names = []
        root = self.panel.tree.invisibleRootItem()
        for i in range(root.childCount()):
            group = root.child(i)
            for j in range(group.childCount()):
                names.append(group.child(j).text(0))
        return names

    def _send_labels(self):
        labels = []
        root = self.panel.send_tree.invisibleRootItem()
        for i in range(root.childCount()):
            group = root.child(i)
            for j in range(group.childCount()):
                labels.append(group.child(j).text(0))
        return labels

    def test_basket_groups_by_type(self):
        self.assertTrue(self.manager.add_to_basket(self.treeset.id))
        self.assertTrue(self.manager.add_to_basket(self.tree.id))
        self.panel.refresh()
        self.assertEqual(self._basket_names(), ["Species tree", "Gene trees"])
        self.assertEqual(self._group_labels(), ["Trees (1)", "Tree sets (1)"])

        self.manager.remove_from_basket(self.treeset.id)
        self.panel.refresh()
        self.assertEqual(self._basket_names(), ["Species tree"])
        self.manager.clear_basket()
        self.panel.refresh()
        self.assertEqual(self._basket_names(), [])

    def test_tree_and_treeset_recommends_treerecs(self):
        self.manager.add_to_basket(self.treeset.id)
        self.manager.add_to_basket(self.tree.id)
        self.panel.refresh()
        labels = self._send_labels()
        self.assertIn("TreeRecs 1.2", labels)
        self.assertIn("Robinson-Foulds comparison", labels)

    def test_treeset_recommends_astral(self):
        self.manager.add_to_basket(self.treeset.id)
        self.panel.refresh()
        labels = self._send_labels()
        self.assertIn("ASTRAL-IV", labels)
        self.assertIn("wASTRAL", labels)
        self.assertIn("Robinson-Foulds comparison", labels)
        self.assertNotIn("TreeRecs 1.2", labels)
        self.assertNotIn("LSD2", labels)

    def test_tree_recommends_lsd2_not_astral(self):
        self.manager.add_to_basket(self.tree.id)
        self.panel.refresh()
        labels = self._send_labels()
        self.assertIn("LSD2", labels)
        self.assertNotIn("ASTRAL-IV", labels)
        self.assertNotIn("wASTRAL", labels)
        self.assertNotIn("TreeRecs 1.2", labels)
        self.assertNotIn("Robinson-Foulds comparison", labels)

    def test_alignment_recommends_iqtree(self):
        self.manager.add_to_basket(self.alignment.id)
        self.panel.refresh()
        labels = self._send_labels()
        self.assertIn("IQ-TREE 3", labels)
        self.assertIn("CASTER-site", labels)
        self.assertIn("MAFFT", labels)
        self.assertIn("TrimAl", labels)
        self.assertIn("ModelTest-NG", labels)
        self.assertIn("ML pairwise distance", labels)
        self.assertNotIn("ASTRAL-IV", labels)
        self.assertNotIn("DecentTree", labels)

    def test_sequence_recommends_aligners_only(self):
        sequence = make_item(
            self.manager, self.dataset_id, ITEM_TYPE_SEQUENCE, "Raw locus",
            {"content": ">A\nACGT\n"},
            sequences=[SeqRecord(Seq("ACGT"), id="A")],
            sequence_count=1, is_aligned=False,
        )
        self.manager.add_to_basket(sequence.id)
        self.panel.refresh()
        labels = self._send_labels()
        self.assertIn("MAFFT", labels)
        self.assertIn("Muscle5", labels)
        self.assertNotIn("IQ-TREE 3", labels)
        self.assertNotIn("TrimAl", labels)

    def test_distance_recommends_decenttree(self):
        distance = make_item(
            self.manager, self.dataset_id, ITEM_TYPE_DISTANCE, "Distances",
            {"content": "4\nA 0 1 2 3\nB 1 0 2 3\nC 2 2 0 3\nD 3 3 3 0\n"},
        )
        self.manager.add_to_basket(distance.id)
        self.panel.refresh()
        labels = self._send_labels()
        self.assertIn("DecentTree", labels)
        self.assertNotIn("IQ-TREE 3", labels)
        self.assertNotIn("ASTRAL-IV", labels)

    def test_rf_needs_two_trees(self):
        self.manager.add_to_basket(self.tree.id)
        self.panel.refresh()
        self.assertNotIn("Robinson-Foulds comparison", self._send_labels())
        second = make_item(
            self.manager, self.dataset_id, ITEM_TYPE_PHYLOGENY, "Gene tree",
            {"newick": "(A,(B,C));"},
        )
        self.manager.add_to_basket(second.id)
        self.panel.refresh()
        self.assertIn("Robinson-Foulds comparison", self._send_labels())
        self.assertIn("ASTRAL-IV", self._send_labels())

    def test_send_tree_single_click_emits(self):
        self.manager.add_to_basket(self.treeset.id)
        self.panel.refresh()
        spy = QSignalSpy(self.panel.send_requested)
        root = self.panel.send_tree.invisibleRootItem()
        tool_item = root.child(0).child(0)
        self.panel.send_tree.itemClicked.emit(
            tool_item, self.panel.send_tree.currentColumn()
        )
        self.assertEqual(len(spy), 1)
        self.assertEqual(spy[0][0], "astral4")
        self.assertEqual(list(spy[0][1]), [self.treeset.id])

    def test_selected_button_accepts_item_drop(self):
        class FakeDrop:
            def __init__(self, mime):
                self._mime = mime
                self.accepted = None

            def mimeData(self):
                return self._mime

            def setDropAction(self, action):
                self.action = action

            def accept(self):
                self.accepted = True

            def ignore(self):
                self.accepted = False

        button = SelectedDropButton()
        spy = QSignalSpy(button.item_dropped)

        mime = QMimeData()
        mime.setData(ITEM_MIME_TYPE, self.tree.id.encode("utf-8"))
        event = FakeDrop(mime)
        button.dropEvent(event)
        self.assertEqual(len(spy), 1)
        self.assertEqual(spy[0][0], self.tree.id)
        self.assertTrue(event.accepted)

        empty = FakeDrop(QMimeData())
        button.dropEvent(empty)
        self.assertEqual(len(spy), 1)
        self.assertFalse(empty.accepted)

    def test_basket_persists_in_legacy_state(self):
        self.manager.disable_auto_save = False
        self.manager.add_to_basket(self.treeset.id)
        self.manager._save_state()
        loaded = DatasetSelectionManager()
        loaded.workspace_path = str(self.root)
        loaded._load_state()
        self.assertIn(self.treeset.id, loaded.basket_items)


if __name__ == "__main__":
    unittest.main(verbosity=2)
