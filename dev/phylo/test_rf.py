"""Tests for the pairwise Robinson-Foulds comparison."""

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

from YR_MPE.plugins.rf_plugin import (  # noqa: E402
    RFComparisonThread, classify_pairs, count_newick_trees, normalize_tree_file,
    normalized_rf, pairwise_rows, parse_rfdist,
)

APP = QApplication.instance() or QApplication(sys.argv)
IQTREE = REPO / "YR_MPE" / "bin" / "iqtree3" / "iqtree3.exe"


class RFHelperTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.rfdist = self.root / "rf.rfdist"
        self.rfdist.write_text(
            "4 4\n"
            "Tree0 0 2 0 0\n"
            "Tree1 2 0 2 2\n"
            "Tree2 0 2 0 0\n"
            "Tree3 0 2 0 0\n",
            encoding="utf-8",
        )

    def test_parse_rfdist(self):
        labels, matrix = parse_rfdist(str(self.rfdist))
        self.assertEqual(labels, ["Tree0", "Tree1", "Tree2", "Tree3"])
        self.assertEqual(len(matrix), 4)
        self.assertEqual(matrix[0][1], 2)
        self.assertEqual(matrix[1][2], 2)
        self.assertEqual(matrix[0][0], 0)

    def test_pairwise_rows_use_lower_triangle(self):
        labels, matrix = parse_rfdist(str(self.rfdist))
        rows = pairwise_rows(labels, matrix)
        self.assertEqual(len(rows), 6)
        self.assertIn(("Tree0", "Tree1", 2), rows)
        self.assertIn(("Tree2", "Tree3", 0), rows)

    def test_normalized_rf(self):
        self.assertEqual(normalized_rf(2, 4), 1.0)
        self.assertEqual(normalized_rf(1, 4), 0.5)
        self.assertEqual(normalized_rf(2, 3), 1.0)
        self.assertIsNone(normalized_rf(0, 2))

    def test_malformed_rfdist_raises(self):
        bad = self.root / "bad.rfdist"
        bad.write_text("3 3\nTree0 0 1\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            parse_rfdist(str(bad))

    def test_count_newick_trees(self):
        trees = self.root / "trees.nwk"
        trees.write_text(
            "# comment\n(A,B);\n\n(A,C);\n", encoding="utf-8"
        )
        self.assertEqual(count_newick_trees(str(trees)), 2)

    def test_classify_pairs(self):
        groups = [
            {"group_id": "a", "group_name": "A", "kind": "treeset"},
            {"group_id": "a", "group_name": "A", "kind": "treeset"},
            {"group_id": "b", "group_name": "B", "kind": "treeset"},
            {"group_id": "t", "group_name": "T", "kind": "tree"},
        ]
        self.assertEqual(
            classify_pairs(groups),
            [
                "Within tree set", "Between tree sets", "Between objects",
                "Between tree sets", "Between objects", "Between objects",
            ],
        )
        self.assertEqual(
            classify_pairs([None, None]), ["All pairs"]
        )

    def test_selected_tree_payload_preserves_tree_set_groups(self):
        from types import SimpleNamespace

        from YR_MPE.platforms.methods.dataset_models import (
            ITEM_TYPE_PHYLOGENY, ITEM_TYPE_TREESET,
        )
        from YR_MPE.platforms.mpea import YR_MPEA_Widget

        tree_set = SimpleNamespace(
            item_type=ITEM_TYPE_TREESET,
            id="set-1",
            loci_name="Set A",
            data={"data": [
                {"content": "(A,B,(C,D));"},
                {"content": "(A,(B,C),D);"},
            ]},
        )
        tree = SimpleNamespace(
            item_type=ITEM_TYPE_PHYLOGENY,
            id="tree-1",
            loci_name="Species tree",
            data={"content": "(A,B,(C,D));"},
        )
        payload = YR_MPEA_Widget._selected_trees_payload(
            object(), [tree_set, tree]
        )
        self.assertIsNotNone(payload)
        self.assertEqual(
            payload["labels"],
            ["Set A #1", "Set A #2", "Species tree"],
        )
        self.assertEqual(
            [group["kind"] for group in payload["groups"]],
            ["treeset", "treeset", "tree"],
        )
        self.assertEqual(payload["groups"][0]["group_id"], "treeset:set-1")
        Path(payload["tree_file"]).unlink(missing_ok=True)

    def test_normalize_tree_file_unroots_rooted_tree(self):
        source = self.root / "rooted.nwk"
        target = self.root / "unrooted.nwk"
        source.write_text(
            "(A,(B,(C,D)));\n(A,B,(C,D));\n", encoding="utf-8"
        )
        summary = normalize_tree_file(str(source), str(target))
        self.assertEqual(summary["parsed"], 2)
        self.assertEqual(summary["unrooted"], 1)
        self.assertEqual(summary["failed"], 0)
        from io import StringIO

        from Bio import Phylo

        with open(target, "r", encoding="utf-8") as handle:
            first = Phylo.read(StringIO(handle.readline()), "newick")
        self.assertEqual(len(first.root.clades), 3)


@unittest.skipUnless(IQTREE.is_file(), "IQ-TREE 3 runtime not available")
class RFIntegrationTests(unittest.TestCase):
    def test_iqtree_rfdist_roundtrip(self):
        root = Path(tempfile.mkdtemp())
        trees = root / "trees.nwk"
        trees.write_text(
            "(A,B,(C,D));\n"
            "(A,(B,C),D);\n"
            "((A,B),C,D);\n"
            "(A,B,(C,D));\n"
            "(A,(B,(C,D)));\n",
            encoding="utf-8",
        )
        groups = [
            {"group_id": "a", "group_name": "A", "kind": "treeset"},
            {"group_id": "a", "group_name": "A", "kind": "treeset"},
            {"group_id": "b", "group_name": "B", "kind": "treeset"},
            {"group_id": "b", "group_name": "B", "kind": "treeset"},
            {"group_id": "t", "group_name": "T", "kind": "tree"},
        ]
        payloads = []
        thread = RFComparisonThread(
            str(IQTREE), str(trees), ["T1", "T2", "T3", "T4", "T5"],
            str(root / "out"), threads=1, groups=groups,
        )
        thread.result_ready.connect(payloads.append)
        thread.execute_commands()
        self.assertEqual(len(payloads), 1)
        payload = payloads[0]
        self.assertEqual(payload["labels"], ["T1", "T2", "T3", "T4", "T5"])
        self.assertEqual(len(payload["matrix"]), 5)
        self.assertEqual(payload["matrix"][0][1], 2)
        self.assertEqual(payload["matrix"][0][2], 0)
        # The rooted fifth tree is the same unrooted topology as tree 1.
        self.assertEqual(payload["matrix"][0][4], 0)
        self.assertEqual(payload["ntaxa"], 4)
        self.assertEqual(len(payload["groups"]), 5)
        self.assertTrue(Path(payload["rfdist"]).is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
