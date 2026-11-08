"""Tests for the uncorrected p-distance plugin migration."""

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

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.plugins.uncorrected_distance_plugin import (  # noqa: E402
    UncorrectedDistancePlugin, UncorrectedDistanceThread,
)

APP = QApplication.instance() or QApplication(sys.argv)


def records(*ids):
    return [SeqRecord(Seq("ATGCAA"), id=name) for name in ids]


class UncorrectedDistanceTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = UncorrectedDistancePlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "Built-in Python Implementation")
        self.assertIn("nei2000molecular", self.plugin.citation_keys)
        self.assertEqual(self.plugin.run_button.text(), "Calculate p-distance")

    def test_import_once(self):
        self.plugin.close()
        self.plugin = UncorrectedDistancePlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)

    def test_result_spec_from_memory(self):
        self.plugin.calculated_distance_results = [{
            "filename": "locus.dist",
            "content": "2\nA 0 0.1\nB 0.1 0\n",
            "matrix": [[0.0, 0.1], [0.1, 0.0]],
            "sequence_names": ["A", "B"],
        }]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "p_distance")
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "distance")
        self.assertEqual(item.data["labels"], ["A", "B"])
        self.assertTrue(Path(item.file_path).is_file())


class UncorrectedDistanceThreadTests(unittest.TestCase):
    def test_matrix_pairwise(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">A\nATGCAAATGCAA\n>B\nATGCAAATGTAA\n>C\nATGCAAATGCAT\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = UncorrectedDistanceThread(
            [str(source)],
            {"gap_treatment": "pairwise", "substitution_filter": None},
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        result = payloads[0][0]
        self.assertEqual(result["sequence_names"], ["A", "B", "C"])
        self.assertAlmostEqual(result["matrix"][0][1], 1 / 12, places=6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
