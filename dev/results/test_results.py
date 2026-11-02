"""Tests for the plugin result contract and the workspace importer.

Run from the repository root:

    python dev/results/test_results.py

PyQt runs with the offscreen platform; no display is required.
"""

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

from YR_MPE.platforms.methods.dataset_selection_manager import DatasetSelectionManager  # noqa: E402
from YR_MPE.platforms.mpea import SingleGeneWorkspace  # noqa: E402
from YR_MPE.runtime.results import (  # noqa: E402
    DatasetItemSpec,
    DatasetSpec,
    ResultSpec,
)

APP = QApplication.instance() or QApplication(sys.argv)


def _records(*ids):
    return [SeqRecord(Seq("ATGCAA"), id=name) for name in ids]


class ResultModelTests(unittest.TestCase):
    def test_item_defaults(self):
        item = DatasetItemSpec(loci_name="OG1", sequences=_records("a", "b"))
        self.assertEqual(item.sequence_count, 2)
        self.assertEqual(item.length, 6)
        self.assertEqual(item.validate(), [])

    def test_validation_rejects_incomplete_items(self):
        result = ResultSpec(source="test")
        self.assertTrue(result.validate())
        dataset = DatasetSpec(name="empty")
        self.assertTrue(dataset.validate())
        item = DatasetItemSpec(loci_name="OG1")
        self.assertTrue(item.validate())
        result = ResultSpec(source="test", datasets=[DatasetSpec(name="d", items=[item])])
        self.assertTrue(result.validate())

    def test_settings_dict_is_jsonable(self):
        item = DatasetItemSpec(
            loci_name="OG1",
            sequences=_records("a"),
            data={"path": Path("x/y"), "nested": (1, 2)},
        )
        payload = item.to_settings_dict()
        self.assertEqual(payload["sequence_count"], 1)
        self.assertEqual(payload["data"]["path"], str(Path("x/y")))
        self.assertEqual(payload["data"]["nested"], [1, 2])


class WorkspaceImportTests(unittest.TestCase):
    def setUp(self):
        self.manager = DatasetSelectionManager(disable_auto_save=True)
        self.workspace = SingleGeneWorkspace(dataset_selection_manager=self.manager)

    def test_import_result_spec(self):
        spec = DatasetSpec(
            name="OrthoFinder_test",
            source="one_step_orthofinder",
            summary={"single_copy_ogs": 2},
            metadata={"input": "gb/"},
            items=[
                DatasetItemSpec(loci_name="OG1", sequences=_records("sp1", "sp2")),
                DatasetItemSpec(loci_name="OG2", sequences=_records("sp1", "sp2")),
            ],
        )
        ids = self.workspace.add_result_spec(ResultSpec(source="test", datasets=[spec]))
        self.assertEqual(len(ids), 1)

        dataset = self.manager.get_dataset(ids[0])
        self.assertEqual(dataset.name, "OrthoFinder_test")
        self.assertEqual(dataset.loci_count, 2)
        self.assertEqual(dataset.taxa_count, 2)
        self.assertEqual(dataset.settings["source"], "one_step_orthofinder")

        items = self.manager.get_items_by_dataset(ids[0])
        self.assertEqual([item.loci_name for item in items], ["OG1", "OG2"])
        self.assertEqual([item.item_type for item in items], ["sequence", "sequence"])
        self.assertEqual(len(dataset.settings["dataset_items"]), 2)
        self.assertIn("2 loci", self.workspace.item_buttons["datasets"][0].toolTip())

    def test_legacy_dataset_adapter(self):
        class LegacyItem:
            def __init__(self, name):
                self.loci_name = name
                self.sequences = _records("x")
                self.length = 6
                self.sequence_count = 1
                self.is_aligned = True
                self.file_path = ""
                self.selected = False
                self.name = name

        class LegacyDataset:
            dataset_name = "Legacy"
            items = [LegacyItem("gene1")]
            partition_scheme = None
            summary = {}

        ids = self.workspace.add_result_spec(LegacyDataset())
        items = self.manager.get_items_by_dataset(ids[0])
        self.assertEqual([(item.loci_name, item.item_type) for item in items], [("gene1", "alignment")])

    def test_unique_dataset_names(self):
        spec = DatasetSpec(name="same", items=[DatasetItemSpec(loci_name="x", sequences=_records("a"))])
        first = self.workspace.add_result_spec(spec)[0]
        second = self.workspace.add_result_spec(spec)[0]
        self.assertNotEqual(first, second)
        self.assertEqual(self.manager.get_dataset(second).name, "same (2)")


class OrthoFinderResultTests(unittest.TestCase):
    def test_build_result_spec_from_output(self):
        from YR_MPE.plugins.orthofinder_plugin import OrthoFinderPlugin

        class TestPlugin(OrthoFinderPlugin):
            def config(self):
                return True

        root = Path(tempfile.mkdtemp())
        single = root / "single_copy_ogs"
        single.mkdir(parents=True)
        (single / "OG0000001.cds.fna").write_text(
            ">sp1_1\nATGCAA\n>sp2_1\nATGCAA\n", encoding="utf-8"
        )
        (single / "OG0000002.cds.fna").write_text(
            ">sp1_2\nATGAAA\n>sp2_2\nATGAAA\n", encoding="utf-8"
        )
        (single / "OG0000001.faa").write_text(">sp1_1\nMQ\n", encoding="utf-8")
        (single / "OG0000001.map.tsv").write_text("orthogroup\tspecies\n", encoding="utf-8")
        (single / "single_copy_og_summary.tsv").write_text(
            "orthogroup\tspecies_count\tprotein_count\tstatus\n", encoding="utf-8"
        )

        plugin = TestPlugin(import_from="YR_MPEA", workdir=str(root))
        plugin.input_edit.setText(str(root / "gb"))
        plugin.output_edit.setText(str(root))
        plugin._last_summary = {
            "species": 2,
            "cds": 4,
            "orthogroups": 2,
            "single_copy_ogs": 2,
            "output_dir": str(root),
            "elapsed": 1.0,
            "run_mode": "full",
            "cds_mapping": True,
        }
        result = plugin._build_result_spec()
        plugin.close()

        self.assertEqual(result.schema, "yrmpe.result.v1")
        self.assertEqual(result.validate(), [])
        dataset = result.datasets[0]
        self.assertEqual(len(dataset.items), 2)
        self.assertEqual(dataset.items[0].loci_name, "OG0000001")
        self.assertEqual(dataset.items[0].item_type, "sequence")
        self.assertEqual(dataset.items[0].sequence_count, 2)
        self.assertTrue(dataset.items[0].file_path.endswith(".cds.fna"))
        self.assertEqual(dataset.items[0].data["orthogroup"], "OG0000001")


if __name__ == "__main__":
    unittest.main(verbosity=2)
