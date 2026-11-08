"""Tests for the MAFFT plugin migration (Phase A + B)."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import SingleGeneWorkspace  # noqa: E402
from YR_MPE.plugins.mafft_plugin import MAFFTPlugin, MAFFTThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
MAFFT = REPO / "YR_MPE" / "bin" / "mafft" / "mafft.bat"


def records(*ids):
    return [SeqRecord(Seq("ATGCAAATGCAA"), id=name) for name in ids]


class MAFFTPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = MAFFTPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata_and_parameters(self):
        self.assertEqual(self.plugin.tool_name, "MAFFT")
        self.assertIn("katoh2013mafft", self.plugin.citation_keys)
        self.assertTrue(self.plugin.tool_path)
        self.assertIn("--auto", self.plugin.get_parameters())
        self.assertEqual(self.plugin.mode_combo.currentIndex(), 0)

    def test_analysis_finished_does_not_touch_removed_report_combo(self):
        """The native preview removed report_combo; finishing must not abort."""
        from YR_MPE.plugins import mafft_plugin as module
        with mock.patch.object(module.QMessageBox, "information", return_value=None):
            self.plugin.analysis_finished([], [])
        self.assertEqual(self.plugin.alignment_output_files, [])

    def test_import_list_is_handled_once(self):
        self.plugin.close()
        self.plugin = MAFFTPlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)
        self.assertFalse(self.plugin.batch_mode)
        text = Path(self.plugin.imported_files[0]).read_text(encoding="utf-8")
        self.assertIn(">a", text)
        self.assertIn(">b", text)

    def test_import_batch_list(self):
        self.plugin.close()
        self.plugin = MAFFTPlugin(
            import_from="YR_MPEA",
            import_data=[records("a", "b"), records("c", "d")],
            workdir=str(self.root),
        )
        self.assertTrue(self.plugin.batch_mode)
        self.assertEqual(len(self.plugin.imported_files), 2)

    def test_import_dataset_items(self):
        self.plugin.close()
        payload = {
            "dataset_items": [
                SimpleNamespace(file_path="", sequences=records("a", "b")),
                SimpleNamespace(file_path="", sequences=records("c", "d")),
            ]
        }
        self.plugin = MAFFTPlugin(
            import_from="DATASET_MANAGER", import_data=payload,
            workdir=str(self.root),
        )
        self.assertTrue(self.plugin.batch_mode)
        self.assertEqual(len(self.plugin.imported_files), 2)
        self.assertIn("2 files", self.plugin.file_path_edit.text())

    def test_result_spec_from_output(self):
        output = self.root / "locus_mafft.fas"
        output.write_text(
            ">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8"
        )
        self.plugin.alignment_output_files = [str(output)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "mafft")
        dataset = result.datasets[0]
        self.assertEqual(dataset.source, "mafft")
        self.assertEqual(dataset.summary["alignments"], 1)
        self.assertEqual(dataset.summary["sequences"], 2)
        item = dataset.items[0]
        self.assertEqual(item.item_type, "alignment")
        self.assertEqual(item.file_path, str(output))
        self.assertTrue(item.is_aligned)
        self.assertEqual(item.sequence_count, 2)

    def test_reversed_direction_prefix_is_removed(self):
        output = self.root / "reversed.fas"
        output.write_text(
            ">_R_a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8"
        )
        records, _ = self.plugin._read_alignment(output)
        self.assertEqual([record.id for record in records], ["a", "b"])
        self.assertEqual([record.description for record in records], ["a", "b"])

    def test_result_spec_rejects_ragged_alignment(self):
        output = self.root / "bad.fas"
        output.write_text(
            ">a\nATGCAA\n>b\nATGC\n", encoding="utf-8"
        )
        self.plugin.alignment_output_files = [str(output)]
        with self.assertRaises(ValueError):
            self.plugin._build_result_spec()

    def test_alignment_result_updates_source_dataset(self):
        """MAFFT alignments replace their source items, not a new dataset."""
        from YR_MPE.platforms.methods.dataset_models import DatasetItem
        from YR_MPE.platforms.methods.dataset_selection_manager import (
            DatasetSelectionManager,
        )
        from YR_MPE.platforms.mpea import SingleGeneWorkspace

        manager = DatasetSelectionManager(disable_auto_save=True)
        workspace = SingleGeneWorkspace(dataset_selection_manager=manager)
        dataset_id = manager.create_dataset("OrthoFinder single-copy OGs")

        inputs = []
        for name in ("OG1", "OG2"):
            source = self.root / ("%s.cds.fna" % name)
            source.write_text(">sp1\nATGCAA\n>sp2\nATGCAA\n", encoding="utf-8")
            inputs.append(source)
            item = DatasetItem(item_type="sequence")
            item.loci_name = name
            item.name = name
            item.file_path = str(source)
            item.format = "fasta"
            item.sequences = records("sp1", "sp2")
            manager.add_item(item, dataset_id)

        outputs = []
        for source in inputs:
            output = self.root / (source.name + "_mafft.fas")
            output.write_text(">sp1\nATGCAA\n>sp2\nATGCAA\n", encoding="utf-8")
            outputs.append(str(output))

        plugin = MAFFTPlugin(
            import_from="DATASET_MANAGER",
            import_data={
                "dataset_items": manager.get_items_by_dataset(dataset_id),
                "dataset_config": {"dataset_id": dataset_id},
            },
            workdir=str(self.root),
        )
        plugin.alignment_output_files = outputs
        try:
            result = plugin._build_result_spec()
        finally:
            plugin.close()

        self.assertEqual([item.loci_name for item in result.datasets[0].items],
                         ["OG1", "OG2"])
        self.assertEqual(result.datasets[0].dataset_id, dataset_id)
        self.assertEqual(workspace.add_result_spec(result), [dataset_id])

        items = manager.get_items_by_dataset(dataset_id)
        self.assertEqual(len(items), 2)
        self.assertEqual([item.loci_name for item in items], ["OG1", "OG2"])
        for item in items:
            self.assertEqual(item.item_type, "alignment")
            self.assertTrue(item.is_aligned)
            self.assertTrue(item.file_path.endswith("_mafft.fas"))

    def test_host_imports_result_spec(self):
        output = self.root / "host.fas"
        output.write_text(
            ">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8"
        )
        self.plugin.alignment_output_files = [str(output)]
        result = self.plugin._build_result_spec()

        manager = DatasetSelectionManager(disable_auto_save=True)
        workspace = SingleGeneWorkspace(dataset_selection_manager=manager)
        dataset_ids = workspace.add_result_spec(result)
        self.assertEqual(len(dataset_ids), 1)
        items = manager.get_items_by_dataset(dataset_ids[0])
        self.assertEqual(len(items), 1)
        self.assertEqual(items[0].item_type, "alignment")
        self.assertEqual(items[0].sequence_count, 2)

    def test_file_tag_remove(self):
        path = self.root / "input.fas"
        path.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        self.plugin.add_file_tag(str(path))
        self.assertEqual(self.plugin.imported_files, [str(path)])
        tag_widget = self.plugin.file_tags[0][1]
        self.plugin.remove_file_tag(str(path), tag_widget)
        self.assertEqual(self.plugin.imported_files, [])
        self.assertFalse(self.plugin.file_tags)


@unittest.skipUnless(MAFFT.is_file(), "MAFFT runtime not available")
class MAFFTThreadIntegrationTests(unittest.TestCase):
    def test_mafft_writes_aligned_fasta(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGCAAATGCAA\n>b\nATGCAAATGCAA\n>c\nATGCAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = MAFFTThread(
            str(MAFFT), [str(source)], ["--auto"],
            workdir=str(root), suffix="_mafft",
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()

        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        output_files = payloads[0]["output_files"]
        self.assertEqual(len(output_files), 1)
        output = Path(output_files[0])
        self.assertTrue(output.is_file())
        self.assertGreater(output.stat().st_size, 0)
        self.assertIn(">a", output.read_text(encoding="utf-8").lower())


if __name__ == "__main__":
    unittest.main(verbosity=2)
