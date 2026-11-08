"""Regression tests for Dataset Manager -> workspace alignment handoff."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.platforms.methods.dataset_manager import (  # noqa: E402
    DatasetItem as LocalItem, DatasetManager,
)
from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, SELECTION_STATE_GREEN, SELECTION_STATE_NONE,
)
from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import YR_MPEA_Widget, SingleGeneWorkspace  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


def records(*ids):
    return [SeqRecord(Seq("ATGCAA"), id=name) for name in ids]


class DatasetManagerSyncTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.manager = DatasetSelectionManager(disable_auto_save=True)
        self.workspace = SingleGeneWorkspace(
            dataset_selection_manager=self.manager
        )
        self.dataset_id = self.manager.create_dataset("OrthoFinder")
        source = self.root / "OG1.cds.fna"
        source.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        item = DatasetItem(item_type="sequence")
        item.loci_name = "OG1"
        item.name = "OG1"
        item.file_path = str(source)
        item.format = "fasta"
        item.sequences = records("a", "b")
        item.sequence_count = 2
        item.length = 6
        self.manager.add_item(item, self.dataset_id)
        self.workspace.current_dataset_id = self.dataset_id
        self.manager.get_dataset(self.dataset_id).selection_state = (
            SELECTION_STATE_GREEN
        )
        self.dialog = DatasetManager(workspace=self.workspace)

    def tearDown(self):
        self.dialog.close()

    def test_batch_alignment_updates_selection_manager_item(self):
        local = LocalItem()
        local.loci_name = "OG1"
        local.name = "OG1"
        local.file_path = str(self.root / "OG1.cds.fna")
        local.format = "fasta"
        local.sequences = records("a", "b")
        local.sequence_count = 2
        local.length = 6
        local.selected = True
        self.dialog.dataset_items = [local]
        self.dialog.add_dataset_to_table(local)

        output = self.root / "OG1.cds.fna_mafft.fas"
        output.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")

        class FakePlugin:
            alignment_output_files = [str(output)]

        self.dialog._active_alignment_plugin = FakePlugin()
        with mock.patch(
            "YR_MPE.platforms.methods.dataset_manager.QMessageBox.information",
            return_value=None,
        ):
            self.dialog.handle_batch_alignment_results([records("a", "b")])

        target = self.manager.get_items_by_dataset(self.dataset_id)[0]
        self.assertEqual(target.item_type, "alignment")
        self.assertTrue(target.is_aligned)
        self.assertEqual(target.file_path, str(output))
        self.assertEqual(target.sequence_count, 2)

        # The method lives on the host widget; bind it to a tiny stub.
        class Host:
            pass

        host = Host()
        host.dataset_selection_manager = self.manager
        host.workspace = self.workspace
        prepare = YR_MPEA_Widget._prepare_import_data.__get__(host)
        import_from, payload = prepare(workspace_item_types=["alignments"])
        self.assertEqual(import_from, "DATASET_MANAGER")
        self.assertEqual(
            [item.item_type for item in payload["dataset_items"]],
            ["alignment"],
        )

    def test_current_dataset_is_used_without_green_state(self):
        self.manager.get_dataset(self.dataset_id).selection_state = (
            SELECTION_STATE_NONE
        )
        self.workspace.current_dataset_id = self.dataset_id

        class Host:
            pass

        host = Host()
        host.dataset_selection_manager = self.manager
        host.workspace = self.workspace
        prepare = YR_MPEA_Widget._prepare_import_data.__get__(host)
        import_from, payload = prepare(workspace_item_types=["alignments"])
        self.assertEqual(import_from, "DATASET_MANAGER")
        self.assertEqual(payload["dataset_config"]["dataset_id"], self.dataset_id)


if __name__ == "__main__":
    unittest.main(verbosity=2)
