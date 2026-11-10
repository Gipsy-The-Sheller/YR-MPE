"""Reopen a workspace and check that every item comes back.

This exercises the same path the application uses:

    workspace.json -> WorkspaceStore -> DatasetSelectionManager._load_from_store
                   -> SingleGeneWorkspace.rebuild_from_manager()

Run from the repository root:

    python dev/workspace/test_workspace_reopen.py
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMessageBox

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.platforms.methods.dataset_models import (  # noqa: E402
    DatasetItem, ITEM_TYPE_ALIGNMENT, ITEM_TYPE_CHAIN, ITEM_TYPE_CLOCK,
    ITEM_TYPE_COALESCENT, ITEM_TYPE_DISTANCE, ITEM_TYPE_MODEL,
    ITEM_TYPE_PHYLOGENY, ITEM_TYPE_SEQUENCE, ITEM_TYPE_TABLE,
    ITEM_TYPE_TREESET, ITEM_TYPE_VARIANT, SELECTION_STATE_GREEN,
)
from YR_MPE.platforms.methods.dataset_selection_manager import (  # noqa: E402
    DatasetSelectionManager,
)
from YR_MPE.platforms.mpea import SingleGeneWorkspace  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)

FASTA_A = ">a\nATGCAA\n>b\nATGCAA\n>c\nATGAAA\n"
FASTA_ALN = ">a\nATGC-AAA\n>b\nATGC-AAA\n>c\nATG--AAA\n"

STORE_MODULE = REPO / "YR_MPE" / "runtime" / "workspace_store.py"


def records(text: str):
    from io import StringIO
    from Bio import SeqIO

    return list(SeqIO.parse(StringIO(text), "fasta"))


def make_item(item_type: str, name: str, dataset_id: str, **kwargs) -> DatasetItem:
    item = DatasetItem(item_type=item_type)
    item.dataset_id = dataset_id
    item.loci_name = name
    item.name = name
    for key, value in kwargs.items():
        setattr(item, key, value)
    return item


class WorkspaceReopenTests(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="yrmpe-reopen-")).resolve()
        (self.root / "inputs").mkdir()
        self.fasta = self.root / "inputs" / "OG0001.fna"
        self.fasta.write_text(FASTA_A, encoding="utf-8")
        self.alignment = self.root / "inputs" / "OG0001.aln.fna"
        self.alignment.write_text(FASTA_ALN, encoding="utf-8")
        self.treefile = self.root / "inputs" / "OG0001.treefile"
        self.treefile.write_text("(a:0.1,b:0.2,c:0.3);\n", encoding="utf-8")

        self.build_workspace()
        self.reopened = DatasetSelectionManager()
        self.reopened.set_workspace_path(str(self.root))

    # ------------------------------------------------------------------
    # fixture
    # ------------------------------------------------------------------
    def build_workspace(self):
        manager = DatasetSelectionManager()
        manager.set_workspace_path(str(self.root))

        self.ds_ortho = manager.create_dataset(
            "OrthoFinder", description="One-Step OrthoFinder import"
        )
        manager.get_dataset(self.ds_ortho).settings.update(
            {"source": "one_step_orthofinder", "genetic_code": 5}
        )
        self.ds_bayes = manager.create_dataset("MrBayes")

        self.import_record = manager.create_import(
            [str(self.fasta)],
            "one_step_orthofinder",
            parameters={"genetic_code": 5, "min_sequences": 4},
            dataset_id=self.ds_ortho,
        )
        self.run_record = manager.create_run(
            "mafft", "MAFFT", [], name="OG0001 MAFFT",
            tool_version="7.520", parameters={"auto": True, "adjustdirection": True},
        )

        a = self.ds_ortho
        self.seq = make_item(
            ITEM_TYPE_SEQUENCE, "OG0001", a,
            file_path=str(self.fasta), format="fasta",
            length=6, sequence_count=3, sequences=records(FASTA_A),
            metadata={"input": str(self.fasta), "sha256": "sha256:deadbeef"},
        )
        manager.add_item(self.seq, a)

        self.aln = make_item(
            ITEM_TYPE_ALIGNMENT, "OG0001 alignment", a,
            file_path=str(self.alignment), format="fasta",
            length=8, sequence_count=3, is_aligned=True,
            sequences=records(FASTA_ALN),
            parent_id=self.seq.id,
            import_id=self.import_record.id,
            run_id=self.run_record.id,
            metadata={"tool": "MAFFT", "adjustdirection": True},
            data={"nucleotide": "DNA"},
        )
        self.aln.ownership_uuid = "own-alignment"
        manager.add_item(self.aln, a)

        self.model = make_item(
            ITEM_TYPE_MODEL, "GTR+G", a,
            parent_id=self.aln.id,
            data={
                "model": "GTR+G", "type": "single", "partitioned": False,
                "params": {"gamma": 0.5123, "pinv": 0.0},
            },
            metadata={"engine": "ModelFinder"},
        )
        manager.add_item(self.model, a)

        self.part_model = make_item(
            ITEM_TYPE_MODEL, "Partitioned model", a,
            parent_id=self.aln.id,
            model_sub_type="partitioned",
            partition_mode="EL",
            data={
                "type": "partitioned", "partitioned": True,
                "partition_mode": "EL",
                "partitions": [{"name": "p1", "model": "GTR+G"}],
            },
        )
        manager.add_item(self.part_model, a)

        self.distance = make_item(
            ITEM_TYPE_DISTANCE, "OG0001 distance", a,
            parent_id=self.model.id,
            data={"labels": ["a", "b", "c"],
                  "matrix": [[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]],
                  "metric": "JC69"},
        )
        manager.add_item(self.distance, a)

        self.tree = make_item(
            ITEM_TYPE_PHYLOGENY, "OG0001 tree", a,
            parent_id=self.aln.id,
            file_path=str(self.treefile), format="newick",
            data={"newick": "(a:0.1,b:0.2,c:0.3);", "log_likelihood": -1234.5},
            metadata={"support": "ultrafast", "model": "GTR+G"},
        )
        manager.add_item(self.tree, a)

        self.treeset = make_item(
            ITEM_TYPE_TREESET, "OG0001 treeset", a,
            parent_id=self.tree.id,
            data={"count": 150, "trees": ["(a,b);", "(a,c);"]},
        )
        manager.add_item(self.treeset, a)

        self.chain = make_item(
            ITEM_TYPE_CHAIN, "mrbayes_run2", a,
            data={"file_paths": [str(self.root / "inputs" / "run2.p")],
                  "run_number": 2, "chain_count": 2, "tool": "mrbayes",
                  "generation_count": 100000},
        )
        manager.add_item(self.chain, a)

        self.table = make_item(
            ITEM_TYPE_TABLE, "saturation table", a,
            data={"columns": ["taxon", "value"], "rows": [["a", 1], ["b", 2]]},
        )
        manager.add_item(self.table, a)

        self.variant = make_item(
            ITEM_TYPE_VARIANT, "variants", a,
            data={"sites": [{"position": 12, "ref": "A", "alt": "G"}]},
        )
        manager.add_item(self.variant, a)

        self.coalescent = make_item(
            ITEM_TYPE_COALESCENT, "species tree", a,
            data={"species_tree": "((a,b),c);"},
        )
        manager.add_item(self.coalescent, a)

        self.clock = make_item(
            ITEM_TYPE_CLOCK, "clock", a,
            data={"rate": 0.0123, "prior": "lognormal"},
        )
        manager.add_item(self.clock, a)

        self.chain_b = make_item(
            ITEM_TYPE_CHAIN, "phylobayes_run1", self.ds_bayes,
            data={"file_paths": [], "run_number": 1, "chain_count": 1,
                  "tool": "phylobayes"},
        )
        manager.add_item(self.chain_b, self.ds_bayes)

        manager.finish_run(
            self.run_record.id,
            output_ids=[self.aln.id],
            artifacts=[str(self.alignment)],
        )
        self.stream_record = manager.create_stream(
            "OG0001 pipeline", self.ds_ortho, [self.run_record.id]
        )
        manager.select_item(self.aln.id)
        manager.add_to_basket(self.model.id)
        manager.get_dataset(self.ds_ortho).selection_state = SELECTION_STATE_GREEN
        manager._save_state()
        self.manager = manager

    # ------------------------------------------------------------------
    # core restore
    # ------------------------------------------------------------------
    def test_store_load_is_lossless(self):
        """WorkspaceIndex.from_dict -> to_dict must reproduce the file."""
        on_disk = json.loads((self.root / "workspace.json").read_text(encoding="utf-8"))
        script = (
            "import importlib.util, json, sys\n"
            "spec = importlib.util.spec_from_file_location('ws', %r)\n"
            "mod = importlib.util.module_from_spec(spec)\n"
            "sys.modules[spec.name] = mod\n"
            "spec.loader.exec_module(mod)\n"
            "print(json.dumps(mod.WorkspaceIndex.from_dict(json.loads(sys.stdin.read())).to_dict()))\n"
            % str(STORE_MODULE)
        )
        completed = subprocess.run(
            [sys.executable, "-I", "-c", script],
            input=json.dumps(on_disk), capture_output=True, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(completed.stdout), on_disk)

    def test_all_item_types_restore(self):
        expected = {
            self.seq.id: self.seq,
            self.aln.id: self.aln,
            self.model.id: self.model,
            self.part_model.id: self.part_model,
            self.distance.id: self.distance,
            self.tree.id: self.tree,
            self.treeset.id: self.treeset,
            self.chain.id: self.chain,
            self.table.id: self.table,
            self.variant.id: self.variant,
            self.coalescent.id: self.coalescent,
            self.clock.id: self.clock,
            self.chain_b.id: self.chain_b,
        }
        self.assertEqual(set(self.reopened.items), set(expected))
        for item_id, original in expected.items():
            loaded = self.reopened.items[item_id]
            self.assertEqual(loaded.item_type, original.item_type, item_id)
            self.assertEqual(loaded.loci_name, original.loci_name, item_id)
            self.assertEqual(loaded.dataset_id, original.dataset_id, item_id)
            self.assertEqual(loaded.format, original.format, item_id)
            self.assertEqual(loaded.length, original.length, item_id)
            self.assertEqual(loaded.sequence_count, original.sequence_count, item_id)
            self.assertEqual(loaded.is_aligned, original.is_aligned, item_id)
            self.assertEqual(loaded.parent_id, original.parent_id, item_id)
            self.assertEqual(loaded.import_id, original.import_id, item_id)
            self.assertEqual(loaded.run_id, original.run_id, item_id)
            self.assertEqual(loaded.metadata, original.metadata, item_id)
            self.assertEqual(loaded.data, original.data, item_id)

        partitioned = self.reopened.get_item(self.part_model.id)
        self.assertEqual(partitioned.model_sub_type, "partitioned")
        self.assertEqual(partitioned.partition_mode, "EL")
        self.assertEqual(self.reopened.get_item(self.model.id).model_sub_type, "single")

    def test_selection_state_restores(self):
        self.assertEqual(self.reopened.basket_items, {self.model.id})
        self.assertIn(self.model.id, self.reopened.selected_items)
        self.assertIn(self.aln.id, self.reopened.selected_items)

        model = self.reopened.get_item(self.model.id)
        self.assertEqual(model.selection_state, SELECTION_STATE_GREEN)
        self.assertEqual(model.selection_reason, "Selected for workflow")
        aligned = self.reopened.get_item(self.aln.id)
        self.assertEqual(aligned.selection_state, SELECTION_STATE_GREEN)
        self.assertEqual(aligned.selection_reason, "Direct selection")
        self.assertEqual(aligned.ownership_uuid, "own-alignment")
        self.assertEqual(
            self.reopened.get_dataset(self.ds_ortho).selection_state,
            SELECTION_STATE_GREEN,
        )

    def test_sequence_data_hydrates_from_file(self):
        loaded = self.reopened.get_item(self.seq.id)
        self.assertEqual([r.id for r in loaded.sequences], ["a", "b", "c"])
        self.assertEqual([str(r.seq) for r in loaded.sequences],
                         [str(r.seq) for r in self.seq.sequences])
        self.assertTrue(loaded.file_path)
        self.assertTrue(os.path.exists(loaded.file_path))

        loaded_aln = self.reopened.get_item(self.aln.id)
        self.assertEqual(len(loaded_aln.sequences), 3)
        self.assertEqual(str(loaded_aln.sequences[0].seq), "ATGC-AAA")

    def test_datasets_and_order_restore(self):
        self.assertEqual(set(self.reopened.datasets), {self.ds_ortho, self.ds_bayes})
        ortho = self.reopened.get_dataset(self.ds_ortho)
        self.assertEqual(ortho.name, "OrthoFinder")
        self.assertEqual(ortho.description, "One-Step OrthoFinder import")
        self.assertEqual(ortho.settings.get("genetic_code"), 5)
        self.assertEqual(
            ortho.items,
            [self.seq.id, self.aln.id, self.model.id, self.part_model.id,
             self.distance.id, self.tree.id, self.treeset.id, self.chain.id,
             self.table.id, self.variant.id, self.coalescent.id, self.clock.id],
        )

    def test_provenance_restore(self):
        imported = self.reopened.store.get_import(self.import_record.id)
        self.assertEqual(imported.importer, "one_step_orthofinder")
        self.assertEqual(imported.parameters["genetic_code"], 5)
        self.assertIn(self.fasta.name, " ".join(imported.source_hashes))

        run = self.reopened.store.get_run(self.run_record.id)
        self.assertEqual(run.name, "OG0001 MAFFT")
        self.assertEqual(run.tool_version, "7.520")
        self.assertEqual(run.status, "completed")
        self.assertEqual(run.output_ids, [self.aln.id])

        stream = self.reopened.store.get_stream(self.stream_record.id)
        self.assertEqual(stream.run_ids, [self.run_record.id])

    # ------------------------------------------------------------------
    # UI rebuild
    # ------------------------------------------------------------------
    def test_workspace_grid_rebuilds_item_buttons(self):
        view = SingleGeneWorkspace(dataset_selection_manager=self.reopened)
        view.current_dataset_id = self.ds_ortho
        view.rebuild_from_manager()
        counts = {key: len(buttons) for key, buttons in view.item_buttons.items()}
        self.assertEqual(counts["alignments"], 2)
        self.assertEqual(counts["models"], 2)
        self.assertEqual(counts["distances"], 1)
        self.assertEqual(counts["phylogenies"], 1)
        self.assertEqual(counts["treesets"], 1)
        self.assertEqual(counts["chains"], 1)
        self.assertEqual(counts["datasets"], 2)
        names = {button.item.loci_name for button in view.item_buttons["models"]}
        self.assertEqual(names, {"GTR+G", "Partitioned model"})

    def test_dataset_button_switches_canvas(self):
        from PyQt5.QtTest import QTest

        view = SingleGeneWorkspace(dataset_selection_manager=self.reopened)
        view.current_dataset_id = self.ds_ortho
        view.rebuild_from_manager()
        view._on_dataset_click(self.ds_bayes)
        QTest.qWait(20)
        self.assertEqual(view.current_dataset_id, self.ds_bayes)
        self.assertEqual(len(view.item_buttons["chains"]), 1)
        self.assertEqual(len(view.item_buttons["models"]), 0)
        self.assertEqual(len(view.item_buttons["datasets"]), 2)

    def test_full_application_reopen(self):
        """Open the real host widget twice on the same folder."""
        from unittest import mock
        from YR_MPE.platforms.methods.workspace_manager import WorkspaceManager
        from YR_MPE.platforms.mpea import YR_MPEA_Widget

        with mock.patch.object(
            QMessageBox, "information", return_value=None
        ), mock.patch.object(
            QMessageBox, "warning", return_value=None
        ), mock.patch.object(
            QMessageBox, "question", return_value=QMessageBox.Yes
        ), mock.patch.object(
            WorkspaceManager, "_save_workspace_history", return_value=None
        ):
            first = YR_MPEA_Widget()
            try:
                first._switch_workspace(str(self.root))
                self.assertEqual(
                    set(first.dataset_selection_manager.items), set(self.manager.items)
                )
            finally:
                first.close()
                APP.processEvents()

            second = YR_MPEA_Widget()
            try:
                second._switch_workspace(str(self.root))
                manager = second.dataset_selection_manager
                self.assertEqual(set(manager.items), set(self.manager.items))
                self.assertEqual(set(manager.datasets), {self.ds_ortho, self.ds_bayes})
                self.assertEqual(
                    [button.item.loci_name
                     for button in second.workspace.item_buttons["models"]],
                    ["GTR+G", "Partitioned model"],
                )
            finally:
                second.close()
                APP.processEvents()


if __name__ == "__main__":
    unittest.main(verbosity=2)
