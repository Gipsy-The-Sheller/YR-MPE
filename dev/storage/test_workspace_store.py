"""Tests for the workspace store (stdlib only).

Run from the repository root:

    python dev/storage/test_workspace_store.py
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "YR_MPE" / "runtime" / "workspace_store.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("yrmpe_workspace_store", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


store_module = _load_module()


class WorkspaceStoreTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.source_dir = Path(tempfile.mkdtemp())
        self.fasta = self.source_dir / "OG1.cds.fna"
        self.fasta.write_text(">a\nATGC\n>b\nATGA\n", encoding="utf-8")

    def test_import_item_run_stream_roundtrip(self):
        store = store_module.WorkspaceStore(str(self.root), name="test")
        dataset = store.add_dataset("OrthoFinder", source="one_step_orthofinder")

        imported = store.new_import(
            [str(self.fasta)],
            importer="one_step_orthofinder",
            parameters={"genetic_code": 5},
            dataset_id=dataset.id,
        )
        item = store.add_item(
            dataset.id,
            "sequence",
            "OG1",
            source_path=str(self.fasta),
            copy=True,
            file_format="fasta",
            import_id=imported.id,
            length=4,
            sequence_count=2,
        )
        run = store.new_run(
            "mafft",
            "MAFFT",
            [item.id],
            name="MAFFT run 1",
            parameters={"auto": True},
        )
        aligned = store.add_item(
            dataset.id,
            "alignment",
            "OG1 alignment",
            source_path=str(self.fasta),
            copy=True,
            file_format="fasta",
            run_id=run.id,
            is_aligned=True,
            length=4,
            sequence_count=2,
        )
        store.finish_run(run.id, output_ids=[aligned.id], artifacts=[str(self.fasta)])
        stream = store.add_stream("OG1 pipeline", dataset.id, [run.id])
        run.stream_id = stream.id
        store.rename_run(run.id, "OG1 MAFFT")
        store.save()

        self.assertTrue((self.root / "workspace.json").exists())
        self.assertTrue((self.root / item.file).exists())
        self.assertTrue((self.root / imported.manifest).exists())
        self.assertTrue((self.root / run.manifest).exists())
        self.assertTrue((self.root / run.log).exists())

        reloaded = store_module.WorkspaceStore(str(self.root))
        self.assertEqual(reloaded.index.name, "test")
        self.assertEqual(reloaded.get_dataset(dataset.id).name, "OrthoFinder")
        loaded_item = reloaded.get_item(item.id)
        self.assertEqual(loaded_item.import_id, imported.id)
        self.assertEqual(loaded_item.type, "sequence")
        self.assertEqual(loaded_item.sequence_count, 2)
        loaded_run = reloaded.get_run(run.id)
        self.assertEqual(loaded_run.name, "OG1 MAFFT")
        self.assertEqual(loaded_run.output_ids, [aligned.id])
        self.assertEqual(loaded_run.status, "completed")
        self.assertEqual(reloaded.get_stream(stream.id).run_ids, [run.id])
        self.assertEqual(reloaded.get_import(imported.id).parameters["genetic_code"], 5)

        with (self.root / imported.manifest).open(encoding="utf-8") as handle:
            import_manifest = json.load(handle)
        self.assertEqual(import_manifest["items"][0]["id"], item.id)

        with (self.root / run.manifest).open(encoding="utf-8") as handle:
            run_manifest = json.load(handle)
        self.assertEqual(run_manifest["inputs"], [item.id])
        self.assertEqual(run_manifest["outputs"], [aligned.id])

    def test_import_hashes_and_missing_source(self):
        store = store_module.WorkspaceStore(str(self.root))
        dataset = store.add_dataset("D")
        imported = store.new_import([str(self.fasta), str(self.source_dir / "nope")], "test")
        self.assertIn(str(self.fasta), imported.source_hashes)
        self.assertTrue(imported.source_hashes[str(self.fasta)].startswith("sha256:"))
        self.assertEqual(len(imported.source_paths), 2)
        self.assertTrue((self.root / "imports" / imported.id / "files" / "OG1.cds.fna").exists())

        item = store.add_item(
            dataset.id, "sequence", "missing", source_path=str(self.source_dir / "nope")
        )
        self.assertEqual(item.file, "")
        self.assertIn("missing_source", item.metadata)

    def test_migrate_legacy_selection_state(self):
        legacy = self.root / "dataset_selection_state.json"
        legacy.write_text(
            json.dumps(
                {
                    "version": "2.0",
                    "datasets": {
                        "ds1": {
                            "id": "ds1",
                            "name": "Legacy",
                            "description": "",
                            "settings": {},
                        }
                    },
                    "items": {
                        "item1": {
                            "id": "item1",
                            "dataset_id": "ds1",
                            "item_type": "sequence",
                            "loci_name": "OG1",
                            "file_path": str(self.fasta),
                            "length": 4,
                            "sequence_count": 2,
                            "selection_state": "green",
                            "selection_reason": "Selected for workflow",
                        }
                    },
                    "selected_items": ["item1"],
                    "basket_items": ["item1"],
                }
            ),
            encoding="utf-8",
        )
        store = store_module.WorkspaceStore(str(self.root))
        migrated = store.migrate_legacy_selection_state(str(legacy))
        self.assertEqual(migrated, 1)
        self.assertTrue((self.root / "dataset_selection_state.json.bak").exists())
        self.assertEqual(store.get_item("item1").dataset_id, "ds1")
        self.assertEqual(store.get_dataset("ds1").item_ids, ["item1"])
        self.assertEqual(store.get_item("item1").metadata["migrated_from"], legacy.name)
        selection = store.index.settings["selection"]
        self.assertEqual(selection["basket"], ["item1"])
        self.assertEqual(selection["selected"], ["item1"])
        self.assertEqual(selection["items"]["item1"]["state"], "green")

        reloaded = store_module.WorkspaceStore(str(self.root))
        self.assertEqual(reloaded.get_item("item1").name, "OG1")
        self.assertEqual(
            reloaded.index.settings["selection"]["items"]["item1"]["reason"],
            "Selected for workflow",
        )

    def test_rename_manifest_and_stream_links(self):
        store = store_module.WorkspaceStore(str(self.root))
        dataset = store.add_dataset("D")
        run = store.new_run("iqtree", "IQ-TREE", [])
        stream = store.add_stream("S", dataset.id)
        store.add_run_to_stream(stream.id, run.id)
        store.rename_run(run.id, "Tree run A")
        store.save()

        reloaded = store_module.WorkspaceStore(str(self.root))
        self.assertEqual(reloaded.get_run(run.id).stream_id, stream.id)
        self.assertEqual(reloaded.get_stream(stream.id).run_ids, [run.id])
        with (self.root / run.manifest).open(encoding="utf-8") as handle:
            self.assertEqual(json.load(handle)["name"], "Tree run A")


if __name__ == "__main__":
    unittest.main(verbosity=2)
