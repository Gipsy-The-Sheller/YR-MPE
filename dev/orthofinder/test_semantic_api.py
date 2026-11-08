"""Tests for the OrthoFinder plugin's agent-facing semantic API."""

from __future__ import annotations

import json
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

from YR_MPE.plugins.orthofinder_plugin import OrthoFinderPlugin  # noqa: E402
from YR_MPE.platforms.methods.dataset_manager import sequence_file_stem  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class TestPlugin(OrthoFinderPlugin):
    def config(self):
        return True


class OrthoFinderApiTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = TestPlugin(import_from="YR_MPEA", workdir=str(self.root))
        self.plugin.result_ready.connect(lambda result: None)

    def tearDown(self):
        self.plugin.is_running = False
        self.plugin.close()

    def test_inspect_options(self):
        payload = self.plugin.inspect_options()
        self.assertEqual(payload["schema"], "yrmpe.orthofinder.options.v1")
        self.assertEqual(payload["options"]["source_type"], "gbff")
        self.assertIsNone(payload["options"]["genetic_code"])
        self.assertTrue(payload["options"]["skip_pseudogenes"])
        self.assertEqual(payload["options"]["threads"], 2)

        as_json = self.plugin.inspect_options(as_json=True)
        self.assertEqual(json.loads(as_json)["tool"], "OrthoFinder")

    def test_modify_options_atomic(self):
        result = self.plugin.modify_options({
            "source_type": "proteome",
            "input_path": str(self.root / "proteomes"),
            "output_dir": str(self.root / "out"),
            "genetic_code": 5,
            "skip_pseudogenes": False,
            "internal_stop": "warn",
            "keep_terminal_stop": True,
            "run_mode": "parse_only",
            "threads": 4,
        })
        self.assertTrue(result["ok"])
        self.assertIn("source_type", result["applied"])
        options = result["options"]
        self.assertEqual(options["source_type"], "proteome")
        self.assertEqual(options["genetic_code"], 5)
        self.assertFalse(options["skip_pseudogenes"])
        self.assertTrue(options["keep_terminal_stop"])
        self.assertEqual(options["internal_stop"], "warn")
        self.assertEqual(options["run_mode"], "parse_only")
        self.assertEqual(options["threads"], 4)

        before = self.plugin._collect_options()
        result = self.plugin.modify_options({"threads": 8, "unknown": 1})
        self.assertFalse(result["ok"])
        self.assertIn("unknown", result["rejected"])
        self.assertEqual(self.plugin._collect_options(), before)

        result = self.plugin.modify_options('{"threads": 3}')
        self.assertTrue(result["ok"])
        self.assertEqual(self.plugin._collect_options()["threads"], 3)

    def test_modify_rejects_invalid_values(self):
        result = self.plugin.modify_options({"threads": 99})
        self.assertFalse(result["ok"])
        self.assertIn("threads", result["rejected"])

        result = self.plugin.modify_options({"genetic_code": 3})
        self.assertFalse(result["ok"])
        self.assertIn("genetic_code", result["rejected"])

    def test_status_and_actions(self):
        status = self.plugin.inspect_status()
        self.assertEqual(status["state"], "idle")
        self.assertTrue(status["can_run"])
        self.assertTrue(status["can_close"])
        self.assertFalse(status["can_export"])

        actions = self.plugin.inspect_actions()
        ids = {action["id"] for action in actions["actions"]}
        self.assertIn("orthofinder.run", ids)
        self.assertIn("orthofinder.close", ids)
        self.assertIn("orthofinder.export_to_workspace", ids)
        self.assertEqual(actions["next_import"]["action"], "orthofinder.export_to_workspace")

        self.plugin._status_state = "completed"
        self.plugin._last_summary = {"output_dir": str(self.root), "single_copy_ogs": 2}
        status = self.plugin.inspect_status()
        self.assertTrue(status["can_export"])

    def test_console_inspection(self):
        self.plugin.console_output = []
        self.plugin.add_console_message("first", "info")
        self.plugin.add_console_message("second", "error")
        result = self.plugin.inspect_console(0, 10)
        self.assertEqual(result["total"], 2)
        self.assertEqual([line["message"] for line in result["lines"]], ["first", "second"])
        self.assertEqual(result["next"], 2)

        result = self.plugin.inspect_console(1, 10)
        self.assertEqual([line["message"] for line in result["lines"]], ["second"])

        result = self.plugin.inspect_console(0, 10, level="error")
        self.assertEqual([line["message"] for line in result["lines"]], ["second"])
        self.assertEqual(result["lines"][0]["index"], 1)
        self.assertEqual(result["next"], 2)

    def test_preview_results(self):
        out = self.root / "out"
        (out / "proteomes").mkdir(parents=True)
        (out / "cds").mkdir()
        (out / "proteomes" / "sp1.faa").write_text(">a\nMQ\n", encoding="utf-8")
        (out / "cds" / "sp1.fna").write_text(">a\nATG\n", encoding="utf-8")
        (out / "cds_map.tsv").write_text("gene\n", encoding="utf-8")
        single = out / "single_copy_ogs"
        single.mkdir()
        (single / "OG1.cds.fna").write_text(">a\nATG\n", encoding="utf-8")
        (single / "OG1.faa").write_text(">a\nM\n", encoding="utf-8")
        (single / "OG1.map.tsv").write_text("orthogroup\n", encoding="utf-8")
        (single / "single_copy_og_summary.tsv").write_text("orthogroup\n", encoding="utf-8")
        og_dir = out / "orthofinder" / "Results_1"
        og_dir.mkdir(parents=True)
        (og_dir / "Orthogroups.tsv").write_text("Orthogroup\tsp1\nOG1\ta\n", encoding="utf-8")

        self.plugin._last_summary = {"output_dir": str(out), "single_copy_ogs": 1}
        result = self.plugin.preview_results()
        self.assertTrue(result["ready"])
        self.assertEqual(result["single_copy_ogs"][0]["orthogroup"], "OG1")
        self.assertTrue(result["orthogroups_tsv"].endswith("Orthogroups.tsv"))
        roles = {item["role"] for item in result["files"]}
        self.assertIn("cds_map", roles)
        self.assertIn("proteome", roles)
        self.assertIn("cds", roles)

    def test_export_to_workspace_api(self):
        out = self.root / "out"
        single = out / "single_copy_ogs"
        single.mkdir(parents=True)
        (single / "OG1.cds.fna").write_text(">a\nATG\n>b\nATG\n", encoding="utf-8")
        (single / "single_copy_og_summary.tsv").write_text("orthogroup\n", encoding="utf-8")
        self.plugin._last_summary = {
            "output_dir": str(out), "species": 2, "cds": 2,
            "orthogroups": 1, "single_copy_ogs": 1,
        }
        result = self.plugin.agent_invoke("orthofinder.export_to_workspace")
        self.assertTrue(result["ok"])
        self.assertEqual(result["loci"], 1)

    def test_sequence_file_stem_removes_compound_suffixes(self):
        self.assertEqual(sequence_file_stem("OG0000008.cds.fna"), "OG0000008")
        self.assertEqual(sequence_file_stem("OG1.cds.fasta"), "OG1")
        self.assertEqual(sequence_file_stem("gene.fas"), "gene")

    def test_export_normalizes_og_names_and_species_ids(self):
        out = self.root / "out"
        single = out / "single_copy_ogs"
        single.mkdir(parents=True)
        (single / "OG1.cds.fna").write_text(
            ">sp1_protein1\nATG\n>sp2_protein2\nATG\n", encoding="utf-8"
        )
        (single / "OG1.map.tsv").write_text(
            "orthogroup\tspecies\tprotein_id\n"
            "OG1\tsp1\tsp1_protein1\n"
            "OG1\tsp2\tsp2_protein2\n",
            encoding="utf-8",
        )
        self.plugin._last_summary = {
            "output_dir": str(out), "species": 2, "cds": 2,
            "orthogroups": 1, "single_copy_ogs": 1,
        }
        result = self.plugin._build_result_spec()
        item = result.datasets[0].items[0]
        self.assertEqual(item.loci_name, "OG1")
        self.assertEqual([record.id for record in item.sequences], ["sp1", "sp2"])
        self.assertEqual(
            [record.description for record in item.sequences], ["sp1", "sp2"]
        )

    def test_agent_contract_and_invoke(self):
        contract = self.plugin.agent_contract()
        self.assertEqual(contract["schema"], "yrmpe.ui.orthofinder.v1")
        self.assertIn("options", contract)
        self.assertIn("actions", contract)

        result = self.plugin.agent_invoke("orthofinder.inspect_status")
        self.assertTrue(result["ok"])

        result = self.plugin.agent_invoke(
            "orthofinder.modify_options", {"options": {"threads": 5}}
        )
        self.assertTrue(result["ok"])
        self.assertEqual(self.plugin._collect_options()["threads"], 5)

        result = self.plugin.agent_invoke("orthofinder.preview_results")
        self.assertTrue(result["ok"])

    def test_run_parse_only_and_status(self):
        input_dir = REPO / "dev" / "gbff_to_orthofinder" / "testdata" / "gb"
        self.plugin.agent_invoke("orthofinder.modify_options", {"options": {
            "source_type": "gbff",
            "input_path": str(input_dir),
            "output_dir": str(self.root / "out"),
            "run_mode": "parse_only",
            "threads": 1,
        }})
        result = self.plugin.agent_invoke("orthofinder.run", {"show_dialogs": False})
        self.assertTrue(result["ok"])
        if self.plugin.analysis_thread:
            self.plugin.analysis_thread.wait(60000)
            for _ in range(20):
                APP.processEvents()

        status = self.plugin.inspect_status()
        self.assertEqual(status["state"], "completed")
        self.assertGreater(status["summary"].get("cds", 0), 0)
        results = self.plugin.preview_results()
        self.assertTrue(results["ready"])
        self.assertGreater(len(results["files"]), 0)

    def test_validate_parse_only(self):
        input_dir = self.root / "gb"
        input_dir.mkdir()
        (input_dir / "a.gb").write_text("LOCUS       test\n", encoding="utf-8")
        self.plugin.modify_options({
            "source_type": "gbff",
            "input_path": str(input_dir),
            "output_dir": str(self.root / "out"),
            "run_mode": "parse_only",
        })
        result = self.plugin.validate_configuration()
        self.assertTrue(result["input"])
        self.assertTrue(result["runtime"])
        self.assertTrue(result["ok"])

    def test_running_blocks_modification_and_close(self):
        self.plugin.is_running = True
        self.plugin._status_state = "running"
        status = self.plugin.inspect_status()
        self.assertFalse(status["can_run"])
        self.assertFalse(status["can_close"])
        self.assertFalse(self.plugin.can_close()["can_close"])
        result = self.plugin.modify_options({"threads": 3})
        self.assertFalse(result["ok"])
        self.assertIn("running", result["errors"][0])
        self.assertFalse(self.plugin.start_pipeline(show_dialogs=False))


if __name__ == "__main__":
    unittest.main(verbosity=2)
