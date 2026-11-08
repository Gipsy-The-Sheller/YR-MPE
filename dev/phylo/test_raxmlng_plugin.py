"""Tests for the RAxML-NG plugin."""

from __future__ import annotations

import os
import subprocess
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

from YR_MPE.factories import PluginFactory  # noqa: E402
from YR_MPE.plugins.raxmlng_plugin import (  # noqa: E402
    RAxMLNGPlugin, RAxMLNGPluginEntry, RAxMLNGThread,
    normalize_model_name,
)

APP = QApplication.instance() or QApplication(sys.argv)
FIXTURE = REPO / "dev" / "model" / "testdata" / "four_taxa.fas"
TOOL = REPO / "YR_MPE" / "bin" / "raxml-ng" / "raxml-ng.exe"


class RAxMLNGModelTests(unittest.TestCase):
    def test_normalizes_iqtree_equal_frequencies(self):
        self.assertEqual(normalize_model_name("GTR+FQ+G4"), "GTR+FE+G4")
        self.assertEqual(normalize_model_name("GTR+F+G4"), "GTR+F+G4")

    def test_help_doc_maps_every_cli_control_to_help(self):
        plugin = RAxMLNGPlugin()
        try:
            doc = plugin.help_doc()
            self.assertEqual(doc["id"], "raxml-ng")
            self.assertTrue(doc["markdown"])
            self.assertEqual(doc["warnings"], [])
            controls = {row["id"]: row for row in doc["controls"]}
            self.assertEqual(plugin.bootstrap_spinbox.singleStep(), 100)
            self.assertEqual(controls["bootstrap"]["step"], 100)
            self.assertEqual(controls["model"]["cli_args"], ["--model", "GTR+G4"])
            self.assertEqual(controls["start_trees"]["cli_args"], [])
            help_text = subprocess.run(
                [str(TOOL), "--help"], capture_output=True, text=True
            ).stdout
            for row in doc["controls"]:
                for flag in row.get("cli", []):
                    if str(flag).startswith("-"):
                        self.assertIn(flag, help_text)
        finally:
            plugin.close()

    def test_plugin_metadata_and_factory(self):
        plugin = RAxMLNGPlugin()
        try:
            self.assertEqual(plugin.tool_name, "RAxML-NG")
            self.assertTrue(plugin.tool_path)
            self.assertEqual(plugin.get_citation_keys(), ["kozlov2019raxmlng"])
            plugin.set_model("GTR+FQ+G4")
            self.assertEqual(plugin.model_combo.currentText(), "GTR+FE+G4")
            factory = PluginFactory()
            self.assertIsInstance(
                factory.get_raxmlng_plugin(), RAxMLNGPluginEntry
            )
        finally:
            plugin.close()


@unittest.skipUnless(TOOL.is_file(), "RAxML-NG runtime not available")
class RAxMLNGThreadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(tempfile.mkdtemp())
        cls.payloads = []
        cls.errors = []
        cls.thread = RAxMLNGThread(
            str(TOOL), [str(FIXTURE)], "GTR+G4", threads=2,
            bootstrap=10, auto_mre=False, start_trees="",
            workdir=str(cls.root),
        )
        cls.thread.result_ready.connect(cls.payloads.append)
        cls.thread.error.connect(cls.errors.append)
        cls.thread.run()

    def test_run_succeeds(self):
        self.assertEqual(self.errors, [])
        self.assertEqual(len(self.payloads), 1)
        output = self.payloads[0]["outputs"][0]
        self.assertTrue(Path(output["best_tree"]).is_file())
        self.assertTrue(Path(output["support_tree"]).is_file())
        self.assertTrue(Path(output["bootstraps"]).is_file())
        self.assertIsInstance(output["log_likelihood"], float)

    def test_result_spec_has_tree_and_treeset(self):
        plugin = RAxMLNGPlugin(workdir=str(self.root))
        try:
            plugin._on_result(self.payloads[0])
            result = plugin._build_result_spec()
            self.assertEqual(result.validate(), [])
            types = {item.item_type for item in result.datasets[0].items}
            self.assertEqual(types, {"phylogeny", "treeset"})
            tree = next(
                item for item in result.datasets[0].items
                if item.item_type == "phylogeny"
            )
            self.assertEqual(tree.format, "newick")
            self.assertIn("GTR+G4", tree.data["model"])
        finally:
            plugin.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
