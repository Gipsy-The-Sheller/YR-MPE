"""Tests for the TreeRecs plugin migration."""

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

from YR_MPE.plugins.treerecs_plugin import TreeRecsPlugin, TreeRecsThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
TREERECS = REPO / "YR_MPE" / "bin" / "treerecs" / "treerecs.exe"


class TreeRecsPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = TreeRecsPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "TreeRecs")
        self.assertIn("penel2020treerecs", self.plugin.citation_keys)

    def test_result_spec_reconciled_tree(self):
        tree = self.root / "gene.nwk_recs.nwk"
        tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree)]
        result = self.plugin._build_result_spec()

        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "treerecs")
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "phylogeny")
        self.assertEqual(item.format, "newick")
        self.assertIn("A", item.data["newick"])

    @unittest.skipUnless(TREERECS.is_file(), "TreeRecs runtime not available")
    def test_run_analysis_wires_result_spec(self):
        gene = self.root / "gene.nwk"
        species = self.root / "species.nwk"
        content = "((A:1,B:1):1,(C:1,D:1):1);\n"
        gene.write_text(content, encoding="utf-8")
        species.write_text(content, encoding="utf-8")
        self.plugin.gene_tree_edit.setText(str(gene))
        self.plugin.species_tree_edit.setText(str(species))

        from YR_MPE.plugins import treerecs_plugin as module
        with mock.patch.object(module.QMessageBox, "information", return_value=None), \
                mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.run_analysis()
            self.assertTrue(self.plugin.analysis_thread.wait(30000))
            APP.processEvents()

        self.assertTrue(self.plugin.current_output_files)
        names = [Path(path).name for path in self.plugin.current_output_files]
        self.assertTrue(any(name.endswith("_recs.nwk") for name in names))

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import treerecs_plugin as module
        gene = self.root / "gene.nwk"
        species = self.root / "species.nwk"
        gene.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
        species.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
        self.plugin.gene_tree_edit.setText(str(gene))
        self.plugin.species_tree_edit.setText(str(species))
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(TREERECS.is_file(), "TreeRecs runtime not available")
class TreeRecsThreadTests(unittest.TestCase):
    def test_reconciled_tree_written(self):
        root = Path(tempfile.mkdtemp())
        gene = root / "gene.nwk"
        species = root / "species.nwk"
        content = "((A:1,B:1):1,(C:1,D:1):1);\n"
        gene.write_text(content, encoding="utf-8")
        species.write_text(content, encoding="utf-8")
        outdir = root / "out"

        payloads = []
        errors = []
        thread = TreeRecsThread(
            str(TREERECS), str(gene), str(species), str(outdir),
            ["-d", "2", "-l", "1"],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()

        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        names = [Path(path).name for path in payloads[0]["output_files"]]
        self.assertIn("gene.nwk_recs.nwk", names)
        self.assertIn("gene.nwk_recs.svg", names)
        reconciled = outdir / "gene.nwk_recs.nwk"
        self.assertTrue(reconciled.is_file())
        self.assertIn("(", reconciled.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
