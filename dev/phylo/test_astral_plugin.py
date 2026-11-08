"""Tests for the ASTRAL plugin migration."""

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

from YR_MPE.plugins.astral_plugin import (  # noqa: E402
    AstralPlugin, AstralThread, normalize_weighted_support_labels,
)

APP = QApplication.instance() or QApplication(sys.argv)
ASTRAL = REPO / "YR_MPE" / "bin" / "aster" / "astral.bat"

GENES = ""\
    "((A:1,B:1):1,(C:1,D:1):1);\n"\
    "((A:1,C:1):1,(B:1,D:1):1);\n"\
    "((A:1,D:1):1,(B:1,C:1):1);\n"


class AstralPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = AstralPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata_and_versions(self):
        self.assertEqual(self.plugin.tool_name, "astral")
        self.assertIn("zhang2025aster", self.plugin.citation_keys)
        self.assertNotIn("astral-pro3", self.plugin.ASTRAL_VERSIONS)
        texts = [self.plugin.version_combo.itemText(i)
                 for i in range(self.plugin.version_combo.count())]
        self.assertFalse(any("Pro3" in text for text in texts))

    def test_composite_support_labels_are_reduced_to_one_value(self):
        source = self.root / "composite.trees"
        target = self.root / "normalized.trees"
        source.write_text(
            "((A:1,B:1)/1/100:1,(C:1,D:1)/0.998/93:1);\n"
            "((A:1,C:1)/1/100:1,(B:1,D:1)/0.998/93:1);\n",
            encoding="utf-8",
        )
        count, changed = normalize_weighted_support_labels(
            str(source), str(target)
        )
        self.assertEqual(count, 2)
        self.assertGreater(changed, 0)
        text = target.read_text(encoding="utf-8")
        self.assertNotIn("/", text)
        self.assertIn(")1:1", text)
        self.assertIn(")0.998:1", text)

    def test_preimport_path_fills_the_input_box(self):
        path = self.root / "preimport.trees"
        path.write_text("((A:1,B:1):1,(C:1,D:1):1);\n", encoding="utf-8")
        plugin = AstralPlugin(import_from="YR_MPEA", import_data=str(path),
                              workdir=str(self.root))
        try:
            self.assertEqual(plugin.file_path_edit.text(), str(path))
            self.assertEqual(plugin.imported_files, [str(path)])
            self.assertEqual(plugin.import_file, str(path))
        finally:
            plugin.close()

    def test_selected_phylogenies_fill_the_input_box(self):
        from YR_MPE.platforms.methods.dataset_models import (
            DatasetItem, ITEM_TYPE_PHYLOGENY,
        )

        class Manager:
            def get_selected_items(self):
                items = []
                for name, newick in (
                    ("t1", "((A:1,B:1):1,(C:1,D:1):1);"),
                    ("t2", "((A:1,C:1):1,(B:1,D:1):1);"),
                ):
                    item = DatasetItem(item_type=ITEM_TYPE_PHYLOGENY)
                    item.loci_name = name
                    item.data = {"newick": newick}
                    items.append(item)
                return items

        self.plugin._dataset_selection_manager = Manager()
        self.plugin.import_from_selected_treeset()
        self.assertTrue(self.plugin.file_path_edit.text())
        self.assertEqual(len(self.plugin.imported_files), 1)
        lines = Path(self.plugin.imported_files[0]).read_text(
            encoding="utf-8"
        ).splitlines()
        self.assertEqual(len(lines), 2)

    def test_result_spec_tree(self):
        tree = self.root / "species.tre"
        tree.write_text("(A:0.1,B:0.1,(C:0.1,D:0.1)0.9:0.1);\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "astral")
        self.assertEqual(result.datasets[0].items[0].item_type, "phylogeny")


@unittest.skipUnless(ASTRAL.is_file(), "ASTRAL runtime not available")
class AstralThreadTests(unittest.TestCase):
    def test_species_tree_written(self):
        root = Path(tempfile.mkdtemp())
        genes = root / "genes.tre"
        genes.write_text(GENES, encoding="utf-8")
        payloads = []
        errors = []
        thread = AstralThread(
            str(ASTRAL), [str(genes)],
            ["astral4", "-t", "1", "-u", "1", str(genes)],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        tree = Path(payloads[0]["output_files"][0])
        self.assertTrue(tree.is_file())
        self.assertIn("(", tree.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
