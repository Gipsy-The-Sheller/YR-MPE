"""Tests for the IQ-TREE plugin migration."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.plugins.components.alignment_common import (  # noqa: E402
    collect_partition_sequences,
)
from YR_MPE.plugins.iqtree_plugin import IQTreePlugin, IQTreeThread  # noqa: E402
from YR_MPE.plugins.iqtree_partition_thread import (  # noqa: E402
    IQTreePartitionThread,
)
from YR_MPE.plugins.partition_mode import (  # noqa: E402
    PartitionDefinition, PartitionMode,
)

APP = QApplication.instance() or QApplication(sys.argv)
IQTREE = REPO / "YR_MPE" / "bin" / "iqtree3" / "iqtree3.exe"


def records(*ids):
    return [SeqRecord(Seq("ATGCAAATGCAA"), id=name) for name in ids]


class IQTreePluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = IQTreePlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "IQ-TREE 3")
        self.assertIn("wong2026iqtree", self.plugin.citation_keys)
        self.assertEqual(self.plugin.run_button.text(), "Run IQ-TREE 3")

    def test_import_once(self):
        self.plugin.close()
        self.plugin = IQTreePlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)

    def test_help_doc_lists_models_rh_and_partition_options(self):
        self.assertEqual(self.plugin.ufboot_spinbox.singleStep(), 100)
        doc = self.plugin.help_doc()
        controls = {row["id"]: row for row in doc["controls"]}
        model = controls["model"]
        self.assertGreaterEqual(len(model["choices"]), 60)
        self.assertEqual(model["choices"]["GTR"], "GTR")
        self.assertEqual(model["choices"]["LG"], "LG")
        self.assertEqual(model["choices"]["auto"], "MFP")
        self.assertEqual(controls["ufboot"]["step"], 100)
        self.assertEqual(controls["partition_bootstrap"]["step"], 100)
        rate = controls["rate_heterogeneity"]["choices"]
        self.assertEqual(rate["Equal"], [])
        self.assertEqual(rate["Gamma Distribution (+G)"], ["+G"])
        self.assertEqual(rate["FreeRate Model (+R)"], ["+R"])
        self.assertEqual(controls["invariant_sites"]["cli"], ["+I"])
        freq = controls["state_freq"]["choices"]
        self.assertEqual(freq["Empirical (+F)"], "+F")
        self.assertEqual(freq["Equal (+FQ)"], "+FQ")
        partition = controls["partition_mode"]["choices"]
        self.assertEqual(partition["Edge-linked Proportional (-p)"], "-p")
        self.assertEqual(partition["Separate Tree (-S)"], "-S")

    def test_free_rate_uses_selected_category_count(self):
        self.plugin.model_combo.setCurrentText("GTR")
        self.plugin.rate_combo.setCurrentText("FreeRate Model (+R)")
        self.assertTrue(self.plugin.rate_categories_spinbox.isEnabled())
        self.plugin.rate_categories_spinbox.setValue(8)
        self.plugin.invar_checkbox.setChecked(False)
        self.plugin.state_freq_combo.setCurrentText("Estimated")
        self.assertEqual(self.plugin.build_model_string(), "GTR+R8")

        self.plugin.rate_categories_spinbox.setValue(4)
        self.assertEqual(self.plugin.build_model_string(), "GTR+R4")

        self.plugin.rate_combo.setCurrentText("Gamma Distribution (+G)")
        self.plugin.rate_categories_spinbox.setValue(8)
        self.assertEqual(self.plugin.build_model_string(), "GTR+G8")

        self.plugin.rate_combo.setCurrentText("Equal")
        self.assertFalse(self.plugin.rate_categories_spinbox.isEnabled())
        self.assertEqual(self.plugin.build_model_string(), "GTR")

    def test_single_alignment_without_file_uses_sequences(self):
        from YR_MPE.platforms.methods.dataset_models import DatasetItem

        item = DatasetItem(item_type="alignment")
        item.loci_name = "OG1"
        item.is_aligned = True
        item.sequences = records("a", "b")
        item.sequence_count = 2
        payload = {"dataset_items": [item], "dataset_config": {}}
        self.plugin.handle_import_data(payload)
        self.assertTrue(self.plugin.import_file)
        self.assertTrue(Path(self.plugin.import_file).is_file())
        self.assertIn(">a", Path(self.plugin.import_file).read_text(encoding="utf-8"))

    def test_partition_taxa_are_the_union_not_the_first_locus(self):
        items = [
            SimpleNamespace(sequences=records("D1", "D2"), file_path="",
                            length=12),
            SimpleNamespace(sequences=records("D2", "D3"), file_path="",
                            length=12),
        ]
        taxa, by_item = collect_partition_sequences(items)
        self.assertEqual(taxa, ["D1", "D2", "D3"])
        self.assertEqual([len(records) for records in by_item], [2, 2])

    def test_multilocus_import_keeps_taxa_missing_from_first_locus(self):
        from YR_MPE.platforms.methods.dataset_models import DatasetItem

        items = []
        for names in (("D1", "D2"), ("D2", "D3")):
            item = DatasetItem(item_type="alignment")
            item.loci_name = "locus_%d" % (len(items) + 1)
            item.is_aligned = True
            item.sequences = records(*names)
            item.sequence_count = len(item.sequences)
            item.length = 12
            items.append(item)

        self.plugin.handle_import_data({
            "dataset_items": items,
            "dataset_config": {"topo_linked": False, "edge_linked": False},
        })
        self.assertTrue(self.plugin.import_file)
        imported = list(SeqIO.parse(self.plugin.import_file, "fasta"))
        self.assertEqual(sorted(record.id for record in imported),
                         ["D1", "D2", "D3"])
        for record in imported:
            self.assertEqual(len(str(record.seq)), 24)

    def test_partition_merge_can_use_rcluster(self):
        thread = IQTreePartitionThread(
            str(IQTREE),
            str(self.root / "super.fas"),
            [PartitionDefinition("p1", "", "DNA", "1-10")],
            PartitionMode.EL,
            {"model": "MFP+MERGE", "threads": 1},
            bootstrap_enabled=False,
        )
        cmd = thread.build_iqtree_command("parts.nex")
        self.assertEqual(cmd[cmd.index("-m") + 1], "MFP+MERGE")
        self.assertNotIn("--rcluster", cmd)
        retry = thread.build_iqtree_command("parts.nex", rcluster_percent=10)
        self.assertEqual(retry[retry.index("--rcluster") + 1], "10")

    def test_result_spec_tree_and_model(self):
        tree = self.root / "input.fas.treefile"
        report = self.root / "input.fas.iqtree"
        tree.write_text("(a:0.1,b:0.1,c:0.1);\n", encoding="utf-8")
        report.write_text("Model of substitution: JC\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree), str(report)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "iqtree")
        types = [item.item_type for item in result.datasets[0].items]
        self.assertIn("phylogeny", types)
        self.assertIn("model", types)

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import iqtree_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(IQTREE.is_file(), "IQ-TREE runtime not available")
class IQTreeThreadTests(unittest.TestCase):
    def test_tree_and_report_written(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGCAAATGCAA\n>b\nATGCAAATGCAA\n>c\nATGCAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = IQTreeThread(
            str(IQTREE), [str(source)],
            ["-m", "JC", "-nt", "1", "-redo"],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        files = [Path(path).name for path in payloads[0]["output_files"]]
        self.assertTrue(any(name.endswith(".treefile") for name in files))
        self.assertTrue(any(name.endswith(".iqtree") for name in files))


if __name__ == "__main__":
    unittest.main(verbosity=2)
