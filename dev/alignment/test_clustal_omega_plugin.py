"""Tests for the Clustal Omega plugin migration."""

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
from unittest import mock  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

from YR_MPE.plugins.clustal_omega_plugin import (  # noqa: E402
    ClustalOmegaPlugin, ClustalOmegaThread,
)

APP = QApplication.instance() or QApplication(sys.argv)
CLUSTAL = REPO / "YR_MPE" / "bin" / "clustalo" / "clustalo.exe"


def records(*ids):
    return [SeqRecord(Seq("ATGCAAATGCAA"), id=name) for name in ids]


class ClustalOmegaPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = ClustalOmegaPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "Clustal Omega")
        self.assertIn("sievers2011fast", self.plugin.citation_keys)
        self.assertTrue(self.plugin.tool_path)
        self.assertIn("--auto", self.plugin.get_parameters())
        self.assertEqual(self.plugin.run_button.text(), "Run Clustal Omega")

    def test_import_list_handled_once(self):
        self.plugin.close()
        self.plugin = ClustalOmegaPlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)
        self.assertFalse(self.plugin.batch_mode)
        self.assertIn(">a", Path(self.plugin.imported_files[0]).read_text(encoding="utf-8"))

    def test_import_batch(self):
        self.plugin.close()
        self.plugin = ClustalOmegaPlugin(
            import_from="YR_MPEA",
            import_data=[records("a", "b"), records("c", "d")],
            workdir=str(self.root),
        )
        self.assertTrue(self.plugin.batch_mode)
        self.assertEqual(len(self.plugin.imported_files), 2)

    def test_result_spec(self):
        output = self.root / "locus_clustal.fas"
        output.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        self.plugin.alignment_output_files = [str(output)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "clustal_omega")
        self.assertEqual(result.datasets[0].source, "clustal_omega")
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "alignment")
        self.assertTrue(item.is_aligned)
        self.assertEqual(item.sequence_count, 2)

    def test_run_without_input_warns(self):
        from YR_MPE.plugins import clustal_omega_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(CLUSTAL.is_file(), "Clustal Omega runtime not available")
class ClustalOmegaThreadTests(unittest.TestCase):
    def test_clustal_writes_aligned_fasta(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGCAAATGCAA\n>b\nATGCAAATGCAA\n>c\nATGCAAATGCAA\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = ClustalOmegaThread(
            str(CLUSTAL), [str(source)], ["--auto"],
            workdir=str(root), suffix="_clustal",
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        output = Path(payloads[0]["output_files"][0])
        self.assertTrue(output.is_file())
        self.assertIn(">a", output.read_text(encoding="utf-8").lower())


if __name__ == "__main__":
    unittest.main(verbosity=2)
