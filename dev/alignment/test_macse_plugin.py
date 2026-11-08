"""Tests for the MACSE2 plugin migration."""

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

from YR_MPE.plugins.macse_plugin import (  # noqa: E402
    MACSEPlugin, MACSEThread,
)

APP = QApplication.instance() or QApplication(sys.argv)
MACSE = REPO / "YR_MPE" / "bin" / "macse" / "macse.bat"


def records(*ids):
    return [SeqRecord(Seq("ATGCAAATGCAA"), id=name) for name in ids]


class MACSEPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = MACSEPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "MACSE2")
        self.assertIn("ranwez2018macse", self.plugin.citation_keys)
        self.assertIn("-gc_def", self.plugin.get_parameters())
        self.assertEqual(self.plugin.run_button.text(), "Run MACSE2")

    def test_import_list_once(self):
        self.plugin.close()
        self.plugin = MACSEPlugin(
            import_from="YR_MPEA", import_data=records("a", "b"),
            workdir=str(self.root),
        )
        self.assertEqual(len(self.plugin.imported_files), 1)
        self.assertFalse(self.plugin.batch_mode)

    def test_frameshift_filter_is_codon_aware(self):
        nt = [
            SeqRecord(Seq("ATG---AAA"), id="a"),
            SeqRecord(Seq("ATG---AAA"), id="b"),
        ]
        aa = [
            SeqRecord(Seq("M-K"), id="a"),
            SeqRecord(Seq("M-K"), id="b"),
        ]
        nt_out, aa_out = MACSEThread.remove_frameshift_gaps(nt, aa)
        self.assertEqual(str(nt_out[0].seq), "ATGAAA")
        self.assertEqual(str(aa_out[0].seq), "MK")

    def test_result_spec_has_nt_and_aa(self):
        nt = self.root / "locus_NT.fasta"
        aa = self.root / "locus_AA.fasta"
        nt.write_text(">a\nATGCAA\n>b\nATGCAA\n", encoding="utf-8")
        aa.write_text(">a\nMK\n>b\nMK\n", encoding="utf-8")
        self.plugin.alignment_output_files = [str(nt), str(aa)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "macse2")
        self.assertEqual(len(result.datasets[0].items), 2)
        self.assertEqual(result.datasets[0].summary["alignments"], 1)
        self.assertTrue(result.datasets[0].items[0].selected)
        self.assertFalse(result.datasets[0].items[1].selected)

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import macse_plugin as module
        with mock.patch.object(module.QMessageBox, "critical", return_value=None):
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(MACSE.is_file(), "MACSE2 runtime not available")
class MACSEThreadTests(unittest.TestCase):
    def test_macse_writes_nt_and_aa(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">a\nATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n"
            ">b\nATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n"
            ">c\nATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = MACSEThread(
            str(MACSE), [str(source)],
            ["-gc_def", "1", "-alphabet_AA", "SE_B_8",
             "-score_matrix", "BLOSUM62", "-optim", "2"],
            "align", workdir=str(root), suffix="_macse",
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()
        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        self.assertEqual(len(payloads[0]["output_files"]), 2)
        for path in payloads[0]["output_files"]:
            self.assertTrue(Path(path).is_file())
            self.assertGreater(Path(path).stat().st_size, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
