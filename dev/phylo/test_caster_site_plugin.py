"""Tests for the CASTER-site plugin migration."""

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

from YR_MPE.plugins.caster_site_plugin import (  # noqa: E402
    CasterSitePlugin, CasterSiteThread,
)

APP = QApplication.instance() or QApplication(sys.argv)
CASTER = REPO / "YR_MPE" / "bin" / "aster" / "caster-site.exe"

ALIGNMENT = (
    ">A\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCAA\n"
    ">B\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGTAA\n"
    ">C\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCAT\n"
    ">D\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCTT\n"
)


class CasterSitePluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = CasterSitePlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "caster-site")
        self.assertIn("zhang2025caster", self.plugin.citation_keys)

    def test_result_spec_tree(self):
        tree = self.root / "species.tre"
        tree.write_text("(A:0.1,B:0.1,(C:0.1,D:0.1)0.9:0.1);\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree)]
        result = self.plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "caster_site")
        self.assertEqual(result.datasets[0].items[0].item_type, "phylogeny")


@unittest.skipUnless(CASTER.is_file(), "CASTER-site runtime not available")
class CasterSiteThreadTests(unittest.TestCase):
    def test_species_tree_written(self):
        root = Path(tempfile.mkdtemp())
        aln = root / "aln.fas"
        aln.write_text(ALIGNMENT, encoding="utf-8")
        payloads = []
        errors = []
        thread = CasterSiteThread(
            str(CASTER), [str(aln)],
            ["--format", "fasta", "--thread", "1", "--support", "1", str(aln)],
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
