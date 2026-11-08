"""Tests for the MrBayes plugin migration."""

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

from YR_MPE.plugins.mrbayes_plugin import MrBayesPlugin, MrBayesThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)
MB = REPO / "YR_MPE" / "bin" / "mb" / "mb.exe"

CON_TRE = """#NEXUS
begin taxa;
\tdimensions ntax=4;
\ttaxlabels
\t\tA
\t\tB
\t\tC
\t\tD
\t\t;
end;
begin trees;
\ttranslate
\t\t1\tA,
\t\t2\tB,
\t\t3\tC,
\t\t4\tD
\t\t;
   tree con_50_majrule = [&U] (1[&prob=1.00]:0.03,2[&prob=1.00]:0.02,(3[&prob=1.00]:0.01,4[&prob=1.00]:0.02)[&prob=0.67]:0.003);
end;
"""


def short_parameters():
    return [{
        'datatype': 'DNA',
        'nst': 2,
        'statefreq': 'fixed(equal)',
        'rates': 'equal',
        'ngammacat': 4,
        'use_beagle': False,
        'beagle_device': 'cpu',
        'beagle_precision': 'double',
        'beagle_scaling': 'dynamic',
        'ngen': 1000,
        'samplefreq': 100,
        'nchains': 2,
        'nruns': 1,
        'burnin_as_fraction': True,
        'burnin_frac': 0.25,
        'contype': 'majorityrule',
        'conformat': 'figtree',
        'enable_dating': False,
        'clock_model': 'strict',
        'clock_params': {},
        'calibrations': [],
        'outgroup': '',
    }]


class MrBayesPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = MrBayesPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "MrBayes-MPI-BEAGLE")
        self.assertIn("ronquist2012mrbayes", self.plugin.citation_keys)

    def test_result_spec_tree_and_chains(self):
        tree = self.root / "input.nex.con.tre"
        tree.write_text(CON_TRE, encoding="utf-8")
        run1 = self.root / "input.nex.run1.p"
        run1.write_text("Gen\tlnLike\n", encoding="utf-8")
        run2 = self.root / "input.nex.run2.p"
        run2.write_text("Gen\tlnLike\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree), str(run1), str(run2)]
        result = self.plugin._build_result_spec()

        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "mrbayes")
        items = result.datasets[0].items
        phylogeny = [item for item in items if item.item_type == "phylogeny"]
        chains = [item for item in items if item.item_type == "chain"]
        self.assertEqual(len(phylogeny), 1)
        self.assertEqual(len(chains), 2)
        newick = phylogeny[0].data["newick"]
        for name in ("A", "B", "C", "D"):
            self.assertIn(name, newick)
        self.assertNotIn("prob=", newick)
        self.assertIn("0.67", newick)  # posterior support kept as a node label
        self.assertEqual(chains[0].data["file_paths"], [str(run1)])

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import mrbayes_plugin as module
        source = self.root / "input.fas"
        source.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.import_file = str(source)
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(MB.is_file(), "MrBayes runtime not available")
class MrBayesThreadTests(unittest.TestCase):
    def test_consensus_tree_and_chains_written(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.fas"
        source.write_text(
            ">A\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCAA\n"
            ">B\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGTAA\n"
            ">C\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCAT\n"
            ">D\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCTT\n",
            encoding="utf-8",
        )
        payloads = []
        errors = []
        thread = MrBayesThread(
            str(MB), None, [str(source)], short_parameters(),
            run_data_block_checked=True, workdir=str(root), use_mpi=False,
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()

        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        files = payloads[0]["output_files"]
        names = [Path(path).name for path in files]
        self.assertTrue(any(name.endswith(".con.tre") for name in names))
        self.assertTrue(any(name.endswith(".p") for name in names))
        tree = next(Path(path) for path in files if str(path).endswith(".con.tre"))
        self.assertIn("begin trees", tree.read_text(encoding="utf-8").lower())


if __name__ == "__main__":
    unittest.main(verbosity=2)
