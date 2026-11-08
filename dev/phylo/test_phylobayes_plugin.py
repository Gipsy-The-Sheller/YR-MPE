"""Tests for the PhyloBayes plugin migration."""

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

from YR_MPE.plugins.phylobayes_plugin import (  # noqa: E402
    BpcompThread, PhyloBayesPlugin, PhyloBayesThread,
)

APP = QApplication.instance() or QApplication(sys.argv)
PB = REPO / "YR_MPE" / "bin" / "pbmpi-1.9" / "pb_mpi.exe"
BPCOMP = REPO / "YR_MPE" / "bin" / "pbmpi-1.9" / "bpcomp.exe"
MPIEXEC = REPO / "YR_MPE" / "bin" / "msmpi" / "mpiexec.exe"

ALIGNMENT = (
    ">A\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCAA\n"
    ">B\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGTAA\n"
    ">C\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCAT\n"
    ">D\nATGCAAATGCAAATGCAAATGCAAATGCAAATGCAAATGCTT\n"
)


def write_phylip(path: Path) -> None:
    records = []
    for line in ALIGNMENT.splitlines():
        if line.startswith(">"):
            records.append([line[1:], ""])
        else:
            records[-1][1] += line.strip()
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("%d %d\n" % (len(records), len(records[0][1])))
        for name, sequence in records:
            handle.write("%s\t%s\n" % (name, sequence))


class PhyloBayesPluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.plugin = PhyloBayesPlugin(workdir=str(self.root))

    def tearDown(self):
        self.plugin.close()

    def test_metadata(self):
        self.assertEqual(self.plugin.tool_name, "PhyloBayes-MPI")
        self.assertIn("lartillot2013phylobayes", self.plugin.citation_keys)

    def test_result_spec_tree_and_trace(self):
        tree = self.root / "output.con.tre"
        tree.write_text("(A:0.1,(B:0.1,C:0.1)1:0.1,D:0.1);\n", encoding="utf-8")
        trace = self.root / "A_chain.trace"
        trace.write_text("Gen\tlnL\n1\t-1.0\n", encoding="utf-8")
        self.plugin.current_output_files = [str(tree), str(trace)]
        result = self.plugin._build_result_spec()

        self.assertEqual(result.validate(), [])
        self.assertEqual(result.source, "phylobayes")
        items = result.datasets[0].items
        phylogeny = [item for item in items if item.item_type == "phylogeny"]
        chains = [item for item in items if item.item_type == "chain"]
        self.assertEqual(len(phylogeny), 1)
        self.assertEqual(len(chains), 1)
        for name in ("A", "B", "C", "D"):
            self.assertIn(name, phylogeny[0].data["newick"])
        self.assertEqual(chains[0].data["file_paths"], [str(trace)])

    def test_missing_tool_does_not_start(self):
        from YR_MPE.plugins import phylobayes_plugin as module
        source = self.root / "input.phy"
        write_phylip(source)
        with mock.patch.object(module.QMessageBox, "critical", return_value=None), \
                mock.patch.object(module.QMessageBox, "warning", return_value=None):
            self.plugin.imported_files = [str(source)]
            self.plugin.tool_path = ""
            self.plugin.run_analysis()
        self.assertFalse(self.plugin.is_running)


@unittest.skipUnless(
    PB.is_file() and BPCOMP.is_file() and MPIEXEC.is_file(),
    "PhyloBayes runtime not available",
)
class PhyloBayesThreadTests(unittest.TestCase):
    def test_chain_and_consensus_written(self):
        root = Path(tempfile.mkdtemp())
        source = root / "input.phy"
        write_phylip(source)

        payloads = []
        errors = []
        thread = PhyloBayesThread(
            str(PB), str(MPIEXEC), [str(source)],
            ["2", "-ncat", "1", "-poisson", "-x", "1", "3"],
            workdir=str(root),
        )
        thread.result_ready.connect(payloads.append)
        thread.error.connect(errors.append)
        thread.run()

        self.assertEqual(errors, [])
        self.assertEqual(len(payloads), 1)
        names = [Path(path).name for path in payloads[0]["output_files"]]
        self.assertTrue(any(name.endswith(".trace") for name in names))
        self.assertTrue(any(name.endswith(".treelist") for name in names))

        chain_prefix = str(root / "input_chain")
        bp_payloads = []
        bp_errors = []
        bpcomp = BpcompThread(
            str(PB), [chain_prefix],
            {"cutoff": 0.5, "burnin": 1, "every": 1, "until": 3},
            workdir=str(root),
        )
        bpcomp.result_ready.connect(bp_payloads.append)
        bpcomp.error.connect(bp_errors.append)
        bpcomp.run()

        self.assertEqual(bp_errors, [])
        self.assertEqual(len(bp_payloads), 1)
        tree = Path(bp_payloads[0]["output_files"][0])
        self.assertTrue(tree.is_file())
        self.assertIn("(", tree.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
