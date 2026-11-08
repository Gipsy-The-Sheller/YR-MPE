"""Tests for the NG86 / YN00 / Z-test / IQ-TREE omega utilities."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "YR_MPE" / "plugins" / "components" / "methods" / "selection_utils.py"


def _load():
    spec = importlib.util.spec_from_file_location("yrmpe_selection_utils", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


selection = _load()
TABLE = CodonTable.unambiguous_dna_by_id[1]
PAML_EXAMPLE = Path("E:/Codes/paml/examples/abglobin.nuc")


def read_paml_nuc(path: Path):
    """Read PAML's sequential format; '.' means 'same as first sequence'."""
    lines = path.read_text().splitlines()
    nseq, _ = map(int, lines[0].split()[:2])
    names = []
    index = 1
    while len(names) < nseq:
        if lines[index].strip():
            names.append(lines[index].strip())
        index += 1
    chunks = [[] for _ in range(nseq)]
    block = []
    for line in lines[index:]:
        text = line.strip()
        if not text:
            continue
        parts = text.split()
        if parts and all(part.isdigit() for part in parts):
            continue  # position marker
        if not any(ch in "ACGTN.-?U" for ch in text.upper()):
            continue
        block.append("".join(parts))
        if len(block) == nseq:
            for row, chunk in enumerate(block):
                chunks[row].append(chunk)
            block = []
    raw = ["".join(parts) for parts in chunks]
    reference = raw[0]
    expanded = [
        "".join(reference[i] if base == "." else base for i, base in enumerate(sequence))
        for sequence in raw
    ]
    return [SeqRecord(Seq(sequence), id=name)
            for name, sequence in zip(names, expanded)]


class CodonCountingTests(unittest.TestCase):
    def test_synonymous_sites(self):
        self.assertAlmostEqual(selection.synonymous_sites("ATG", TABLE), 0.0)
        self.assertAlmostEqual(selection.synonymous_sites("AAA", TABLE), 1.0 / 3.0)
        # CTG is fourfold synonymous at the third position.
        self.assertAlmostEqual(selection.synonymous_sites("CTG", TABLE), 4.0 / 3.0)
        # TAT -> TAC is synonymous; TAA and TAG are stops and are excluded.
        self.assertAlmostEqual(selection.synonymous_sites("TAT", TABLE), 1.0)

    def test_one_step_differences(self):
        self.assertEqual(selection.count_codon_differences("AAA", "AAG", TABLE), (1.0, 0.0))
        self.assertEqual(selection.count_codon_differences("AAA", "GAA", TABLE), (0.0, 1.0))
        self.assertEqual(selection.count_codon_differences("AAA", "AAA", TABLE), (0.0, 0.0))

    def test_multi_step_differences(self):
        syn, non = selection.count_codon_differences("TTT", "TTC", TABLE)
        self.assertEqual((syn, non), (1.0, 0.0))
        syn, non = selection.count_codon_differences("TTT", "GCC", TABLE)
        self.assertGreaterEqual(syn + non, 1.0)

    def test_pairwise_and_z(self):
        seq1 = "ATG" + "AAA" * 10 + "CTG" * 10
        seq2 = "ATG" + "AAG" * 10 + "CTG" * 10
        result = selection.pairwise_ng86(seq1, seq2, TABLE, name1="a", name2="b")
        self.assertIsNotNone(result)
        self.assertGreater(result.synonymous_differences, 0)
        self.assertEqual(result.nonsynonymous_differences, 0)
        self.assertIsNotNone(result.dS)
        self.assertIsNotNone(result.dN)
        self.assertIsNotNone(result.z)
        self.assertLess(result.z, 0)  # dN < dS

    def test_z_variance_uses_site_counts(self):
        seq1 = "TAT" * 50
        seq2 = "TAC" * 5 + "TAT" * 45
        result = selection.pairwise_ng86(seq1, seq2, TABLE)
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result.synonymous_sites, 50.0)
        self.assertAlmostEqual(result.nonsynonymous_sites, 100.0)
        # Var = Var(dS) + Var(dN), with pS(1-pS)/S as the binomial variance.
        pS = 5.0 / 50.0
        var_dS = (pS * (1.0 - pS) / 50.0) / (1.0 - 4.0 * pS / 3.0) ** 2
        expected_z = -selection._jc_distance(pS) / var_dS ** 0.5
        self.assertAlmostEqual(result.z, expected_z)

    def test_alignment_summary(self):
        sequences = [
            SeqRecord(Seq("ATG" + "AAA" * 10 + "CTG" * 10), id="a"),
            SeqRecord(Seq("ATG" + "AAG" * 10 + "CTG" * 10), id="b"),
            SeqRecord(Seq("ATG" + "AAA" * 10 + "CTT" * 10), id="c"),
            SeqRecord(Seq("ATG" + "AAG" * 10 + "CTC" * 10), id="d"),
        ]
        summary = selection.analyze_alignment(sequences, code_id=1, alignment="test")
        self.assertEqual(summary.sequences, 4)
        self.assertEqual(len(summary.results), 6)
        self.assertGreater(summary.usable_codons, 0)
        for result in summary.results:
            self.assertIn(result.conclusion, {
                "positive selection", "purifying selection",
                "not significant", "saturated/not estimable", "not estimable",
            })

    def test_stop_and_ambiguous_codons_skipped(self):
        seq1 = "ATG" + "TAA" + "AAA" + "NNN"
        seq2 = "ATG" + "TAA" + "AAG" + "AAA"
        result = selection.pairwise_ng86(seq1, seq2, TABLE)
        self.assertIsNotNone(result)
        self.assertEqual(result.sites, 2)  # ATG/ATG and AAA/AAG remain

    def test_csv_and_iqtree_parser(self):
        seq1 = "CTG" * 50
        seq2 = "CTA" * 5 + "CTG" * 45
        result = selection.pairwise_ng86(seq1, seq2, TABLE)
        self.assertIsNotNone(result)
        csv = selection.format_results_csv([result])
        self.assertIn("alignment\tseq1\tseq2", csv)
        header = csv.splitlines()[0].split("\t")
        self.assertIn("method", header)
        self.assertIn("SE_dS", header)
        self.assertIn("SE_dN", header)
        self.assertLess(result.dN, result.dS)
        self.assertIn(result.conclusion, {"purifying selection", "not significant"})

        text = """
        Optimal log-likelihood: -15.988
        Nonsynonymous/synonymous ratio (omega): 0.381
        """
        omega, logl = selection.parse_iqtree_omega(text)
        self.assertAlmostEqual(omega, 0.381)
        self.assertAlmostEqual(logl, -15.988)


class YN00Tests(unittest.TestCase):
    def test_yn00_synthetic(self):
        records = [
            SeqRecord(Seq("ATG" + "AAA" * 10 + "CTG" * 10), id="a"),
            SeqRecord(Seq("ATG" + "AAG" * 10 + "CTG" * 10), id="b"),
            SeqRecord(Seq("ATG" + "AAA" * 10 + "CTT" * 10), id="c"),
        ]
        summary = selection.analyze_alignment_yn00(records, code_id=1, alignment="test")
        self.assertEqual(summary.method, "YN00")
        self.assertEqual(len(summary.results), 3)
        self.assertTrue(any(result.se_dS is not None for result in summary.results))

    @unittest.skipUnless(PAML_EXAMPLE.exists(), "local PAML example not found")
    def test_yn00_matches_paml_abglobin(self):
        records = read_paml_nuc(PAML_EXAMPLE)
        summary = selection.analyze_alignment_yn00(records, code_id=1,
                                                   alignment="abglobin")
        result = summary.results[0]
        self.assertEqual((result.seq1, result.seq2), ("human", "goat-cow"))
        self.assertAlmostEqual(result.synonymous_sites, 183.7, places=1)
        self.assertAlmostEqual(result.nonsynonymous_sites, 671.3, places=1)
        self.assertAlmostEqual(result.dN, 0.0818, places=4)
        self.assertAlmostEqual(result.dS, 0.5031, places=4)
        self.assertAlmostEqual(result.omega, 0.1625, places=4)
        self.assertAlmostEqual(result.se_dN, 0.0115, places=4)
        self.assertAlmostEqual(result.se_dS, 0.0930, places=4)

        result = next(r for r in summary.results
                      if {r.seq1, r.seq2} == {"rat", "marsupial"})
        self.assertAlmostEqual(result.synonymous_sites, 194.5, places=1)
        self.assertAlmostEqual(result.nonsynonymous_sites, 660.5, places=1)
        self.assertAlmostEqual(result.dN, 0.2354, places=4)
        self.assertAlmostEqual(result.dS, 1.4797, places=4)
        self.assertAlmostEqual(result.omega, 0.1591, places=4)
        self.assertAlmostEqual(result.se_dN, 0.0215, places=4)
        self.assertAlmostEqual(result.se_dS, 0.5173, places=4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
