"""Tests for the BEST / MrBayes MSC configuration block."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.plugins.best_plugin import (  # noqa: E402
    BestConfig, build_mrbayes_block, parse_nexus_charsets, parse_nexus_taxa,
)

FINCH = Path("E:/Codes/YR-MPE-Windows-x86_64-softwares/"
             "MrBayes-MPI-BEAGLE3/examples/finch.nex")


class BestConfigTests(unittest.TestCase):
    def test_block_contains_msc_commands(self):
        config = BestConfig(
            partitions=[("gene1", "1-100", 1), ("gene2", "101-200", 1),
                        ("gene3", "201-300", 2)],
            species_partition="SpA: 1-2, SpB: 3-4",
            ploidy="haploid",
            species_prior="yule",
            ngen=5000,
        )
        block = build_mrbayes_block(config)
        self.assertIn("partition genes = 3:gene1,gene2,gene3;", block)
        self.assertIn("link revmat=(1,2)", block)
        self.assertNotIn("link revmat=(1,2,3)", block)
        self.assertIn("speciespartition species = SpA: 1-2, SpB: 3-4;", block)
        self.assertIn("lset ploidy=haploid;", block)
        self.assertIn("prset topologypr=speciestree;", block)
        self.assertIn("prset brlenspr=clock:speciestree;", block)
        self.assertIn("mcmc ngen=5000", block)

    def test_clock_and_gene_rate_layers_are_distinct(self):
        config = BestConfig(
            clock_model="wn",
            clock_rate_prior="lognormal(-0.7,0.3)",
            partition_rate_prior="variable",
            gene_rate_prior="variable",
        )
        block = build_mrbayes_block(config)
        self.assertIn("prset clockvarpr=wn;", block)
        self.assertIn("prset clockratepr=lognormal(-0.7,0.3);", block)
        self.assertIn("prset ratepr=variable;", block)
        self.assertIn("prset generatepr=variable;", block)

    def test_birth_death_uses_parameters(self):
        config = BestConfig(species_prior="birth-death",
                            extinction_prior="beta(1,2)", sample_prob="0.8")
        block = build_mrbayes_block(config)
        self.assertIn("prset extinctionpr=beta(1,2);", block)
        self.assertIn("prset sampleprob=0.8;", block)

    def test_invalid_yule_priors_are_normalized(self):
        config = BestConfig(
            partitions=[("p1", "1-10", 1)],
            species_partition="Sp1: 1, Sp2: 2",
            species_prior="yule",
            speciation_prior="gamma(1,10)",
            extinction_prior="fixed(0)",
        )
        block = build_mrbayes_block(config)
        self.assertIn("prset speciationpr=exponential(10.0);", block)
        self.assertNotIn("prset extinctionpr", block)

    @unittest.skipUnless(FINCH.exists(), "local MrBayes finch example not found")
    def test_parse_finch_charsets(self):
        partitions = parse_nexus_charsets(str(FINCH))
        self.assertEqual(len(partitions), 30)
        self.assertEqual(partitions[0][0].lower(), "locus097")
        self.assertEqual(len(parse_nexus_taxa(str(FINCH))), 4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
