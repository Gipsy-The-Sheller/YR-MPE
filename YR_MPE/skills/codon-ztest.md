# Codon-based Z-test

## Purpose

Test whether synonymous and nonsynonymous substitution rates differ
pairwise among codon-aligned sequences. The plugin implements NG86
(Nei-Gojobori counting) and YN00 (Yang-Nielsen, F84 correction) with an
analytic delta-method Z-test.

## Input

- Codon-aligned nucleotide FASTA.
- All sequences must be in frame and have the same aligned length.
- Stop codons and frameshifts should be masked before running.

## GUI workflow

1. Select a codon-aligned FASTA.
2. Choose NG86 or YN00.
3. Choose the genetic code.
4. Run.
5. Inspect the pairwise table and the omega matrix; export TSV if needed.

## Statistics

- Z-test is two-sided by default; one-sided p-values are also reported.
- The upper triangle of the omega matrix reports pairwise omega values; the
  lower triangle reports the corresponding p-values.
- Significant pairs are highlighted.
- The default variance follows the MEGA-style NG86 implementation; YN00 uses
  the F84 correction.

## Output

- Pairwise table: sequence pair, synonymous/nonsynonymous rates, omega, Z,
  p-values.
- Lower/upper triangular omega/p matrices.
- TSV export.

## Pitfalls

- The test is pairwise; it does not fit a single tree-wide dN/dS model.
- NG86 can be biased when transition/transversion rates are unequal or
  divergence is high; YN00 is the F84-corrected alternative.
- Alignment gaps and ambiguous codons reduce the usable codon count.
- A non-significant Z-test does not prove neutrality; it only fails to reject
  equal rates.
- Genetic code must match the organism; the default is the standard code.
