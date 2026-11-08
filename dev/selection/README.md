# SELECTION module

Two plugins behind the `SELECTION` toolbar menu:

| Plugin | Input | Method |
| --- | --- | --- |
| Codon-based Z-test (NG86 / YN00) | one in-frame codon FASTA | pairwise NG86 or YN00 + analytic Z-test |
| Batch dN/dS (NG86 / YN00 / IQ-TREE) | folder / workspace alignments | NG86 or YN00 counting, or IQ-TREE 3 GY/MG MLE |

## NG86 details

- Synonymous sites per codon: at each position, synonymous alternatives are
  counted after excluding changes that create a stop codon. The valid
  alternatives form the denominator (MEGA convention).
- Multi-step codon differences: all shortest mutation paths are enumerated;
  paths through stop codons are discarded and the remaining paths averaged.
- `pS = Sd / S`, `pN = Nd / N`, then Jukes-Cantor correction.
- Analytic Z-test: `Z = (dN - dS) / sqrt(Var(dN) + Var(dS))`, with
  `Var(dS) = pS(1-pS) / (1 - 4pS/3)^2 / S` (likewise for `N`).
  This is MEGA's large-sample Z-test, not a codon-count binomial.
- Two-sided `p` is shown in the matrix upper triangle; omega (`dN/dS`) in the
  lower triangle. Significant `p < 0.05` cells are bold and tinted
  (red `dN > dS`, green `dN < dS`).

## YN00 details

- Codon frequencies: average 3x4 base composition of the two sequences;
  `pi(codon) = f1(c1) f2(c2) f3(c3)`, normalised over sense codons.
- Kappa: estimated once per alignment from non-degenerate and fourfold
  degenerate sites, using the same F84 routine as PAML `yn00`.
- Sites: `CountSites` weights each one-step mutation by the target codon
  frequency and by kappa for transitions; S+N is scaled to 3 x codons.
- Differences: equal weighting over shortest paths; stop paths are dropped.
  This is PAML `yn00` with `weighting=0` (the default).
- dS/dN: F84 correction with K80/JC fallbacks; SEs from the F84 variance.
- Z: `(dN - dS) / sqrt(SE(dN)^2 + SE(dS)^2)`; the two-sided p-value is shown
  in the matrix upper triangle.
- Cross-check: `test_yn00_matches_paml_abglobin` compares against
  `PAML/examples/abglobin.nuc` (S=183.7, N=671.3, dN=0.0818, dS=0.5031,
  omega=0.1625, SE(dN)=0.0115, SE(dS)=0.0930 for human vs goat-cow).

## Batch MLE

```
iqtree3 -s ALN --seqtype CODON -m GY|MG -T N --prefix OUT --redo
```

The plugin parses the global omega and log-likelihood from stdout/`.iqtree`.
Counting mode does not need IQ-TREE.

## Tests and screenshots

```
python -m unittest dev.selection.test_selection_utils
python -m unittest dev.selection.test_selection_plugins
python dev/selection/render_plugins.py
```

Screenshots: `dev/selection/screenshots/01–06*.png`.

Reference for the counting/Z-test conventions (GPL-3, not vendored):
<https://github.com/KumarMEGALab/MEGA-source-code>, especially
`MEGA12-source/Compute/msynnonsyndist.pas` and `mcodons.pas`.
