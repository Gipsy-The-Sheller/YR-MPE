# Case-study example data

These files are the data used by `docs/cases/cases.tex`. They are small
on purpose: each case should run on a laptop in a few minutes.

| Directory | Contents | Source |
|---|---|---|
| `single_gene/` | ITS2 raw/aligned/trimmed FASTA, IQ-TREE tree and ML distance | Wu & Liu (2019), ZooKeys 870:33-50; GenBank accessions listed in the chapter |
| `sanger/` | 18S, ITS and 28S amplicon FASTA and a reference tree | `benchmark/open_tasks/task02` |
| `organelle/` | Four annotated chloroplast GenBank files | `dev/gbff_to_orthofinder/testdata/gb` (Eucalyptus, Sesamum, Andrographis, Camellia) |
| `mrbayes/` | `locus001.nex` and two MrBayes `.p` chains | `dev/fig3/mrbayes` |
| `infcat/` | `prmA.phy` (55 taxa, 411 aa) | BMGE example alignment |
| `distance/` | ITS2 trimmed alignment, ML distance, square subset | derived from `single_gene/` |
| `species_tree/` | 10-taxon species tree, 150 ILS gene trees, gCF/sCF/MSC tables, duplicated gene tree | `dev/fig3/ils` and `dev/fig3/treerecs` |
| `dating/` | Simulated 6-taxon tree, date file, LSD2 report and dated NEXUS | generated for this manual |
| `selection/` | Four-sequence, 26-codon alignment | `dev/selection/testdata` |

## Notes

- The ITS2 accessions were retrieved from NCBI with `efetch`; the raw file is
  retained so the alignment can be reproduced without network access.
- `distance/its2_subset_square.dist` is an eight-taxon square matrix for the
  display figure; `its2.mldist` is the full IQ-TREE lower-triangular matrix.
- `dating/simulated6.dates` sets all tips to 0 and the root to -60. The
  negative sign is intentional: these are years before present. LSD2 writes
  the dated tree to `simulated6_dated.nexus` and the report to
  `simulated6_lsd2.log`.
- The selection alignment is a synthetic teaching example; it is not a
  biological data set.
- The dating example was produced with:

```text
lsd2 -i simulated6.nwk -s 1000 -d simulated6.dates -r a \
     -o simulated6_dated.nwk
```

  `simulated6_dated.nexus` is the NEXUS tree with `[&date="..."]`
  annotations; `simulated6_dated.nwk` is the LSD2 report.

