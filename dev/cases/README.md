# Case-study UI assets

Screenshots and composite figures for the three case studies. All UI images
were captured from the real PyQt widgets; no mock-ups.

| Figure | Panels | Data |
|---|---|---|
| `figures/case_sanger.png` | A SeqMatrix, B MAFFT advanced dialog, C concatenated ML tree | `benchmark/open_tasks/task02` (12 taxa; 18S/ITS/28S; one missing ITS amplicon) |
| `figures/case_organelle.png` | A SeqDBG gene statistics + Sankey, B organelle tree | Four annotated chloroplast genomes in `dev/gbff_to_orthofinder/testdata/gb` |
| `figures/case_bayes.png` | A PhyloBayes infinite CAT setup, B MiniTracer trace/density/ESS, C MiniTracer two-chain comparison with PSRF | `prmA.phy` (55-taxon protein alignment) for the setup; MrBayes `locus001` run1/run2 chains for diagnostics |

## Reproduce

```bash
cd E:/YR-MPE_tutorial/YR-MPE_dev
python dev/cases/render_cases.py
python dev/cases/build_case_figures.py
```

`render_cases.py` writes the raw panel assets to `dev/cases/assets/`;
`build_case_figures.py` assembles the three composite figures in
`dev/cases/figures/`.

## Captions

**Sanger case.** (A) SeqMatrix matrix for 12 taxa and three amplicons (18S,
ITS, 28S); the blank ITS cell for B3421 is retained rather than imputed.
(B) MAFFT advanced parameters with `--adjustdirection` selected, used to
reverse-complement reads whose sequencing orientation is unknown. (C) The
concatenated maximum-likelihood tree inferred with IQ-TREE 3 and displayed in
IcyTree.

**Organelle case.** (A) SeqDBG gene-statistics table and adjacency Sankey for
four annotated chloroplast genomes; the selected gene (`ycf15`) and its
left/right neighbours are highlighted. (B) The species tree recovered from
the linked organelle CDS workflow, displayed in IcyTree.

**Bayesian case.** (A) PhyloBayes-MPI configured with GTR, four discrete
gamma categories and the infinite CAT mixture; two independent chains are
requested. (B) MiniTracer single-chain view: post-burn-in trace, marginal
density, ESS and 95% HPD for `lnLike`. (C) MiniTracer two-chain comparison
with the Gelman--Rubin PSRF (R-hat); the current implementation reports
PSRF in both the table and the plot title.

Note: the MiniTracer panels use the two MrBayes chains from the same case
workspace because they provide a converged pair for the PSRF demonstration;
MiniTracer reads MrBayes `.p`, PhyloBayes `.trace` and BEAST `.log` files
through the same import path.
