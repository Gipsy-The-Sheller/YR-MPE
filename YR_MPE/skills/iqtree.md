# IQ-TREE 3 skill

## Purpose

Maximum-likelihood tree inference from an alignment, with optional partition
mode, model selection and branch support.

## GUI to CLI: basic controls

| GUI control | CLI option | Notes |
|---|---|---|
| File input | `-s` | One IQ-TREE process per alignment |
| Sequence type | `-st DNA` / `-st AA` | `AUTO` omits the flag |
| Threads | `-nt` | Omitted when 1/AUTO |
| UFBoot | `-bb N` | Only when the UFBoot checkbox is enabled |
| aBayes | `-abayes` | Boolean flag |
| SH-aLRT | `-alrt N` | Only when the SH-aLRT checkbox is enabled |
| (fixed) | `-redo` | Plugin always asks IQ-TREE to overwrite |

## Substitution model

The GUI model combo contains the following base models. The live
`workbench.help_doc` returns the same catalog under
`controls[id=model].choices`, mapping each GUI label to the exact `-m` token.

- **Automatic:** `auto` -> `-m MFP` (ModelFinder plus tree inference).
- **DNA:** JC69/JC, F81, K2P/K80, HKY85/HKY, TNe, TN93/TN, K3P/K81,
  K81u/K3Pu, TPM2, TPM2u, TPM3, TPM3u, TIM, TIMe, TIM2, TIM2e, TIM3,
  TIM3e, TVM, TVMe, SYM, GTR.
- **Protein:** Blosum62, cpREV, Dayhoff, DCMut, EAL, ELM, FLAVI, FLU,
  GTR20, HIVb, HIVw, JTT, JTTDCMut, LG, mtART, mtMAM, mtREV, mtZOA,
  mtMet, mtVer, mtInv, NQ.bird, NQ.insect, NQ.mammal, NQ.pfam,
  NQ.plant, NQ.yeast, Poisson, PMB, Q.bird, Q.insect, Q.mammal,
  Q.pfam, Q.plant, Q.yeast, rtREV, VT, WAG.

The base model is extended into one `-m` string:

| GUI control | Model suffix | Notes |
|---|---|---|
| Rate Heterogeneity = Equal | none | No rate-heterogeneity operation |
| Rate Heterogeneity = Gamma Distribution (+G) | `+G{n}` | `n` from the Rate categories spin box |
| Rate Heterogeneity = FreeRate Model (+R) | `+R{n}` | `n` from the Rate categories spin box (default 4) |
| Invariable Sites | `+I` | Boolean |
| Ascertainment Bias | `+ASC` | Boolean |
| State Frequencies = Estimated | none | ML estimates |
| State Frequencies = Empirical (+F) | `+F` | Empirical |
| State Frequencies = ML-optimized (+FO) | `+FO` | ML-optimized |
| State Frequencies = Equal (+FQ) | `+FQ` | Equal |

Example: `GTR` + Gamma 4 + Invariant + Empirical = `-m GTR+G4+I+F`.
Example: `GTR` + FreeRate 8 = `-m GTR+R8`.

## Partition mode

Partition mode is configured in **Configure Partitions** and uses a partition
file generated from the selected dataset:

| GUI partition option | CLI flag | Notes |
|---|---|---|
| Edge-linked Proportional | `-p partition.nex` | Proportional branch lengths |
| Edge-linked Equal | `-q partition.nex` | Equal branch lengths |
| Edge-unlinked | `-Q partition.nex` | Unlinked branch lengths |
| Separate Tree | `-S partition.nex` | Separate tree per partition |
| Enable Bootstrap + replicates | `-bb N` | Partition thread also adds `-abayes` |
| Sampling: SITE / GENE / GENESITE | (reserved) | Current command builder records the choice but does not yet pass a sampling flag |
| RCluster + percentage | `--rcluster N` | Accelerated partition merging; passed when the merge path is used |

The partition thread also passes `-m MFP+MERGE` when no per-partition model is
set; if a `charpartition` model is present, models are read from the partition
file instead.

Taxa for the supermatrix are the **union across all selected partitions**;
a taxon absent from a locus is padded with `?` for that locus length. Do not
take the first partition as the taxon list when the dataset has missing data.

### Partition-merger workaround

IQ-TREE 3.0.1 can abort partition merging with:

```text
Fixing branch lengths not supported under specified site rate model
```

This happens inside the partition-merger evaluation, not in the final tree
search. RCluster takes a different pair-search path and avoids it. The plugin
passes the RCluster percentage when configured; if the first attempt still
fails with this message, it retries once with `--rcluster 10`.


## Flow

1. Select an alignment; optionally select a model item first.
2. Open `PHYLOGENY -> Maximum Likelihood (ML) -> IQ-Tree 3`.
3. Choose model/rate/state-frequency/support options and run.
4. Inspect `cli.preview` in `workbench.help_doc`; it shows the exact command
   assembled from the current controls.
5. Import the resulting phylogeny or tree set.
