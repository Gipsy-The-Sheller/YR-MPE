# Robinson-Foulds Comparison

## Purpose

Compute pairwise unrooted Robinson-Foulds (RF) distances among trees or tree
sets with IQ-TREE's `--tree-dist-all`. The plugin exposes the RF matrix,
per-pair table, normalized RF and an RF distribution plot.

## GUI to CLI

| GUI control | CLI option | Notes |
|---|---|---|
| Tree file / TreeSet | `-t PATH` | One Newick tree per line |
| Threads | `-T N` | Default: half the CPU count, capped at 4 |
| Run | `--tree-dist-all` | Always added by the runner |
| Output prefix | `--prefix DIR/rf` | Plugin-managed |
| Rerun | `--redo` | Always added |

Rooted input trees are converted to unrooted topologies in a temporary file
before the IQ-TREE call. This normalization is GUI-only and does not appear as
a CLI flag.

## Comparison groups

When the input comes from selected TreeSets, the plugin records which trees
belong to the same set. Pairwise rows are labelled `within` or `between` so the
box plot can separate within-set and between-set RF distributions. Plain
single-tree input produces no group labels.

## Output

- Lower-triangular RF matrix (TSV and figure).
- Pairwise table: tree 1, tree 2, comparison type, RF, normalized RF.
- Normalized RF uses `RF / (2(n-3))` for `n >= 4` taxa.
- RF distribution box plot.

## Pitfalls

- RF is unrooted; do not interpret it as a rooted-tree distance.
- The normalized denominator depends on the taxon count; compare normalized
  values only for the same taxon set.
- Two trees are the minimum input; duplicate splits are counted, not deduplicated.
- IQ-TREE must be configured before Run, even though the dialog opens without it.
