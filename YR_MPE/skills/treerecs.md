# TreeRecs

## Purpose

Reconcile a gene tree with a species tree, report duplication/loss events and
optionally write orthology/paralogy events.

## GUI to CLI

| GUI control | CLI option | Notes |
|---|---|---|
| Gene tree | `-g PATH` | Newick |
| Species tree | `-s PATH` | Newick |
| Output directory | `-o DIR` | Plugin-managed |
| Fix root / quiet | `-f -q` | Always added |
| Duplication cost | `-d N` | Default 2 |
| Loss cost | `-l N` | Default 1 |
| Best root by reconciliation | `-r` | Boolean |
| Contraction threshold | `-t X` | Omitted for `none` |
| Species prefix mode | `-p Y` | Only `Y` is emitted |
| Species separator | `-c CHAR` | Omitted for `_` |
| Orthology/paralogy events | `--fevent` | |
| Gene-to-species mapping | `-S PATH` | |

## Flow

1. Import or select a gene tree and a species tree.
2. Set duplication/loss costs and optional mapping.
3. Run; inspect the reconciliation and import the event tables/tree.

## Pitfalls

- TreeRecs 1.2 aborts with `std::stod` when given `-t none`; the plugin leaves
  the flag out unless a numeric threshold is set.
- Passing `-p N` or `-c _` explicitly switches separator mapping on and breaks
  leaf names such as `A/B/C`; both flags are only emitted when changed from
  their defaults.
- `-p` uses a species-prefix convention; for a separate mapping file use
  `-S` instead.
- The output directory is always created by the plugin; `-f -q` makes the run
  non-interactive.
