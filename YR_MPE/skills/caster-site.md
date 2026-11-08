# CASTER-site

## Purpose

Infer a species tree directly from a multiple-sequence alignment by sampling
quartets from alignment sites. CASTER-site does not use gene trees, so it is
an alternative to ASTRAL when site-level information should drive the
estimate.

## GUI to CLI

`control_bindings()` is the source of truth.

| GUI control | CLI option | Notes |
|---|---|---|
| Input alignment | positional | FASTA, PHYLIP, or a file list |
| Format | `--format fasta\|list\|phylip` | `auto` is omitted |
| Threads | `--thread N` | |
| Output tree | `--output PATH` | Optional |
| Support option | `--support 0..2` | 0 none, 1 length+support, 2 detailed |
| Subsample rounds | `--subsample N` | Advanced; default 4 |
| Constraint tree | `--constraint PATH` | Advanced |
| Guide trees | `--guide PATH` | Advanced |
| Rounds | `--round N` | Advanced; default 4 |
| Gene-to-taxon map | `--mapping PATH` | Advanced |
| Root taxon | `--root NAME` | Advanced |
| Subsample proportion | `--proportion X` | Advanced; default 0.25 |
| Seed | `--seed N` | Advanced; default 233 |
| Chunk size | `--chunk N` | Advanced; default 10000 |
| Ambiguity handling | `--ambiguity 0\|1` | 0=N, 1=diploid unphased |
| Branch lengths | `--length SULength\|CULength` | Advanced; default SULength |
| Scoring only | `--scoring` | |
| More rounds | `--moreround` | |

## Flow

1. Import one alignment, or a dataset that concatenates to one supermatrix.
2. Set format (usually `auto` or `fasta`), threads and support option.
3. Run and send the tree to the workspace.

## Pitfalls

- The input is an alignment, not a gene-tree file.
- `--format list` expects a file containing paths to per-locus alignments,
  not a FASTA alignment.
- CASTER-site can be slow on large alignments; increase threads before
  increasing rounds.
- The output is a Newick tree with branch lengths; terminal lengths may be
  missing and are added by the plugin for viewers.
