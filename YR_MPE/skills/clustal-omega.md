# Clustal Omega

## Purpose

Multiple-sequence alignment with Clustal Omega 1.2.2.

## GUI to CLI

The runner calls:

```text
clustalo -i <input> -o <output> --outfmt=fa [--auto] [--seqtype T] [--threads N] [--outorder=tree-order]
```

| GUI control | CLI option | Notes |
|---|---|---|
| Input FASTA | `-i PATH` | |
| Auto options | `--auto` | Default on |
| Sequence type | `--seqtype Protein\|RNA\|DNA` | Auto is omitted |
| Threads | `--threads N` | Omitted when 1 |
| Output order = tree-order | `--outorder=tree-order` | Input order adds nothing |
| Output | `-o PATH` | Plugin-managed |
| Output format | `--outfmt=fa` | Always FASTA |

## Flow

1. Import one or more FASTA files.
2. Leave Auto on unless you need explicit guide-tree/HMM control.
3. Set sequence type and threads.
4. Run; send the aligned FASTA to the workspace.

## Pitfalls

- `--auto` overrides manual guide-tree/HMM settings; turn it off before
  expecting those options to matter.
- `--seqtype` is case-sensitive in the GUI mapping; Auto means Clustal detects
  the type itself.
- The output is always FASTA (`--outfmt=fa`); the GUI suffix only changes the
  filename.
- Batch mode runs one process per input file.
