# TrimAl

## Purpose

Trim gap-rich and ambiguous columns from a multiple-sequence alignment with
TrimAl 1.2rev57.

## GUI to CLI

The runner builds `trimal -in <input> -out <output> [mode] [format]`.

| GUI control | CLI option | Notes |
|---|---|---|
| Input alignment | `-in PATH` | One run per file in batch mode |
| Trimming mode | `-gappyout`, `-strict`, `-strictplus`, `-automated1` | Default: automated1 |
| Output format | `-clustal`, `-nexus`, `-phylip`, `-phylip3.2`, `-mega`, `-nbrf` | FASTA adds no flag |
| Output file | `-out PATH` | Plugin-managed |
| Output suffix | none | GUI-only; output is beside the input |

Advanced parameters (applied only after accepting the dialog):

| GUI control | CLI option |
|---|---|
| Gap threshold | `-gt X` |
| Similarity threshold | `-st X` |
| Consistency threshold | `-ct X` |
| Minimum positions | `-cons N` |
| Window size | `-w N` |
| Gap window | `-gw N` |
| Similarity window | `-sw N` |
| Consistency window | `-cw N` |
| Residue overlap | `-resoverlap X` |
| Sequence overlap | `-seqoverlap N` |
| Remove all-gap columns | `-nogaps` |
| Remove columns composed only of gaps | `-noallgaps` |
| Complementary alignment | `-complementary` |
| Column numbering | `-colnumbering` |

Thresholds are only emitted when greater than zero; the default mode supplies
its own internal thresholds.

## Flow

1. Import an alignment.
2. Choose a trimming mode and output format.
3. Optionally set advanced thresholds/windows.
4. Run; send the trimmed alignment to the workspace.

## Pitfalls

- `-strict` and `-automated1` are alternatives; only one mode flag is emitted.
- Advanced threshold values of zero are omitted, not passed as `0`.
- `-nogaps` removes any column containing a gap; `-noallgaps` only removes
  columns made entirely of gaps. They are not interchangeable.
- FASTA output has no format flag; the other formats add one.
