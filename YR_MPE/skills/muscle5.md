# MUSCLE5

## Purpose

Fast multiple-sequence alignment with MUSCLE5.

## GUI to CLI

The runner calls `muscle -align|-super5 <input> -output <output>`.

| GUI control | CLI option | Notes |
|---|---|---|
| Input FASTA | positional | |
| Mode = Standard alignment | `-align` | Default |
| Mode = Large dataset | `-super5` | For many sequences |
| Output | `-output PATH` | Plugin-managed |
| Output suffix | none | GUI-only; default `_m5.fas` beside the input |

## Flow

1. Import one or more FASTA files.
2. Choose standard or large-dataset mode.
3. Run; the alignment is written beside the input (or to a temporary file)
   and can be sent to the workspace.

## Pitfalls

- `-align` and `-super5` are mutually exclusive modes, not tuning flags.
- MUSCLE5 writes through `-output`; the GUI suffix is not a MUSCLE option.
- Batch mode runs one process per input file.
- `-super5` is intended for large datasets and can differ slightly from
  `-align`; do not mix modes within one comparison.
