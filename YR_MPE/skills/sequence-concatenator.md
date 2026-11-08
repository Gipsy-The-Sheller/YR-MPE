# Sequence Concatenator

## Purpose

Concatenate aligned FASTA partitions into one supermatrix and generate a NEXUS
charset file describing the partition coordinates.

## GUI workflow

1. **Add Files** selects the aligned partition FASTA files.
2. Each file must be internally rectangular: all sequences in that file have
   the same aligned length.
3. Set the output FASTA path.
4. Choose whether the NEXUS partition text is shown in a dialog or written to a
   file.
5. **Concatenate Sequences** writes the supermatrix and charset definitions.

There is no external command-line tool; concatenation is performed in Python.

## Matching and missing data

- Sequences are matched by their **full FASTA description**, not by a separate
  taxon ID field.
- Descriptions present in any partition are included in the supermatrix.
- A taxon absent from one partition is filled with `?` for that partition's
  full length.
- Partition lengths come from each input file; a file with ragged sequence
  lengths is rejected.

## Output contract

- Output FASTA contains one record per unique description, with concatenated
  sequence.
- The NEXUS text contains one `charset` per input file, with 1-based
  coordinates covering the concatenated length.
- File names are sanitised into NEXUS labels (punctuation replaced by `_`).

## Pitfalls

- Matching by full description means headers must be identical across
  partitions; `taxonA` and `taxonA description` are different taxa.
- Every input must be aligned; the concatenator does not align or trim.
- Missing data is `?`, not `-` or `N`.
- Partition order follows the file table, so reorder files before running if
  the charset order matters.
- The NEXUS output describes partitions only; it does not contain the sequence
  matrix.
