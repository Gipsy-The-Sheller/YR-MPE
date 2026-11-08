# SeqMatrix

## Purpose

Build and edit multi-locus datasets in a spreadsheet-style matrix, including
sequence import, NCBI downloads, partition definitions and export.

## GUI features

- Import from sequences, NEXUS, partitioned NEXUS, or Excel.
- Download sequences from NCBI by accession.
- Append/delete rows and columns.
- Edit sequence names and partition metadata.
- Export partitions as FASTA files or as partitioned NEXUS.
- **Import as a dataset** sends the matrix to the workspace as typed items.

There is no command-line interface; SeqMatrix is a data-management GUI.

## Workflow

1. Create or import a matrix.
2. Add taxa (rows) and loci/partitions (columns), or download sequences.
3. Edit names/partitions so rows and columns are consistent.
4. Export FASTA/NEXUS, or import as a dataset for alignment/tree plugins.
5. Continue with MAFFT/TrimAl/IQ-TREE from the dataset.

## Pitfalls

- Rows are taxa and columns are partitions; mixing the two produces a
  transposed dataset.
- NCBI downloads require network access and valid accessions; downloaded
  sequences are grouped separately until placed.
- Partition definitions must match the exported FASTA order.
- Importing as a dataset creates workspace items; exporting to files does not.
- Sequence names must be unique within a partition for downstream matching.
