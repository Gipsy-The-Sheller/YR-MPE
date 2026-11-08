# SeqDBG

## Purpose

Compare gene order and annotation across organelle genomes, reconcile gene
names, inspect adjacency relationships and export a normalized dataset.

## GUI features

- Load multiple annotated genomes.
- Build a gene-adjacency graph.
- Normalize gene names with synonym rules.
- View left/right adjacency counts per gene.
- Edit and export/import normalization rules.
- Export selected loci as FASTA or as a workspace dataset.

There is no external command-line tool; SeqDBG runs in Python.

## Workflow

1. Load annotated genomes.
2. Review the gene list and adjacency relationships.
3. Add normalization rules where annotations use different names for the same
   gene.
4. Refresh the graph/table so adjacency counts use the normalized names.
5. Select loci and export as a dataset for alignment/tree inference.

## Pitfalls

- Adjacency is computed within each genome; a gene missing from some genomes
  appears with lower coverage.
- Normalization rules change labels, not the underlying sequences.
- Gene-order comparisons are not orthology calls by themselves; verify
  one-to-one correspondence before concatenation.
- Exporting to the workspace creates typed items; exporting FASTA only writes
  files.
- Strand/orientation is not normalized by renaming; check reverse-complemented
  genes separately.
