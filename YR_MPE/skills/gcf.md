# Gene Concordance Factors (gCF)

## Purpose

Annotate a fixed species tree with IQ-TREE 3 gene concordance factors. Input
is one species tree plus one file of gene trees (one Newick tree per line).
The plugin writes a TSV table with one row per internal bipartition. Each row
names the split by its **smaller taxon side**, so the representation is stable
across programs.

## Flow

1. Open `SPECIES TREE -> Gene Concordance (gCF)`.
2. Select the species tree and the gene-tree file.
3. Choose an output directory and thread count.
4. Run. The output tab has a Table tab and an Overview tab.
5. Send the table to the workspace; it appears as a `table` item.

## Columns

`split`, `split_size`, `branch_id`, `gCF`, `gDF1`, `gDF2`, `gDFP`, `gN`,
`branch_length`.

`gCF + gDF1 + gDF2 + gDFP = 100%`; `gN` is the number of decisive gene trees.
`gDFP` is polyphyly/other, not a clean NNI alternative.

## Pitfalls

- gCF is gene-tree support, not a site-level statistic.
- A gCF at or below 1/3 cannot be inverted to a positive CU length.
- Composite or missing gene-tree support labels are fine for gCF, but not for
  the wASTRAL part of the MSC diagnostics plugin.
