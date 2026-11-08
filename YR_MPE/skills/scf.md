# Site Concordance Factors (sCF)

## Purpose

Annotate a fixed species tree with IQ-TREE 3 site concordance factors. Input
is one species tree plus one alignment; gene trees are not required. The
plugin uses likelihood-based sCF (`--scfl`) and writes a TSV table with one
row per internal bipartition, named by the **smaller taxon side**.

## Flow

1. Open `SPECIES TREE -> Site Concordance (sCF)`.
2. Select the species tree and alignment.
3. Set the number of quartets and threads.
4. Run, inspect Table/Overview, then send the table to the workspace.

## Columns

`split`, `split_size`, `branch_id`, `sCF`, `sDF1`, `sDF2`, `sN`,
`branch_length`.

## Pitfalls

- sCF is site-level support and should not be substituted for gCF in the
  gCF-derived coalescent-length formula.
- sCF depends on the alignment and substitution model; an uninformative or
  misaligned partition can return zero informative sites.
