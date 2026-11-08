# MSC Branch Diagnostics (gCF + ASTER CU)

## Purpose

Compare three views of the same fixed species tree and gene-tree set:

- IQ-TREE 3 gCF/gDF categories;
- ASTRAL-IV quartet frequencies and CU branch lengths (fixed-tree scoring);
- wASTRAL weighted quartet frequencies and CU branch lengths.

The result is one TSV table with one row per internal bipartition, named by
the **smaller taxon side**. The plugin does not infer a new species tree.

## Flow

1. Open `SPECIES TREE -> MSC Diagnostics (gCF + ASTER CU)`.
2. Select the fixed species tree and gene trees.
3. Run. IQ-TREE computes gCF; ASTER scores the fixed tree with ASTRAL-IV and,
   if the gene-tree labels are valid, wASTRAL.
4. Inspect the Table and Overview tabs; send the table to the workspace.

## Columns

`split`, `split_size`, `branch_id`, `gCF`, `gDF1`, `gDF2`, `gDFP`, `gN`,
`qCF`, `qDF1`, `qDF2`, `CU_astral`, `CU_wastral`, `x_gcf`, `gDFP_flag`,
`wastral_status`.

## wASTRAL label check

wASTRAL needs a finite numeric branch label on every internal gene-tree
branch. Missing labels, composite labels such as `/1/100`, and non-numeric
labels cause the plugin to reject wASTRAL automatically; ASTRAL-IV and the
gCF table are still produced. The reason is written to `wastral_status` and
the console.

## Pitfalls

- `x_gcf` is blank when gCF <= 1/3 or gCF = 100%; it is not a CU estimate at
  those boundary values.
- `gDFP > 0` is flagged; the gCF inversion assumes the three-way MSC
  resolution and is invalid under polyphyly/gene-tree error.
- ASTRAL-IV, wASTRAL and gCF are not three independent truth sources; their
  disagreement is diagnostic evidence.
- `CU_astral` is ASTER's estimate under its model, not ground truth.
