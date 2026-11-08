# Site Concordance Diagnostics (sCF + CASTER q1-q3)

## Purpose

Compare two site-level views of the same fixed species tree and alignment:

- IQ-TREE 3 likelihood site concordance factors (`--scfl`) with
  `sCF/sDF1/sDF2`;
- CASTER-site fixed-tree quartet scores (`-C -c ... -u 2`) with
  `q1/q2/q3`.

The result is one TSV table with one row per internal bipartition, named by
the **smaller taxon side**. The plugin does not infer a new species tree and
does not estimate MSC parameters.

## Flow

1. Open `SPECIES TREE -> Site Diagnostics (sCF + CASTER)`.
2. Select the fixed species tree and alignment.
3. Set the number of sCF quartets and threads.
4. Run. IQ-TREE computes sCF; CASTER-site scores the same fixed tree.
5. Inspect the Table and Overview tabs; send the table to the workspace.

## Columns

`split`, `split_size`, `branch_id`, `sCF`, `sDF1`, `sDF2`, `sN`,
`q1_caster`, `q2_caster`, `q3_caster`, `q_negative`.

## Interpretation

- `sCF/sDF1/sDF2` are percentages from IQ-TREE; `q1/q2/q3` are normalized
  CASTER-site quartet scores (not percentages).
- CASTER-site scores can be negative because of negative quartet weights;
  `q_negative` marks rows containing a negative q value.
- sCF is sensitive to multiple hits. Neither sCF nor CASTER q1-q3 should be
  inverted to coalescent units, theta, or any other MSC parameter. They are
  descriptive concordance diagnostics only.
- q1-q3 and sCF are not independent truth sources; their disagreement is
  evidence about model fit or data quality, not a validation target.

## Pitfalls

- Use the same alignment, taxon set and **binary** species tree for both
  programs; `-C -c` fixes the CASTER topology.
- CASTER `--chunk` can affect its quartet scores; keep it fixed across
  comparisons.
- Blank sCF or q cells mean the split was absent from that program.
- Negative q values are retained rather than clipped; clipping would hide
  the negative-weight behaviour.
