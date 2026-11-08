# Substitution saturation (BaCoCa C-value)

## Purpose

Assess substitution saturation in multiple-sequence alignments with the
BaCoCa saturation index C. The plugin replaces the earlier simulated
saturation analysis with the BaCoCa C-value.

The C-value is used instead of Xia's Iss because Iss is not a clean saturation
decision rule (see Xu 2026 preprints). The C-value was introduced by Struck et
al. (2008) and implemented in BaCoCa; it has no intrinsic exclusion threshold.
The plugin reports it descriptively, with Ti/Tv versus p-distance scatter
plots, which is useful for direct inspection when only a few partitions are
available.

## Method

- The plugin computes the BaCoCa saturation index
  `C = std(Ti/Tv) / std(p)` from pairwise Ti/Tv counts and p-distances.
- Gap treatment can be pairwise deletion or complete deletion.
- Each alignment gets a saturation plot and summary statistics.
- A high C indicates that substitutions are approaching saturation; low C
  indicates that the alignment still carries phylogenetic signal.

This is an internal Python analysis; no external executable is invoked.

## GUI workflow

1. Import one or more alignments.
2. Choose gap treatment.
3. Run.
4. Select a file in the result tab to view its plot and statistics.
5. Export the current plot if needed.

## Output

- Saturation plot per alignment.
- Summary statistics including the C index.
- Per-file selector for multi-file runs.

## Pitfalls

- The C index depends on the alignment and gap treatment; compare only runs
  with the same settings.
- Saturation is not the same as model selection: a saturated alignment may
  still be usable after trimming or with a site-heterogeneous model.
- Pairwise and complete deletion can change the ranking of alignments.
- The plot is descriptive; it is not a formal hypothesis test.
