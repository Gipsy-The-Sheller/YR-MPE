# Uncorrected p-distance

## Purpose

Compute an uncorrected p-distance matrix from an alignment without a
substitution model. This is the observed proportion of differing sites, not an
evolutionary distance corrected for multiple hits.

## GUI controls

| Control | Effect |
|---|---|
| Input alignment | FASTA alignment |
| Sequence type | AUTO, DNA or Protein |
| Gap treatment | Pairwise deletion or complete deletion |
| DNA rule | All substitutions, transitions only, transversions only, or advanced rules |

The calculation is implemented in Python; there is no external command-line
tool.

## Flow

1. Import an alignment.
2. Choose sequence type and gap treatment.
3. For DNA, optionally restrict the counted substitutions.
4. Run; inspect the matrix and send it to the workspace.

## Output

- Square pairwise p-distance matrix.
- Values range from 0 to 1 (proportion of differing sites).
- The matrix can feed distance-based tree methods, but it is uncorrected.

## Pitfalls

- p-distance underestimates evolutionary distance when multiple substitutions
  occur at the same site.
- Pairwise deletion changes the compared length per pair; complete deletion
  uses one common set of sites.
- Advanced DNA rules are not the same as a transition/transversion ratio
  model.
- Protein and DNA rules are not interchangeable; set the sequence type first.
- Do not compare p-distances computed with different gap treatments.
