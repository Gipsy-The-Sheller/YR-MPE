# DecentTree

## Purpose

Build a distance-based tree from an ML distance matrix with DecentTree's
scalable neighbour-joining variants.

## GUI to CLI

| GUI control | CLI option | Notes |
|---|---|---|
| Input matrix | `-in PATH` | `.mldist` or `.dist` |
| Algorithm | `-t BIONJ-R\|BIONJ-V\|NJ-R-D\|ONJ-R\|ONJ-R-V\|AUCTION` | |
| Output tree | `-out PATH` | Plugin-managed: `<input>.nwk` |
| Banner | `-no-banner` | Always added |

The current Windows build avoids the classic BIONJ/NJ paths, which crash; the
exposed algorithms are the verified variants. `BIONJ-R` is the default.

## Flow

1. Compute or import an ML distance matrix (for example from the ML distance
   plugin).
2. Choose an algorithm.
3. Run; send the resulting Newick tree to the workspace.

## Pitfalls

- The input must be a distance matrix, not an alignment.
- DecentTree writes the tree next to the input as `<input>.nwk`; the GUI does
  not expose a separate output picker.
- Batch mode runs one command per matrix file.
- The algorithm names are case-sensitive CLI tokens; the GUI strips the
  parenthetical description before passing them.
