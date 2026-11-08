# PD-Guide

## Purpose

Set up molecular-dating calibrations and run a first dating analysis. PD-Guide
is a guidance layer over the LSD2 controls, with method slots for other dating
engines.

## Method choices

| Method | Status |
|---|---|
| LSD2 | Implemented |
| IQ-TREE tree inference | Implemented as the tree-preparation step |
| PhyloBayes | Not yet implemented |
| MrBayes | Not yet implemented |
| MCMCTree | Not yet implemented |

## GUI workflow

1. Load a tree, or infer one with IQ-TREE from an alignment.
2. Configure tip dates and node calibrations in the shared dating panel.
3. Optionally export the calibration table as a text file.
4. Select LSD2 and run; the dated tree and `.result` files are returned.
5. Inspect the console and send the dated tree to the workspace.

## Calibration inputs

- Tip dates: one date per sampled sequence.
- Node calibrations: minimum/maximum or distribution-based bounds on internal
  nodes.
- Calibration tables are written to a temporary file and passed to LSD2 with
  `-d`.

## Pitfalls

- Only LSD2 is currently wired end-to-end; the other method entries report
  "not implemented" rather than silently running.
- The tree must be rooted before dating; PD-Guide can use a chosen root or ask
  LSD2 to estimate it.
- Calibrations must be compatible with the tree taxa; unmatched labels are not
  applied.
- Dating requires branch lengths in the input tree.
- Export the calibration table if you need to reuse it with MrBayes/MCMCTree
  later.
