# Model Parameter Estimation

## Purpose

Estimate substitution-model parameters with IQ-TREE and export the fitted
model (Q matrix, base frequencies, rate parameters) as a workspace model item.
It can also parse existing `.iqtree` reports without rerunning IQ-TREE.

## GUI to CLI

`control_bindings()` is the source of truth. The plugin builds a single `-m`
string, appending rate suffixes to the base model.

| GUI control | CLI option | Notes |
|---|---|---|
| Input alignment | `-s PATH` | One IQ-TREE process per file in batch mode |
| Substitution model | `-m MODEL` | `auto` becomes `MFP` |
| Discrete Gamma | `+G` | Appended inside `-m` |
| Gamma categories | `+G{n}` | Appended only when +G is checked |
| Invariant sites | `+I` | Appended inside `-m` |
| Empirical frequencies | `+F` | Appended inside `-m` |
| Threads | `-nt N` | Omitted when 1 |
| Load existing `.iqtree` | none | GUI-only parsing path |

The runner always appends `-redo`.

## Flow

1. Import an alignment or a set of alignments.
2. Choose a model, or keep `auto` (MFP).
3. Set Gamma/+I/+F if needed.
4. Run; inspect the parameter table and send the model to the workspace.

## Pitfalls

- `+G`, `+I` and `+F` are part of the same `-m` token, not separate options.
- `auto` is a GUI label; the command uses `MFP`.
- Batch mode runs IQ-TREE once per input file, so a shared model is not fitted
  across loci.
- Loading an existing `.iqtree` file skips estimation and does not probe the
  sequence data.
