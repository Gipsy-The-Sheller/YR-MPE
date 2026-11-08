# ML distance calculation

## Purpose

Compute a maximum-likelihood pairwise distance matrix with IQ-TREE and import
it as a distance item. The matrix can feed DecentTree or other distance-based
tree builders.

## GUI to CLI

The runner calls:

```text
iqtree3 -s alignment [-st DNA|AA] [-m MODEL] [-nt N] -redo
```

| GUI control | CLI option | Notes |
|---|---|---|
| Input alignment | `-s PATH` | One IQ-TREE run per file in batch mode |
| Sequence type | `-st DNA\|AA` | AUTO is omitted |
| Substitution model | `-m MODEL` | `auto` omits `-m` |
| Discrete Gamma | `+G{n}` | Appended to `-m` |
| Invariable sites | `+I` | Appended to `-m` |
| FreeRate | `+R` | Appended to `-m`; category count is not exposed here |
| Ascertainment bias | `+ASC` | Appended to `-m` |
| State frequencies | `+F`, `+FO`, `+FQ` | Estimated adds nothing |
| Threads | `-nt N` | Omitted when 1 |
| Rerun | `-redo` | Always added |

IQ-TREE writes `<alignment>.mldist`; the plugin parses it into a distance
matrix item.

## Flow

1. Import one or more alignments.
2. Choose sequence type and model.
3. Run; inspect the distance table and send it to the workspace.

## Pitfalls

- `auto` does not mean ModelFinder here: the plugin simply omits `-m` and lets
  IQ-TREE use its distance default.
- `+R` in this plugin does not take a category count; use the IQ-TREE plugin
  when you need `+R{n}`.
- The output is a distance matrix, not a tree.
- Batch mode writes one `.mldist` per input.
- Uncorrected p-distances are a separate plugin and do not use IQ-TREE.
