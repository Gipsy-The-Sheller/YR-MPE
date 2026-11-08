# RAxML-NG skill

## Purpose

Infer a maximum-likelihood tree from an aligned FASTA/PHYLIP file. The plugin
can run a search only, or an all-in-one search plus bootstrap analysis, and
imports the best tree and bootstrap tree set.

## Preconditions

- An aligned file with at least four sequences.
- A RAxML-NG model string, for example `GTR+G4` or `LG+G4`.
- Optionally, a green-selected model item from ModelFinder or ModelTest-NG.

## GUI to CLI

The plugin's `control_bindings()` is the source of truth. The expected mapping
is:

| GUI control | CLI option | Notes |
|---|---|---|
| File input | `--msa` | One RAxML-NG process per selected file |
| Model | `--model` | RAxML-NG syntax; `+FQ` is normalized to `+FE` |
| Threads | `--threads` | |
| Bootstrap replicates | `--bs-trees` | Only with `--all`; `0` means no bootstrap unless autoMRE is enabled |
| autoMRE convergence | `--bs-trees autoMRE{N}` | `N` is the replicate limit; `1000` when the spin box is 0 |
| Starting trees | `--tree` | Empty/default omits the flag; `pars{10}`, `rand{10}`, `rand{1}` are accepted |
| Outgroup | `--outgroup` | Comma-separated taxa; drawing/annotation option in RAxML-NG |
| (fixed) | `--search` or `--all` | `--all` when bootstrap or autoMRE is active |
| (fixed) | `--prefix`, `--redo` | Plugin-generated output prefix |

Verify options with:

```text
raxml-ng --help
```

## Handoff from model selection

1. Run ModelFinder or ModelTest-NG and import its model report.
2. Click the model item in the workspace so it is green (or add it to Selected).
3. Open `PHYLOGENY -> Maximum Likelihood (ML) -> RAxML-NG`.
4. The plugin calls `set_model(best_model)`; for example `TrN+I` or `TPM2uf`.

## Outputs

- `*.raxml.bestTree` - best ML tree.
- `*.raxml.support` - tree with bootstrap support when `--all` is used.
- `*.raxml.bootstraps` - bootstrap replicate trees, imported as a `treeset`.
- `*.raxml.bestModel`, `*.raxml.log` - retained as run artifacts.

## Pitfalls

- The Windows MSYS build strips braces in `--tree`; the plugin escapes them
  internally. Do not pass already-escaped braces through the GUI.
- RAxML-NG treats `--outgroup` as a drawing option, not a constraint on the
  inference itself.
- Bootstrap runs are much slower; start with a search-only run when checking a
  new alignment.
