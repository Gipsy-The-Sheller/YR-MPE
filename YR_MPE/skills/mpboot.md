# MPBoot

## Purpose

Maximum-parsimony tree inference with MPBoot's ultrafast bootstrap and
consensus reconstruction.

## GUI to CLI

`control_bindings()` covers the top-level controls. Advanced Parameters are
collected in a dialog; their flags are listed below.

| GUI control | CLI option | Notes |
|---|---|---|
| Input alignment | `-s PATH` | One MPBoot run per file in batch mode |
| Sequence type | `-st DNA\|AA\|CODON\|BIN\|MORPH` | AUTO is omitted |
| Threads | `-nt N` | Omitted when 1 |
| Ultrafast bootstrap | `-bb N` | Enabled by the UFBoot checkbox |
| UFBoot max iterations | `-nm N` | |
| UFBoot stopping rule | `-nstep N` | |
| UFBoot min correlation | `-bcor X` | |
| RELL epsilon | `-beps X` | |
| Consensus threshold | `-t X` | Omitted at 0 (extended consensus); 0.5 = majority-rule |
| Consensus burn-in | `-bi N` | Omitted at 0 |
| Consensus tree | `-con` | Always appended by the plugin |

Advanced Parameters (dialog values are applied when the dialog is accepted):

| GUI control | CLI option |
|---|---|
| Initial trees | `-numpars N` |
| Best trees | `-toppars N` |
| Candidate trees | `-numcand N` |
| Iterations = fixed | `-n N` |
| Perturbation strength | `-pers X` |
| Method = IQP | `-iqp` |
| Method = IQPNNI (old) | `-iqpnni` |
| Ratchet iterations | `-ratchet_iter N` |
| Ratchet weight | `-ratchet_wgt N` |
| Ratchet percent | `-ratchet_percent N` |
| Turn off Ratchet | `-ratchet_off` |
| SPR radius | `-spr_rad N` |
| Candidate cutoff | `-cand_cutoff N` |
| Turn off tree refinement | `-opt_btree_off` |
| Use NNI instead of SPR | `-nni_pars` |

## Flow

1. Import an alignment.
2. Choose sequence type and threads.
3. Enable UFBoot and set the replicate/tuning values if needed.
4. Set the consensus threshold; run.
5. Send the `.treefile`/`.contree` result to the workspace.

## Pitfalls

- `-con` is always emitted, so a consensus tree is produced even when the
  threshold is 0.
- `-bb` should be at least 1000; the GUI exposes the tuning options separately.
- Advanced dialog values only take effect after accepting the dialog.
- MPBoot uses the input alignment name for output; batch mode therefore writes
  one `.treefile`/`.contree` per input.
