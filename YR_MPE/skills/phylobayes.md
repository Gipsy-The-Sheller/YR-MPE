# PhyloBayes-MPI

## Purpose

Bayesian phylogenetic inference with PhyloBayes-MPI, including CAT mixture
models, discrete Gamma rates and site-heterogeneous analyses.

## GUI to CLI

The runner builds an `mpirun` command:

```text
mpirun -np N pb_mpi -d alignment [model options] chain_name
```

| GUI control | CLI option | Notes |
|---|---|---|
| Alignment | `-d PATH` | PHYLIP input |
| MPI processes | `mpirun -np N` | |
| Parallel chains | none | GUI-only: launches N independent chains |
| CAT mixture | `-cat` / `-ncat N` | Infinite or finite CAT |
| nCAT components | `-ncat N` | Enabled by the CAT checkbox |
| Discrete Gamma | `-dgam N` | |
| Substitution model | `-gtr`, `-poisson`, `-jtt`, `-wag`, `-lg` | |
| Sampling frequency + generations | `-x FREQ GEN` | Emitted together |
| Chain name | positional | Plugin-generated from the input name |

Chain post-processing uses `bpcomp`:

| GUI control | bpcomp option |
|---|---|
| Consensus cutoff | `-c X` (0.5 default) |
| Burn-in | `-x BURNIN` |
| Sample every | `-x BURNIN EVERY` |
| Stop at | `-x BURNIN EVERY UNTIL` |

## Flow

1. Import a PHYLIP alignment.
2. Set MPI processes, CAT/Gamma and substitution model.
3. Run one or more chains.
4. Select the chain files and run `bpcomp` to obtain the consensus tree.

## Pitfalls

- CAT is enabled by default (`-cat`); switch to finite `-ncat` explicitly.
- Gamma and CAT can be combined; both checkboxes are on by default.
- The MCMC generation control is combined with the sampling frequency into a
  single `-x` pair, not two separate flags.
- `-1` generations is converted to a large finite value for the CLI.
- Chain analysis uses `bpcomp`, not `pb_mpi`; its cutoff/burn-in options are
  separate from the MCMC flags.
- PhyloBayes requires PHYLIP input; convert FASTA first.
