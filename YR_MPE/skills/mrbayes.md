# MrBayes

## Purpose

Bayesian phylogenetic inference with MrBayes 3.2.7a, including partitioned
models, mixed protein models, MCMC diagnostics and optional dating priors.
The Windows MPI+BEAGLE build is exposed as a separate plugin entry.

## How the plugin runs

The plugin writes a `begin mrbayes; ... end;` block into a NEXUS file and runs:

```text
mb -i input.nex
```

With MPI enabled it runs `mpirun -np N mb -i input.nex`. If the input already
contains a MRBAYES block and **Run MrBayes Data Block directly if present** is
checked, that block is used instead of generating a new one.

## GUI to NEXUS

| GUI control | NEXUS directive |
|---|---|
| Data type = DNA | `lset nucmodel=4by4;` |
| Data type = Protein | `lset nucmodel=protein;` |
| DNA substitution number | `lset nst=N;` |
| DNA state frequencies | `prset statefreqpr=dirichlet(1,1,1,1);`, `fixed(equal)`, or `fixed(empirical)` |
| Protein model | `prset aamodelpr=fixed(model);` or `prset aamodelpr=mixed;` |
| Protein state frequencies | `prset statefreqpr=...;` |
| Rate heterogeneity | `lset rates=equal\|gamma\|invgamma\|propinv\|lnorm\|adgamma;` |
| Gamma categories | `ngammacat=N` in the `lset rates=...` line |
| MCMC settings | `mcmcp ngen=N samplefreq=N nchains=N nruns=N printfreq=N savebrlens=yes checkpoint=yes checkfreq=5000;` then `mcmc;` |
| Consensus type | `sumt contype=majorityrule\|halfcompat\|allcompat;` |
| Consensus format | `conformat=simple\|figtree` |
| Burn-in as fraction | `sump relburnin=yes burninfrac=X;` and `sumt relburnin=yes burninfrac=X;` |
| Burn-in as states | `sump burnin=N;` and `sumt burnin=N;` |
| Partition mode | `charset`/`partition` definitions in NEXUS, then `prset applyto=(all)` as appropriate |
| Dating | `prset brlenspr=clock:uniform;` plus clock and calibration priors |
| BEAGLE | `set usebeagle=yes beagledevice=... beagleprecision=... beaglescaling=...;` |

## Flow

1. Import a NEXUS/FASTA alignment, or a NEXUS file with an existing MRBAYES
   block.
2. Set data type, model, rate heterogeneity, MCMC length and consensus
   options.
3. For partitioned data, enable partition mode and configure partitions before
   running.
4. Run; inspect the `.con.tre`, `.p` and `.mcmc` outputs; import the consensus
   tree.

## Pitfalls

- The generated block must contain `mcmc;`; without it MrBayes only reads the
  data.
- `sump` and `sumt` must run after `mcmc`; the plugin appends them
  automatically.
- Burn-in is either a fraction (`relburnin=yes burninfrac=`) or a state count
  (`burnin=`), never both.
- Partition definitions must live inside the NEXUS file; GUI partition mode
  only writes them when configured.
- The supermatrix uses the union of taxa across partitions; missing taxa are
  padded with `?` for that partition, not removed.
- An existing MRBAYES block is not silently replaced: uncheck the direct-run
  option if the GUI settings should regenerate it.
- The MPI entry needs `mpirun` in addition to `mb`; the standard entry runs the
  serial executable.
