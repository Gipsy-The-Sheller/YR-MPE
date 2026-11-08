# MiniTracer skill

## Purpose

Inspect MCMC chain files without an external Tracer installation: trace plots,
ESS, 95% HPD, harmonic means, skyline plots and two-chain comparisons.
MiniTracer is a Python implementation of Tracer with a redesigned interface
layout, not a pixel-for-pixel port of the Java application. Cite Tracer as
`rambaut2018posterior` when reporting trace diagnostics.

## Inputs

- MrBayes `.p` files (tab-delimited chain parameters).
- PhyloBayes `.trace` files.
- Plain tab-delimited trace files with a header row.

## GUI workflow

1. Open `PHYLOGENY -> Bayesian Inference (BI) -> MCMC Diagnostics (MiniTracer)`.
2. Add one or more chain files; set burn-in.
3. Select parameters in the statistics table (Ctrl/Shift for multiple).
4. Inspect the summary form, box plot or skyline plot.
5. Open **Compare 2 Chains** to select two traces; the scatter panel shows
   chain 1 against chain 2 with the y = x reference line, and the statistics
   table reports sample counts, means, the correlation coefficient and the
   Gelman--Rubin PSRF (R-hat). The harmonic-mean interval plot sits below the
   right-hand statistics table.

## Notes for agents

MiniTracer is a GUI analysis tool, not a command-line binary. It consumes chain
items from the workspace and produces no new file artifact unless the caller
exports one. Use `workbench.workspace` to discover `chain` items and
`workbench.selection` to see whether a chain can be sent to MiniTracer.

## Interpretation

- ESS should be checked before trusting posterior summaries. The statistics
  table colours rows green at ESS >= 200, yellow at 100--200 and red below 100.
- The two-chain PSRF is a per-parameter Gelman--Rubin R-hat computed from the
  two post-burn-in traces after truncating to the shorter length; values close
  to 1 indicate agreement. It is not a multivariate PSRF over the whole model.
- The harmonic-mean panel is descriptive; its approximate delta-method test is
  not a substitute for a proper marginal-likelihood comparison.
- Skyline plots connect the selected parameters across posterior samples; the
  dashed bounds are 95% HPD.
