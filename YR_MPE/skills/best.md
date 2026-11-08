# BEST / MrBayes Multi-Species Coalescent

BEST is the full-likelihood multispecies-coalescent analysis implemented by
MrBayes. YPhylo Workbench exposes it as a block-oriented agent surface:

- `best.inspect_options` returns the current structured options and the
  canonical MRBAYES block.
- `best.apply_block` accepts a complete MRBAYES block. The form is updated
  from the parsed values and the response returns the canonical block plus
  any lines that were not parsed.
- `best.modify_options` sets individual controls directly.
- `best.run` starts MrBayes through the task queue. `best.stop` cancels it.
- `best.export_to_workspace` imports the consensus tree and chain files.

The input must be a NEXUS data file. Partitions may be loaded from charset
definitions; the species partition uses MrBayes `speciespartition` syntax.
A round-trip response should be checked with `unparsed_lines` before running.
