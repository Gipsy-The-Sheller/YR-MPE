# MACSE2

## Purpose

Codon-aware alignment of coding sequences with MACSE2 2.07. It preserves
reading frames, handles frameshifts/stops, and writes both nucleotide and
amino-acid alignments.

## GUI to CLI

The runner builds:

```text
macse -prog alignSequences -seq <input> [options] -out_NT <nt> -out_AA <aa>
macse -prog refineAlignment -align <input> [options] -out_NT <nt> -out_AA <aa>
```

| GUI control | CLI option | Notes |
|---|---|---|
| Mode = Align sequences | `-prog alignSequences -seq PATH` | |
| Mode = Refine existing alignment | `-prog refineAlignment -align PATH` | |
| Genetic code | `-gc_def N` | |
| AA alphabet | `-alphabet_AA NAME` | Default SE_B_8 |
| Score matrix | `-score_matrix NAME` | BLOSUM62, VTML200_BIS, VTML240 |
| Optimization cycles | `-optim N` | Default 2 |
| NT output | `-out_NT PATH` | Plugin-managed |
| AA output | `-out_AA PATH` | Plugin-managed |
| Preferred result | none | GUI-only: selects NT or AA for import |

Advanced parameters are free text and are appended verbatim to the MACSE
command after shell-style splitting.

## Post-processing (GUI-only)

- **Mask exceptional characters** rewrites `!` to `?`.
- **Mask terminal stop codon** replaces a terminal stop with `---`.
- **Drop all-gap codon columns** removes columns where every sequence has a
  gap.

These are applied by the plugin after MACSE, not by MACSE flags.

## Flow

1. Import nucleotide FASTA files.
2. Choose align or refine mode.
3. Set genetic code, alphabet, matrix and optimization.
4. Run; send the NT or AA alignment to the workspace.

## Pitfalls

- Align mode expects unaligned CDS; refine mode expects an existing codon
  alignment.
- The genetic code is passed by its NCBI number, not the display label.
- MACSE writes both NT and AA files; the preferred-result combo only controls
  which one is selected for import.
- Free-text advanced parameters are passed through without validation.
