# MAFFT skill

## Purpose

Align nucleotide or protein FASTA files and import each alignment into the
workspace.

## GUI to CLI

`control_bindings()` is the source of truth. The expected mapping is:

| GUI control | CLI option | Notes |
|---|---|---|
| File input | positional | Appended after all options |
| Alignment mode | `--auto`, `--localpair`, `--genafpair`, `--globalpair`, `--retree` | High-accuracy modes also add `--maxiterate 1000` |
| Threads | `--thread` | `-1`/automatic is omitted |
| Direction adjustment | `--adjustdirection` | Detects and reverse-complements sequences that appear reversed; marks them with `_R_` |
| Direction adjustment (accurate) | `--adjustdirectionaccurately` | Slower, more reliable for divergent or long reads; also marks reversed records with `_R_` |
| Advanced penalties | `--op`, `--ep`, `--maxiterate`, `--clustalout`, `--reorder`, ... | Set in Advanced Parameters |
| Output suffix | none | GUI-only; MAFFT writes to stdout and the plugin saves it beside the input |

Verify options with:

```text
mafft --help
```

## Flow

1. Select an alignment/sequence or open `ALIGN -> Align by... -> MAFFT`.
2. Choose mode and threads; run.
3. Import the resulting alignment. When MAFFT was launched from a workspace
   dataset, each alignment updates the matching item in that same dataset
   (the dataset itself is not duplicated).
4. Continue with TrimAl/GBlocks, ModelTest-NG or ModelFinder, then RAxML-NG/
   IQ-TREE.

## Pitfalls

- MAFFT writes to stdout; the plugin, not MAFFT's `-o`, creates the output
  file.
- Batch mode runs one MAFFT process per input file.
- Advanced parameters are only included when they differ from MAFFT defaults.
- With direction adjustment enabled, MAFFT prefixes reverse-complemented
  records with `_R_`. The plugin strips this prefix from output headers so
  downstream taxon matching (OrthoFinder/concatenation) is not broken; the
  sequences themselves remain reverse-complemented. Do not re-add `_R_` when
  constructing taxon labels.
- Direction adjustment is off by default. Enable it for Sanger amplicons or
  other data of unknown orientation.
