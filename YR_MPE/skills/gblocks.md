# GBlocks

## Purpose

Select conserved alignment blocks with GBlocks 0.91b.

## GUI to CLI

The runner calls `GBlocks <input> [options] -p=n -q=n -r=n -e=-gb`.

| GUI control | CLI option | Notes |
|---|---|---|
| Input alignment | positional | First argument |
| Sequence type | `-t=p`, `-t=d`, `-t=c` | Auto adds nothing |
| Minimum block length | `-b1=N` | Omitted at 0 |
| Minimum conserved positions | `-b2=N` | Omitted at 0 |
| Max contiguous non-conserved | `-b3=N` | Omitted at 8 |
| Minimum flank positions | `-b4=N` | Omitted at 10 |
| Gap positions | `-b5=n\|h\|a` | None is omitted |
| Similarity matrix | `-b6=y\|n` | Always emitted |
| Output file | `-e=-gb` | Plugin-managed extension |
| Batch/report flags | `-p=n -q=n -r=n` | Always added |

## Flow

1. Import an alignment.
2. Choose sequence type and block parameters.
3. Run; send the `.gb` alignment to the workspace.

## Pitfalls

- `-b1`/`-b2` are omitted at 0, not passed as zero; GBlocks then uses its own
  defaults.
- `-b3` and `-b4` have non-zero defaults (8 and 10); only changed values are
  emitted.
- `-b6` is always emitted, including `-b6=n` when similarity matrices are off.
- The plugin suppresses interactive prompts with `-p=n -q=n -r=n`.
- GBlocks output is a block-selected alignment, not a tree.
