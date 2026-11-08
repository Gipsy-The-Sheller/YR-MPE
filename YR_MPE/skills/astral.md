# ASTRAL-IV / wASTRAL

## Purpose

Infer a species tree from a set of unrooted gene trees. YR-MPE wraps the
ASTER executables:

- **ASTRAL-IV** (`astral4`) — standard quartet-based species-tree inference
  with CASTLES-II branch lengths.
- **wASTRAL** (`wastral`) — weighted ASTRAL; weights quartets by gene-tree
  support/length and reports coalescent-unit branch lengths.
- **ASTRAL-Hybrid** (`astral-hybrid`) — a synonym/preset of wASTRAL using
  support values and branch lengths; not a separate method.

## Citation

- ASTER package: `zhang2025aster`.
- wASTRAL weighting: `zhang2022wastral` (Zhang & Mirarab, 2022).

## GUI to CLI

The first argv token selects the executable; the remaining options are passed
to it. `control_bindings()` is the source of truth.

| GUI control | CLI option | Notes |
|---|---|---|
| Input gene trees | positional | Appended last |
| ASTRAL version | `astral4`, `wastral`, `astral-hybrid` | First argv token |
| Threads | `-t N` | |
| Output tree | `-o PATH` | Optional; a temporary `.tre` is used if empty |
| Support option | `-u 0..3` | 0 no support, 1 length+support, 2 detailed, 3 freqQuad.csv |
| Constraint tree | `-c PATH` | Advanced |
| Guide trees | `-g PATH` | Advanced |
| Rounds | `-r N` | Advanced; default 4 |
| Subsample rounds | `-s N` | Advanced; default 4 |
| Gene-to-taxon map | `-a PATH` | Advanced |
| Root taxon | `--root NAME` | Advanced |
| Subsample proportion | `--proportion X` | Advanced; default 0.25 |
| Seed | `--seed N` | Advanced; default 233 |
| Yule lambda | `-l X` | Advanced; default 0.5 |
| Downweight repeats | `-w X` | Advanced; default 1.0 |
| Verbose level | `-v 1\|2` | Advanced; default 1 |
| Scoring only | `-C` | |
| More rounds | `-R` | |
| Branch lengths (ASTRAL-IV) | `--length SULength\|CULength` | ASTRAL-IV only |
| Gene length (ASTRAL-IV) | `--genelength X` | ASTRAL-IV only |
| Hybrid preset | `-B`, `-S`, `-L` | ASTRAL-Hybrid only |
| Tree weights | `--treeweights PATH` | ASTRAL-Hybrid only |
| Weighting mode | `--mode 1..4` | ASTRAL-Hybrid only |
| Min/max/default support | `-n`, `-x`, `-d` | ASTRAL-Hybrid only |

## Flow

1. Provide one Newick file with one gene tree per line. Opening ASTRAL from
   Selected with Phylogeny/TreeSet objects preloads their trees; the
   **Import from Selected TreeSet** button is the fallback and also accepts
   plain phylogenies.
2. Choose ASTRAL-IV, wASTRAL or ASTRAL-Hybrid.
3. Set support option, threads and output path.
4. Run; send the resulting tree to the workspace.

## Weighted-mode support labels

ASTER's weighted parsers (`wastral`, `astral-hybrid`) call `std::stod` on
node labels. IQ-TREE can write composite labels such as `/1/100`
(SH-aLRT/aBayes/UFBoot). The leading slash makes `std::stod` throw
`std::invalid_argument` and the process aborts.

YR-MPE normalizes composite labels before weighted runs: it keeps the first
numeric field (aBayes when present, otherwise aLRT/UFBoot) and drops the
others. A single support type such as aBayes-only also avoids the crash.
ASTRAL-IV accepts composite labels directly.

## Pitfalls

- ASTER writes the tree to `-o`; stdout contains only logs.
- wASTRAL leaves terminal branch lengths empty in coalescent units; the
  plugin adds display-only `1e-5` terminal lengths for viewers.
- `--length` and `--genelength` are ASTRAL-IV options; do not pass them to
  wASTRAL or ASTRAL-Hybrid.
- Constraint and guide trees must be binary Newick trees.
