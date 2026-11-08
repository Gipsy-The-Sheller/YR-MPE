# IQ-TREE ModelFinder skill

## Purpose

Run IQ-TREE ModelFinder (`-m TESTONLY`) and import the best substitution model
as a workspace `model` item.

## GUI to CLI

| GUI control | CLI option | Notes |
|---|---|---|
| File input | `-s` | One IQ-TREE process per alignment (or a supermatrix in partition mode) |
| Sequence type | `-st DNA` or `-st AA` | `AUTO` omits `-st` |
| Criterion | default BIC; `-AICc` / `-AIC` | BIC is IQ-TREE's default |
| Threads | `-nt` | Omitted when 1/AUTO |
| Model set | `-mset ALL` | The plugin searches all models |
| Partition mode | `-p` / `-q` / `-m` scheme | Built by the partition thread |

The core command is `iqtree3 -s <alignment> -m TESTONLY <options>`.

## Flow

1. Select an alignment or open `MODEL -> Find Best Substitution Model (ML)`.
2. Run and inspect the model table.
3. Import the model item.
4. Select the model and open IQ-TREE or RAxML-NG; the host applies the best
   model when the naming is compatible.

## Pitfalls

- ModelFinder's model names use IQ-TREE syntax (`GTR+F+G4`); RAxML-NG accepts
  most of it but `+FQ` is normalized to `+FE`.
- Partition mode uses a different thread and can take much longer than a
  single-alignment search.
- In partition mode the supermatrix taxa are the union across all selected
  loci; taxa missing from a locus are padded with `?`, not dropped.
