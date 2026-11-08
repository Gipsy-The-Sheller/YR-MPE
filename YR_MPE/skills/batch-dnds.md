# Batch dN/dS

## Purpose

Estimate dN/dS across many codon alignments in one run, using either fast
counting methods or IQ-TREE codon-model MLE.

## Methods

| GUI method | Implementation | Output |
|---|---|---|
| Counting — NG86 | Nei-Gojobori counting | Pairwise omega table |
| Counting — YN00 | Yang-Nielsen F84 correction, one kappa per alignment | Pairwise omega table |
| MLE — IQ-TREE GY | `iqtree -s ALN --seqtype CODON -m GY` | Global omega + logL |
| MLE — IQ-TREE MG | `iqtree -s ALN --seqtype CODON -m MG` | Global omega + logL |

## GUI workflow

1. Select an input folder of codon-aligned FASTA files.
2. Choose a method and genetic code.
3. For MLE, choose an output folder and threads.
4. Run; inspect the batch table and export TSV.

## MLE command

```text
iqtree -s alignment.fasta --seqtype CODON -m GY|MG -T N --prefix OUT/stem_gy|mg --redo
```

The plugin parses omega and log-likelihood from the IQ-TREE report/stdout.

## Output

- One row per alignment: method, omega, logL, status, artifact path.
- Failed alignments are marked and kept in the table.
- TSV export of the batch table.

## Pitfalls

- Counting methods skip stop/ambiguous codons pairwise; MLE requires clean,
  in-frame codon alignments.
- NG86/YN00 report pairwise omega; GY/MG report one global omega per
  alignment.
- The MLE output folder is required; picking a file instead of a folder is not
  supported.
- IQ-TREE must be configured before MLE runs.
- Genetic code applies to both counting and MLE paths.
