# Tree Merger

## Purpose

Collect trees from multiple Newick files into one tree file, usually to prepare
a gene-tree set for ASTRAL/wASTRAL or another supertree method.

## GUI workflow

1. **Add Files** selects one or more Newick files.
2. **Remove Selected** / **Clear All** edit the file list.
3. Set the output file.
4. **Merge Trees** writes every non-empty line from every input file, one tree
   per line, to the output.

There is no external command-line tool; the merge is performed in Python.

## Output contract

- Output is plain Newick text, one tree per line.
- Files are processed in the order shown in the file table.
- If a file contains several trees, each non-empty line is appended.
- The merger does not reroot, deduplicate, validate, or relabel trees.

## Pitfalls

- Lines that are not Newick are copied verbatim; check the output before
  running ASTRAL.
- A NEXUS tree block is not parsed; use plain Newick files.
- Taxon labels must already be consistent across files.
- The output order follows the table order, so reorder files if the downstream
  tool depends on locus order.
