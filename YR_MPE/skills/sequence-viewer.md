# Sequence Viewer

## Purpose

View and edit a multiple-sequence alignment in a spreadsheet-style window, and
launch alignment/trimming plugins on the current data.

## GUI features

- Open FASTA/alignment files.
- Insert sequences from another file.
- Save changes, Save As, and export the current alignment.
- Switch nucleotide/protein colouring and font.
- Launch MAFFT, Clustal Omega, MUSCLE5, TrimAl or GBlocks from the viewer.
- Navigate by sequence and column position.

There is no command-line interface; the viewer also serves as a lightweight
editor for alignments produced by the plugins.

## Workflow

1. Open an alignment (or receive one from a plugin).
2. Inspect/align columns, insert missing sequences if needed.
3. Optionally launch an aligner or trimmer from the toolbar.
4. Save changes or Save As before importing elsewhere.

## Pitfalls

- Editing is in-memory until Save/Save As; closing without saving loses edits.
- Inserting sequences can change the alignment width; check that all rows stay
  rectangular.
- The viewer is not a phylogeny tool; use the tree viewers for Newick/NEXUS.
- Launching a plugin from the viewer passes the current alignment, but the
  returned result is a new item; it does not overwrite the viewer silently.
- Colour schemes are display-only and do not affect the exported sequences.
