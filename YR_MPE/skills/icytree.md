# IcyTree

## Purpose

Interactive tree viewing and annotation in an embedded IcyTree web view.

## Inputs

- Newick (`.nwk`, `.newick`, `.tre`, `.tree`)
- NEXUS (`.nex`, `.nexus`) with a trees block

NEXUS is preferred when the tree carries rich annotations, because IcyTree can
display support values, branch lengths and metadata that plain Newick may lose.

## GUI workflow

1. Open IcyTree from the Phylogeny menu or send a tree/TreeSet to it.
2. The tree is loaded into the embedded viewer.
3. Use IcyTree's own controls to reroot, change layout, colour by support and
   export.
4. Downloads are offered as SVG, PNG, JPEG, Newick, NEXUS or XML.

There is no command-line interface; IcyTree runs entirely in the embedded web
view.

## Pitfalls

- IcyTree is a viewer, not an inference tool; it does not modify the workspace
  tree unless you export and re-import it.
- Very large trees may be slow in the embedded browser.
- Export format matters: SVG/PNG for figures, Newick/NEXUS for downstream
  software.
- Rooting in the viewer is display-side unless the exported tree is re-imported.
- The viewer needs the bundled IcyTree HTML assets; if they are missing, the
  plugin shows an error page instead of a tree.
