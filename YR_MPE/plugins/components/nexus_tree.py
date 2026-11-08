"""Read trees from NEXUS files written by MrBayes / PhyloBayes.

Biopython's NEXUS reader rejects the ``[&prob=...]`` annotations MrBayes puts
on every node, so this helper strips bracket comments, applies the file's
``translate`` table and returns plain Newick strings that the rest of YR-MPE
can consume.
"""

from __future__ import annotations

import io
import re
from pathlib import Path
from typing import Dict, List

from Bio import Phylo

_COMMENT_RE = re.compile(r"\[[^\]]*\]")
_INTERNAL_PROB_RE = re.compile(r"\)\[&prob=([0-9eE+\-.]+)[^\]]*\]")
_TRANSLATE_RE = re.compile(r"\btranslate\b(.*?);", re.IGNORECASE | re.DOTALL)
_TREE_RE = re.compile(
    r"\btree\s+\S+\s*=\s*(?:\[[^\]]*\]\s*)?(\(.*\))\s*;",
    re.IGNORECASE | re.DOTALL,
)


def _parse_translate(text: str) -> Dict[str, str]:
    match = _TRANSLATE_RE.search(text)
    if not match:
        return {}
    table: Dict[str, str] = {}
    for entry in match.group(1).split(","):
        fields = entry.split(None, 1)
        if len(fields) != 2:
            continue
        key, value = fields[0].strip(), fields[1].strip().strip("'\"")
        if key:
            table[key] = value
    return table


def read_nexus_trees(path: str) -> List[str]:
    """Return every tree in *path* as a Newick string."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    return _read_nexus_text(text)


def read_tree_file(path: str) -> List[str]:
    """Read trees from a NEXUS file or a plain-Newick file."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    stripped = text.lstrip()
    if stripped.upper().startswith("#NEXUS") or "begin trees" in text.lower():
        return _read_nexus_text(text)
    trees: List[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            tree = Phylo.read(io.StringIO(line), "newick")
        except Exception:
            continue
        handle = io.StringIO()
        Phylo.write(tree, handle, "newick")
        trees.append(handle.getvalue().strip())
    return trees


def _read_nexus_text(text: str) -> List[str]:
    """Parse all tree statements out of NEXUS *text*."""
    translate = _parse_translate(text)
    trees: List[str] = []
    for match in _TREE_RE.finditer(text):
        newick = match.group(1)
        # MrBayes stores posterior probabilities in comments appended to
        # internal nodes. Promote them to node labels so support survives the
        # conversion to plain Newick; leaf probabilities are discarded.
        newick = _INTERNAL_PROB_RE.sub(
            lambda m: ")" + ("%.6g" % float(m.group(1))), newick
        )
        newick = _COMMENT_RE.sub("", newick).strip()
        if not newick:
            continue
        try:
            tree = Phylo.read(io.StringIO(newick), "newick")
        except Exception:
            continue
        if translate:
            for terminal in tree.get_terminals():
                if terminal.name in translate:
                    terminal.name = translate[terminal.name]
        handle = io.StringIO()
        Phylo.write(tree, handle, "newick")
        trees.append(handle.getvalue().strip())
    return trees
