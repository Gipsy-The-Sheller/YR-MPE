#!/usr/bin/env python3
"""Generate the human-readable tool reference from the agent skill files.

The skill Markdown is deliberately regular: headings, paragraphs, bullet or
numbered lists, fenced code, and pipe tables. This converter handles that
subset and writes ``docs/manual/tool-reference.tex``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SKILLS = REPO / "YR_MPE" / "skills"
OUT = REPO / "docs" / "manual" / "tool-reference.tex"

GROUPS = [
    ("Alignment and trimming", [
        "mafft", "macse", "muscle5", "clustal-omega", "trimal", "gblocks",
    ]),
    ("Sequence and dataset management", [
        "orthofinder", "sequence-viewer", "seqmatrix", "seqdbg",
        "sequence-concatenator", "tree-merger",
    ]),
    ("Model selection", [
        "modelfinder", "modeltest-ng", "model-parameters",
    ]),
    ("Distance and saturation", [
        "ml-distance", "uncorrected-distance", "saturation", "decenttree",
    ]),
    ("Phylogenetic inference", [
        "iqtree", "raxmlng", "mrbayes", "phylobayes", "mpboot",
    ]),
    ("Species trees and reconciliation", [
        "astral", "caster-site", "treerecs", "best-mrbayes",
    ]),
    ("Dating and selection", [
        "lsd2", "pdguide", "codon-ztest", "batch-dnds",
    ]),
    ("Visualization and diagnostics", [
        "icytree", "minitracer", "rf", "console",
    ]),
]

UNICODE = {
    "→": r"$\rightarrow$",
    "←": r"$\leftarrow$",
    "≥": r"$\geq$",
    "≤": r"$\leq$",
    "×": r"$\times$",
    "–": "--",
    "—": "---",
    "…": r"\ldots{}",
    "±": r"$\pm$",
    "µ": r"$\mu$",
    "°": r"$^\circ$",
    "·": r"$\cdot$",
    "≈": r"$\approx$",
    "\u00a0": "~",
}

SPECIAL = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}


def latex_escape(text: str) -> str:
    text = "".join(UNICODE.get(ch, ch) for ch in text)
    return "".join(SPECIAL.get(ch, ch) for ch in text)


def inline(text: str) -> str:
    """Convert inline Markdown, protecting code spans and links first."""
    placeholders = []

    def hold(value: str) -> str:
        token = "\uE100%d\uE101" % len(placeholders)
        placeholders.append(value)
        return token

    text = text.replace("\uE000", "|")

    def code(match):
        return hold(r"\texttt{%s}" % latex_escape(match.group(1)))

    text = re.sub(r"`([^`]+)`", code, text)

    def link(match):
        label, url = match.group(1), match.group(2)
        return hold(r"\href{%s}{%s}" % (url, inline(label)))

    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
    text = latex_escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", text)
    text = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"\\emph{\1}", text)
    for index, value in enumerate(placeholders):
        text = text.replace("\uE100%d\uE101" % index, value)
    return text


def table(lines: list[str]) -> str:
    rows = []
    for line in lines:
        line = line.strip().replace("\\|", "\uE000")
        if line.startswith("|"):
            line = line[1:]
        if line.endswith("|"):
            line = line[:-1]
        cells = [cell.strip() for cell in line.split("|")]
        rows.append(cells)
    if len(rows) >= 2 and all(
        re.fullmatch(r":?-{2,}:?", cell.replace(" ", "")) for cell in rows[1]
    ):
        rows = [rows[0]] + rows[2:]
    if not rows:
        return ""
    ncols = max(len(row) for row in rows)
    width = 0.94 / max(1, ncols)
    spec = "@{}" + "".join(
        r">{\raggedright\arraybackslash}p{%.2f\linewidth}" % width
        for _ in range(ncols)
    ) + "@{}"
    out = [r"\begin{longtable}{%s}" % spec, r"\toprule"]
    header = rows[0] + [""] * (ncols - len(rows[0]))
    out.append(" & ".join(inline(cell) for cell in header) + r" \\")
    out.extend([r"\midrule", r"\endhead"])
    for row in rows[1:]:
        row = row + [""] * (ncols - len(row))
        out.append(" & ".join(inline(cell) for cell in row) + r" \\")
    out.extend([r"\bottomrule", r"\end{longtable}", ""])
    return "\n".join(out)


def convert(markdown: str, level_offset: int = 0) -> str:
    lines = markdown.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if stripped.startswith("```"):
            i += 1
            code_lines = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code_lines.append(lines[i])
                i += 1
            i += 1
            out.extend([
                r"\begin{lstlisting}",
                "\n".join(code_lines),
                r"\end{lstlisting}",
                "",
            ])
            continue

        if stripped.startswith("|"):
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i])
                i += 1
            out.append(table(block))
            continue

        heading = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if heading:
            depth = len(heading.group(1)) + level_offset
            command = {1: "section", 2: "subsection",
                       3: "subsubsection"}.get(depth, "paragraph")
            title = inline(heading.group(2))
            if command == "section":
                out.extend([r"\section{%s}" % title, ""])
            elif command == "subsection":
                out.extend([r"\subsection{%s}" % title, ""])
            elif command == "subsubsection":
                out.extend([r"\subsubsection{%s}" % title, ""])
            else:
                out.extend([r"\paragraph{%s}" % title, ""])
            i += 1
            continue

        if stripped.startswith("- "):
            items = []
            while i < len(lines):
                current = lines[i].strip()
                if current.startswith("- "):
                    items.append(current[2:].strip())
                elif current and items and not re.match(
                        r"^(#|\||```|\d+\.\s)", current):
                    items[-1] += " " + current
                else:
                    break
                i += 1
            out.append(r"\begin{itemize}")
            for item in items:
                out.append(r"  \item %s" % inline(item))
            out.extend([r"\end{itemize}", ""])
            continue

        numbered = re.match(r"^\d+\.\s+(.*)$", stripped)
        if numbered:
            items = []
            while i < len(lines):
                current = lines[i].strip()
                match = re.match(r"^\d+\.\s+(.*)$", current)
                if match:
                    items.append(match.group(1))
                elif current and items and not re.match(
                        r"^(#|\||```|-\s)", current):
                    items[-1] += " " + current
                else:
                    break
                i += 1
            out.append(r"\begin{enumerate}")
            for item in items:
                out.append(r"  \item %s" % inline(item))
            out.extend([r"\end{enumerate}", ""])
            continue

        if not stripped:
            i += 1
            continue

        # Paragraph: consume until a blank line or another block.
        para = [stripped]
        i += 1
        while i < len(lines):
            current = lines[i].strip()
            if not current or current.startswith(("#", "|", "```", "- ")):
                break
            if re.match(r"^\d+\.\s", current):
                break
            para.append(current)
            i += 1
        out.extend([inline(" ".join(para)), ""])
    return "\n".join(out)


def main():
    manifest = json.loads((SKILLS / "manifest.json").read_text(encoding="utf-8"))
    entries = {entry["id"]: entry for entry in manifest["skills"]}
    lines = [
        r"\documentclass[11pt]{report}",
        r"\input{shared/preamble}",
        r"\usepackage{lmodern}",
        r"\title{YPhylo Workbench Tool Reference}",
        r"\author{Zhi-Jie Xu}",
        r"\date{Version 1.1.0 --- \today}",
        r"\begin{document}",
        r"\maketitle",
        r"\setcounter{tocdepth}{1}",
        r"\tableofcontents",
        r"\chapter*{Preface}",
        r"\addcontentsline{toc}{chapter}{Preface}",
        "This reference is generated from the agent-facing skill files in "
        r"\texttt{YR\_MPE/skills}. Each tool section states what the tool does, "
        "how the GUI controls map to command-line options or NEXUS directives, "
        "the recommended workflow, and the pitfalls that are easy to hit. "
        "Tools that are implemented entirely in Python or in the GUI say so "
        "explicitly instead of inventing a command line.",
        "",
        "The skill files remain the source of truth. Regenerate this document "
        "with \\texttt{python docs/tools/build\\_tool\\_reference.py} from the "
        "repository root, then build it with "
        "\\texttt{latexmk -lualatex manual/tool-reference.tex} from "
        "\\texttt{docs/}.",
        "",
    ]

    for chapter, ids in GROUPS:
        lines.append(r"\chapter{%s}" % latex_escape(chapter))
        for skill_id in ids:
            entry = entries.get(skill_id)
            if entry is None:
                continue
            path = SKILLS / entry["file"]
            if not path.is_file():
                continue
            text = path.read_text(encoding="utf-8")
            title = entry.get("title", skill_id)
            lines.append(r"\section{%s}" % inline(title))
            menu = entry.get("menu_path") or []
            if menu:
                menu_text = r" $\rightarrow$ ".join(inline(x) for x in menu)
                lines.append(r"\noindent\textit{Menu:} %s\\" % menu_text)
            lines.append(r"\noindent\textit{Skill id:} \texttt{%s}" %
                         latex_escape(skill_id))
            lines.append("")
            body = re.sub(r"^#\s+.*\n+", "", text, count=1)
            lines.append(convert(body, level_offset=0))
            lines.append("")

    lines.extend([
        r"\appendix",
        r"\chapter{Tool and menu index}",
        r"\begin{longtable}{@{}p{0.28\linewidth}p{0.48\linewidth}p{0.18\linewidth}@{}}",
        r"\toprule",
        r"Tool & Menu & Skill id \\",
        r"\midrule",
        r"\endhead",
    ])
    for entry in manifest["skills"]:
        menu = " $\\rightarrow$ ".join(inline(x) for x in entry.get("menu_path", []))
        lines.append("%s & %s & \\texttt{%s} \\\\" % (
            inline(entry.get("title", entry["id"])), menu,
            latex_escape(entry["id"]),
        ))
    lines.extend([
        r"\bottomrule",
        r"\end{longtable}",
        r"\end{document}",
        "",
    ])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
