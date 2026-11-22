"""Render a YGBlocks masked-alignment report for the GBlocks Output tab.

YGBlocks 1.0.0 writes ``<input>-gb.txt``: the first line is
``<ntaxa> <nchar>``, followed by interleaved PHYLIP blocks of 60 columns.
Kept positions are uppercase, removed positions lowercase, and gaps are
written as ``.``.  Older builds wrote ``<input>-gb.html`` instead; both forms
are accepted and converted to the small HTML subset QTextBrowser supports.
"""

from __future__ import annotations

import html
import re
from pathlib import Path

MAX_COLUMNS = 20000

_PRE = re.compile(r"<pre>(.*?)</pre>", re.S | re.I)
_LINE = re.compile(r"^([^\s<]+)(\s+)(.*)$")


def parse_report(path):
    """Return ``(names, rows, nchar)``.

    ``rows`` is a list of ``[(char, kept)]`` entries.  A ``.txt`` report is
    parsed directly; a ``.html`` report is parsed from its ``<pre>`` block.
    """
    path = Path(path)
    if path.suffix.lower() in (".html", ".htm"):
        return _parse_html(path)
    return _parse_text(path)


def _parse_text(path):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise ValueError("empty report")
    header = lines[0].split()
    if len(header) < 2:
        raise ValueError("report header is not '<ntaxa> <nchar>'")
    nchar = int(header[1])

    blocks = []
    current = []
    for line in lines[1:]:
        if line.strip():
            current.append(line)
        elif current:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    if not blocks:
        raise ValueError("report contains no alignment blocks")

    names = []
    rows = []
    for index, block in enumerate(blocks):
        for row, line in enumerate(block):
            if index == 0:
                parts = line.strip().split(None, 1)
                if len(parts) != 2:
                    raise ValueError("cannot parse the first report block")
                names.append(parts[0])
                rows.append([])
                text = parts[1]
            else:
                if row >= len(rows):
                    raise ValueError("report block has more rows than the first")
                text = line
            for char in "".join(text.split()):
                rows[row].append((char, char.isupper()))

    lengths = {len(row) for row in rows}
    if len(lengths) != 1 or lengths.pop() != nchar:
        raise ValueError("report rows do not match the declared alignment length")
    return names, rows, nchar


def _parse_html(path):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    match = _PRE.search(text)
    if not match:
        raise ValueError("no <pre> alignment in the YGBlocks report")

    names = []
    rows = []
    for raw in match.group(1).splitlines():
        if "<span" not in raw and not re.search(r"[A-Za-z]\.?[A-Za-z]", raw):
            continue
        parsed = _LINE.match(raw)
        if not parsed:
            continue
        name, _spaces, fragment = parsed.groups()
        if not name.strip() or name.startswith("="):
            continue
        chars = _line_chars(fragment)
        if not chars:
            continue
        names.append(name)
        rows.append(chars)

    if not rows:
        raise ValueError("no sequence rows in the YGBlocks report")
    lengths = {len(row) for row in rows}
    if len(lengths) != 1:
        raise ValueError("YGBlocks report rows have different lengths")
    return names, rows, lengths.pop()


def _line_chars(fragment: str):
    chars = []
    selected = False
    index = 0
    while index < len(fragment):
        if fragment.startswith("<span", index):
            end = fragment.find(">", index)
            if end < 0:
                break
            selected = True
            index = end + 1
            continue
        if fragment.startswith("</span>", index):
            selected = False
            index += len("</span>")
            continue
        char = fragment[index]
        if not char.isspace():
            chars.append((html.unescape(char), selected))
        index += 1
    return chars


def render_report_html(path, title="YGBlocks masked alignment") -> str:
    """Return a small HTML document suitable for ``QTextBrowser.setHtml``."""
    names, rows, nchar = parse_report(path)
    shown = min(nchar, MAX_COLUMNS)
    width = max((len(name) for name in names), default=1) + 2
    kept = sum(1 for _char, selected in rows[0] if selected)
    kept_pct = (100.0 * kept / nchar) if nchar else 0.0

    lines = []
    for name, row in zip(names, rows):
        cells = []
        for char, selected in row[:shown]:
            escaped = html.escape(char)
            if selected:
                cells.append(
                    '<span style="color:#116329;background-color:#dafbe1">'
                    "%s</span>" % escaped)
            elif char == ".":
                cells.append(
                    '<span style="color:#a40e26;background-color:#ffebe9">'
                    "%s</span>" % escaped)
            else:
                cells.append(
                    '<span style="color:#8c959f;background-color:#f6f8fa">'
                    "%s</span>" % escaped)
        lines.append(
            '<div><span style="color:#57606a">%s</span>%s</div>'
            % (html.escape(name.ljust(width)), "".join(cells))
        )

    truncation = ""
    if shown < nchar:
        truncation = (
            '<p style="color:#9a6700">Preview limited to the first %d of %d '
            "columns.</p>" % (shown, nchar)
        )

    return (
        "<html><body>"
        '<h3 style="margin:0 0 2px 0">%s</h3>'
        '<p style="color:#57606a;margin:0 0 6px 0">%d taxa &middot; %d columns '
        "&middot; %d kept (%.1f%%)</p>"
        '<p style="margin:0 0 8px 0">'
        '<span style="color:#116329;background-color:#dafbe1">kept</span> '
        '<span style="color:#8c959f;background-color:#f6f8fa">removed</span> '
        '<span style="color:#a40e26;background-color:#ffebe9">gap</span></p>'
        "%s"
        '<pre style="font-family:Consolas,monospace;font-size:11px">%s</pre>'
        "</body></html>"
        % (html.escape(title), len(names), nchar, kept, kept_pct,
           truncation, "".join(lines))
    )
