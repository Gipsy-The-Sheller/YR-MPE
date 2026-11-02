"""Software citations as data, not code.

Citations live in ``YR_MPE/citations/*.bib`` (or ``*.ris``). ``config.json``
maps each tool to one or more citation keys; a plugin never stores the
formatted citation itself. The Citation tab parses only the keys it needs,
at display time, and re-parses automatically when a citation file changes.

The module is stdlib-only so it can be used by headless tests and CLI tools.
"""

from __future__ import annotations

import glob
import html
import os
import re
import threading
import unicodedata
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# LaTeX-lite / BibTeX helpers
# ---------------------------------------------------------------------------

_ENTRY_RE = re.compile(r"@([A-Za-z][A-Za-z0-9_-]*)\s*([{(])")
_FIELD_NAME_RE = re.compile(r"[A-Za-z][A-Za-z0-9_:-]*")

_ACCENTS = {
    "'": "\u0301",  # acute
    "`": "\u0300",  # grave
    "^": "\u0302",  # circumflex
    '"': "\u0308",  # diaeresis
    "~": "\u0303",  # tilde
    "=": "\u0304",  # macron
    ".": "\u0307",  # dot above
    "c": "\u0327",  # cedilla
    "v": "\u030c",  # caron
    "u": "\u0306",  # breve
    "H": "\u030b",  # double acute
    "r": "\u030a",  # ring
}

_TEXT_COMMAND_RE = re.compile(
    r"\\(?:textit|emph|textbf|textrm|textsc|mathrm|mathbf|text)\s*\{([^{}]*)\}"
)
_ACCENT_RE = re.compile(r"\\(['`^\"~=.]|[cvuHr])\s*\{?([A-Za-z])\}?")
_SPECIAL_LETTERS = {
    "ss": "\u00df",
    "ae": "\u00e6",
    "AE": "\u00c6",
    "oe": "\u0153",
    "OE": "\u0152",
    "o": "\u00f8",
    "O": "\u00d8",
    "aa": "\u00e5",
    "AA": "\u00c5",
    "l": "\u0142",
    "L": "\u0141",
}


def _accent_repl(match: "re.Match[str]") -> str:
    accent, char = match.group(1), match.group(2)
    combining = _ACCENTS.get(accent)
    if not combining:
        return char
    return unicodedata.normalize("NFC", char + combining)


def clean_latex(value: str) -> str:
    """Turn a small, safe subset of LaTeX into plain Unicode text.

    This is deliberately not a TeX engine: it handles accents, escaped
    characters, text commands, dashes and grouping braces, which covers
    normal bibliography metadata.
    """
    if not value:
        return ""

    text = str(value)
    text = _TEXT_COMMAND_RE.sub(r"\1", text)

    for name, replacement in _SPECIAL_LETTERS.items():
        text = re.sub(r"\\%s(?![A-Za-z])" % re.escape(name), replacement, text)

    text = re.sub(r"\\([&%$#_{}])", r"\1", text)
    text = _ACCENT_RE.sub(_accent_repl, text)
    text = text.replace("---", "\u2014").replace("--", "\u2013")
    text = text.replace("{", "").replace("}", "")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def _find_matching(text: str, start: int, open_ch: str, close_ch: str) -> int:
    """Return the index of the delimiter matching the one before ``start``."""
    depth = 1
    i = start
    while i < len(text):
        ch = text[i]
        if ch == "\\":
            i += 2
            continue
        if ch == open_ch:
            depth += 1
        elif ch == close_ch:
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


def _read_bib_value(text: str, index: int) -> Tuple[str, int]:
    pieces: List[str] = []
    i = index
    while True:
        while i < len(text) and text[i].isspace():
            i += 1
        if i >= len(text):
            break

        ch = text[i]
        if ch == "{":
            end = _find_matching(text, i + 1, "{", "}")
            if end < 0:
                pieces.append(text[i + 1 :])
                i = len(text)
                break
            pieces.append(text[i + 1 : end])
            i = end + 1
        elif ch == '"':
            j = i + 1
            depth = 0
            buf: List[str] = []
            while j < len(text):
                c = text[j]
                if c == "\\" and j + 1 < len(text):
                    buf.append(text[j : j + 2])
                    j += 2
                    continue
                if c == "{":
                    depth += 1
                elif c == "}":
                    depth = max(0, depth - 1)
                elif c == '"' and depth == 0:
                    break
                buf.append(c)
                j += 1
            pieces.append("".join(buf))
            i = j + 1
        else:
            j = i
            while j < len(text) and text[j] not in ",#\n\r\t ":
                j += 1
            pieces.append(text[i:j])
            i = j

        while i < len(text) and text[i].isspace():
            i += 1
        if i < len(text) and text[i] == "#":
            i += 1
            continue
        break

    return "".join(pieces).strip(), i


def _parse_entry_body(body: str) -> Tuple[str, Dict[str, str]]:
    """Split ``key, field = {value}, ...`` into a key and a field dict."""
    i = 0
    n = len(body)
    depth = 0
    while i < n:
        ch = body[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch == "," and depth == 0:
            break
        i += 1

    key = body[:i].strip()
    i += 1
    fields: Dict[str, str] = {}

    while i < n:
        while i < n and body[i] in " \t\r\n,":
            i += 1
        if i >= n:
            break

        match = _FIELD_NAME_RE.match(body, i)
        if not match:
            # Unknown token: skip to the next top-level comma.
            depth = 0
            while i < n:
                if body[i] == "{":
                    depth += 1
                elif body[i] == "}":
                    depth = max(0, depth - 1)
                elif body[i] == "," and depth == 0:
                    break
                i += 1
            continue

        name = match.group(0).lower()
        i = match.end()
        while i < n and body[i].isspace():
            i += 1
        if i >= n or body[i] != "=":
            continue
        i += 1
        value, i = _read_bib_value(body, i)
        fields[name] = value

    return key, fields


def parse_bibtex(text: str, source: str = "") -> Tuple[List["Citation"], List[str]]:
    """Parse a BibTeX string. Returns ``(entries, errors)``."""
    entries: List["Citation"] = []
    errors: List[str] = []
    pos = 0

    while True:
        match = _ENTRY_RE.search(text, pos)
        if not match:
            break

        entry_type = match.group(1).lower()
        open_ch = match.group(2)
        close_ch = "}" if open_ch == "{" else ")"
        end = _find_matching(text, match.end(), open_ch, close_ch)
        if end < 0:
            errors.append("%s: unterminated @%s entry" % (source, entry_type))
            break

        body = text[match.end() : end]
        pos = end + 1

        if entry_type in ("comment", "string", "preamble"):
            continue

        key, fields = _parse_entry_body(body)
        if not key:
            errors.append("%s: @%s entry without a key" % (source, entry_type))
            continue
        entries.append(
            Citation(key=key, entry_type=entry_type, fields=fields, source=source)
        )

    return entries, errors


# ---------------------------------------------------------------------------
# RIS
# ---------------------------------------------------------------------------

_RIS_TYPE_MAP = {
    "JOUR": "article",
    "JFULL": "article",
    "MGZN": "article",
    "BOOK": "book",
    "CHAP": "incollection",
    "CONF": "inproceedings",
    "CPAPER": "inproceedings",
    "THES": "phdthesis",
    "RPRT": "techreport",
    "GEN": "misc",
    "COMP": "software",
    "ELEC": "online",
}

_RIS_TAG_RE = re.compile(r"^([A-Z][A-Z0-9])  - ?(.*)$")


def _first(record: Dict[str, List[str]], *tags: str) -> str:
    for tag in tags:
        values = record.get(tag)
        if values:
            return values[0]
    return ""


def parse_ris(text: str, source: str = "") -> Tuple[List["Citation"], List[str]]:
    """Parse RIS records into citations."""
    records: List[Dict[str, List[str]]] = []
    current: Dict[str, List[str]] = {}
    last_tag: Optional[str] = None
    errors: List[str] = []

    for line in text.splitlines():
        if line.strip() == "ER  -":
            if current:
                records.append(current)
            current = {}
            last_tag = None
            continue

        match = _RIS_TAG_RE.match(line)
        if match:
            tag, value = match.group(1), match.group(2).strip()
            current.setdefault(tag, []).append(value)
            last_tag = tag
        elif line.startswith(" ") and last_tag and line.strip():
            # Continuation of the previous tag.
            current[last_tag][-1] = (current[last_tag][-1] + " " + line.strip()).strip()

    if current:
        records.append(current)

    entries: List["Citation"] = []
    for index, record in enumerate(records):
        ris_type = _first(record, "TY").upper()
        entry_type = _RIS_TYPE_MAP.get(ris_type, "misc")
        key = _first(record, "ID", "N1") or "ris_%d" % (index + 1)

        fields: Dict[str, str] = {}
        authors = record.get("AU") or record.get("A1") or []
        if authors:
            fields["author"] = " and ".join(authors)
        editors = record.get("ED") or record.get("A2") or []
        if editors:
            fields["editor"] = " and ".join(editors)

        title = _first(record, "TI", "T1", "CT")
        if title:
            fields["title"] = title

        journal = _first(record, "JO", "JF", "T2", "JA", "BT")
        if journal:
            fields["journal" if entry_type == "article" else "booktitle"] = journal

        year = _first(record, "PY", "Y1", "DA")
        if year:
            fields["year"] = year

        volume = _first(record, "VL")
        if volume:
            fields["volume"] = volume
        issue = _first(record, "IS")
        if issue:
            fields["number"] = issue

        start = _first(record, "SP")
        end = _first(record, "EP")
        if start:
            fields["pages"] = start + ("-" + end if end else "")

        doi = _first(record, "DO")
        if doi:
            fields["doi"] = doi
        url = _first(record, "UR")
        if url:
            fields["url"] = url

        publisher = _first(record, "PB")
        if publisher:
            fields["publisher"] = publisher
        address = _first(record, "CY")
        if address:
            fields["address"] = address
        school = _first(record, "PB") if entry_type.endswith("thesis") else ""
        if school:
            fields["school"] = school
        version = _first(record, "VN")
        if version:
            fields["version"] = version

        if not fields.get("title"):
            errors.append("%s: RIS record %s has no title" % (source, key))
            continue

        entries.append(
            Citation(key=key, entry_type=entry_type, fields=fields, source=source)
        )

    return entries, errors


# ---------------------------------------------------------------------------
# Author / rendering helpers
# ---------------------------------------------------------------------------

def _split_authors(raw: str) -> List[str]:
    if not raw:
        return []

    parts: List[str] = []
    buf: List[str] = []
    depth = 0
    i = 0
    while i < len(raw):
        ch = raw[i]
        if ch == "{":
            depth += 1
            buf.append(ch)
            i += 1
        elif ch == "}":
            depth = max(0, depth - 1)
            buf.append(ch)
            i += 1
        elif (
            depth == 0
            and raw[i : i + 3].lower() == "and"
            and (i == 0 or raw[i - 1].isspace())
            and (i + 3 == len(raw) or raw[i + 3].isspace())
        ):
            parts.append("".join(buf).strip())
            buf = []
            i += 3
        else:
            buf.append(ch)
            i += 1

    if buf:
        parts.append("".join(buf).strip())
    return [part for part in parts if part]


def _initials(given: str) -> str:
    initials: List[str] = []
    for token in given.split():
        if not token:
            continue
        if "-" in token:
            pieces = [piece for piece in token.split("-") if piece]
            initials.append("-".join(piece[0].upper() + "." for piece in pieces))
        else:
            initials.append(token[0].upper() + ".")
    return " ".join(initials)


def _format_author(name: str) -> str:
    name = clean_latex(name).strip()
    if not name:
        return ""

    if "," in name:
        parts = [part.strip() for part in name.split(",")]
        last = parts[0]
        first = parts[-1] if len(parts) >= 2 else ""
    else:
        tokens = name.split()
        if len(tokens) == 1:
            return tokens[0]
        last, first = tokens[-1], " ".join(tokens[:-1])

    initials = _initials(first)
    return "%s, %s" % (last, initials) if initials else last


def format_authors(raw: str) -> str:
    names = [_format_author(name) for name in _split_authors(raw)]
    names = [name for name in names if name]
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + " & " + names[-1]


def _sentence(text: str) -> str:
    text = text.strip()
    if not text:
        return ""
    if text[-1] in ".?!":
        return text
    return text + "."


# ---------------------------------------------------------------------------
# Citation model
# ---------------------------------------------------------------------------

@dataclass
class Citation:
    key: str
    entry_type: str
    fields: Dict[str, str] = field(default_factory=dict)
    source: str = ""

    def get(self, *names: str) -> str:
        for name in names:
            value = self.fields.get(name)
            if value:
                return value
        return ""

    @property
    def year(self) -> str:
        match = re.search(r"\d{4}", self.get("year", "date"))
        return match.group(0) if match else ""

    @property
    def authors(self) -> str:
        return format_authors(self.get("author", "editor"))

    @property
    def doi(self) -> str:
        return self.get("doi").strip()

    @property
    def link(self) -> str:
        doi = self.doi
        if doi:
            doi = re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", doi, flags=re.I)
            return "https://doi.org/" + doi
        return self.get("url").strip()

    def _venue(self) -> str:
        return clean_latex(self.get("journal", "booktitle"))

    def _volume(self) -> str:
        return clean_latex(self.get("volume"))

    def _issue(self) -> str:
        return clean_latex(self.get("number", "issue"))

    def _pages(self) -> str:
        pages = clean_latex(self.get("pages"))
        return pages.replace("-", "\u2013") if re.fullmatch(r"\d+-\d+", pages) else pages

    def _plain_title(self) -> str:
        return clean_latex(self.get("title"))

    def format_text(self) -> str:
        """Render a human-readable, copyable citation without HTML."""
        authors = self.authors
        year = self.year
        title = self._plain_title()
        venue = self._venue()
        volume = self._volume()
        issue = self._issue()
        pages = self._pages()
        publisher = clean_latex(self.get("publisher"))
        school = clean_latex(self.get("school"))
        version = clean_latex(self.get("version"))
        howpublished = clean_latex(self.get("howpublished"))

        head = ""
        if authors:
            head += authors
        if year:
            head += (" " if head else "") + "(%s)." % year

        entry_type = self.entry_type.lower()

        if entry_type == "book":
            body = _sentence(title)
            if publisher:
                body += " " + _sentence(publisher)
            text = _join(head, body)
        elif entry_type in ("phdthesis", "mastersthesis"):
            kind = "PhD thesis" if entry_type == "phdthesis" else "Master's thesis"
            body = _sentence(title)
            tail = ", ".join(part for part in (kind, school) if part)
            if tail:
                body += " " + _sentence(tail)
            text = _join(head, body)
        elif entry_type in ("techreport",):
            body = _sentence(title)
            if publisher or school:
                body += " " + _sentence(publisher or school)
            text = _join(head, body)
        else:
            body = _sentence(title)
            if venue:
                body += " " + venue
                if volume:
                    body += " " + volume
                    if issue:
                        body += "(%s)" % issue
                if pages:
                    body += ": " + pages
                body = _sentence(body)
            elif publisher:
                body += " " + _sentence(publisher)
            if version:
                body += " Version %s." % version
            if howpublished:
                body += " " + _sentence(howpublished)
            text = _join(head, body)

        if self.link:
            text += " " + self.link
        return re.sub(r"\s+", " ", text).strip()

    def format_html(self) -> str:
        """Render the citation for the Qt Citation tab (HTML subset)."""
        authors = html.escape(self.authors, quote=False)
        year = html.escape(self.year, quote=False)
        title = html.escape(self._plain_title(), quote=False)
        venue = html.escape(self._venue(), quote=False)
        volume = html.escape(self._volume(), quote=False)
        issue = html.escape(self._issue(), quote=False)
        pages = html.escape(self._pages(), quote=False)
        publisher = html.escape(clean_latex(self.get("publisher")), quote=False)
        school = html.escape(clean_latex(self.get("school")), quote=False)
        version = html.escape(clean_latex(self.get("version")), quote=False)
        howpublished = html.escape(clean_latex(self.get("howpublished")), quote=False)

        head = ""
        if authors:
            head += authors
        if year:
            head += (" " if head else "") + "(%s)." % year

        entry_type = self.entry_type.lower()

        if entry_type == "book":
            body = _sentence("<i>%s</i>" % title)
            if publisher:
                body += " " + _sentence(publisher)
            text = _join(head, body)
        elif entry_type in ("phdthesis", "mastersthesis"):
            kind = "PhD thesis" if entry_type == "phdthesis" else "Master's thesis"
            body = _sentence("<i>%s</i>" % title)
            tail = ", ".join(part for part in (kind, school) if part)
            if tail:
                body += " " + _sentence(tail)
            text = _join(head, body)
        elif entry_type in ("techreport",):
            body = _sentence("<i>%s</i>" % title)
            if publisher or school:
                body += " " + _sentence(publisher or school)
            text = _join(head, body)
        else:
            body = _sentence(title)
            if venue:
                body += " <i>%s</i>" % venue
                if volume:
                    body += " <b>%s</b>" % volume
                    if issue:
                        body += "(%s)" % issue
                if pages:
                    body += ": " + pages
                body = _sentence(body)
            elif publisher:
                body += " " + _sentence(publisher)
            if version:
                body += " Version %s." % version
            if howpublished:
                body += " " + _sentence(howpublished)
            text = _join(head, body)

        if self.link:
            link = html.escape(self.link, quote=True)
            text += ' <a href="%s">%s</a>' % (link, html.escape(self.link, quote=False))
        return re.sub(r"\s+", " ", text).strip()


def _join(head: str, body: str) -> str:
    if head and body:
        return head + " " + body
    return head or body


# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------

class CitationLibrary:
    """Lazily parsed citation files in one directory.

    ``reload_if_changed`` compares mtimes and sizes, so editing a ``.bib``
    file is reflected without restarting the application.
    """

    def __init__(self, directory: str):
        self.directory = os.path.abspath(directory)
        self.entries: Dict[str, Citation] = {}
        self.errors: List[str] = []
        self.duplicates: List[str] = []
        self._signature: Optional[Tuple[Tuple[str, int, int], ...]] = None
        self.reload_if_changed()

    def _files(self) -> List[str]:
        patterns = ("*.bib", "*.ris")
        files: List[str] = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(self.directory, pattern)))
        # Prefer BibTeX over RIS when the same key appears in both.
        return sorted(files, key=lambda path: (not path.lower().endswith(".bib"), path.lower()))

    def _current_signature(self) -> Tuple[Tuple[str, int, int], ...]:
        signature = []
        for path in self._files():
            try:
                stat = os.stat(path)
            except OSError:
                continue
            signature.append((path, stat.st_mtime_ns, stat.st_size))
        return tuple(signature)

    def reload_if_changed(self) -> bool:
        signature = self._current_signature()
        if signature == self._signature:
            return False

        self._signature = signature
        self.entries = {}
        self.errors = []
        self.duplicates = []

        for path in self._files():
            try:
                with open(path, "r", encoding="utf-8-sig") as handle:
                    text = handle.read()
            except OSError as exc:
                self.errors.append("%s: %s" % (path, exc))
                continue

            if path.lower().endswith(".ris"):
                entries, errors = parse_ris(text, source=path)
            else:
                entries, errors = parse_bibtex(text, source=path)
            self.errors.extend(errors)

            for entry in entries:
                if entry.key in self.entries:
                    self.duplicates.append(
                        "%s: duplicate key %r (kept first)" % (path, entry.key)
                    )
                    continue
                self.entries[entry.key] = entry

        return True

    def get(self, key: str) -> Optional[Citation]:
        self.reload_if_changed()
        return self.entries.get(key)

    def get_many(self, keys: Iterable[str]) -> Tuple[List[Citation], List[str]]:
        self.reload_if_changed()
        found: List[Citation] = []
        missing: List[str] = []
        seen = set()
        for key in keys:
            key = (key or "").strip()
            if not key or key in seen:
                continue
            seen.add(key)
            entry = self.entries.get(key)
            if entry is None:
                missing.append(key)
            else:
                found.append(entry)
        return found, missing

    def keys(self) -> List[str]:
        self.reload_if_changed()
        return sorted(self.entries)


_LIBRARY_CACHE: Dict[str, CitationLibrary] = {}
_LIBRARY_LOCK = threading.Lock()


def get_citation_library(directory: str) -> CitationLibrary:
    """Return a shared, mtime-aware library for ``directory``."""
    path = os.path.abspath(directory)
    with _LIBRARY_LOCK:
        library = _LIBRARY_CACHE.get(path)
        if library is None:
            library = CitationLibrary(path)
            _LIBRARY_CACHE[path] = library
    library.reload_if_changed()
    return library
