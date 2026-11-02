# Citation contract (v1)

**Rule: citation text is data, not code.**

All software citations live in `YR_MPE/citations/*.bib` (or `*.ris`).
`YR_MPE/config.json` maps a tool to one or more citation keys. The plugin
stores neither the formatted string nor a copy of the metadata; it asks the
citation library to parse the keys it needs when the Citation tab is built.

This also means a binary update and its citation update are one edit:

```json
{
    "name": "OrthoFinder",
    "path": "/bin/orthofinder/orthofinder.exe",
    "citations": ["emms2019orthofinder", "camacho2009blast"]
}
```

Replace the binary path, add/remove a key, restart or reopen the panel. The
citation page follows the config immediately; no plugin code changes.

## Layout

```
YR_MPE/
  citations/
    software.bib        # canonical, BibTeX
    *.ris               # optional, parsed too
  runtime/
    citations.py        # stdlib-only parser + renderer
  config.json           # tool -> citation keys
```

BibTeX is preferred because it is easier to review and merge. RIS is accepted
because publisher sites and reference managers export it directly. If the same
key appears in both a `.bib` and a `.ris` file, the `.bib` entry wins and the
duplicate is reported by the checker.

## Config schema

Each tool entry may contain:

| Field | Type | Meaning |
|---|---|---|
| `citations` | list[str] | Ordered keys; a tool may need several papers |
| `citation` | str | Convenience for a single key; `citations` wins if both exist |

Unknown keys are not silently ignored: the Citation tab shows a warning and
the console records the missing key. `python dev/citations/check_citations.py`
fails the build/check when config, plugins and files disagree.

Plugins with no 1:1 binary (composite wizards, built-in Python methods) may add
`self.citation_keys = [...]` in `init_plugin_info`. The keys still point into
`citations/`; only the mapping lives in code. Config keys come first, then
plugin keys, deduplicated.

## Rendering

`Citation.format_html()` produces the Qt rich-text line:

> Authors (Year). Title. *Journal* **Volume**(Issue): pages.
> https://doi.org/...

`Citation.format_text()` produces the same citation without tags for the
clipboard. The parser handles UTF-8 and a small LaTeX subset (accents, escaped
characters, `--`/`---`, simple text commands); it is not a TeX engine.

## Refresh

`CitationLibrary` records `(path, mtime_ns, size)`. Every lookup calls
`reload_if_changed()`, so editing `software.bib` in a text editor is visible
the next time the Citation tab is opened or a plugin asks for a citation.
No restart is required.

## Adding a citation

1. Add a `@article`/`@book`/`@software` entry to `software.bib` with a stable
   key (`authorYYYYtopic`).
2. Add the key to the tool's `citations` array in `config.json`.
3. Run `python dev/citations/check_citations.py`.
4. Run `python dev/citations/test_citations.py`.

## Current coverage

- 28 entries, 16 configurable tools.
- Multiple citations per tool: IQ-TREE (2), PhyloBayes (4), MrBayes/BEAGLE
  (2), caster-site (2), MACSE (2), OrthoFinder (4).
- Built-in/composite: BaCoCa method (`kuck2014bacoca`), PD-Guide (LSD2 from
  config + IQ-TREE + MCMCTree keys).
