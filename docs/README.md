# YR-MPE documentation (LaTeX)

**Status: Phase 0 skeleton (2026-10-30).** The manual and case chapters contain
headings and `\todo{}` markers; content will be imported from the wiki and the
read-only Word archive in later phases.

## Build

From this `docs/` directory:

```bash
latexmk -lualatex manual/manual.tex
latexmk -lualatex cases/cases.tex

# Human-readable tool reference generated from YR_MPE/skills/*.md
python docs/tools/build_tool_reference.py
latexmk -lualatex -outdir=manual manual/tool-reference.tex
```

Outputs: `manual.pdf`, `cases.pdf`, and `manual/tool-reference.pdf`.
`latexmkrc` sets `lualatex` and `biber`. The tool reference is regenerated
from the skill Markdown; edit the skill file, not the generated `.tex`.
Build artifacts are ignored by `docs/.gitignore`.

## Structure

```text
docs/
├── shared/
│   ├── preamble.tex      # packages, macros, tip/warning boxes
│   └── software.bib      # @software citations for wrapped tools
├── manual/
│   ├── manual.tex        # user manual skeleton
│   └── tool-reference.tex  # generated from YR_MPE/skills/*.md
├── tools/
│   └── build_tool_reference.py  # skill Markdown -> tool-reference.tex
├── cases/
│   └── cases.tex         # 9 case-study chapters
├── plugin-development-standard.md  # developer contract for new/refactored plugins
├── tools_inventory.md    # verified local tool versions/paths
└── latexmkrc
```

## Plugin development

`plugin-development-standard.md` is the maintained contract for new plugins
and the old-plugin refactor: metadata, `BasePlugin`/`BaseProcessThread`,
hidden subprocesses, citations, `ResultSpec` export, agent actions, host
registration, testing, and the migration mapping from legacy signals.

## Content sources

- `YR-MPE.wiki/Usage.md` — per-tool manual text (English).
- `YR-MPE.wiki/User-Guide.md` — installation, interface, basic workflow.
- `YR-MPE.wiki/workflow/*.md` — existing tutorials (single-gene, multi-gene,
  MrBayes, PhyloBayes, ASTER).
- `../user manual/User Manual-final.docx` — clean manual text.
- `../user manual/revision-judgments.md` — factual corrections; apply the
  accepted (✅/⚠️) items, reject the fabricated (❌) ones.

## Rules

- English-only prose, PhyloSuite appendix style.
- Quote exact UI strings from `YR_MPE/platforms/mpea.py` and the plugin code.
- One `@software` citation per wrapped tool; verify version/DOI.
- Do not use `\checkmark` (no glyph in Latin Modern); use `\cmark`/`\xmark`.
- Do not load `amssymb` with `unicode-math`.
- No Docker instructions; YR-MPE is distributed as a local YRTools plugin.
