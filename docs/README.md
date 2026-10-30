# YR-MPE documentation (LaTeX)

**Status: Phase 0 skeleton (2026-10-30).** The manual and case chapters contain
headings and `\todo{}` markers; content will be imported from the wiki and the
read-only Word archive in later phases.

## Build

From this `docs/` directory:

```bash
latexmk -lualatex manual/manual.tex
latexmk -lualatex cases/cases.tex
```

Outputs: `manual.pdf`, `cases.pdf` in this `docs/` directory (latexmk writes
them next to the build aux files). `latexmkrc` sets `lualatex` and `biber`.
Build artifacts are ignored by `docs/.gitignore`.

## Structure

```text
docs/
├── shared/
│   ├── preamble.tex      # packages, macros, tip/warning boxes
│   └── software.bib      # @software citations for wrapped tools
├── manual/
│   └── manual.tex        # user manual skeleton
├── cases/
│   └── cases.tex         # 9 case-study chapters
├── tools_inventory.md    # verified local tool versions/paths
└── latexmkrc
```

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
