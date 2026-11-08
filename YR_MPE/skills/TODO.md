# Skill coverage TODO

Status date: 2026-11-05.

The catalog has 30 entries but only 8 long-form documents. Twenty-two
manifest entries have no `file`, and most plugins do not yet implement
`control_bindings()`, so `workbench.help_doc` cannot show an accurate
GUI-to-CLI table. This file tracks the repair in small batches.

Contract for every completed skill:

- `YR_MPE/skills/<id>.md` exists and is referenced by `manifest.json`.
- `control_bindings()` lists every CLI-capable control with the exact flag.
- `widget_found` is true for every binding (dotted paths allowed).
- GUI-only controls use `"cli": []` with a note.
- The options are checked against the bundled executable, not memory.
- `dev` tests still pass.

## Phase 0 — infrastructure

- [x] Support dotted widget paths in `runtime/skills.py`
      (`advanced_dialog.round_spin`, `control_panel.methods_combo`).
- [x] Add a manifest audit script:
      every entry has a valid `file`; every plugin in `settings.json` maps to
      a skill entry. Run `python dev/skills/audit_skills.py`.
- [~] Add a help-doc audit test that instantiates each plugin offscreen and
      asserts every binding resolves to a real widget. Covered for ASTRAL,
      CASTER-site, Model Parameter and RF; extend as skills land.

## Phase 1 — species tree / coalescent

- [x] `astral.md` + bindings. Cover ASTRAL-IV, wASTRAL and ASTRAL-Hybrid;
      document the composite-support normalization for weighted modes.
- [x] `caster-site.md` + bindings.
- [x] `mrbayes.md` + bindings. The non-MPI and MPI-BEAGLE entries share it.
- [x] `phylobayes.md` + bindings.
- [x] `mpboot.md` + bindings.
- [x] `decenttree.md` + bindings.
- [x] `treerecs.md` + bindings.

## Phase 2 — alignment and trimming

- [x] `macse.md` + bindings (NT/AA output, frameshift and stop options).
- [x] `trimal.md` + bindings.
- [x] `gblocks.md` + bindings.
- [x] `muscle5.md` + bindings.
- [x] `clustal-omega.md` + bindings.

## Phase 3 — model, dating and comparison

- [x] `model-parameters.md` + bindings.
- [x] `lsd2.md` + bindings. Controls live in `control_panel`; requires
      dotted paths.
- [x] `rf.md` + bindings.
- [ ] Retrofit `orthofinder.md` with `control_bindings()`.
- [ ] Retrofit `minitracer.md` with `control_bindings()` where a CLI exists;
      mark the rest GUI-only.

## Phase 4 — datasets and selection

- [x] `seqmatrix.md` (GUI/format workflow; mostly GUI-only).
- [x] `seqdbg.md` (annotation reconciliation; GUI-only).
- [x] `sequence-concatenator.md` (GUI-only + partition output contract).
- [x] `tree-merger.md` (GUI-only + Newick output contract).
- [x] `codon-ztest.md` (GUI-only; no external CLI).
- [x] `batch-dnds.md` (IQ-TREE MLE command documented).
- [x] `icytree.md` (viewer; GUI-only, document accepted formats).

## Phase 5 — tools missing from the manifest

- [x] `ml-distance` — `distance calculation (ML)`.
- [x] `uncorrected-distance` — `Uncorrected p-distance`.
- [x] `saturation` — Substitution Saturation (BaCoCa C-value).
- [x] `pdguide` — molecular-dating calibration setup.
- [x] `sequence-viewer` — viewer/editor; GUI-only.
- [x] `console` — bundled-tool console; GUI-only.

## Explicitly out of scope

- `Example Plugin` (`hello_plugin.py`) — template, not a user tool.
- `About YR-MPE`, `YR-MPEA`, `YR-MPEA [Multigene devel]` — application
  shells, not analysis plugins.

## Batch log

- 2026-11-05: audit is clean: 36 catalog entries cover all 35 configured
  analysis plugins, every entry has a document. Completed all Phase 1–5
  docs. Remaining: add/retrofit control_bindings for OrthoFinder and
  MiniTracer, and add CLI bindings for GUI-only dataset/selection tools where
  a real CLI exists.
