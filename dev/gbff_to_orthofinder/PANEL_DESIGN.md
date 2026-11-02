# OrthoFinder plugin panel — design

## Goal

A three-step panel for **GBFF → CDS/AA → OrthoFinder → single-copy OG → CDS**.
It should make the pipeline obvious, validate inputs before running, and never
block the platform during a long OrthoFinder run.

The panel follows the existing `BasePlugin` four-tab contract:

1. **Input & Parameters** — step wizard (Input → OrthoFinder → Review)
2. **Output Preview** — summary cards + single-copy OG table
3. **Console** — streamed logs (base)
4. **Citation** — OrthoFinder / BLAST+ / MCL / Biopython

## Layout

```text
┌────────────────────────────────────────────────────────────────────┐
│ Genomic Dataset Formation (OrthoFinder)                            │
│ GBFF → CDS/AA → OrthoFinder → single-copy orthogroups → CDS        │
├────────────────────────────────────────────────────────────────────┤
│  ● 1 Input          ○ 2 OrthoFinder          ○ 3 Review            │
├────────────────────────────────────────────────────────────────────┤
│ Input source                                                       │
│   ( ) GenBank (GBFF/GB/GBK)     ( ) Proteome FASTA                 │
│   Folder: [____________________________________] [Browse]           │
│   4 files · 4 species · CDS counted on run                          │
│                                                                    │
│ CDS parsing                                                        │
│   Genetic code: [Auto (from /transl_table) ▾]                      │
│   [x] Skip pseudogenes     Internal stop: [Skip ▾]                 │
│   [ ] Keep terminal stop codon                                     │
│                                                                    │
│ Output                                                             │
│   Folder: [____________________________________] [Browse]           │
│                                                                    │
│                                      [Next: OrthoFinder →]         │
├────────────────────────────────────────────────────────────────────┤
│ OrthoFinder runtime                                                │
│   Executable: [______________________________] [Browse]             │
│   BLAST+ dir: [______________________________] [Browse]             │
│   Threads: [2 ▲▼]     [Check runtime]  ✓ OrthoFinder, blastp, ...  │
│                                                                    │
│ Run mode                                                           │
│   (•) Parse + OrthoFinder + single-copy OGs                        │
│   ( ) Parse only                                                   │
│   ( ) Use existing OrthoFinder results                              │
│                                                                    │
│                         [← Back]      [Next: Review →]             │
├────────────────────────────────────────────────────────────────────┤
│ Review                                                             │
│   Input: 4 GBFF files · 4 species                                  │
│   Output: .../orthofinder_20261102_205500                          │
│   Runtime: OrthoFinder 2.5.5 · BLAST+ 2.12 · threads 2             │
│   [Run pipeline]  (also available as the bottom Run button)         │
└────────────────────────────────────────────────────────────────────┘
```

## Output tab

```text
┌────────────────────────────────────────────────────────────────────┐
│ Pipeline complete                                                   │
│ ┌──────────┐ ┌──────────┐ ┌──────────────┐ ┌───────────────────┐   │
│ │ Species  │ │ CDS      │ │ Orthogroups  │ │ Single-copy OGs   │   │
│ │    4     │ │   335    │ │     80       │ │       62          │   │
│ └──────────┘ └──────────┘ └──────────────┘ └───────────────────┘   │
│ Single-copy orthogroups                                            │
│ ┌────────────┬──────────┬───────────┬────────┐                     │
│ │ OG         │ Species  │ Proteins  │ Status │                     │
│ │ OG0000008  │    4     │     4     │ ok     │                     │
│ │ ...        │          │           │        │                     │
│ └────────────┴──────────┴───────────┴────────┘                     │
│ Output: .../single_copy_ogs                                         │
│ [Open output folder]  [Export summary CSV]                         │
└────────────────────────────────────────────────────────────────────┘
```

## Validation

- **Input folder** must exist and contain at least one `.gb/.gbff/.gbk`
  (GBFF mode) or `.fa/.faa/.fasta` (proteome mode).
- **Output folder** must be writable; default to
  `<workdir>/orthofinder_<timestamp>` or the system temp dir.
- **Genetic code** Auto uses `/transl_table` when present; fallback to the
  selected table.
- **Full/Existing mode** requires OrthoFinder executable and, for full mode,
  `blastp.exe`, `makeblastdb.exe`, `mcl.exe` in the BLAST/runtime dir.
- **Parse only** does not require OrthoFinder/BLAST.
- **Existing results** requires a folder containing `Orthogroups.tsv`.
- The Review page shows a compact summary; Run is disabled until page 3.

## States

- **Empty**: no folder selected → helper text + disabled Next.
- **Scanning**: "Scanning…" label; do not block the UI for more than 400 ms;
  for large folders, count in a worker or show indeterminate progress.
- **Valid**: green check + text "Ready"; not color alone.
- **Warning**: amber warning icon + text ("3 files, 2 species; check naming").
- **Error**: red error icon + text; inline message under the field.
- **Running**: indeterminate progress bar + Console tab; Stop enabled.
- **Done**: summary cards + table; "Open output folder".
- **No single-copy OGs**: empty state "No single-copy orthogroups found;
  check species overlap or OrthoFinder log."

## Visual rules (from qt-ui-design)

- 8 px spacing grid; section gaps 16 px; body 10 pt, section heading 13 pt,
  title 18 pt; system font (Segoe UI on Windows).
- Use layout managers; no absolute positions.
- Status uses **icon + text + color**; never color alone.
- Primary action per page: one button; Back/Next secondary.
- Keyboard: Tab order follows visual order; Enter = Next; Esc = close;
  Ctrl+O = browse input; F5 = Run.
- No QWebEngine in the output tab; native cards/table are lighter.
- Long tasks are non-modal; Console streams logs; Stop uses
  `BaseProcessThread.stop()`.

## Implementation notes

- `gbff.py` and `orthogroups.py` live in
  `YR_MPE/plugins/components/methods/`.
- `OrthoFinderThread(BaseProcessThread)` parses, runs OrthoFinder with
  `-S blast -og -X`, and post-processes; it emits `result_summary(dict)`.
- The plugin overrides `setup_input_tab`, `setup_output_tab`, `run_analysis`,
  `stop_analysis`, `analysis_finished`, `analysis_error`.
- Runtime default: `self.tool_path` from `config.json` (`OrthoFinder`);
  BLAST dir defaults to the same directory; user can override.
- Add the menu action under WORKSPACE → Dataset:
  "Genomic Dataset Formation (OrthoFinder)".
- Placeholder icon: `YR_MPE/icons/file/dataset.svg` (user will replace the
  final vector art).
