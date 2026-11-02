# GBFF -> OrthoFinder -> single-copy OG -> CDS (prototype)

Prototype for the YR-MPE "Genomic Dataset Formation (OrthoFinder)" feature.
It takes organelle GenBank files (GBFF/GB/GBK), extracts CDS features, writes
protein FASTA + in-frame CDS FASTA + a mapping table, runs OrthoFinder, and
maps single-copy orthogroups back to CDS.

## What it does

```text
GBFF/GBK/GB folder
  -> parse CDS (join/complement, codon_start, transl_table, /translation)
  -> proteomes/<species>.faa
  -> cds/<species>.fna
  -> cds_map.tsv
  -> OrthoFinder (-S blast -og -X)
  -> single_copy_ogs/<OG>.faa, <OG>.cds.fna, <OG>.map.tsv
```

## Usage

```bash
python cli.py \
  --input testdata/gb \
  --out testdata/work \
  --default-table 11 \
  --run-orthofinder \
  --orthofinder E:/YR-MPE_tutorial/YR-MPE_dev/dev/bin/orthofinder/orthofinder.exe \
  --blast-dir E:/YR-MPE_tutorial/YR-MPE_dev/dev/bin/orthofinder \
  --threads 2 \
  --postprocess
```

Parse only (no OrthoFinder):

```bash
python cli.py --input testdata/gb --out testdata/work --default-table 11
```

Tests:

```bash
python -m unittest discover -s tests -v
```

## Test result (4 chloroplast genomes)

Input: `testdata/gb/` — `HM347959.gb`, `JN637766.gb`, `KF150644.gb`,
`KJ806277.gb` (public GenBank records).

- Parsed: **335 CDS** (74 + 87 + 87 + 87), all CDS lengths multiples of 3.
- OrthoFinder 2.5.5 (BLAST+ 2.12.0, 2 threads):
  - 335 genes assigned to **80 orthogroups**
  - 68 orthogroups with all species present
  - **62 single-copy orthogroups**
- Back-mapping: **62 OG × 4 species = 248 proteins**; every OG has exactly one
  CDS per species.
- Translation check: 238/248 exact; **10/248 differ only at the first residue**
  (`/translation` starts M; CDS translates I/V/T because of alternative start
  codons). No real mismatches. This is expected for chloroplast genomes.

## Dependencies / local findings

- Biopython 1.84 (available in `D:\PyMOL\python.exe`).
- OrthoFinder 2.5.5 at `E:/Codes/gf_mito_b/orthofinder.exe`.
  - This build only accepts **`-S blast`**; its help mentions diamond, but
    `-S diamond` returns "Valid options are: {blast}".
  - It needs `blastp.exe`/`makeblastdb.exe` and `_internal/python38.dll` next
    to the executable. A self-contained dev runtime was assembled at
    `dev/bin/orthofinder/` (264 MB, gitignored).
  - DIAMOND 2.2.2 exists at `E:/Bsim/diamond.exe`, but this OrthoFinder build
    does not accept it. Keep it for other tools.
- OrthoFinder expects a relative `..\temp\WorkingDirectory` path; create
  `dev/temp/WorkingDirectory` (or set `cwd` accordingly) before running.
- On Windows, `subprocess` resolves `blastp` using the parent process PATH;
  the CLI temporarily prepends the blast directory to `os.environ["PATH"]`.

## Plugin integration (2026-11-02)

The prototype is now wired into the dev copy:

- Parser: `YR_MPE/plugins/components/methods/gbff_utils.py`
- Post-processing: `YR_MPE/plugins/components/methods/orthogroup_utils.py`
- Plugin/panel: `YR_MPE/plugins/orthofinder_plugin.py`
- Menu: `WORKSPACE → Dataset → Genomic Dataset Formation (OrthoFinder)`
- Factory getter: `factories/plugin_factory.py::get_orthofinder_plugin`
- Config: `YR_MPE/config.json` entry `OrthoFinder`
- Settings entry: `settings.json` (working tree only; the file already had
  pre-existing uncommitted changes, so it was not committed)
- Panel design: `dev/gbff_to_orthofinder/PANEL_DESIGN.md`

Headless checks passed:

- Plugin instantiates with `QtWebEngineWidgets` attribute set before
  `QApplication`.
- Parse-only thread: 4 species / 335 CDS.
- Full thread: 4 species / 335 CDS / 80 orthogroups / 62 single-copy OGs in
  19.2 s with 1 thread.
- Output tab: 62 table rows; cards show 4 / 335 / 80 / 62.

## Known issues / next steps

- GFF3 path is not implemented; GBFF is the priority (per 2026-11-02 decision).
- Alternative start codons are reported, not "corrected" (never rewrite the
  biological sequence). For codon analyses, use the CDS FASTA; for search,
  use the protein FASTA.
- Add `@software` citations (OrthoFinder, BLAST+, MCL, Biopython) — the panel
  Citation tab already lists them.
- Package OrthoFinder + BLAST + MCL via the shared runtime/manifest (T20).
- Add a real plugin test with 3–5 organelle GBFF files and a manifest entry.
- Decide how to import the single-copy OG set into the YR-MPE workspace
  (dataset item / ResultSpec); currently the plugin writes files and shows
  the summary table.
