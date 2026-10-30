# Tool inventory (verified 2026-10-29)

Human-readable snapshot for the binary manager and documentation. The eventual
machine-readable manifest should be `YR_MPE/tools/manifest.json`; this file is
the source of truth until that exists.

## Wired in `YR_MPE/bin/` and `config.json`

| Tool | Path | Version | Notes |
|---|---|---|---|
| IQ-TREE 3 | `YR_MPE/bin/iqtree3/iqtree3.exe` | 3.0.1 (2025-05-05) | `settings.json` says 3.0.6; fix |
| ASTER/ASTRAL | `YR_MPE/bin/aster/astral4.exe`, `wastral.exe`, `astral-pro3.exe`, `astral.bat` | ASTRAL-IV 1.24.4.8 | `astral.bat` takes a version token |
| CASTER-site | `YR_MPE/bin/aster/caster-site.exe` | 1.24.2.6 | — |
| LSD2 | `YR_MPE/bin/lsd2/lsd2.exe` | 2.4.4 | no standalone `settings.json` entry |
| MrBayes | `YR_MPE/bin/mb/mb.exe` | 3.2.7a | BEAGLE/MS-MPI DLLs alongside |
| PhyloBayes | `YR_MPE/bin/pbmpi-1.9/` | 1.9 | also 1.8c |
| MPBoot | `YR_MPE/bin/mpboot/mpboot.exe` | 1.1.1 | folder says 1.2; fix |
| DecentTree | `YR_MPE/bin/decenttree/decenttree.exe` | 1.0.0 | no standalone `settings.json` entry |
| MAFFT | `YR_MPE/bin/mafft/mafft.bat` | 7.505 | — |
| MACSE | `YR_MPE/bin/macse/macse_v2.07.jar` + JRE | 2.07 | — |
| TrimAl | `YR_MPE/bin/trimAl/bin/trimal.exe` | 1.2rev57 | — |
| Gblocks | `YR_MPE/bin/Gblocks/Gblocks.exe` | 1.0 | manifest says 0.91b; fix |
| MUSCLE | `YR_MPE/bin/muscle5/muscle.exe` | 5.3 | — |
| Clustal Omega | `YR_MPE/bin/clustalo/clustalo.exe` | 1.2.2 | — |
| FastTree | `YR_MPE/bin/FastTree/FastTree.exe` | 2.1.11 No SSE3 | **not wired** |

## Available locally, not wired

| Tool | Path | Version | Notes |
|---|---|---|---|
| RAxML-NG | `E:\Codes\raxml-ng\build\bin\raxml-ng.exe` | 1.2.2-master | bundle `raxml-ng/bin` lacks MSYS DLLs (exit 127); package `build/bin` + 3 DLLs |
| OrthoFinder | `E:\Codes\gf_mito_b\orthofinder.exe` | 2.5.5 | BLAST+/MCL in the same folder |
| MCMCTree | `E:\Codes\PYML\release\mcmctree_{scalar,SSE,AVX2,AVX512}.exe` | PAML 4.10.10 | default to `scalar`; SIMD builds lack node scaling |
| baseml | `E:\Codes\gf_mito_b\baseml.exe` | PAML 4.10.7 | — |
| FastME | `E:\Codes\fastme\src\fastme.exe` | — | not wired |
| BMGE | `E:\Codes\BMGE\src` | — | source only |

## Not available locally

| Tool | Status | Impact |
|---|---|---|
| TreeRecs | no exe/JAR/source | P1 reconciliation blocked; need upstream + JRE |
| IQ2MC | no converter | MCMCTree preprocessing needs a custom converter |
| HyPhy | no exe/source | Selection HyPhy entry is a planned placeholder |
| PAML codeml | source only (`E:\Codes\paml\src\codeml.c`, 4.10.10) | compile before CodeML entry |
| PAML yn00 | source only (`E:\Codes\paml\src\yn00.c`, 4.10.10) | compile before YN00, or implement NG86/YN00 |

## Packaging hazards

- `Initializer.finish_install()` rewrites `YR_MPE/config.json` from the 6-tool
  `config.template.json`, dropping 10 of the 16 current entries. Fix before
  shipping a new Initializer.
- Version mismatches: IQ-TREE 3.0.6 vs 3.0.1; Gblocks 0.91b vs 1.0; MPBoot 1.2
  vs 1.1.1; LSD2 2.4.1 vs 2.4.4; `__version__` 0.2.0 vs manifest 0.1.0.
