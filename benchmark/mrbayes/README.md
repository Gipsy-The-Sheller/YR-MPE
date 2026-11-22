# MrBayes acceleration benchmark

Replicates Supplementary Note S2: official stock MrBayes 3.2.9 (serial, no
SSE/AVX/FMA/BEAGLE/MPI) versus the CYAN MrBayes 3.2.8 build under MS-MPI and
BEAGLE 4 (eight processes, two CPU threads each).

- `bench_mrbayes.py` — simulates the four 20/50-OTU × 100/1000-bp datasets,
  runs both builds under identical priors and MCMC settings, and records
  wall-clock time and RF distance.
- `mrbayes_beagle.json` — the run reported in the manuscript
  (`ngen=200000`, `nruns=2`, `nchains=4`, three replicates).
- `summary.txt` — per-case means, ranges and RF counts.
- `summarize_mrbayes_bench.py` — regenerates `summary.txt`.
- `redraw_figS2.py` — regenerates Fig. S2 from the JSON.

Reproduce:

```bash
gh release download v2026.11.22 --repo Gipsy-The-Sheller/CYAN \
  --pattern "mrbayes-3.2.8-Windows-x64-*" --dir _mb
unzip _mb/mrbayes-3.2.8-Windows-x64-*.zip -d _mb/pkg
gh release download v3.2.9 -R NBISweden/MrBayes \
  --pattern "mb_3.2.9_win_x64.exe" --dir _mb

python bench_mrbayes.py --mb-dir _mb/pkg/mb \
  --serial-mb _mb/mb_3.2.9_win_x64.exe \
  --ngen 200000 --reps 3 --nruns 2 --mpi 8 --threads 2 \
  --out mrbayes_beagle.json
```
