"""Print the numbers used in Supplementary Note S2."""

from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

path = Path(sys.argv[1] if len(sys.argv) > 1 else
            "tools/bench_2026-11-22/mrbayes_beagle.json")
rows = json.loads(path.read_text(encoding="utf-8"))

ratios, speedups = [], []
rf20, rf50 = [], []
for row in rows:
    pairs = list(zip(row["serial_times"], row["mpi_times"]))
    case_ratios = [mpi / serial for serial, mpi in pairs]
    case_speedups = [serial / mpi for serial, mpi in pairs]
    ratios.extend(case_ratios)
    speedups.extend(case_speedups)
    print("%-14s ratio %.3f  speedup %.2f  rf %s"
          % (row["case"], statistics.mean(case_ratios),
             statistics.mean(case_speedups), row["rf"]))
    (rf20 if row["ntaxa"] == 20 else rf50).extend(row["rf"])

print("ratio range %.3f-%.3f" % (min(ratios), max(ratios)))
print("speedup range %.2f-%.2f" % (min(speedups), max(speedups)))
print("RF 20:", {value: rf20.count(value) for value in sorted(set(rf20))})
print("RF 50:", {value: rf50.count(value) for value in sorted(set(rf50))})
print("RF total:", len(rf20) + len(rf50), "zero:",
      sum(1 for value in rf20 + rf50 if value == 0))
