"""Redraw Fig. S2 (MrBayes acceleration) from the benchmark JSON.

The input is ``tools/bench_2026-11-22/mrbayes_beagle.json``, written by
``tools/bench_mrbayes.py`` for MrBayes 3.2.8 + BEAGLE 4.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy.stats import gaussian_kde

matplotlib.use("Agg")
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.linewidth"] = 0.7

RED = "#ba3e45"
BLUE = "#4e699a"
TEXT = "#1f2937"
GRID = "#d1d5db"

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tools" / "bench_2026-11-22" / "mrbayes_beagle.json"
OUT_PDF = ROOT / "figures" / "figS2_mrbayes_restyled.pdf"
OUT_PNG = ROOT / "figures" / "figS2_mrbayes_restyled.png"


def load_rows():
    rows = json.loads(DATA.read_text(encoding="utf-8"))
    for row in rows:
        pairs = list(zip(row["serial_times"], row["mpi_times"]))
        row["speedups"] = [serial / mpi for serial, mpi in pairs]
    return rows


def case_label(row):
    return "%d OTU\n%d bp" % (row["ntaxa"], row["length"])


def style_axis(ax, title, ylabel):
    ax.set_title(title, fontsize=14, fontweight="bold", pad=9, color=TEXT)
    ax.set_ylabel(ylabel, fontsize=12, color=TEXT)
    ax.tick_params(axis="both", labelsize=10, colors=TEXT, width=0.6,
                   length=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.7)
    ax.spines["bottom"].set_linewidth(0.7)
    ax.yaxis.grid(True, linestyle="--", alpha=0.35, linewidth=0.6)
    ax.set_axisbelow(True)


def panel_a(ax):
    rows = load_rows()
    labels = [case_label(row) for row in rows] + ["Stock"]
    ratios = [
        np.array(row["mpi_times"], dtype=float)
        / np.array(row["serial_times"], dtype=float)
        for row in rows
    ]
    means = np.array([values.mean() for values in ratios] + [1.0])
    errors = np.array([
        values.std(ddof=1) if len(values) > 1 else 0.0 for values in ratios
    ] + [0.0])
    colors = [BLUE] * len(rows) + [RED]
    x = np.arange(len(labels))
    ax.bar(x, means, width=0.64, color=colors, edgecolor="white",
           linewidth=0.8, yerr=errors, capsize=4,
           error_kw={"elinewidth": 0.8, "capthick": 0.8,
                     "ecolor": "#374151"})
    for i, (mean, error) in enumerate(zip(means, errors)):
        label = "(1.0)" if i == len(labels) - 1 else (
            "%.3f±%.3f" % (mean, error)
        )
        ax.text(i, mean + error + 0.035, label, ha="center", va="bottom",
                fontsize=9, color=TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, 1.30)
    ax.set_yticks(np.arange(0, 1.3, 0.2))
    style_axis(ax, "Normalized Time Ratio\n(MPI+BEAGLE / Stock)",
               "Time Ratio")
    ax.legend(
        handles=[
            Patch(facecolor=BLUE, edgecolor="white", label="MPI+BEAGLE"),
            Patch(facecolor=RED, edgecolor="white", label="Stock"),
        ],
        loc="upper left", frameon=False, fontsize=9,
    )


def panel_b(ax):
    rows = load_rows()
    labels = [case_label(row) for row in rows] + ["Stock"]
    means = np.array([
        np.mean(row["speedups"]) for row in rows
    ] + [1.0])
    errors = np.array([
        np.std(row["speedups"], ddof=1) if len(row["speedups"]) > 1 else 0.0
        for row in rows
    ] + [0.0])
    colors = [BLUE] * len(rows) + [RED]
    x = np.arange(len(labels))
    ax.bar(x, means, width=0.64, color=colors, edgecolor="white",
           linewidth=0.8, yerr=errors, capsize=4,
           error_kw={"elinewidth": 0.8, "capthick": 0.8,
                     "ecolor": "#374151"})
    for i, (mean, error) in enumerate(zip(means, errors)):
        label = "(1.0)" if i == len(labels) - 1 else (
            "%.2f±%.2f" % (mean, error)
        )
        ax.text(i, mean + error + 0.04, label, ha="center", va="bottom",
                fontsize=9, color=TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    top = max(1.2, float((means + errors).max()) * 1.18)
    ax.set_ylim(0, top)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    style_axis(ax, "Speedup Factor\n(Stock / MPI+BEAGLE)", "Speedup")


def panel_c(ax):
    rows = load_rows()
    rf20 = [rf for row in rows if row["ntaxa"] == 20 for rf in row["rf"]]
    rf50 = [rf for row in rows if row["ntaxa"] == 50 for rf in row["rf"]]
    max_rf = max(rf20 + rf50) if (rf20 or rf50) else 0
    positions = np.arange(max_rf + 1)
    blue_counts = np.array([rf20.count(int(rf)) for rf in positions])
    red_counts = np.array([rf50.count(int(rf)) for rf in positions])
    totals = red_counts + blue_counts
    percentages = totals / totals.sum() * 100 if totals.sum() else totals
    ax.bar(positions, red_counts, width=0.62, color=RED,
           edgecolor="white", linewidth=0.8, label="50 OTU")
    ax.bar(positions, blue_counts, bottom=red_counts, width=0.62,
           color=BLUE, edgecolor="white", linewidth=0.8, label="20 OTU")
    for position, red, blue, total, percentage in zip(
        positions, red_counts, blue_counts, totals, percentages
    ):
        if red:
            ax.text(position, red / 2, str(red), ha="center", va="center",
                    color="white", fontsize=9, fontweight="bold")
        if blue:
            ax.text(position, red + blue / 2, str(blue), ha="center",
                    va="center", color="white", fontsize=9,
                    fontweight="bold")
        ax.text(position, total + 0.7, "%.1f%%" % percentage,
                ha="center", va="bottom", fontsize=9, color=TEXT)
    ax.set_xticks(positions)
    ax.set_xticklabels(["%.1f" % rf for rf in positions], fontsize=9)
    ax.set_xlabel("Robinson–Foulds Distance", fontsize=12, color=TEXT)
    ax.set_ylim(0, max(5, totals.max() * 1.35))
    ax.set_yticks(np.linspace(0, max(5, totals.max() * 1.35), 5))
    style_axis(
        ax,
        "RF Distance Frequency Distribution\n(All datasets, n=%d)"
        % totals.sum(),
        "Frequency")
    ax.legend(loc="upper right", frameon=False, fontsize=9)


def panel_d(ax):
    rows = load_rows()
    rf20 = np.array(
        [rf for row in rows if row["ntaxa"] == 20 for rf in row["rf"]],
        dtype=float)
    rf50 = np.array(
        [rf for row in rows if row["ntaxa"] == 50 for rf in row["rf"]],
        dtype=float)
    x = np.linspace(0, 100, 1200)
    peak = 0.0
    for values, color, label in [
        (rf20, BLUE, "20 OTU"),
        (rf50, RED, "50 OTU"),
    ]:
        if len(values) < 2 or np.std(values) == 0:
            # A single observed value cannot support a KDE; mark the point.
            ax.axvline(values[0] if len(values) else 0, color=color,
                       linewidth=1.8, alpha=0.75, label=label)
            peak = max(peak, 0.2)
            continue
        density = gaussian_kde(values, bw_method=0.55)(x)
        peak = max(peak, float(density.max()))
        ax.plot(x, density, color=color, linewidth=1.8, label=label)
    ax.axvline(34, color=BLUE, linestyle="--", linewidth=1.1, alpha=0.8)
    ax.axvline(94, color=RED, linestyle="--", linewidth=1.1, alpha=0.8)
    ax.text(20, peak * 1.12, "20 OTU\ntheoretical max RF = 34", color=BLUE,
            ha="center", va="center", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": BLUE, "alpha": 0.8,
                  "pad": 2.5})
    ax.text(70, peak * 0.95, "50 OTU\ntheoretical max RF = 94", color=RED,
            ha="center", va="center", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": RED, "alpha": 0.8,
                  "pad": 2.5})
    ax.set_xlim(0, 100)
    ax.set_ylim(0, max(0.25, peak * 1.35))
    ax.set_xticks(np.arange(0, 101, 20))
    ax.set_xlabel("Robinson–Foulds Distance", fontsize=12, color=TEXT)
    style_axis(ax, "RF Distance Distribution Overview (KDE)\n"
                   "(All datasets, n=%d)" % (len(rf20) + len(rf50)),
               "Density")
    ax.legend(loc="upper right", frameon=False, fontsize=9)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    panel_a(axes[0])
    panel_b(axes[1])
    for axis, letter in zip(axes.flat, "AB"):
        axis.text(-0.10, 1.04, letter, transform=axis.transAxes,
                  fontsize=12, fontweight="bold", va="bottom", color=TEXT)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.16,
                        wspace=0.24)
    fig.savefig(OUT_PDF, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT_PNG, dpi=180, bbox_inches="tight", pad_inches=0.02)
    print(OUT_PDF)
    print(OUT_PNG)


if __name__ == "__main__":
    main()
