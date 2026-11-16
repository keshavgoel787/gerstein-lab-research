#!/usr/bin/env python3
"""Plots of the KING results for the meeting write-up.

  results/plots/king_distribution.png   genome-wide KING score, pedigree-related vs not
  results/plots/chr22_vs_genome.png     unrelated pairs: chr22 alone vs 20 autosomes
  results/plots/heatmap_CHS.png         all CHS pairs, samples grouped by pedigree family

Usage:
    python3 plots.py
"""
import csv
import glob
import gzip

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from pedigree_check import read_pedigree, relationship

PED = "data/20130606_g1k_3202_samples_ped_population.txt"
OUT = "results/plots"
CUTOFFS = [(0.0442, "3rd"), (0.0884, "2nd"), (0.177, "1st")]

# reference palette (dataviz skill): categorical slots 1-2, blue sequential ramp
UNRELATED, RELATED = "#2a78d6", "#eb6834"
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUES = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 10, "text.color": INK,
    "axes.edgecolor": GRID, "axes.labelcolor": INK_2, "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE, "xtick.color": INK_2, "ytick.color": INK_2,
    "axes.spines.top": False, "axes.spines.right": False,
})


def load_genome(ped):
    rows = []
    for path in sorted(glob.glob("results/genome/king_*.kin0.gz")):
        with gzip.open(path, "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                rows.append((r["#IID1"], r["IID2"], float(r["KINSHIP"])))
    return rows


def cutoff_lines(ax, ymax_label):
    for x, name in CUTOFFS:
        ax.axvline(x, color=INK_2, lw=0.8, ls=(0, (2, 3)), zorder=0)
        ax.text(x + 0.004, ymax_label, name, color=INK_2, fontsize=8, va="top")


def distribution(rows, ped):
    rel = np.array([relationship(a, b, ped) != "none in pedigree" for a, b, _ in rows])
    phi = np.array([p for _, _, p in rows])
    bins = np.arange(-0.15, 0.40, 0.005)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist([phi[~rel], phi[rel]], bins=bins, stacked=True, color=[UNRELATED, RELATED],
            edgecolor=SURFACE, linewidth=0.5,
            label=[f"Not related in pedigree ({(~rel).sum():,} pairs)",
                   f"Related in pedigree ({rel.sum():,} pairs)"])
    ax.set_xlim(-0.16, 0.32)
    ax.set_yscale("log")
    ax.set_ylim(0.8, 3e5)
    cutoff_lines(ax, 2.5e5)
    ax.set_xlabel("KING kinship (between-family, plink2)", labelpad=2)
    ax.set_ylabel("Pairs (log scale)")
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.annotate("unrelated pairs\ncentred on 0", (0.0, 4e4), (-0.14, 2e4), color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.6))
    ax.annotate("parent–child and\nsiblings at 0.25", (0.245, 300), (0.12, 2000), color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.6))
    ax.annotate("20 unlisted 1st-degree pairs\n(the blue inside the orange peak)", (0.252, 5), (0.19, 2.2e4),
                color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.6))
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=2)
    ax.set_title("Genome-wide KING score for every within-population pair (20 autosomes)",
                 loc="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"{OUT}/king_distribution.png", dpi=200)


def chr22_vs_genome(rows, ped):
    chr22 = []
    with gzip.open("results/pairs_chr22_maf05.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["pedigree"] == "none in pedigree":
                chr22.append(float(r["king_between"]))
    genome = [p for a, b, p in rows if relationship(a, b, ped) == "none in pedigree"]
    bins = np.arange(-0.3, 0.3, 0.004)
    fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    for ax, vals, name in ((axes[0], np.array(chr22), "chr22 alone"),
                           (axes[1], np.array(genome), "20 autosomes")):
        frac_3rd = (vals > 0.0442).mean()
        ax.hist(vals, bins=bins, color=UNRELATED, edgecolor=SURFACE, linewidth=0.4)
        ax.axvline(0.0442, color=INK_2, lw=0.8, ls=(0, (2, 3)))
        ax.text(0.048, ax.get_ylim()[1] * 0.9, "3rd-degree cutoff", color=INK_2, fontsize=8, va="top")
        ax.text(-0.3, ax.get_ylim()[1] * 0.9,
                f"{name} · sd {vals.std():.3f}\n{frac_3rd:.1%} called 3rd degree or closer",
                fontsize=9, va="top")
        ax.set_ylabel("Pairs")
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    axes[1].set_xlabel("KING kinship, pairs the pedigree lists as unrelated (between-family estimator)")
    axes[0].set_title("More chromosomes shrink the noise for unrelated pairs",
                      loc="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"{OUT}/chr22_vs_genome.png", dpi=200)


def heatmap(rows, ped, pop="CHS"):
    samples = sorted((s for s, r in ped.items() if r["Population"] == pop),
                     key=lambda s: (ped[s]["FamilyID"], s))
    idx = {s: i for i, s in enumerate(samples)}
    m = np.full((len(samples), len(samples)), np.nan)
    for a, b, p in rows:
        if a in idx and b in idx:
            m[idx[a], idx[b]] = m[idx[b], idx[a]] = p
    cmap = LinearSegmentedColormap.from_list("blues", BLUES)
    cmap.set_bad(SURFACE)
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    im = ax.imshow(np.clip(m, 0, 0.3), cmap=cmap, vmin=0, vmax=0.3, interpolation="nearest")
    # outline each pedigree family with more than one sequenced member
    fams = [ped[s]["FamilyID"] for s in samples]
    start = 0
    for i in range(1, len(fams) + 1):
        if i == len(fams) or fams[i] != fams[start]:
            if i - start > 1:
                ax.add_patch(plt.Rectangle((start - 0.5, start - 0.5), i - start, i - start,
                                           fill=False, ec=RELATED, lw=0.9))
            start = i
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(f"{len(samples)} {pop} samples, grouped by pedigree family")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, ticks=[0, 0.0625, 0.125, 0.25])
    cb.ax.set_yticklabels(["0 unrelated", "0.0625 3rd", "0.125 2nd", "0.25 1st"])
    cb.outline.set_visible(False)
    ax.set_title(f"{pop}: KING kinship for every pair (orange boxes = pedigree families)",
                 loc="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(f"{OUT}/heatmap_{pop}.png", dpi=200)


def main():
    import os
    os.makedirs(OUT, exist_ok=True)
    ped = read_pedigree(PED)
    rows = load_genome(ped)
    distribution(rows, ped)
    chr22_vs_genome(rows, ped)
    heatmap(rows, ped)


if __name__ == "__main__":
    main()
