#!/usr/bin/env python3
"""Run my own KING code genome-wide and compare it with plink2.

For each population, exports the same sites plink2 scored (MAF >= 5% within
the population, from data/genome/autosomes) to a PLINK 1 .bed, decodes the
genotypes with numpy, and computes both KING estimators with the counting
code from pedigree_check.py. Then it compares my between-family scores with
plink2's (results/genome/king_*.kin0.gz) and checks both estimators against
the pedigree.

Usage:
    python3 genome_check.py
"""
import csv
import gzip
import os
import subprocess
from collections import Counter, defaultdict

import numpy as np

from pedigree_check import EXPECTED_DEGREE, degree, king_counts, read_pedigree, relationship

PED = "data/20130606_g1k_3202_samples_ped_population.txt"
G = "data/genome"
TMP = "data/genome/bed_tmp"

# PLINK 1 .bed: 2 bits per sample, 4 samples per byte, low bits first.
# 00 = homozygous A1, 01 = missing, 10 = heterozygous, 11 = homozygous A2.
CODE = np.array([2, -1, 1, 0], dtype=np.int8)
LOOKUP = np.array([[CODE[(byte >> (2 * k)) & 3] for k in range(4)] for byte in range(256)], dtype=np.int8)


def read_bed(prefix, n_samples):
    raw = np.fromfile(prefix + ".bed", dtype=np.uint8)
    assert raw[:3].tolist() == [0x6C, 0x1B, 0x01], "not a variant-major .bed"
    per_variant = (n_samples + 3) // 4
    raw = raw[3:].reshape(-1, per_variant)
    return LOOKUP[raw].reshape(raw.shape[0], -1)[:, :n_samples]   # [sites x samples]


def plink_scores(pop):
    out = {}
    with gzip.open(f"results/genome/king_{pop}.kin0.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[frozenset((r["#IID1"], r["IID2"]))] = float(r["KINSHIP"])
    return out


def main():
    ped = read_pedigree(PED)
    pops = sorted({r["Population"] for r in ped.values()})
    mine_between, theirs, tally = [], [], {"within": defaultdict(Counter), "between": defaultdict(Counter)}
    for pop in pops:
        keep = f"{G}/keep_{pop}.txt"
        subprocess.run(["plink2", "--pfile", f"{G}/autosomes", "--keep", keep, "--maf", "0.05",
                        "--make-bed", "--threads", "8", "--out", TMP], check=True, stdout=subprocess.DEVNULL)
        samples = [line.split()[1] for line in open(TMP + ".fam")]
        d = read_bed(TMP, len(samples))
        if (d < 0).any():
            raise SystemExit(f"{pop}: missing genotypes in the .bed")
        het_het, ibs0, het = king_counts(d)
        ref = plink_scores(pop)
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                num = het_het[i, j] - 2 * ibs0[i, j]
                lo = min(het[i], het[j])
                within = num / (het[i] + het[j])
                between = num / (2 * lo) + 0.5 - (het[i] + het[j]) / (4 * lo)
                pair = frozenset((samples[i], samples[j]))
                mine_between.append(between)
                theirs.append(ref[pair])
                rel = relationship(samples[i], samples[j], ped)
                tally["within"][rel][degree(within)] += 1
                tally["between"][rel][degree(between)] += 1
        print(f"{pop}: {d.shape[0]:,} sites, {len(samples)} samples", flush=True)
    for ext in (".bed", ".bim", ".fam", ".log"):
        if os.path.exists(TMP + ext):
            os.remove(TMP + ext)

    a, b = np.array(mine_between), np.array(theirs)
    cut = [0.0442, 0.0884, 0.177, 0.354]
    print(f"\nMy between-family KING vs plink2 --make-king-table, {len(a):,} within-population pairs:")
    print(f"  correlation {np.corrcoef(a, b)[0, 1]:.6f}, mean |diff| {np.abs(a - b).mean():.2e}, "
          f"max |diff| {np.abs(a - b).max():.2e}, same degree call {(np.digitize(a, cut) == np.digitize(b, cut)).mean():.4%}")
    for est, t in tally.items():
        right = sum(t[r][EXPECTED_DEGREE[r]] for r in t)
        total = sum(sum(c.values()) for c in t.values())
        po = t["parent-child"]
        print(f"  my {est}-family estimator: matches pedigree for {right:,} of {total:,} pairs ({right / total:.2%}); "
              f"parent-child called 1st degree {po[1]:,} of {sum(po.values()):,}")


if __name__ == "__main__":
    main()
