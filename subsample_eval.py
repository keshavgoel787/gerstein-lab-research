#!/usr/bin/env python3
"""Summarise subsample.sh: how KING holds up as the number of variants shrinks.

For each variant set (common sites, or common + LD-pruned) and subset size,
pooled over all populations and averaged over the 5 random draws:

  match      inferred degree agrees with the pedigree (all pairs)
  PO found   pedigree parent-child pairs called 1st degree
  unrel FP   pedigree-unrelated pairs called 3rd degree or closer
  unrel sd   spread of the score for pedigree-unrelated pairs
  RMSE       root-mean-square difference from the full genome-wide score

Also prints the three original CHS pairs at each size.

Usage:
    python3 subsample_eval.py
"""
import csv
import glob
import gzip
import os
import re
from collections import defaultdict

import numpy as np

from pedigree_check import EXPECTED_DEGREE, degree, read_pedigree, relationship

PED = "data/20130606_g1k_3202_samples_ped_population.txt"
TRACKED = [("HG00405", "HG00403", "daughter–father"), ("HG00403", "HG00404", "the two parents"),
           ("HG00405", "HG00406", "different families")]


def read_kin0(path, opener=open):
    with opener(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            yield frozenset((r["#IID1"], r["IID2"])), float(r["KINSHIP"]), int(r["NSNP"])


def main():
    ped = read_pedigree(PED)
    full = {}
    for path in glob.glob("results/genome/king_*.kin0.gz"):
        for pair, phi, _ in read_kin0(path, gzip.open):
            full[pair] = phi
    rel = {pair: relationship(*sorted(pair), ped) for pair in full}

    # (set, n) -> seed -> list of (pair, phi)
    runs = defaultdict(lambda: defaultdict(list))
    name = re.compile(r"(\w+?)_(common|pruned)_(\w+)_s(\d+)\.kin0$")
    for path in glob.glob("data/sub/king/*.kin0"):
        pop, vset, n, seed = name.search(os.path.basename(path)).groups()
        for pair, phi, _ in read_kin0(path):
            runs[(vset, n)][seed].append((pair, phi))

    def size_key(k):
        return (k[0], float("inf") if k[1] == "all" else int(k[1]))

    rows = []
    print(f"{'variant set':<14} {'variants':>9} {'match':>8} {'PO found':>9} {'unrel FP':>9} "
          f"{'unrel sd':>9} {'RMSE':>7}")
    for key in sorted(runs, key=size_key):
        stats = []
        for seed, pairs in runs[key].items():
            ok = po = po_ok = un = un_fp = 0
            un_phi, err = [], []
            for pair, phi in pairs:
                r, d = rel[pair], degree(phi)
                ok += d == EXPECTED_DEGREE[r]
                if r == "parent-child":
                    po += 1
                    po_ok += d == 1
                if r == "none in pedigree":
                    un += 1
                    un_fp += d is not None
                    un_phi.append(phi)
                err.append(phi - full[pair])
            stats.append((ok / len(pairs), po_ok / po, un_fp / un, np.std(un_phi),
                          np.sqrt(np.mean(np.square(err)))))
        m = np.mean(stats, axis=0)
        vset, n = key
        label = "pruned (all)" if n == "all" else vset
        n_txt = "all" if n == "all" else f"{int(n):,}"
        print(f"{label:<14} {n_txt:>9} {m[0]:>8.2%} {m[1]:>9.1%} {m[2]:>9.2%} {m[3]:>9.4f} {m[4]:>7.4f}")
        rows.append((label, n_txt, *m))

    print("\nOriginal CHS pairs (mean ± sd over the 5 draws; full genome-wide score in brackets):")
    for a, b, desc in TRACKED:
        pair = frozenset((a, b))
        cells = []
        for key in sorted(runs, key=size_key):
            vals = [phi for pairs in runs[key].values() for p, phi in pairs if p == pair]
            n_txt = "all" if key[1] == "all" else f"{int(key[1]):,}"
            cells.append(f"{key[0][0]}{n_txt}: {np.mean(vals):.3f}±{np.std(vals):.3f}")
        print(f"  {desc} [{full[pair]:.3f}]\n    " + "  ".join(cells))

    with open("results/subsample/summary.tsv", "w") as f:
        f.write("variant_set\tvariants\tmatch\tpo_found\tunrelated_fp\tunrelated_sd\trmse_vs_full\n")
        for r in rows:
            f.write("\t".join([r[0], r[1]] + [f"{x:.5f}" for x in r[2:]]) + "\n")


if __name__ == "__main__":
    main()
