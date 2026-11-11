#!/usr/bin/env python3
"""Check plink2 --make-king-table output against the pedigree.

Labels each pair with its pedigree relationship (same rules as
pedigree_check.py) and reports how often the inferred KING degree matches.
plink2 reports KING's between-family estimator.

Usage:
    python3 king_vs_pedigree.py data/20130606_g1k_3202_samples_ped_population.txt results/genome/king_*.kin0
"""
import argparse
import csv
from collections import Counter, defaultdict

import numpy as np

from pedigree_check import DEGREE_NAME, EXPECTED_DEGREE, degree, read_pedigree, relationship


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pedigree")
    p.add_argument("kin0", nargs="+")
    args = p.parse_args()

    ped = read_pedigree(args.pedigree)
    tally = defaultdict(Counter)
    phis = defaultdict(list)
    sites = []
    conflicts = []
    for path in args.kin0:
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                a, b, phi = r["#IID1"], r["IID2"], float(r["KINSHIP"])
                rel = relationship(a, b, ped)
                d = degree(phi)
                tally[rel][d] += 1
                phis[rel].append(phi)
                sites.append(int(r["NSNP"]))
                if d != EXPECTED_DEGREE[rel] and (rel != "none in pedigree" or d in (0, 1)):
                    conflicts.append((ped[a]["Population"], a, b, rel, int(r["IBS0"]), phi, DEGREE_NAME[d]))

    print(f"{min(sites):,}-{max(sites):,} sites per population (MAF >= 5% within the population)")
    print(f"{'pedigree relationship':<18} {'pairs':>8} {'φ (mean ± sd)':>18}   inferred degree")
    for rel in EXPECTED_DEGREE:
        if rel not in tally:
            continue
        v = np.array(phis[rel])
        calls = ", ".join(f"{DEGREE_NAME[k]} {n}" for k, n in
                          sorted(tally[rel].items(), key=lambda kv: (kv[0] is None, kv[0])))
        print(f"{rel:<18} {len(v):>8,} {v.mean():>8.4f} ± {v.std():.4f}   {calls}")
    right = sum(tally[r][EXPECTED_DEGREE[r]] for r in tally)
    total = sum(sum(c.values()) for c in tally.values())
    print(f"degree matches pedigree for {right:,} of {total:,} pairs ({right / total:.2%})")

    print("\nConflicts (related pairs called wrong, or pedigree-unrelated pairs called 1st degree or closer):")
    for c in sorted(conflicts):
        print("  %s %s %s  pedigree: %s  IBS0 %d  φ %.4f  inferred: %s" % c)


if __name__ == "__main__":
    main()
