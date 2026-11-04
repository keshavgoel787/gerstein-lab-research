#!/usr/bin/env python3
"""Score every within-population pair with KING-robust and check it against the pedigree.

Reads the matrix written by extract.py. For each population, keeps sites whose
minor allele frequency *in that population* is at least --min-maf, computes
both KING estimators for every pair, labels each pair with its pedigree
relationship, and reports how often the inferred degree matches.

  within-family  phi = (N_het,het - 2 N_ibs0) / (N_het,i + N_het,j)       (kinship.py)
  between-family phi = (N_het,het - 2 N_ibs0) / (2 min(N_het,i, N_het,j))
                       + 1/2 - (N_het,i + N_het,j) / (4 min(N_het,i, N_het,j))   (plink2 --make-king)

Usage:
    python3 pedigree_check.py data/chr22_snv.npz data/20130606_g1k_3202_samples_ped_population.txt \
        --min-maf 0 0.01 0.05 --pairs-out results/pairs_maf05.tsv
"""
import argparse
import csv
from collections import Counter, defaultdict

import numpy as np

DEGREES = [(0.354, 0), (0.177, 1), (0.0884, 2), (0.0442, 3)]  # cutoff -> degree
DEGREE_NAME = {0: "duplicate", 1: "1st", 2: "2nd", 3: "3rd", None: "unrelated"}
EXPECTED_DEGREE = {"parent-child": 1, "full sibling": 1, "half sibling": 2,
                   "grandparent": 2, "avuncular": 2, "none in pedigree": None}


def degree(phi):
    for cutoff, d in DEGREES:
        if phi > cutoff:
            return d
    return None


def read_pedigree(path):
    ped = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter=" "):
            ped[row["SampleID"]] = row
    return ped


def relationship(a, b, ped):
    """Pedigree relationship between two samples, from the father/mother columns."""
    def parents(x):
        r = ped[x]
        return {p for p in (r["FatherID"], r["MotherID"]) if p != "0"}
    pa, pb = parents(a), parents(b)
    if a in pb or b in pa:
        return "parent-child"
    shared = pa & pb
    if len(shared) == 2:
        return "full sibling"
    if len(shared) == 1:
        return "half sibling"
    if any(a in parents(p) for p in pb if p in ped) or any(b in parents(p) for p in pa if p in ped):
        return "grandparent"
    # a is a sibling of one of b's parents (or vice versa)
    for x, py in ((a, pb), (b, pa)):
        px = parents(x)
        for p in py:
            if p in ped and len(px) == 2 and px == parents(p):
                return "avuncular"
    return "none in pedigree"


def king_counts(d):
    """Pairwise het/het, IBS0 and per-sample het counts for a [sites x samples] dosage block."""
    k = d.shape[1]
    het_het = np.zeros((k, k))
    ibs0 = np.zeros((k, k))
    het = np.zeros(k)
    for s in range(0, d.shape[0], 100_000):
        c = d[s:s + 100_000]
        h = (c == 1).astype(np.float32)
        h0 = (c == 0).astype(np.float32)
        h2 = (c == 2).astype(np.float32)
        het_het += h.T @ h
        x = h0.T @ h2
        ibs0 += x + x.T
        het += h.sum(0)
    return het_het, ibs0, het


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("matrix")
    p.add_argument("pedigree")
    p.add_argument("--min-maf", type=float, nargs="+", default=[0, 0.05])
    p.add_argument("--pairs-out", help="write every pair at the last --min-maf to this TSV")
    args = p.parse_args()

    z = np.load(args.matrix)
    dosage, samples = z["dosage"], list(z["samples"])
    if (dosage < 0).any():
        raise SystemExit("matrix has missing genotypes; this script assumes a fully called panel")
    ped = read_pedigree(args.pedigree)
    by_pop = defaultdict(list)
    for i, s in enumerate(samples):
        by_pop[ped[s]["Population"]].append(i)

    for maf in args.min_maf:
        tally = defaultdict(Counter)   # relationship -> Counter(inferred degree)
        phis = defaultdict(lambda: defaultdict(list))  # estimator -> relationship -> [phi]
        rows = []
        n_sites = []
        for pop, idx in sorted(by_pop.items()):
            d = dosage[:, idx]
            alt = d.sum(1) / (2 * len(idx))
            d = d[np.minimum(alt, 1 - alt) >= maf] if maf > 0 else d[(alt > 0) & (alt < 1)]
            n_sites.append(d.shape[0])
            het_het, ibs0, het = king_counts(d)
            for i in range(len(idx)):
                for j in range(i + 1, len(idx)):
                    a, b = samples[idx[i]], samples[idx[j]]
                    num = het_het[i, j] - 2 * ibs0[i, j]
                    lo = min(het[i], het[j])
                    within = num / (het[i] + het[j])
                    between = num / (2 * lo) + 0.5 - (het[i] + het[j]) / (4 * lo)
                    rel = relationship(a, b, ped)
                    tally[rel][degree(within)] += 1
                    phis["within"][rel].append(within)
                    phis["between"][rel].append(between)
                    rows.append((pop, a, b, rel, d.shape[0], int(ibs0[i, j]), within, between))

        print(f"\n=== min MAF {maf} (per-population; {min(n_sites):,}-{max(n_sites):,} sites per population) ===")
        print(f"{'pedigree relationship':<18} {'pairs':>8} {'within φ':>18} {'between φ':>18}   inferred degree (within)")
        for rel in EXPECTED_DEGREE:
            if rel not in tally:
                continue
            w, bt = np.array(phis["within"][rel]), np.array(phis["between"][rel])
            calls = ", ".join(f"{DEGREE_NAME[k]} {v}" for k, v in sorted(tally[rel].items(), key=lambda kv: (kv[0] is None, kv[0])))
            print(f"{rel:<18} {len(w):>8,} {w.mean():>8.4f} ± {w.std():.4f} {bt.mean():>8.4f} ± {bt.std():.4f}   {calls}")
        right = sum(tally[r][EXPECTED_DEGREE[r]] for r in tally)
        total = sum(sum(c.values()) for c in tally.values())
        print(f"degree matches pedigree for {right:,} of {total:,} pairs ({right / total:.2%})")

    if args.pairs_out:
        with open(args.pairs_out, "w") as f:
            f.write("pop\tsample_a\tsample_b\tpedigree\tsites\tibs0\tking_within\tking_between\tinferred_within\n")
            for pop, a, b, rel, n, i0, w, bt in rows:
                f.write(f"{pop}\t{a}\t{b}\t{rel}\t{n}\t{i0}\t{w:.4f}\t{bt:.4f}\t{DEGREE_NAME[degree(w)]}\n")


if __name__ == "__main__":
    main()
