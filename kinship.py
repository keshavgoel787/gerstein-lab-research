#!/usr/bin/env python3
"""Estimate relatedness between two samples in a VCF.

Streams a (bgzipped) VCF, keeps biallelic SNVs where both samples have a
called genotype, and reports:

  - IBS0: fraction of sites where the two are opposite homozygotes (0/0 vs 1/1).
    Parent-offspring pairs share an allele at every site, so IBS0 is ~0.
  - KING-robust kinship (Manichaikul et al. 2010):
        phi = (N_het_het - 2 * N_ibs0) / (N_het_i + N_het_j)
    Expected ~0.5 for duplicates/MZ twins, ~0.25 for first-degree,
    ~0.125 for second-degree, ~0 for unrelated.

Usage:
    python3 kinship.py chr22.vcf.gz HG00405 HG00403
"""
import argparse
import gzip

# KING degree cutoffs: midpoints between 1/2^(k+1.5) expected values
DEGREES = [
    (0.354, "duplicate / MZ twin"),
    (0.177, "1st degree (parent-offspring or full sibling)"),
    (0.0884, "2nd degree"),
    (0.0442, "3rd degree"),
]


def dosage(gt):
    """'0|1' -> 1, '1|1' -> 2, missing or multi-allelic -> None."""
    # strip: the last sample column still carries the line's newline
    gt = gt.split(":", 1)[0].strip().replace("|", "/")
    if gt not in ("0/0", "0/1", "1/0", "1/1"):
        return None
    return gt.count("1")


def compare(vcf_path, a, b):
    n = het_het = ibs0 = het_a = het_b = 0
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                try:
                    ia, ib = header.index(a), header.index(b)
                except ValueError as e:
                    raise SystemExit(f"sample not found in VCF: {e}")
                continue
            cols = line.split("\t", max(ia, ib) + 1)
            ref, alt = cols[3], cols[4]
            # biallelic SNVs only
            if len(ref) != 1 or len(alt) != 1 or alt == ".":
                continue
            ga, gb = dosage(cols[ia]), dosage(cols[ib])
            if ga is None or gb is None:
                continue
            n += 1
            het_a += ga == 1
            het_b += gb == 1
            het_het += ga == 1 and gb == 1
            ibs0 += abs(ga - gb) == 2
    phi = (het_het - 2 * ibs0) / (het_a + het_b) if het_a + het_b else float("nan")
    return n, ibs0, phi


def classify(phi):
    if phi != phi:  # NaN: no usable sites
        return "insufficient data"
    for cutoff, label in DEGREES:
        if phi > cutoff:
            return label
    return "unrelated"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf")
    p.add_argument("sample_a")
    p.add_argument("sample_b")
    args = p.parse_args()

    n, ibs0, phi = compare(args.vcf, args.sample_a, args.sample_b)
    print(f"{args.sample_a} vs {args.sample_b}")
    print(f"  sites compared : {n:,}")
    print(f"  IBS0 fraction  : {ibs0 / n:.5f}" if n else "  IBS0 fraction  : n/a")
    print(f"  KING kinship   : {phi:.4f}")
    print(f"  inferred       : {classify(phi)}")


if __name__ == "__main__":
    main()
