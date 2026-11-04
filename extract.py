#!/usr/bin/env python3
"""Turn a phased VCF into a genotype matrix for all-pairs kinship.

One pass over the VCF. Keeps biallelic SNVs (dropping positions that appear
more than once, i.e. split multiallelics) whose overall alt-allele frequency
is at least --min-af, and writes:

  <out>.npz  dosage   int8 [sites x samples], 0/1/2 alt copies, -1 missing
             af       float32 [sites], overall AF from INFO
             pos      int64 [sites]
             samples  str [samples]

Usage:
    python3 extract.py data/chr22.vcf.gz data/chr22_snv
"""
import argparse
import subprocess

import numpy as np


def parse_af(info):
    for field in info.split(";"):
        if field.startswith("AF="):
            return float(field[3:])
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf")
    p.add_argument("out")
    p.add_argument("--min-af", type=float, default=0.001,
                   help="drop sites with overall AF below this (default 0.001)")
    args = p.parse_args()

    # decompress in a separate process so it runs in parallel with parsing
    gz = subprocess.Popen(["gzip", "-cd", args.vcf], stdout=subprocess.PIPE, bufsize=1 << 20)
    rows, afs, poss = [], [], []
    seen_pos = {}
    for line in gz.stdout:
        if line.startswith(b"##"):
            continue
        if line.startswith(b"#CHROM"):
            samples = line.rstrip(b"\n").decode().split("\t")[9:]
            n = len(samples)
            continue
        fixed = line.split(b"\t", 9)
        pos, ref, alt, info = int(fixed[1]), fixed[3], fixed[4], fixed[7].decode()
        seen_pos[pos] = seen_pos.get(pos, 0) + 1
        if len(ref) != 1 or len(alt) != 1 or alt == b".":
            continue
        af = parse_af(info)
        if af is None or af < args.min_af:
            continue
        gts = fixed[9].rstrip(b"\n")
        if len(gts) == 4 * n - 1:
            # every genotype is "a|b": fixed width, read alleles at offsets 0 and 2
            buf = np.frombuffer(gts + b"\t", dtype=np.uint8).reshape(n, 4)
            a, b = buf[:, 0], buf[:, 2]
            if not (np.isin(a, (48, 49)).all() and np.isin(b, (48, 49)).all()):
                continue  # something other than 0/1 alleles
            d = (a - 48) + (b - 48)
            d = d.astype(np.int8)
        else:
            d = np.full(n, -1, dtype=np.int8)
            for i, g in enumerate(gts.split(b"\t")):
                g = g.split(b":", 1)[0].replace(b"|", b"/")
                if g in (b"0/0", b"0/1", b"1/0", b"1/1"):
                    d[i] = g.count(b"1")
        rows.append(d)
        afs.append(af)
        poss.append(pos)
    gz.wait()

    pos = np.array(poss, dtype=np.int64)
    keep = np.array([seen_pos[x] == 1 for x in poss], dtype=bool)
    dosage = np.vstack(rows)[keep]
    print(f"kept {keep.sum():,} biallelic SNVs at unique positions "
          f"({(~keep).sum():,} dropped as duplicate positions), {n} samples")
    np.savez(args.out, dosage=dosage, af=np.array(afs, dtype=np.float32)[keep],
             pos=pos[keep], samples=np.array(samples))


if __name__ == "__main__":
    main()
