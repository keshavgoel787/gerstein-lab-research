#!/usr/bin/env bash
# How few variants does KING need? For each population: keep sites with
# MAF >= 5% in that population, LD-prune them (plink2 --indep-pairwise,
# 500 kb window, r^2 < 0.2), then score every pair on random subsets of
# 100 to 100,000 variants, drawn from both the pruned and the unpruned
# common sites, 5 random draws each. subsample_eval.py summarises the
# results against the pedigree and the full genome-wide scores.
#
# Needs data/genome/autosomes.* from run_genome.sh.
# Usage: ./subsample.sh
set -euo pipefail
cd "$(dirname "$0")"

PED=data/20130606_g1k_3202_samples_ped_population.txt
G=data/genome
OUT=data/sub
mkdir -p "$OUT/king" results/subsample

for pop in $(awk 'NR > 1 {print $6}' "$PED" | sort -u); do
  keep=$G/keep_$pop.txt
  [ -f "$keep" ] || awk -v p="$pop" 'NR > 1 && $6 == p {print $2}' "$PED" > "$keep"
  if [ ! -f "$OUT/${pop}_pruned.pgen" ]; then
    plink2 --pfile $G/autosomes --keep "$keep" --maf 0.05 --make-pgen --threads 8 --out "$OUT/${pop}_common" >/dev/null
    plink2 --pfile "$OUT/${pop}_common" --indep-pairwise 500kb 0.2 --threads 8 --out "$OUT/${pop}_prune" >/dev/null
    plink2 --pfile "$OUT/${pop}_common" --extract "$OUT/${pop}_prune.prune.in" --make-pgen --threads 8 \
      --out "$OUT/${pop}_pruned" >/dev/null
  fi
  # all pruned sites, no subsampling
  plink2 --pfile "$OUT/${pop}_pruned" --make-king-table --threads 8 --out "$OUT/king/${pop}_pruned_all_s0" >/dev/null
  for set in common pruned; do
    for n in 100 1000 10000 100000; do
      if [ "$set" = pruned ] && [ "$n" -gt "$(wc -l < "$OUT/${pop}_prune.prune.in")" ]; then continue; fi
      for seed in 1 2 3 4 5; do
        plink2 --pfile "$OUT/${pop}_${set}" --thin-count $n --seed $seed --make-king-table --threads 8 \
          --out "$OUT/king/${pop}_${set}_${n}_s${seed}" >/dev/null
      done
    done
  done
  rm -f "$OUT/${pop}_common".*   # several hundred MB per population; the pruned set is kept
  echo "$(date +%T) $pop done ($(wc -l < "$OUT/${pop}_prune.prune.in" | tr -d ' ') pruned sites)"
done

python3 subsample_eval.py | tee results/subsample/summary.txt
