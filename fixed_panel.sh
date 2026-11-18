#!/usr/bin/env bash
# A fixed site panel, as a zero-knowledge circuit would need: one list of
# sites decided in advance and used for every pair, instead of a fresh
# random draw per population.
#
# Takes sites with MAF >= 5% in every one of the 26 populations, draws one
# fixed random panel of 1,000 and one of 10,000 (seed 1), and scores every
# within-population pair on each panel with plink2. king_vs_pedigree.py
# checks the calls against the pedigree.
#
# Needs data/genome/autosomes.* from run_genome.sh.
# Usage: ./fixed_panel.sh
set -euo pipefail
cd "$(dirname "$0")"

PED=data/20130606_g1k_3202_samples_ped_population.txt
G=data/genome
OUT=data/panel
RES=results/panel
mkdir -p "$OUT" "$RES"
pops=$(awk 'NR > 1 {print $6}' "$PED" | sort -u)

# sites common in each population, then the ones common in all of them
for pop in $pops; do
  [ -f "$OUT/common_$pop.txt" ] && continue
  plink2 --pfile $G/autosomes --keep $G/keep_$pop.txt --maf 0.05 --write-snplist --threads 8 \
    --out "$OUT/common_$pop" >/dev/null
  mv "$OUT/common_$pop.snplist" "$OUT/common_$pop.txt"
done
cat "$OUT"/common_*.txt | sort | uniq -c | awk -v n="$(echo $pops | wc -w)" '$1 == n {print $2}' > "$OUT/common_all.txt"
echo "$(wc -l < "$OUT/common_all.txt" | tr -d ' ') sites have MAF >= 5% in all populations"

for n in 1000 10000; do
  # one fixed panel, the same for everyone
  python3 -c "import random, sys; random.seed(1); s = open(sys.argv[1]).read().split(); print('\n'.join(sorted(random.sample(s, int(sys.argv[2])))))" \
    "$OUT/common_all.txt" $n > "$RES/panel_$n.txt"
  for pop in $pops; do
    plink2 --pfile $G/autosomes --keep $G/keep_$pop.txt --extract "$RES/panel_$n.txt" \
      --make-king-table counts --threads 8 --out "$OUT/king_${n}_$pop" >/dev/null
  done
  echo "=== fixed panel of $n sites ==="
  python3 king_vs_pedigree.py "$PED" "$OUT"/king_${n}_*.kin0 | tee "$RES/pedigree_check_panel_$n.txt" | sed -n '1,10p'
done
