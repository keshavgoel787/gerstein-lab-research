#!/usr/bin/env bash
# Genome-wide KING-robust for every within-population pair, run locally.
#
# For each autosome: download the 1000G high-coverage VCF, convert biallelic
# SNVs to PLINK 2 format, then delete the VCF. All 22 VCFs are ~31 GB, so
# only one is on disk at a time. Then merge the autosomes and score each
# population separately at minor allele frequency >= 5% in that population.
#
# Usage: ./run_genome.sh            (rerunnable; finished chromosomes are skipped)
set -euo pipefail
cd "$(dirname "$0")"

BASE=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
PED=data/20130606_g1k_3202_samples_ped_population.txt
OUT=data/genome
mkdir -p "$OUT" results/genome
[ -f "$OUT/manifest.txt" ] || curl -s "$BASE/20220804_manifest.txt" -o "$OUT/manifest.txt"

for c in $(seq 1 22); do
  [ -f "$OUT/chr$c.pgen" ] && { echo "chr$c already converted"; continue; }
  vcf=1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  local_vcf=data/$vcf
  size=$(curl -sI "$BASE/$vcf" | awk '/Content-Length/ {print $2}' | tr -d '\r')
  for attempt in 1 2 3; do
    [ -f "$local_vcf.done" ] && break
    # EBI drops long transfers; resume until the file is complete
    until [ -f "$local_vcf" ] && [ "$(stat -f %z "$local_vcf")" = "$size" ]; do
      curl -s -C - --speed-limit 50000 --speed-time 60 -o "$local_vcf" "$BASE/$vcf" || sleep 5
    done
    # a resumed transfer can come back the right size but corrupt: check the
    # manifest's MD5 and start over if it doesn't match
    want=$(awk -v f="$vcf" '$1 == f {print $2}' "$OUT/manifest.txt")
    if [ "$(md5 -q "$local_vcf")" = "$want" ]; then
      touch "$local_vcf.done"
    else
      echo "$(date +%T) chr$c failed MD5 check (attempt $attempt), downloading again"
      rm -f "$local_vcf"
    fi
  done
  if [ ! -f "$local_vcf.done" ]; then
    echo "$(date +%T) chr$c: still corrupt after 3 downloads, skipping it"
    continue
  fi
  echo "$(date +%T) chr$c downloaded, converting"
  plink2 --vcf "$local_vcf" --snps-only just-acgt --max-alleles 2 --set-all-var-ids '@:#' \
    --rm-dup exclude-all --make-pgen pvar-cols= --threads 8 --out "$OUT/chr$c" >/dev/null
  # keep chr22, the rest of the repo's scripts use it
  if [ "$c" != 22 ]; then rm -f "$local_vcf" "$local_vcf.done"; fi
  echo "$(date +%T) chr$c done"
done

for c in $(seq 1 22); do [ -f "$OUT/chr$c.pgen" ] && echo "$OUT/chr$c"; done > "$OUT/merge_list.txt"
echo "merging $(wc -l < "$OUT/merge_list.txt" | tr -d ' ') autosomes"
[ -f "$OUT/autosomes.pgen" ] || plink2 --pmerge-list "$OUT/merge_list.txt" pfile --make-pgen --threads 8 --out "$OUT/autosomes" >/dev/null

for pop in $(awk 'NR > 1 {print $6}' "$PED" | sort -u); do
  awk -v p="$pop" 'NR > 1 && $6 == p {print $2}' "$PED" > "$OUT/keep_$pop.txt"
  plink2 --pfile "$OUT/autosomes" --keep "$OUT/keep_$pop.txt" --maf 0.05 \
    --make-king-table counts --threads 8 --out "results/genome/king_$pop" >/dev/null
done

python3 king_vs_pedigree.py "$PED" results/genome/king_*.kin0 | tee results/genome/pedigree_check_genome.txt
