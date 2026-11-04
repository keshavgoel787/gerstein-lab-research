# Gerstein Lab research: kinship from VCFs

Work with Eric Ni: given two people's genomes, decide whether they are related. Built and tested on chromosome 22 of the 1000 Genomes high-coverage release. The math is written up in [KING_MATH.md](KING_MATH.md).

## Scripts

| Script | What it does |
|---|---|
| `kinship.py` | Compare two people. Streams the VCF, keeps biallelic SNVs, reports IBS0 and KING-robust kinship. Standard library only. |
| `extract.py` | One pass over the VCF into a genotype matrix (`.npz`) for all 3,202 samples. Drops positions that appear more than once (split multiallelics). Needs numpy. |
| `pedigree_check.py` | Scores every pair within each population, filtering to common variants per population, and checks the inferred degree against the pedigree. Needs numpy. |

```
python3 kinship.py data/<chr22 vcf.gz> HG00405 HG00403 --min-maf 0.05 --af-field AF_EAS
python3 extract.py data/<chr22 vcf.gz> data/chr22_snv --min-af 0
python3 pedigree_check.py data/chr22_snv.npz data/20130606_g1k_3202_samples_ped_population.txt --min-maf 0 0.01 0.05
```

For `kinship.py`, `--min-maf` keeps sites whose minor allele frequency (read from the INFO field named by `--af-field`) is at least that value. `AF_EAS` is the East Asian frequency, which suits the CHS samples.

## Data

Not committed; download it into `data/`. Use the `ftp://` address: over HTTPS, resumed downloads of the large VCF hang.

- Pedigree (who is related to whom): `20130606_g1k_3202_samples_ped_population.txt`
- Chromosome 22 phased panel (446 MB): `working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`

Both live under `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/`.

```
curl -C - -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

## Results (chr22)

### 1. Every pair vs. the pedigree

206,609 pairs, all within the same population, at minor allele frequency ≥ 5% in that population (76k–126k sites per population). Full output: `results/pedigree_check_chr22.txt`; every pair: `results/pairs_chr22_maf05.tsv.gz`.

| Pedigree says | Pairs | Mean φ ± sd | Inferred correctly |
|---|---|---|---|
| Parent–child | 1,206 | 0.253 ± 0.015 | 1,202 (99.7%) |
| Full siblings | 5 | 0.322 ± 0.034 | 4 |
| 2nd degree (half-sib, grandparent, aunt/uncle) | 11 | 0.13 ± 0.08 | 4 |
| Not related in the pedigree | 205,387 | **0.002 ± 0.042** | 185,012 (90.1%) |

- **Parent–child detection works**: 1,202 of 1,206 pairs land in the 1st-degree band.
- **Unrelated pairs average 0, not +0.045.** Last week's "uniform offset" came from looking at only two unrelated pairs. Across all 205k pairs the mean is 0.002. The real problem is spread: one chromosome gives a standard deviation of about 0.04, roughly the width of the 3rd-degree band (0.044). So about 10% of unrelated pairs get called 3rd degree by chance, and the two pairs from last week were among them.
- **Some "errors" are the pedigree's, not the score's.** 39 pedigree-unrelated pairs score as 1st degree and 5 as duplicates, with near-zero IBS0 (e.g. CHS HG00501/HG00512: φ = 0.487, IBS0 = 2). Some listed relationships also contradict the genotypes. HG00703 is listed as HG00656's grandson but scores 0.018 against him, with 4,724 IBS0 sites.

### 2. Rare vs. common variants

| Min MAF | Sites / population | Parent–child φ | Unrelated φ (mean ± sd) | Unrelated called 3rd degree or closer |
|---|---|---|---|---|
| all | 153k–322k | 0.2513 ± 0.0156 | 0.0016 ± 0.0399 | 7.5% |
| 1% | 111k–252k | 0.2519 ± 0.0155 | 0.0017 ± 0.0407 | 8.1% |
| 5% | 76k–126k | 0.2525 ± 0.0152 | 0.0024 ± 0.0423 | 9.9% |

Dropping rare variants barely moves the scores. Rare sites are almost never heterozygous in both people or opposite homozygotes, so they contribute little to the KING numerator. Filtering removes sites without removing noise.

The original three CHS pairs, rerun with `--min-maf 0.05 --af-field AF_EAS` (82,415 sites):

| Pair | Pedigree | Unfiltered φ | MAF ≥ 5% φ |
|---|---|---|---|
| HG00405 vs HG00403 | daughter and father | 0.297 | 0.301 |
| HG00403 vs HG00404 | the two parents | 0.045 | 0.047 |
| HG00405 vs HG00406 | different families | 0.051 | 0.056 |

### 3. Cross-check with PLINK 2

`plink2 --make-king-table` on all 163 CHS samples at MAF ≥ 5% (`results/plink_chs_maf05.kin0`). PLINK uses KING's between-family form, so it is compared against the same form computed in `pedigree_check.py`:

- 13,203 CHS pairs: correlation 0.99997, mean |difference| 0.0004, max 0.0017
- Same degree call for 99.8% of pairs
- The leftover gap comes from site lists: PLINK keeps 80,634 sites, I keep 79,574 because I also drop SNVs that share a position with an indel.

## Next steps

1. Run more chromosomes. The noise comes from using one short chromosome, and averaging across the genome should shrink it to the point where 3rd-degree calls mean something. This needs the cluster.
2. List the pedigree conflicts (unlisted relatives, contradicted relationships) and check them against the 1000G related-samples list.
3. Look at the 4 parent–child pairs that fall outside the 1st-degree band (e.g. KHV HG02059: φ ≈ 0.14 with IBS0 ≈ 0).
