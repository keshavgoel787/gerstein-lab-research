# Gerstein Lab research: kinship from VCFs

First task for my work with Eric Ni: given two people's genomes, decide whether they are related. Built and tested on chromosome 22 of the 1000 Genomes high-coverage release.

## What the script does

`kinship.py` streams a VCF and compares two samples site by site. It keeps biallelic SNVs where both samples have a genotype call, converts each genotype to 0, 1 or 2 copies of the alternate allele, and reports:

- **IBS0**: the fraction of sites where the two are opposite homozygotes (0/0 vs 1/1). A parent and child always share an allele, so this is close to zero for them.
- **KING-robust kinship** (Manichaikul et al., 2010): about 0.25 for parent–child or full siblings, 0.125 for second-degree relatives, about 0 for unrelated people.

It uses only the Python standard library.

```
python3 kinship.py data/<chr22 vcf.gz> HG00405 HG00403
```

## Data

Not committed; download it into `data/`. Use the `ftp://` address: over HTTPS, resumed downloads of the large VCF hang.

- Pedigree (who is related to whom): `20130606_g1k_3202_samples_ped_population.txt`
- Chromosome 22 phased panel (446 MB): `working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`

Both live under `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/`.

```
curl -C - -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
```

## Results (chr22, 925,730 sites per pair)

All three pairs are Southern Han Chinese (CHS). Raw output is in `results_*.txt`.

| Pair | Pedigree | IBS0 | KING φ | Script's label |
|---|---|---|---|---|
| HG00405 vs HG00403 | daughter and father | 0.00001 | 0.297 | 1st degree ✓ |
| HG00403 vs HG00404 | the two parents | 0.00495 | 0.045 | 3rd degree ✗ |
| HG00405 vs HG00406 | different families | 0.00426 | 0.051 | 3rd degree ✗ |

## Findings

- IBS0 separates the parent–child pair from the unrelated pairs by about 400×.
- Every kinship score sits about 0.045 above its expected value, which pushes the unrelated pairs just past the 3rd-degree cutoff (0.044). The shift is the same for every pair, so it likely comes from the input sites rather than the formula. This is untested.

## Next steps

1. Filter out rare variants and duplicate positions, then re-run the same pairs.
2. Compare the scores with plink2 or KING on the cluster.
3. Score every pair within a population against the pedigree, then run all chromosomes.
