# Gerstein Lab research: kinship from VCFs

Work with Eric Ni: given two people's genomes, decide whether they are related. Built and tested on chromosome 22 of the 1000 Genomes high-coverage release. The math is written up in [KING_MATH.md](KING_MATH.md).

## Scripts

| Script | What it does |
|---|---|
| `kinship.py` | Compare two people. Streams the VCF, keeps biallelic SNVs, reports IBS0 and KING-robust kinship. Standard library only. |
| `extract.py` | One pass over the VCF into a genotype matrix (`.npz`) for all 3,202 samples. Drops positions that appear more than once (split multiallelics). Needs numpy. |
| `pedigree_check.py` | Scores every pair within each population, filtering to common variants per population, and checks the inferred degree against the pedigree. Needs numpy. |
| `run_genome.sh` | Genome-wide run on a laptop: downloads each autosome, checks its MD5, converts SNVs to PLINK 2 format, deletes the VCF, then runs `plink2 --make-king-table` per population. Needs plink2 (`brew install plink-ng`). |
| `king_vs_pedigree.py` | Checks plink2 `.kin0` tables against the pedigree and lists the conflicts. |
| `subsample.sh` / `subsample_eval.py` | How few variants KING needs: scores every pair on random subsets of 100–100,000 common variants, with and without LD pruning, and compares against the pedigree and the full scores. |
| `plots.py` | The figures in `results/plots/`. |
| `genome_check.py` | Runs my own KING code on the same genome-wide sites plink2 used and compares the two. |
| `fixed_panel.sh` | One fixed site panel for everyone, as a zero-knowledge circuit would need. See [ZK_NOTES.md](ZK_NOTES.md). |

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

## Genome-wide results (20 autosomes)

`run_genome.sh` covers chromosomes 1–22 except 6 and 9: those two failed the manifest MD5 check on three separate downloads each. That's 4.8–7.9 million common SNVs per population (MAF ≥ 5% within the population), scored with plink2's between-family KING. Output: `results/genome/pedigree_check_genome.txt`, with per-population tables in `results/genome/king_*.kin0.gz`.

| Pedigree says | Pairs | Mean φ ± sd | Inferred correctly |
|---|---|---|---|
| Parent–child | 1,206 | 0.249 ± 0.004 | 1,206 (100%) |
| Full siblings | 5 | 0.252 ± 0.018 | 5 |
| 2nd degree (half-sib, grandparent, aunt/uncle) | 11 | 0.13 ± 0.02 | 11 |
| Not related in the pedigree | 205,387 | −0.008 ± 0.020 | 205,282 (99.95%) |

Overall, 206,504 of 206,609 pairs (99.95%) match the pedigree, up from 90% on chr22 alone. Only 0.05% of pedigree-unrelated pairs score above the 3rd-degree cutoff, and the 99.9th percentile is 0.025.

The three original CHS pairs, genome-wide:

| Pair | Pedigree | chr22 φ | Genome-wide φ | Call |
|---|---|---|---|---|
| HG00405 vs HG00403 | daughter and father | 0.297 | 0.252 | 1st degree ✓ |
| HG00403 vs HG00404 | the two parents | 0.045 | 0.004 | unrelated ✓ |
| HG00405 vs HG00406 | different families | 0.051 | 0.010 | unrelated ✓ |

**What this settles:**

- The false 3rd-degree calls on chr22 were noise from using one chromosome. With 20 chromosomes, the spread for unrelated pairs halves and the misfires mostly disappear.
- The few remaining conflicts are relatives the pedigree doesn't list. 20 pedigree-unrelated pairs score as 1st degree (φ ≈ 0.25). IBS0 says which kind: a few hundred sites means parent–child (e.g. ASW NA20317/NA20318, IBS0 = 299), and 60k–110k sites means full siblings (e.g. CHS HG00501/HG00512, IBS0 = 80,458).
- Some chr22 "pedigree conflicts" were chr22 noise, not pedigree errors. HG00656 vs his listed grandson HG00703 scored 0.018 on chr22 but 0.114 genome-wide (2nd degree, as listed); a grandchild can easily inherit none of a grandparent's DNA on one short chromosome. KHV HG02059 and her listed father score 0.246 genome-wide, a normal parent–child value.

### My code vs plink2, genome-wide

`genome_check.py` exports the sites plink2 scored to `.bed` files, decodes them with numpy and computes KING with the counting code from `pedigree_check.py` (`results/genome/my_king_vs_plink.txt`). Across all 206,609 pairs my between-family scores match plink2's to within 5 × 10⁻⁷, which is plink2's printing precision, and every degree call is the same. My within-family estimator matches the pedigree for 206,502 pairs (99.95%) and calls all 1,206 parent–child pairs correctly.

## How many variants does KING need? (for the zero-knowledge version)

`subsample.sh` keeps sites with MAF ≥ 5% in each population, optionally LD-prunes them (`plink2 --indep-pairwise 500kb 0.2`, about 100k–380k sites left per population), then scores all 206,609 pairs on random subsets, 5 draws per size. Full table: `results/subsample/summary.txt`; plot: `results/plots/subsample.png`.

| Variant set | Variants | Matches pedigree | Parent–child found | Strangers called related | sd, unrelated | RMSE vs full score |
|---|---|---|---|---|---|---|
| common | 100 | 80.9% | 85.2% | 19.1% | 0.120 | 0.129 |
| common | 1,000 | 97.2% | 100% | 2.84% | 0.038 | 0.035 |
| common | 10,000 | 99.95% | 100% | 0.05% | 0.022 | 0.011 |
| common | 100,000 | 99.95% | 100% | 0.05% | 0.020 | 0.003 |
| LD-pruned | 1,000 | 92.4% | 100% | 7.60% | 0.036 | 0.033 |
| LD-pruned | 10,000 | 99.87% | 100% | 0.13% | 0.019 | 0.018 |
| LD-pruned | all (~160k) | 99.92% | 100% | 0.08% | 0.017 | 0.017 |

- **10,000 random common variants reproduce the full result**: 99.95% agreement with the pedigree (the same as all ~5–8 million sites), and scores within about 0.01 of the full value.
- **1,000 variants still find every parent–child pair**, but 2.8% of strangers get called 3rd degree or closer.
- **100 is too few**: scores swing by ±0.06 between draws.
- **LD pruning didn't help** here and was slightly worse at 1,000–10,000 variants. A random draw of a few thousand sites from 3 billion bases is already spread out, so the sites are mostly uncorrelated anyway. Why pruned subsets are noisier is not yet checked (one guess: the pruned sites have a different allele-frequency mix).
- Daughter–father (full score 0.252) stays at 0.253 ± 0.005 with 10,000 common variants; the two unrelated pairs stay at 0.000 ± 0.007 and 0.012 ± 0.010.

### A fixed panel (what a circuit would use)

`fixed_panel.sh` draws one random panel from the 2,650,251 SNVs with MAF ≥ 5% in all 26 populations and scores every pair on that same panel (`results/panel/`):

| Fixed panel | Matches pedigree | Parent–child found | Strangers called related |
|---|---|---|---|
| 1,000 sites | 97.55% | 1,205 of 1,206 | 2.5% |
| 10,000 sites | 99.95% | 1,206 of 1,206 | 0.05% |

A shared 10,000-site panel is as good as the full genome. [ZK_NOTES.md](ZK_NOTES.md) sketches how the KING test becomes integer-only arithmetic over committed genotype bits.

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
- **Unlisted relatives.** 39 pedigree-unrelated pairs score as 1st degree and 5 as duplicates (e.g. CHS HG00501/HG00512: φ = 0.487, IBS0 = 2). Genome-wide data shows these are real 1st-degree relatives. A "duplicate" call on chr22 alone is a sibling pair that happened to inherit the same two copies of chr22 (HG00501/HG00512: 0.253 genome-wide). The apparent contradictions on chr22 (e.g. HG00703 scoring 0.018 against his listed grandfather) turn out to be chr22 noise; see the genome-wide section.

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

1. Get chr6 and chr9: retry the download, or fetch them on the cluster.
2. Check the 20 unlisted 1st-degree pairs against the 1000G related-samples list.
3. Find out why LD-pruned subsets are noisier than random common ones.
4. Settle the open questions in ZK_NOTES.md (panel choice, commitments, what the proof reveals), then build a first circuit for the 10,000-site panel.
