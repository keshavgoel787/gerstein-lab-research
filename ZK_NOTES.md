# Notes: making KING circuit-friendly

First thoughts on what a zero-knowledge version of the kinship check needs, based on the experiments in this repo. These are not a design.

## 1. A fixed, public site panel

A circuit needs every proof to compare the same positions, so the sites must be chosen in advance and published. `fixed_panel.sh` tests this. It takes the 2,650,251 SNVs with MAF ≥ 5% in all 26 populations, draws one random panel (seed 1), and scores every within-population pair on that same panel:

| Fixed panel | Matches pedigree | Parent–child called 1st degree | Unrelated called 3rd degree or closer |
|---|---|---|---|
| 1,000 sites | 97.55% | 1,205 of 1,206 | 2.5% |
| 10,000 sites | 99.95% | 1,206 of 1,206 | 0.05% |

A 10,000-site panel does as well as the full 5–8 million sites per population (`results/genome/pedigree_check_genome.txt`). Panels: `results/panel/panel_*.txt`.

## 2. Integer-only arithmetic

Each person commits to one genotype per panel site, as three bits: hom-ref `r`, het `h`, hom-alt `a`, with `r + h + a = 1`. For persons *i* and *j*, every KING count is a sum of bit products:

- N_het,het = Σ hᵢ·hⱼ
- N_IBS0 = Σ (rᵢ·aⱼ + aᵢ·rⱼ)
- N_het,i = Σ hᵢ, N_het,j = Σ hⱼ (no products needed)

That is 3 multiplications per site, so about 30,000 for a 10,000-site panel, plus the bit constraints.

The within-family estimator is φ = (N_het,het − 2·N_IBS0) / (N_het,i + N_het,j). The circuit never has to divide. Proving "φ > cutoff" is the same as proving

  D · (N_het,het − 2·N_IBS0) > C · (N_het,i + N_het,j)

with the cutoff written as a fraction C/D (e.g. 3rd degree: 0.0442 ≈ 221/5000). The output can be a single relationship bit, or the degree, without revealing φ or any genotype.

The between-family estimator (the one plink2 uses) also needs min(N_het,i, N_het,j). That is one extra comparison and is still integer-only, but the within-family form is simpler.

## 3. Open questions for Eric

- Who chooses the panel, and does it have to be common in everyone's ancestry? This test used sites common in all 26 populations.
- Commitments: are the genotype bits committed in a Merkle tree, as in the paper, and opened inside the circuit?
- Is revealing the degree acceptable, or only "related / not related"?
- LD pruning didn't improve accuracy here (`results/subsample/`). Does it matter for privacy instead, e.g. fewer correlated sites leaking less?
