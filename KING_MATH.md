# The math behind KING-robust kinship

Source: Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010). Robust relationship inference in genome-wide association studies. *Bioinformatics* 26(22):2867–2873. doi:10.1093/bioinformatics/btq559

## 1. What a kinship coefficient is

Pick one allele at random from person *i* and one from person *j* at the same spot in the genome. The **kinship coefficient φ** is the probability that the two are *identical by descent* (IBD): copies of the same ancestral DNA, inherited through a shared ancestor.

**Parent and child.** The child got exactly one of their two alleles from this parent. Pick one of the child's alleles at random: there is a 1/2 chance it is the one from this parent. Pick one of the parent's two alleles at random: there is a 1/2 chance it is the very one they passed on. So φ = 1/2 × 1/2 = **1/4**.

**Identical twins.** Both twins carry the same two copies. Pick one allele from each, and there is a 1/2 chance you picked the same copy. φ = **1/2**.

**Each meiosis halves φ:** full siblings 1/4, grandparent / half-sibling / aunt–uncle 1/8, first cousins 1/16, unrelated 0.

The IBD-sharing probabilities k0, k1, k2 (sharing 0, 1 or 2 alleles IBD at a site) relate to φ by **φ = k1/4 + k2/2**:

| Relationship | k0 | k1 | k2 | φ |
|---|---|---|---|---|
| Identical twins | 0 | 0 | 1 | 1/2 |
| Parent–child | 0 | 1 | 0 | 1/4 |
| Full siblings | 1/4 | 1/2 | 1/4 | 1/4 |
| 2nd degree | 1/2 | 1/2 | 0 | 1/8 |
| Unrelated | 1 | 0 | 0 | 0 |

Parent–child and full siblings have the same φ but different k0. That is why IBS0 (sites where the two people share no allele) is ~0 for parent–child but not for siblings: a sibling pair has a 1/4 chance of sharing nothing IBD at any site.

## 2. Why the estimator works

At one biallelic site with allele frequencies *p* and *q = 1 − p*, under Hardy–Weinberg equilibrium a random person is heterozygous with probability 2pq.

KING compares two counts that are equal by chance but diverge when people are related:

- **Both heterozygous (Aa, Aa).** For two unrelated people: 2pq × 2pq = **4p²q²**.
- **Opposite homozygotes (AA, aa) in either order.** For two unrelated people: 2 × p² × q² = **2p²q²**.

Hence for unrelated people:

  E[N_Aa,Aa − 2·N_AA,aa] = 4p²q² − 2 × 2p²q² = **0**

Being related **raises** the het/het count (if they share an allele IBD, both being heterozygous is more likely) and **lowers** the opposite-homozygote count (it becomes impossible when they share an allele IBD). Work through the IBD states:

| IBD state (prob.) | P(Aa, Aa) | P(AA, aa or aa, AA) | Contribution to (het,het − 2·IBS0) |
|---|---|---|---|
| share 0 (k0) | 4p²q² | 2p²q² | 0 |
| share 1 (k1) | pq | 0 | pq |
| share 2 (k2) | 2pq | 0 | 2pq |

So E[N_Aa,Aa − 2·N_AA,aa] per site = pq·(k1 + 2·k2) = **4pq · φ**.

The expected heterozygosity of one person is 2pq per site, so the two people together have 4pq. Dividing gives the **within-family estimator** (what `kinship.py` uses):

  φ̂ = (N_Aa,Aa − 2·N_AA,aa) / (N_Aa^(i) + N_Aa^(j))  ≈  4pq·φ / 4pq = φ

**Why "robust":** the allele frequency *p* never appears in the formula. The person's own heterozygosity stands in for 2pq. Older estimators (PLINK 1.9 `--genome`, for example) plug in population allele frequencies, so they go wrong when the sample mixes ancestries.

**Check with parent–child:** k1 = 1, so the numerator is pq per site and the denominator is 2pq + 2pq = 4pq, giving 1/4. ✓

## 3. The between-family estimator (what PLINK reports)

People from different sub-populations, or with different amounts of inbreeding, have different heterozygosity. KING's second form normalises by the *less* heterozygous person and corrects for the gap:

  φ̂ = (N_Aa,Aa − 2·N_AA,aa) / (2·min(N_Aa^(i), N_Aa^(j)))  +  1/2  −  (N_Aa^(i) + N_Aa^(j)) / (4·min(N_Aa^(i), N_Aa^(j)))

If the two people are equally heterozygous, the last two terms cancel and this equals the within-family form. `plink2 --make-king` uses this form for every pair.

## 4. Degree cutoffs

The expected φ halves with each degree, so the cutoffs sit at the geometric midpoints, 2^−(k+1.5):

| Inferred | φ range | Expected |
|---|---|---|
| Duplicate / MZ twin | > 0.354 | 0.5 |
| 1st degree | 0.177 – 0.354 | 0.25 |
| 2nd degree | 0.0884 – 0.177 | 0.125 |
| 3rd degree | 0.0442 – 0.0884 | 0.0625 |
| Unrelated | < 0.0442 | 0 |

## 5. Things the formula assumes

- **Many independent sites.** φ̂ is an average over sites. One chromosome holds only a few dozen independent chunks of inheritance, so a single-chromosome estimate is noisy, especially at the 3rd-degree level.
- **A homogeneous population** for the within-family form. Population structure inflates φ̂ for pairs from the same sub-group.
- **Accurate genotypes.** A genotyping error turns a parent–child pair's "impossible" IBS0 site into an observed one. That is why the father–daughter pair has 7 IBS0 sites rather than 0.
