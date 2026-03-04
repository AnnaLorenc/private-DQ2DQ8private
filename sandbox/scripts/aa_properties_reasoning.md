# Kidera factor analysis – analytical decisions

## Context

Amino-acid IMGT frequency profiles per individual and per TCR beta-chain
position are the starting point. Each individual, at each IMGT position (and
CDR3 length, where applicable), has a frequency distribution over 20 amino
acids. We want to ask: does the **physicochemical character** of the residue
usage differ across HLA genotype groups (hoDQ2, hoDQ8, heDQ2DQ8)?

Testing 20 amino acids × ~30 IMGT positions × multiple lengths × cell types
individually would inflate the multiple-testing burden severely and is hard to
interpret. Kidera factors collapse those 20 dimensions into 10 orthogonal
physicochemical axes, making the comparison both more interpretable and more
powerful.

---

## Script A – Aggregation (`aa_kidera_aggregate.py`)

### Weighted mean for Kidera aggregation

Each amino acid *a* at position *j* for individual *i* has an observed
frequency $f_{ija}$ (the `value` column). The Kidera factors $K_k(a)$ for
$k = 1 \ldots 10$ are fixed scalar properties of each amino acid.

The individual-level Kidera vector at position *j* is:

$$
\hat{K}_k^{(ij)} = \frac{\sum_a f_{ija} \cdot K_k(a)}{\sum_a f_{ija}}
$$

This is the **weighted mean** of factor *k*, with per-amino-acid frequency as
weight.

**Why weighted mean, not equal-weight mean?**

* Amino acids are not equally used; the frequency $f_{ija}$ already encodes
  the individual's usage pattern. Ignoring frequency would give the same result
  for an individual dominated by Gly (neutral, flexible) as for one with only
  a trace of Gly at that position.
* The weighted mean is interpretable as the *expected* Kidera value of a
  randomly drawn residue from that individual's repertoire at that position.
* Normalising by $\sum_a f_{ija}$ means the result is invariant to total
  sequencing depth, which varies across individuals and subsampling runs.

---

## Script B – Testing (`aa_kidera_test.py`)

### Hotelling's T² for joint multivariate test

Each individual contributes a 10-dimensional observation
$\mathbf{x}_i = (\hat{K}_1, \ldots, \hat{K}_{10})$.

The 10 Kidera factors were derived from PCA/PCoA and are approximately
orthogonal across the 20 amino acids, but within a cohort the per-individual
vectors $\mathbf{x}_i$ may still show inter-factor correlations. A
**univariate test applied independently to each factor** would:

1. ignore those correlations (losing power when factors co-vary), and
2. require additional multiple-testing correction across factors.

**Hotelling's T²** is the natural multivariate extension of the two-sample
t-test. It accounts for the covariance structure of the 10-factor vectors and
yields a single p-value per position, controlling the per-position false
positive rate at the intended level.

#### Formula

For groups A ($n_a$ individuals) and B ($n_b$ individuals):

$$
T^2 = \frac{n_a n_b}{n_a + n_b}
      (\bar{\mathbf{x}}_A - \bar{\mathbf{x}}_B)^\top
      \mathbf{S}_\text{pooled}^{-1}
      (\bar{\mathbf{x}}_A - \mathbf{\bar{x}}_B)
$$

$$
F = \frac{n_a + n_b - p - 1}{p(n_a + n_b - 2)} T^2
\;\sim\; F\!\left(p,\; n_a + n_b - p - 1\right)
$$

where $p = 10$ (number of Kidera factors) and $\mathbf{S}_\text{pooled}$ is
the pooled within-group covariance matrix.

#### Effect size – partial η²

$$
\eta^2 = \frac{T^2}{T^2 + n_a + n_b - 2}
$$

This is the multivariate analog of the proportion of variance explained. It
ranges from 0 to 1 and is directly comparable across positions.

#### Limitations and safeguards

* If $n_a + n_b - 2 < p$ (total df < number of factors), the pooled covariance
  is rank-deficient and the test is undefined. Such groups are reported as NaN.
* A condition-number check (threshold $10^{12}$) flags near-singular covariance
  matrices (e.g., when several factors are constant within the cohort).

---

### Post-hoc per-factor Welch t-tests

The Hotelling test answers "is there *any* difference in the Kidera profile?"
but does not identify *which* factors drive it. For each factor independently
we run a **Welch two-sample t-test** (unequal-variance assumption, more
robust for small unequal groups).

**Why Welch rather than Student's?**

Student's t-test assumes equal variances; cohort sizes and variances can differ
between genotype groups. Welch's correction to the degrees of freedom is
conservative but valid without this assumption.

**Effect size – Cohen's d**

$$
d = \frac{\bar{x}_A - \bar{x}_B}{s_\text{pooled}}, \quad
s_\text{pooled} = \sqrt{\frac{(n_a-1)s_A^2 + (n_b-1)s_B^2}{n_a+n_b-2}}
$$

$|d| \approx 0.2$ small, $0.5$ medium, $0.8$ large (Cohen 1988).

The per-factor p-values are **FDR-corrected jointly across all
(position × factor) combinations** (within each comparison), not
per-position. This is consistent with viewing the experiment as a screen over
all positions and all factors simultaneously.

---

### Comparisons

| Tag | Group A        | Group B               | Interpretation            |
|-----|----------------|-----------------------|---------------------------|
| oo  | hoDQ2          | hoDQ8                 | Homozygous DQ2 vs DQ8     |
| eo  | heDQ2DQ8       | hoDQ2 + hoDQ8         | Heterozygous vs homozygous|

Comparisons are performed **within cell type** (N or E cells). Because `cells`
is included among the aggregation columns in Script A, each unique (position ×
cells) combination is automatically tested within one cell type.

---

### FDR correction – Benjamini–Hochberg

Raw p-values from many simultaneous tests are corrected using the
Benjamini–Hochberg (BH) procedure:

* **Hotelling file**: FDR applied across all (position × cell-type) rows,
  separately for oo and eo.
* **Per-factor file**: FDR applied across all (position × cell-type × factor)
  rows, separately for oo and eo.

BH controls the **false discovery rate** (expected proportion of false
positives among rejected tests) at the nominal level. It is preferred over
Bonferroni correction here because we expect many tests to be truly null and
we are more concerned with the rate of false discoveries than with the
probability of any single false positive.

---

## Relationship between Script A output and Script B input

```
raw TSV (sample × AA × IMGT_position × cells × patient × genotype × value)
         │
         ▼  aa_kidera_aggregate.py
aggregated TSV (IMGT_position × cells × patient × genotype × KF1..KF10)
         │
         ▼  aa_kidera_test.py
kidera_hotelling.tsv    (one row per position × cells; T², F, p, q, η²)
kidera_per_factor.tsv   (one row per position × cells × KF; t, d, p, q)
```

The `length` column (CDR3 length) is included in the aggregation by default
when present. For the "WL" (whole-length collapsed) input it is absent, and
the `--index_columns` / `--agg_columns` arguments should be adjusted
accordingly.
