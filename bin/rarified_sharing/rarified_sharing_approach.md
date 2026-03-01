# Rarefaction-Based Clonotype Publicity Scoring — Approach

## Goal

Determine which TCR clonotypes are **"public"** (shared across many individuals)
in a statistically robust way, correcting for variable sequencing depth and
sampling stochasticity.

---

## Input

A matrix of **all clonotypes × all samples**, where:
- Rows are clonotypes, uniquely identified by one or more **index columns**
  (e.g. `aminoAcid`, `vFamilyName`, `jGeneName`).
- Columns are sample-level integer counts (how many reads/templates were
  observed for that clonotype in that sample).

Default input: `results/merged/merged_aminoAcid_vFamilyName_jGeneName_productive.tsv.gz`

---

## Algorithm

### Step 1 — Subset

Select the user-specified samples **S** and filter to clonotypes with count > 0
in at least one sample in S. This defines the working set:
`clonotypes_non0_in_S` (rows) × `|S|` (columns).

### Step 2 — Define sharing groups

Groups G are named subsets of S, e.g.:
- `homo`: samples from HLA-homozygous individuals
- `hete`: samples from HLA-heterozygous individuals

Sharing is counted **within each group independently**, allowing comparison
of publicity patterns across biological categories.

Groups can be specified in two ways:

**A. On the command line** (suitable for a small number of groups/samples):
```bash
--groups homo:F28396E,F28396N,F20790N hete:T04918N,T06306E,T13372E
```
If `--samples` is omitted, S is automatically derived as the union of all
group samples.

**B. From a YAML file** (recommended for many samples or reproducible runs):
```bash
--groups_file groups_E.yml
```

The YAML format is:
```yaml
# Each top-level key is a group name; value is a list of sample names.
homo:
  - F28396E
  - F28396N
  - F20790N
hete:
  - T04918N
  - T06306E
  - T13372E
```

#### Preparing a groups YAML with `make_groups_yml.py`

`make_groups_yml.py` builds the YAML directly from a sample-sheet CSV/TSV by
specifying which column contains the group factor and which contains the sample
name. It also accepts input from **stdin**, so rows can be pre-filtered before
piping in.

```bash
# All samples, grouped by genotype:
python bin/rarified_sharing/make_groups_yml.py \
    --input data/collated_info.csv \
    --group_col genotype_short \
    --sample_col sample_short \
    --rename "heteroDQ2DQ8=hete,homoDQ8=homo" \
    --output groups_all.yml

# Only antigen-expanded (cells=E) samples, filtered via awk:
awk -F',' 'NR==1 || $3=="E"' data/collated_info.csv | \
    python bin/rarified_sharing/make_groups_yml.py \
        --group_col genotype_short \
        --sample_col sample_short \
        --rename "heteroDQ2DQ8=hete,homoDQ8=homo" \
        --output groups_E.yml
```

The `--rename` flag maps verbose factor levels to short group names used in
output file names. The delimiter is auto-detected (CSV, TSV, etc.).
Any tool that can filter rows and print to stdout (awk, csvkit `csvgrep`,
pandas, etc.) can be piped into `make_groups_yml.py`.

### Step 3 — Sparse precomputation (once, before K rounds)

Compute the **union of all group sample sets**: only columns that appear in at
least one group need to be subsampled. Samples in S but absent from every group
are skipped entirely.

For each column in `union(G)`, extract the **non-zero clonotype indices** and
their normalised probabilities. For real data (~6.6M clonotypes, ~330k non-zero
per sample), this reduces the working vector from 6.6M to ~330k entries,
giving a ~6× speedup per round. This structure is computed once and reused
across all K rounds.

### Step 4 — Subsampling (K rounds)

Repeat K times:

1. **For each sample s in union(G)**: draw **M** clonotypes with replacement
   using `np.random.multinomial(M, p)` where `p` is the sparse probability
   vector from Step 3. Convert to **binary presence/absence** (1 if drawn ≥
   once, 0 otherwise). This mirrors the logic in
   `bin/subsample_merged_sequences.py`.

2. **For each group G**: the sharing count for round k is the **sum of presence
   vectors** across the samples in G. For each clonotype, this number is between
   0 (not detected in any sample in G after subsampling) and |G| (detected in
   every sample in G).

### Step 5 — Output

For each group G, write a TSV matrix of shape **(n_clonotypes × K)**, where:
- Rows = clonotypes (index columns prepended)
- Columns `k1..kK` = sharing count in round k

**Filename**: `{output_dir}/{group_name}_sharing.tsv`

---

## Why This Approach

| Concern | How it is addressed |
|---|---|
| Sequencing depth bias | Each sample is subsampled to the same depth M before counting sharing |
| Stochastic variation | K rounds quantify the uncertainty in sharing counts |
| Simple downstream use | Output matrices can directly compute mean sharing, P(sharing ≥ N), etc. |
| Speed | Sparse precomputation on `union(G)` columns only (samples in S but not in any group are skipped); non-zero indices + probs per column computed once; multinomial draws only on the non-zero subset (~330k vs 6.6M entries); matrix sums for group counts |

---

## Parameters

| Parameter | Meaning | Suggested range |
|---|---|---|
| `M` | Subsampling depth per sample per round | ≥ 1 000; match to minimum sample size |
| `K` | Number of rounds | 25–200 for stable estimates |
| `S` | Samples to include | All samples of interest |
| `G` | Sharing groups | Any biologically meaningful subsets of S |



## Performance (benchmarked on real data)

Data: 6.6M clonotypes × 120 samples, ~330k non-zero entries per sample column.

| Approach | Per round | K=200, M=10 000 |
|---|---|---|
| Naïve (multinomial on all 6.6M rows) | ~12 s | ~40 min |
| **Sparse precomputation (implemented)** | **~2 s** | **~7 min** |

If further speed is needed, parallelising the per-sample multinomial draws
with `joblib` or `multiprocessing` across CPU cores would reduce this to
~1–2 min.

---

## Step 2 — Summarise sharing (`summarise_sharing.py`)

After `rare_publicity.py` produces the per-group raw sharing matrices
(`k1..kK`), `summarise_sharing.py` converts them into a compact, interpretable
summary and identifies public clonotypes.

### `summarise_sharing_group(df, index_cols)`

For each clonotype, counts **in how many of the K rounds the sharing value was
≥ N**, for every N from 2 to the maximum observed sharing in that file.

Example: if a clonotype was shared by 3, 5, 7, 0 individuals across 4 rounds,
the output row is:

| aminoAcid | … | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| CAAGGAGGTGELFF | … | 3 | 3 | 2 | 1 | 0 | 0 | … |

Values are counts of rounds, so they are always ≤ K and decrease monotonically
as N increases.

### Multi-group workflow

1. Run `summarise_sharing_group` on each input file independently.
2. **Identify public clonotypes** per group: rows where the N=A column ≥ 1
   (i.e., the clonotype was shared by ≥ A individuals in at least one round).
3. Take the **union** of public clonotypes across all groups.
4. **Subset** every group summary to that union (zero-fill columns that don't
   exist in a group due to lower max sharing).
5. Write one TSV per group (`{name}_summary.tsv`) and one **merged TSV**
   (`merged_summary.tsv`) with a `group` column appended.

### Output columns

The per-group and merged files contain columns **A .. global_max_N** only
(columns 2 .. A−1 are dropped as uninformative given the threshold A).

### Parameters

| Parameter | Meaning |
|---|---|
| `--A` | Minimum sharing threshold; must be ≥ 2 |
| `--input` | One or more `*_sharing.tsv` files from `rare_publicity.py` |
| `--names` | Optional group labels (defaults to file basenames) |
| `--index_cols` | Auto-detected if omitted |

---

## Files

Scripts have been moved to `bin/rarified_sharing/`.

| File | Description |
|---|---|
| `bin/rarified_sharing/rare_publicity.py` | Step 1: rarefaction subsampling → raw sharing matrices |
| `bin/rarified_sharing/summarise_sharing.py` | Step 2: summarise raw matrices → public clonotype tables |
| `bin/rarified_sharing/make_groups_yml.py` | Helper: build a groups YAML from a sample-sheet CSV/TSV |
| `bin/rarified_sharing/make_test_data.py` | Generates `test_data.tsv` (for rare_publicity.py) |
| `bin/rarified_sharing/make_test_sharing.py` | Generates `test_sharing/` (for summarise_sharing.py) |
| `sandbox/rare/test_data.tsv` | 30 clonotypes × 10 samples toy dataset |
| `sandbox/rare/test_sharing/` | Toy sharing matrices (groupA, groupB) |
| `sandbox/rare/test_out/` | Output from a rare_publicity.py test run |
| `sandbox/rare/approach.md` | This file |

---

## Example Runs

### Step 1 — rare_publicity.py

```bash
# Groups on the command line
python bin/rarified_sharing/rare_publicity.py \
    --input results/merged/merged_aminoAcid_vFamilyName_jGeneName_productive.tsv.gz \
    --K 200 --M 10000 \
    --groups homo:F28396E,F28396N hete:T04918N,T06306E \
    --index_cols aminoAcid,vFamilyName,jGeneName \
    --output results/rare_out --seed 42

# Groups from a YAML file (build it first with make_groups_yml.py)
awk -F',' 'NR==1 || $3=="E"' data/collated_info.csv | \
    python bin/rarified_sharing/make_groups_yml.py \
        --group_col genotype_short --sample_col sample_short \
        --rename "heteroDQ2DQ8=hete,homoDQ8=homo" \
        --output groups_E.yml

python bin/rarified_sharing/rare_publicity.py \
    --input results/merged/merged_aminoAcid_vFamilyName_jGeneName_productive.tsv.gz \
    --K 200 --M 10000 \
    --groups_file groups_E.yml \
    --index_cols aminoAcid,vFamilyName,jGeneName \
    --output results/rare_out --seed 42
```

### Step 2 — summarise_sharing.py

```bash
python bin/rarified_sharing/summarise_sharing.py \
    --input results/rare_out/DQ2_sharing.tsv \
            results/rare_out/DQ8_sharing.tsv \
            results/rare_out/hete_sharing.tsv \
    --A 3 \
    --output results/sharing_summary_A3
```

Output files:
- `results/sharing_summary_A3/DQ2_summary.tsv`
- `results/sharing_summary_A3/DQ8_summary.tsv`
- `results/sharing_summary_A3/hete_summary.tsv`
- `results/sharing_summary_A3/merged_summary.tsv`  ← all groups, long format

---

## Downstream Analysis

The raw k-columns from `rare_publicity.py` output can also be used directly:

```python
import polars as pl

df = pl.read_csv("results/rare_out/hete_sharing.tsv", separator="\t")
k_cols = [c for c in df.columns if c.startswith("k")]

# Mean sharing across rounds
df = df.with_columns(
    pl.concat_list(k_cols).list.mean().alias("mean_sharing")
)

# Probability of being shared in >= N samples
N = 3
df = df.with_columns(
    (pl.concat_list(k_cols).list.eval(pl.element() >= N).list.mean())
    .alias(f"P_shared_ge_{N}")
)
```
