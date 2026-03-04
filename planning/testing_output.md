# IMGT testing script (bin/imgtCDR3/test_imgt.py)

## Goal

Implement a fast, CLI-style Python script that:
- takes long-format IMGT AA usage data per (sample, length, AA, IMGT_position, patient, cells),
- encodes HLA genotype into comparison groups,
- and, for each (cells, length, IMGT_position, AA) combination, performs two simple two-group comparisons capturing effect sizes (Cohen's d) and t-test statistics.

This corresponds to the R-style sketch from testing.md using `group_by(cells, length, IMGT_position, AA)` and then fitting `lm(value ~ geno, data = .x)` and computing Cohen's d with a chosen reference group.

## Input format and grouping

Expected input TSV columns (produced by the first script or equivalent):
- `sample`: sample ID / column name from the wide table (kept for traceability, not used in modelling).
- `length`: CDR3 amino-acid length (cast to integer for grouping).
- `AA`: amino acid.
- `IMGT_position`: IMGT position label; kept as string to allow values like `111A`.
- `patient`: short patient ID (e.g. F16018).
- `cells`: E / N etc. (used as the first grouping key, so AgXP/naive are tested separately).
- `value`: the outcome to test (e.g. frequency from the `_subs_rows_freq` table).
- `genotype`: pre-recoded genotype label (e.g. `hoDQ2`, `hoDQ8`, `heDQ2DQ8`).

The script groups by:
- `cells`,
- `length`,
- `IMGT_position`,
- `AA`.

Within each such group we then have multiple samples (rows) with potentially different genotypes.

## Group definitions

Two groupings are implemented, matching the description in testing.md:

1. **homhom**
   - Group A: `hoDQ2`
   - Group B: `hoDQ8`
   - Label recorded in output: `grouping = "homhom"`, plus `label_group_a = "hoDQ2"`, `label_group_b = "hoDQ8"`.

2. **hethom**
   - Group A: `heDQ2DQ8` (the heterozygous group).
   - Group B: all homozygous genotypes in this context: `hoDQ2` or `hoDQ8`.
   - Label recorded in output: `grouping = "hethom"`, plus `label_group_a = "heDQ2DQ8"`, `label_group_b = "hom"`.

Each combination is tested only when both groups have at least two observations inside that `(cells, length, IMGT_position, AA)` group; otherwise, the corresponding result is skipped (no row for that combination).

## Statistics implemented

For each eligible group of rows, we compute, using pandas / NumPy:

- `n_group_a`, `n_group_b`: sample size in each group.
- `mean_group_a`, `mean_group_b`: mean `value` in each group.
- `var_a`, `var_b`: sample variances, using `ddof = 1`.

From these we construct a standard pooled-variance t-test:

1. Degrees of freedom:
   - `df = n_a + n_b - 2`.
   - If `df <= 0`, we drop the comparison.

2. Pooled variance and SD:
   - `s_p^2 = ((n_a - 1) * var_a + (n_b - 1) * var_b) / df`.
   - `s_p = sqrt(s_p^2)`.
   - If `s_p^2 <= 0`, we drop the comparison.

3. Cohen's d (standardised effect size):
   - `d = (mean_a - mean_b) / s_p`.
   - This matches the usual two-sample Cohen's d using the pooled SD.

4. t-statistic for the difference in means:
   - Standard error of the difference: `SE = s_p * sqrt(1/n_a + 1/n_b)`.
   - If `SE <= 0`, we set `t_stat` and `p_value` to NaN.
   - Otherwise, `t_stat = (mean_a - mean_b) / SE`.

5. p-value:
   - If SciPy is available, we compute a two-sided p-value from the t-distribution:
     - `p = 2 * stats.t.sf(abs(t_stat), df)`.
   - If SciPy is not installed, we leave `p_value` as NaN.

The returned row thus contains both a linear-model-like t-statistic and a Cohen's d effect size per position.

## Implementation details

- **Library choice**: For this second script, I used `pandas` and `numpy` (and optionally `scipy`) rather than Polars, because the data volumes at this stage are smaller (already aggregated per AA/position) and pandas makes the statistical calculations and grouping straightforward and explicit.
- **Type handling**:
  - `length` is cast to `int` to ensure consistent grouping (e.g. avoid `8.0` vs `8` splitting groups).
  - `IMGT_position` is kept as a string to gracefully handle labels like `111A`.
- **Missing data**:
  - Any missing values in `value` are dropped within each group before computing statistics.
  - If either group has fewer than 2 non-missing values, we skip that comparison (no output row).
- **Output location**:
  - Both result tables are written to `results_man/imgt_test` (or a user-specified `--output_dir`).
  - Filenames:
    - `imgt_lm_homhom.tsv`
    - `imgt_lm_hethom.tsv`

## Sanity check run

To validate the implementation without relying on the first script, I created a small synthetic input file at `results_man/imgt_test/test_imgt_input.tsv` with:
- 6 rows at the same `(cells = E, length = 8, IMGT_position = 104, AA = C)`,
- 2 samples in each genotype group: `hoDQ2`, `hoDQ8`, `heDQ2DQ8`.

Running:

```bash
python bin/imgtCDR3/test_imgt.py \
  --input results_man/imgt_test/test_imgt_input.tsv \
  --output_dir results_man/imgt_test
```

produced:
- `results_man/imgt_test/imgt_lm_homhom.tsv`
- `results_man/imgt_test/imgt_lm_hethom.tsv`

with one row each for the `(E, 8, 104, C)` combination, confirming that:
- the column expectations are correct,
- grouping and genotype recoding behave as intended,
- and the output files are written in the expected location and format.
