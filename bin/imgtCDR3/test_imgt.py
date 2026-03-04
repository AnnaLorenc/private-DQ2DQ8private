#!/usr/bin/env python3
"""Test IMGT AA usage using simple linear-model-style comparisons.

Input
-----
A tab-separated file with at least the following columns (as produced by
bin/prepare_for_imgt_test.py):

    sample  length  AA  IMGT_position   patient cells   value   genotype

Example row:
    F16018E_subs_rows_freq  8   C   104 F16018  E   1   hoDQ2

The script fits simple two-group comparisons for each combination of
    (cells + index columns)
under two grouping schemes. By default, the index columns are
    (length, IMGT_position, AA).

  1. grouping1 "homhom": hoDQ2 versus hoDQ8
  2. grouping2 "hethom": heDQ2DQ8 versus (any of hoDQ2, hoDQ8)

For each group and each grouping, it computes:
  - n per group,
  - mean per group,
  - pooled standard deviation,
  - Cohen's d,
  - t statistic and (where SciPy is available) a two-sided p-value.

In addition, it performs a paired test between naive (N) and experienced (E)
cells within each individual, for every (length, IMGT_position, AA) index
combination, run:
    - across all individuals, and
    - separately within each genotype group hoDQ2, hoDQ8, heDQ2DQ8.

Outputs
-------
Two TSV files are written to --output_dir (default: results_man/imgt_test):
    - imgt_lm_combined.tsv : homhom ("_oo") and hethom ("_eo") results joined
        side by side by index columns (cells, length, IMGT_position, AA).
    - imgt_lm_combined_cells.tsv : paired N vs E tests per
        (length, IMGT_position, AA), with statistics for all individuals and for
        each genotype group in separate, suffixed columns.

Each row of imgt_lm_combined corresponds to one non-empty
(cells, length, IMGT_position, AA) combination with sufficient observations in
at least one grouping.

Each row of imgt_lm_combined_cells corresponds to one non-empty
(length, IMGT_position, AA) combination with sufficient paired N/E data in at
least one genotype stratum.
"""

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

try:  # Optional; used only for p-values
    from scipy import stats
except ImportError:  # pragma: no cover - optional dependency
    stats = None


GroupingKey = Tuple[str, int, int, str]


def compute_two_group_stats(
    values: pd.Series,
    genotypes: pd.Series,
    samples: pd.Series,
    group_a: List[str],
    group_b: List[str],
    label_a: str,
    label_b: str,
) -> Dict[str, float] | None:
    """Compute summary stats, Cohen's d, t and p-value for two groups.
    

    The input may contain multiple rows per sample; values are first
    aggregated to one value per (sample, group) by taking the mean, and
    n_group_* then reflects the number of unique samples in each group.

    Returns a dictionary with statistics or None if there are insufficient
    observations in either group.
    """

    df_local = pd.DataFrame({
        "sample": samples.astype(str),
        "value": values.astype(float),
        "genotype": genotypes.astype(str),
    })

    # Group A
    sub_a = df_local[df_local["genotype"].isin(group_a)][["sample", "value"]]
    vals_a = (
        sub_a.dropna(subset=["value"])
        .groupby("sample")["value"]
        .mean()
    )

    # Group B
    sub_b = df_local[df_local["genotype"].isin(group_b)][["sample", "value"]]
    vals_b = (
        sub_b.dropna(subset=["value"])
        .groupby("sample")["value"]
        .mean()
    )

    n_a = len(vals_a)
    n_b = len(vals_b)

    if n_a < 2 or n_b < 2:
        return None

    mean_a = float(vals_a.mean())
    mean_b = float(vals_b.mean())

    var_a = float(vals_a.var(ddof=1))
    var_b = float(vals_b.var(ddof=1))

    # Pooled standard deviation
    df = n_a + n_b - 2
    if df <= 0:
        return None

    pooled_var = ((n_a - 1) * var_a + (n_b - 1) * var_b) / df
    if pooled_var <= 0:
        return None
    pooled_sd = math.sqrt(pooled_var)

    diff = mean_a - mean_b
    cohen_d = diff / pooled_sd if pooled_sd > 0 else np.nan

    # t-statistic for difference in means
    se_diff = pooled_sd * math.sqrt(1.0 / n_a + 1.0 / n_b)
    if se_diff <= 0:
        t_stat = np.nan
        p_value = np.nan
    else:
        t_stat = diff / se_diff
        if stats is not None:
            # Two-sided p-value from t distribution
            p_value = float(stats.t.sf(abs(t_stat), df) * 2.0)
        else:
            p_value = np.nan

    return {
        "n_group_a": n_a,
        "n_group_b": n_b,
        "mean_group_a": mean_a,
        "mean_group_b": mean_b,
        "cohens_d": cohen_d,
        "t_stat": float(t_stat),
        "df": float(df),
        "p_value": p_value,
        "label_group_a": label_a,
        "label_group_b": label_b,
    }


def compute_paired_cells_stats(
    df: pd.DataFrame,
    index_cols: List[str],
    genotype_filter: List[str] | None,
) -> pd.DataFrame:
    """Paired N vs E test per index combination.

    For each combination in index_cols (typically [length, IMGT_position, AA]),
    within the specified genotype_filter (or all genotypes if None), this:

      - aggregates values to one mean per (patient, cells),
      - keeps only patients with both N and E measurements,
      - computes paired differences E - N and tests them with a one-sample t
        test against zero.

    Returns a DataFrame with one row per index combination and columns:
      n_pairs, mean_N, mean_E, mean_diff, sd_diff, t_stat, df, p_value.
    """

    if genotype_filter is not None:
        df_use = df[df["genotype"].isin(genotype_filter)].copy()
    else:
        df_use = df.copy()

    # Restrict to N/E cell types
    df_use = df_use[df_use["cells"].isin(["N", "E"])]
    if df_use.empty:
        return pd.DataFrame(columns=index_cols)

    results: List[Dict[str, float]] = []

    for key, sub in df_use.groupby(index_cols):
        # key is a tuple of index values, in the order of index_cols
        # Aggregate to one value per (patient, cells)
        g = (
            sub.groupby(["patient", "cells"], as_index=False)["value"]
            .mean()
        )

        # Pivot to columns N and E
        wide = g.pivot(index="patient", columns="cells", values="value")

        if "N" not in wide.columns or "E" not in wide.columns:
            continue

        paired = wide.dropna(subset=["N", "E"])
        if len(paired) < 2:
            continue

        vals_N = paired["N"].astype(float)
        vals_E = paired["E"].astype(float)
        diff = vals_E - vals_N

        n_pairs = int(len(diff))
        mean_N = float(vals_N.mean())
        mean_E = float(vals_E.mean())
        mean_diff = float(diff.mean())
        sd_diff = float(diff.std(ddof=1)) if n_pairs > 1 else float("nan")

        if n_pairs > 1 and sd_diff > 0:
            df_deg = n_pairs - 1
            se = sd_diff / math.sqrt(n_pairs)
            t_stat = mean_diff / se if se > 0 else float("nan")
            if stats is not None:
                p_value = float(stats.t.sf(abs(t_stat), df_deg) * 2.0)
            else:
                p_value = float("nan")
        else:
            df_deg = float("nan")
            t_stat = float("nan")
            p_value = float("nan")

        row: Dict[str, float] = {}
        for name, val in zip(index_cols, key):
            row[name] = val
        row.update(
            {
                "n_pairs": n_pairs,
                "mean_N": mean_N,
                "mean_E": mean_E,
                "mean_diff": mean_diff,
                "sd_diff": sd_diff,
                "t_stat": float(t_stat),
                "df": float(df_deg) if not math.isnan(df_deg) else float("nan"),
                "p_value": float(p_value),
            }
        )

        results.append(row)

    if not results:
        return pd.DataFrame(columns=index_cols)

    return pd.DataFrame(results)


def run_tests(input_path: Path, output_dir: Path, index_columns: List[str]) -> None:
    print(f"Reading input: {input_path}")
    df = pd.read_csv(input_path, sep="\t")

    # Required columns: configurable index columns plus these fixed ones
    required_cols = ["sample", "cells", "value", "genotype"] + index_columns
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {missing}")

    # Ensure integer-ish length for grouping if present; IMGT_position can
    # contain non-numeric labels (e.g. "111A"), so keep it as-is.
    if "length" in df.columns:
        df["length"] = df["length"].astype(int)

    # --------------------------------------------------
    # Print number of samples per genotype group for each grouping
    # --------------------------------------------------
    unique_samples = df[["sample", "genotype"]].drop_duplicates()

    # Grouping1: homhom (hoDQ2 vs hoDQ8)
    n_homhom_a = (
        unique_samples.loc[unique_samples["genotype"].isin(["hoDQ2"]), "sample"]
        .nunique()
    )
    n_homhom_b = (
        unique_samples.loc[unique_samples["genotype"].isin(["hoDQ8"]), "sample"]
        .nunique()
    )
    print("homhom grouping (hoDQ2 vs hoDQ8):")
    print("  group_a (hoDQ2) n_samples:", int(n_homhom_a))
    print("  group_b (hoDQ8) n_samples:", int(n_homhom_b))

    # Grouping2: hethom (heDQ2DQ8 vs hom = hoDQ2 or hoDQ8)
    n_hethom_a = (
        unique_samples.loc[unique_samples["genotype"].isin(["heDQ2DQ8"]), "sample"]
        .nunique()
    )
    n_hethom_b = (
        unique_samples.loc[unique_samples["genotype"].isin(["hoDQ2", "hoDQ8"]), "sample"]
        .nunique()
    )
    print("hethom grouping (heDQ2DQ8 vs hom):")
    print("  group_a (heDQ2DQ8) n_samples:", int(n_hethom_a))
    print("  group_b (hom) n_samples:", int(n_hethom_b))

    # Ensure requested index columns are present
    missing_idx = [c for c in index_columns if c not in df.columns]
    if missing_idx:
        raise ValueError(
            f"Index columns {missing_idx} not found in input file {input_path}"
        )

    # Grouping keys for genotype comparison tests: cells + index columns
    group_cols = ["cells"] + index_columns

    results_homhom: List[Dict[str, float]] = []
    results_hethom: List[Dict[str, float]] = []

    print("Running per-position tests for each combination of cells + index columns")

    for key, sub in df.groupby(group_cols):
        # key is (cells, *index_columns)
        cells = key[0]
        index_vals = key[1:]

        values = sub["value"]
        genotypes = sub["genotype"].astype(str)
        samples = sub["sample"].astype(str)

        # grouping1: hoDQ2 vs hoDQ8
        stats_homhom = compute_two_group_stats(
            values=values,
            genotypes=genotypes,
            samples=samples,
            group_a=["hoDQ2"],
            group_b=["hoDQ8"],
            label_a="hoDQ2",
            label_b="hoDQ8",
        )
        if stats_homhom is not None:
            row = {"cells": cells}
            row.update({name: val for name, val in zip(index_columns, index_vals)})
            row.update(stats_homhom)
            results_homhom.append(row)

        # grouping2: heDQ2DQ8 vs (hoDQ2, hoDQ8)
        stats_hethom = compute_two_group_stats(
            values=values,
            genotypes=genotypes,
            samples=samples,
            group_a=["heDQ2DQ8"],
            group_b=["hoDQ2", "hoDQ8"],
            label_a="heDQ2DQ8",
            label_b="hom",
        )
        if stats_hethom is not None:
            row = {"cells": cells}
            row.update({name: val for name, val in zip(index_columns, index_vals)})
            row.update(stats_hethom)
            results_hethom.append(row)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build dataframes and suffix statistic columns with grouping tags
    df_homhom = pd.DataFrame(results_homhom) if results_homhom else pd.DataFrame()
    df_hethom = pd.DataFrame(results_hethom) if results_hethom else pd.DataFrame()

    # Index columns for merging genotype-comparison results: cells + index_columns
    index_cols = ["cells"] + index_columns

    if not df_homhom.empty:
        stats_cols_homhom = [c for c in df_homhom.columns if c not in index_cols]
        df_homhom = df_homhom.rename(
            columns={c: f"{c}_oo" for c in stats_cols_homhom}
        )

    if not df_hethom.empty:
        stats_cols_hethom = [c for c in df_hethom.columns if c not in index_cols]
        df_hethom = df_hethom.rename(
            columns={c: f"{c}_eo" for c in stats_cols_hethom}
        )

    if df_homhom.empty and df_hethom.empty:
        print("No valid homhom or hethom comparisons produced any results.")
        return

    if df_homhom.empty:
        df_combined = df_hethom.copy()
    elif df_hethom.empty:
        df_combined = df_homhom.copy()
    else:
        df_combined = pd.merge(
            df_homhom,
            df_hethom,
            on=index_cols,
            how="outer",
            sort=True,
        )

    out_path = output_dir / "imgt_lm_combined.tsv"
    print(f"Writing combined results to: {out_path}")
    df_combined.to_csv(out_path, sep="\t", index=False)

    # --------------------------------------------------
    # Paired N vs E tests per combination of index columns
    # --------------------------------------------------
    cells_index_cols = index_columns

    # All individuals together
    df_cells_all = compute_paired_cells_stats(
        df=df,
        index_cols=cells_index_cols,
        genotype_filter=None,
    )

    # Genotype-specific strata
    df_cells_hoDQ2 = compute_paired_cells_stats(
        df=df,
        index_cols=cells_index_cols,
        genotype_filter=["hoDQ2"],
    )
    df_cells_hoDQ8 = compute_paired_cells_stats(
        df=df,
        index_cols=cells_index_cols,
        genotype_filter=["hoDQ8"],
    )
    df_cells_heDQ2DQ8 = compute_paired_cells_stats(
        df=df,
        index_cols=cells_index_cols,
        genotype_filter=["heDQ2DQ8"],
    )

    # If all strata are empty, nothing to write
    if (
        df_cells_all.empty
        and df_cells_hoDQ2.empty
        and df_cells_hoDQ8.empty
        and df_cells_heDQ2DQ8.empty
    ):
        print("No valid paired N/E cell comparisons produced any results.")
        return

    # Suffix statistic columns by stratum
    def _suffix_stats(df_stratum: pd.DataFrame, suffix: str) -> pd.DataFrame:
        if df_stratum.empty:
            return df_stratum
        stat_cols = [c for c in df_stratum.columns if c not in cells_index_cols]
        return df_stratum.rename(columns={c: f"{c}{suffix}" for c in stat_cols})

    df_cells_all = _suffix_stats(df_cells_all, "_all")
    df_cells_hoDQ2 = _suffix_stats(df_cells_hoDQ2, "_hoDQ2")
    df_cells_hoDQ8 = _suffix_stats(df_cells_hoDQ8, "_hoDQ8")
    df_cells_heDQ2DQ8 = _suffix_stats(df_cells_heDQ2DQ8, "_heDQ2DQ8")

    # Merge all available strata on the index columns
    df_cells_combined = None
    for part in [
        df_cells_all,
        df_cells_hoDQ2,
        df_cells_hoDQ8,
        df_cells_heDQ2DQ8,
    ]:
        if part is None or part.empty:
            continue
        if df_cells_combined is None:
            df_cells_combined = part.copy()
        else:
            df_cells_combined = pd.merge(
                df_cells_combined,
                part,
                on=cells_index_cols,
                how="outer",
                sort=True,
            )

    if df_cells_combined is not None and not df_cells_combined.empty:
        out_cells = output_dir / "imgt_lm_combined_cells.tsv"
        print(f"Writing paired N/E cell results to: {out_cells}")
        df_cells_combined.to_csv(out_cells, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit simple linear-model-style IMGT AA tests for specific "
            "genotype groupings across positions."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input TSV produced by bin/prepare_for_imgt_test.py.",
    )
    parser.add_argument(
        "--output_dir",
        default="results_man/imgt_test",
        help=(
            "Directory for output TSV files. "
            "Default: results_man/imgt_test."
        ),
    )

    parser.add_argument(
        "--index_columns",
        nargs="+",
        default=["length", "AA", "IMGT_position"],
        help=(
            "Index columns defining per-position tests (cells is implicit). "
            "Default: length AA IMGT_position."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)

    run_tests(
        input_path=input_path,
        output_dir=output_dir,
        index_columns=args.index_columns,
    )


if __name__ == "__main__":
    main()
