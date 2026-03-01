#!/usr/bin/env python3
"""Script: aa_properties_test.py

Summary
-------
Performs multivariate (Hotelling's T²) and univariate (per-property Welch
t-test) comparisons of property profiles across genotype groups, for each
combination of aggregation columns (e.g. IMGT_position × cells, or
IMGT_position × length × cells).

The input is the output of imgt_wide_to_properties.py--> properties_add_annotation.py: one row per
(agg_cols × patient × sample × genotype), with property columns specified by --prop_columns.

Two comparisons are performed:
    - oo  :  hoDQ2  vs  hoDQ8
    - eo  :  heDQ2DQ8  vs  (hoDQ2 + hoDQ8)

Both comparisons are run within each cell type (implicit if "cells" is among
the aggregation columns).

Output files
------------
1. <output_dir>/hotelling.tsv
     One row per aggregation-column combination.
     Columns (side-by-side _oo and _eo):
         T2, F, df1, df2, p_value, q_value (BH-FDR), eta_sq, n_a, n_b

2. <output_dir>/per_property.tsv
     One row per (aggregation-column combination × property).
     Columns (side-by-side _oo and _eo):
         mean_a, mean_b, mean_diff, cohens_d, t_stat, df, p_value, q_value (BH-FDR)
     Plus: p_hotelling_oo, q_hotelling_oo, p_hotelling_eo, q_hotelling_eo

Usage
-----
        python aa_properties_test.py \
                --input  sandbox/aggregated.tsv \
                --output_dir sandbox/results \
                [--agg_columns IMGT_position cells] \
                --prop_columns <property1> <property2> ... \
                --propname <property>
"""

import argparse
import math
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.stats import f as f_dist

try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("Warning: statsmodels not found; q-values will be NaN.")


DEFAULT_AGG_COLUMNS = ["IMGT_position", "cells"]



# ---------------------------------------------------------------------------
# Hotelling's T² for two independent groups
# ---------------------------------------------------------------------------

def hotelling_t2(
    X: np.ndarray,
    Y: np.ndarray,
) -> Optional[dict]:
    """Compute Hotelling's T² for two independent multivariate groups.

    Parameters
    ----------
    X : (n_a, p) array for group A
    Y : (n_b, p) array for group B

    Returns
    -------
    dict with keys: T2, F, df1, df2, p_value, eta_sq, n_a, n_b
    or None if the test cannot be computed (insufficient data or singular matrix).

    Notes
    -----
    T² = (n_a·n_b)/(n_a+n_b) · (x̄_a − x̄_b)ᵀ S_pooled⁻¹ (x̄_a − x̄_b)
    F  = (n_a+n_b−p−1) / [p·(n_a+n_b−2)] · T²
    F  ~ F(p, n_a+n_b−p−1)
    η² = T² / (T² + n_a+n_b−2)   [partial eta-squared analog]
    """
    n_a, p = X.shape
    n_b     = Y.shape[0]

    if n_a < 2 or n_b < 2:
        return None
    if n_a + n_b - 2 < p:
        # Degrees of freedom would be non-positive; pooled cov is rank-deficient
        return None

    mean_diff = X.mean(axis=0) - Y.mean(axis=0)

    # Pooled covariance
    S_a = np.cov(X, rowvar=False, ddof=1) if n_a > 1 else np.zeros((p, p))
    S_b = np.cov(Y, rowvar=False, ddof=1) if n_b > 1 else np.zeros((p, p))
    S_pooled = ((n_a - 1) * S_a + (n_b - 1) * S_b) / (n_a + n_b - 2)

    try:
        S_inv = np.linalg.inv(S_pooled)
    except np.linalg.LinAlgError:
        return None

    # Check condition number – if very large, matrix is near-singular
    if np.linalg.cond(S_pooled) > 1e12:
        return None

    T2 = (n_a * n_b / (n_a + n_b)) * float(mean_diff @ S_inv @ mean_diff)

    df1 = p
    df2 = n_a + n_b - p - 1
    if df2 <= 0:
        return None

    F = (n_a + n_b - p - 1) / (p * (n_a + n_b - 2)) * T2
    p_value = float(1.0 - f_dist.cdf(F, df1, df2))
    eta_sq  = T2 / (T2 + n_a + n_b - 2)

    return {
        "T2": T2, "F": F, "df1": df1, "df2": df2,
        "p_value": p_value, "eta_sq": eta_sq,
        "n_a": n_a, "n_b": n_b,
    }


# ---------------------------------------------------------------------------
# Per-factor Welch t-test with Cohen's d
# ---------------------------------------------------------------------------

def welch_t_cohens_d(
    a: np.ndarray,
    b: np.ndarray,
) -> Optional[dict]:
    """Welch two-sample t-test plus Cohen's d for two 1-D arrays.

    Cohen's d uses the pooled SD (Glass's variant would use SD of one group;
    we use pooled to be consistent with the paired two-group framing).

    Returns dict or None if either group has fewer than 2 observations.
    """
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return None

    # Ensure numeric precision and avoid catastrophic cancellation
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    mean_a  = float(a.mean())
    mean_b  = float(b.mean())
    diff    = mean_a - mean_b

    var_a   = float(a.var(ddof=1))
    var_b   = float(b.var(ddof=1))
    n_a, n_b = len(a), len(b)

    # Pooled SD for Cohen's d (use pooled sample variance)
    pooled_var = ((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2) if (n_a + n_b - 2) > 0 else float("nan")
    pooled_sd  = math.sqrt(pooled_var) if pooled_var > 0 else float("nan")
    cohens_d   = diff / pooled_sd if pooled_sd > 0 else float("nan")

    # Manual Welch t-statistic and df calculation to avoid scipy's moment-based
    # routines which can emit precision-loss warnings when data are nearly
    # identical. Short-circuit when variances are zero.
    se_a = var_a / n_a
    se_b = var_b / n_b
    se_sum = se_a + se_b
    if se_sum == 0 or (n_a - 1) <= 0 or (n_b - 1) <= 0:
        return None

    t_stat = diff / math.sqrt(se_sum)

    denom_df = (se_a ** 2) / (n_a - 1) + (se_b ** 2) / (n_b - 1)
    if denom_df <= 0:
        return None
    df_welch = float((se_sum ** 2) / denom_df)

    # two-sided p-value from Student's t survival function
    p_value = float(2.0 * scipy_stats.t.sf(abs(t_stat), df_welch))

    return {
        "mean_a": mean_a, "mean_b": mean_b, "mean_diff": diff,
        "cohens_d": cohens_d,
        "t_stat": float(t_stat), "df": df_welch,
        "p_value": float(p_value),
    }


# ---------------------------------------------------------------------------
# FDR correction helper
# ---------------------------------------------------------------------------

def fdr_correct(p_values: pd.Series) -> pd.Series:
    """Apply Benjamini–Hochberg FDR correction to a Series of p-values.

    NaN p-values are preserved as NaN in the output.
    """
    if not HAS_STATSMODELS:
        return pd.Series([float("nan")] * len(p_values), index=p_values.index)

    mask = p_values.notna()
    q = pd.Series([float("nan")] * len(p_values), index=p_values.index)
    if mask.sum() == 0:
        return q

    _, q_vals, _, _ = multipletests(p_values[mask].values, method="fdr_bh")
    q[mask] = q_vals
    return q


# ---------------------------------------------------------------------------
# Main testing logic
# ---------------------------------------------------------------------------

def run_hotelling(
    df: pd.DataFrame,
    agg_cols: list[str],
    prop_cols: list[str],
    genotypes_a: list[str],
    genotypes_b: list[str],
) -> pd.DataFrame:
    """Run Hotelling's T² per agg_cols group for one comparison.

    Returns a DataFrame with one row per group and Hotelling stats columns.
    """
    rows = []
    for key, sub in df.groupby(agg_cols):
        key_dict = dict(zip(agg_cols, key if isinstance(key, tuple) else (key,)))

        X = sub[sub["genotype"].isin(genotypes_a)][prop_cols].dropna().values
        Y = sub[sub["genotype"].isin(genotypes_b)][prop_cols].dropna().values

        result = hotelling_t2(X, Y)
        row = {**key_dict, "n_a": len(X), "n_b": len(Y)}
        if result:
            result.pop("n_a", None)
            result.pop("n_b", None)
            row.update(result)
        else:
            row.update({k: float("nan") for k in
                        ["T2", "F", "df1", "df2", "p_value", "eta_sq"]})
        rows.append(row)

    out = pd.DataFrame(rows)
    out["q_value"] = fdr_correct(out["p_value"])
    return out


def run_per_factor(
    df: pd.DataFrame,
    agg_cols: list[str],
    prop_cols: list[str],
    genotypes_a: list[str],
    genotypes_b: list[str],
) -> pd.DataFrame:
    """Run per-property Welch t-tests per agg_cols group for one comparison.

    Returns a DataFrame with one row per (group × property).
    """
    rows = []
    for key, sub in df.groupby(agg_cols):
        key_dict = dict(zip(agg_cols, key if isinstance(key, tuple) else (key,)))

        sub_a = sub[sub["genotype"].isin(genotypes_a)]
        sub_b = sub[sub["genotype"].isin(genotypes_b)]

        for prop in prop_cols:
            a = sub_a[prop].dropna().values
            b = sub_b[prop].dropna().values
            result = welch_t_cohens_d(a, b)
            row = {**key_dict, "prop": prop, "n_a": len(a), "n_b": len(b)}
            if result:
                row.update(result)
            else:
                row.update({k: float("nan") for k in
                            ["mean_a", "mean_b", "mean_diff", "cohens_d",
                             "t_stat", "df", "p_value"]})
            rows.append(row)

    out = pd.DataFrame(rows)
    out["q_value"] = fdr_correct(out["p_value"])
    return out


def run_tests(
    df: pd.DataFrame,
    agg_cols: list[str],
    prop_cols: list[str],
    output_dir: Path,
    property: str,
    base_name: str,
) -> None:
    """Orchestrate oo and eo comparisons, write output files."""

    comparisons = {
        "oo": (["hoDQ2"],     ["hoDQ8"],           "hoDQ2", "hoDQ8"),
        "eo": (["heDQ2DQ8"],  ["hoDQ2", "hoDQ8"],  "heDQ2DQ8", "hom"),
    }

    # ---- Hotelling --------------------------------------------------------
    hotelling_frames = {}
    for tag, (ga, gb, _, _) in comparisons.items():
        hotelling_frames[tag] = run_hotelling(df, agg_cols, prop_cols, ga, gb)

    # Merge oo and eo side by side on agg_cols
    h_oo = hotelling_frames["oo"].rename(
        columns={c: f"{c}_oo" for c in hotelling_frames["oo"].columns if c not in agg_cols}
    )
    h_eo = hotelling_frames["eo"].rename(
        columns={c: f"{c}_eo" for c in hotelling_frames["eo"].columns if c not in agg_cols}
    )
    df_hotelling = pd.merge(h_oo, h_eo, on=agg_cols, how="outer").sort_values(agg_cols)

    out1 = output_dir / f"{base_name}_{property}_hotelling.tsv"
    df_hotelling.to_csv(out1, sep="\t", index=False)
    print(f"Written Hotelling results ({len(df_hotelling)} rows) → {out1}")

    # ---- Per-factor -------------------------------------------------------
    pf_frames = {}
    for tag, (ga, gb, _, _) in comparisons.items():
        pf_frames[tag] = run_per_factor(df, agg_cols, prop_cols, ga, gb)

    pf_oo = pf_frames["oo"].rename(
        columns={c: f"{c}_oo" for c in pf_frames["oo"].columns
                 if c not in agg_cols and c != "prop"}
    )
    pf_eo = pf_frames["eo"].rename(
        columns={c: f"{c}_eo" for c in pf_frames["eo"].columns
                 if c not in agg_cols and c != "prop"}
    )
    df_pf = pd.merge(pf_oo, pf_eo, on=agg_cols + ["prop"], how="outer")

    # Attach Hotelling p and q values (per group key, broadcast across KFs)
    hot_cols_oo = hotelling_frames["oo"].rename(
        columns={"p_value": "p_hotelling_oo", "q_value": "q_hotelling_oo"}
    )[agg_cols + ["p_hotelling_oo", "q_hotelling_oo"]]
    hot_cols_eo = hotelling_frames["eo"].rename(
        columns={"p_value": "p_hotelling_eo", "q_value": "q_hotelling_eo"}
    )[agg_cols + ["p_hotelling_eo", "q_hotelling_eo"]]

    df_pf = df_pf.merge(hot_cols_oo, on=agg_cols, how="left")
    df_pf = df_pf.merge(hot_cols_eo, on=agg_cols, how="left")
    df_pf = df_pf.sort_values(agg_cols + ["prop"])

    out2 = output_dir / f"{base_name}_{property}_per_factor.tsv"
    df_pf.to_csv(out2, sep="\t", index=False)
    print(f"Written per-factor results ({len(df_pf)} rows)  → {out2}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Multivariate (Hotelling's T²) and per-factor (Welch t-test) "
            "comparisons of Kidera-factor profiles across genotype groups."
        )
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to the aggregated Kidera TSV produced by aa_properties_aggregate.py.",
    )
    parser.add_argument(
        "--output_dir", default="properties_results",
        help="Directory for output TSV files. Default: <propname>_properties_results.",
    )
    parser.add_argument(
        "--agg_columns", nargs="+",
        default=DEFAULT_AGG_COLUMNS,
        help=(
            "Columns defining each test group (cells is implicit if present). "
            "Default: IMGT_position cells."
        ),
    )
    parser.add_argument(
        "--prop_columns", nargs="+", required=False,
        help="Property column names to test (e.g. physicochemical factors). If not given, all columns starting with --propname will be used.",
    )
    parser.add_argument(
        "--propname", required=True,
        help="Name of the property being tested (used for output file naming). Compulsory."
    )


    args = parser.parse_args()

    print(f"Reading: {args.input}")
    df = pd.read_csv(args.input, sep="\t")

    # Determine property columns
    if args.prop_columns is not None and len(args.prop_columns) > 0:
        prop_columns = args.prop_columns
    else:
        prop_columns = [c for c in df.columns if c.startswith(args.propname)]
        if not prop_columns:
            raise ValueError(f"No columns found starting with '{args.propname}'. Please specify --prop_columns explicitly.")

    needed = args.agg_columns + prop_columns + ["genotype"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing columns: {missing}")

    # Drop rows with any missing property values before testing
    n_before = len(df)
    df = df.dropna(subset=prop_columns)
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"Dropped {n_dropped} rows with missing property values ({len(df)} remaining)")

    od_arg = args.output_dir
    # If the user did not provide --output_dir (left as default), use <propname>_properties_results
    if od_arg is None or od_arg == "properties_results" or str(od_arg).strip() == "":
        output_dir = Path(f"{args.propname}_properties_results")
    else:
        output_dir = Path(od_arg)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Aggregation columns : {args.agg_columns}")
    print(f"Property columns     : {prop_columns}")
    genotypes_present = df["genotype"].unique().tolist()
    print(f"Genotypes present   : {sorted(genotypes_present)}")

    base_name = Path(args.input).stem

    run_tests(
        df=df,
        agg_cols=args.agg_columns,
        prop_cols=prop_columns,
        output_dir=output_dir,
        property=args.propname,
        base_name=base_name,
    )


if __name__ == "__main__":
    main()
