#!/usr/bin/env python3
"""Script: aa_kidera_aggregate.py

Summary
-------
Aggregates amino-acid IMGT frequency data into Kidera-factor weighted-mean
vectors per individual per position group.

For each combination of (index columns minus the AA column), and for each
individual (patient), computes the weighted mean of each of the 10 Kidera
factors using the per-row frequency (field "value") as weights.  This collapses
the amino-acid dimension while retaining per-individual, per-position information
needed for multivariate statistical testing.

Input
-----
A tab-separated test-ready file produced by bin/prepare_for_imgt_test.py, e.g.:
    results/imgt_aa_fortest/comb_aa_imgt_full_subs_rows_freq_5_1_test_ready.tsv

Expected columns (by default):
    sample  length  AA  IMGT_position  cells  patient  value  genotype

Output
------
Tab-separated file with one row per
    (aggregation columns) × patient × sample × genotype
and 10 additional columns KF1 … KF10 containing the Kidera-factor
weighted means.

Usage
-----
    python aa_kidera_aggregate.py \\
        --input  results/imgt_aa_fortest/comb_aa_imgt_full_subs_rows_freq_5_1_test_ready.tsv \\
        --output sandbox/kidera_aggregated.tsv \\
        [--index_columns length AA IMGT_position cells] \\
        [--aa_col AA] \\
        [--value_col value]
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Kidera factor table
# Rows = amino acids, columns = KF1 … KF10
# Source: Kidera et al. (1985) Biochim Biophys Acta 828:357-364
# ---------------------------------------------------------------------------
# Kidera factor column names
KF_COLS = [f"KF{i}" for i in range(1, 11)]


def load_kidera_df(kidera_path: str) -> pd.DataFrame:
    """Load Kidera factor table from a CSV/TSV file.

    Expected format: first column named "AA" (single-letter amino-acid),
    followed by columns KF1..KF10 (names must match). Delimiter may be
    comma or tab; the reader will attempt to auto-detect.
    """
    # try to read with pandas, letting it infer the delimiter
    kf_raw = pd.read_csv(kidera_path, sep=None, engine="python")

    cols = [c.strip() for c in kf_raw.columns]
    kf_raw.columns = cols

    # Case 1: "AA" column present (AA rows, KF columns)
    if "AA" in kf_raw.columns:
        missing = [c for c in KF_COLS if c not in kf_raw.columns]
        if missing:
            raise ValueError(f"Kidera file {kidera_path} is missing columns: {missing}")
        kf = kf_raw.set_index("AA")[KF_COLS]
        return kf

    # Case 2: transposed format with a "prop" (or "property") column
    # where rows are KF1.. and columns are amino-acid single-letter names.
    prop_col = None
    for candidate in ("prop", "property"):
        if candidate in kf_raw.columns:
            prop_col = candidate
            break
    if prop_col is not None:
        kf_t = kf_raw.set_index(prop_col).T
        # Ensure index is stripped and uppercase (amino-acid letters)
        kf_t.index = [str(x).strip().upper() for x in kf_t.index]
        # Clean column names
        kf_t.columns = [str(c).strip() for c in kf_t.columns]
        # Coerce values to numeric where possible and drop all-NaN columns
        kf_t = kf_t.apply(pd.to_numeric, errors="coerce")
        kf_t = kf_t.loc[:, kf_t.notna().any(axis=0)]
        if kf_t.shape[1] == 0:
            raise ValueError(f"Transposed Kidera file {kidera_path} contains no numeric property columns")
        return kf_t

    # Unknown layout
    raise ValueError(
        f"Kidera file {kidera_path} must contain either an 'AA' column or a 'prop'/'property' column (transposed layout)"
    )


def aggregate_kidera(
    df: pd.DataFrame,
    index_columns: list[str],
    aa_col: str,
    value_col: str,
    kf_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compute Kidera-factor weighted means per individual per position group.

    For every combination of (index_columns minus aa_col) × patient × sample ×
    genotype, a weighted mean of each KF is computed:

        KF_weighted_mean_i = Σ(value_j * KF_i_j) / Σ(value_j)

    where j runs over all amino-acid rows belonging to that group.

    Parameters
    ----------
    df            : long-format input DataFrame
    index_columns : column names treated as index; must include aa_col
    aa_col        : name of the amino-acid column within index_columns
    value_col     : name of the frequency / weight column

    Returns
    -------
    DataFrame with one row per group and columns KF1 … KF10.
    """
    # kf_df should be a DataFrame indexed by amino-acid single-letter code
    # with columns KF1..KF10

    # Normalise AA column to uppercase and coerce weights to numeric
    df = df.copy()
    df[aa_col] = df[aa_col].astype(str).str.strip().str.upper()
    # Coerce weights
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    n_before = len(df)
    df = df[~df[value_col].isna()].copy()
    if len(df) < n_before:
        print(f"  Dropped {n_before - len(df)} rows with non-numeric '{value_col}' values")

    # Restrict to amino acids with known factor values
    unknown = set(df[aa_col].unique()) - set(kf_df.index)
    if unknown:
        n_drop = df[aa_col].isin(unknown).sum()
        print(f"  Dropping {n_drop} rows with unknown AAs: {unknown}")
        df = df[~df[aa_col].isin(unknown)].copy()

    # Map AA -> factor values (aligned to df index)
    kf_values = kf_df.loc[df[aa_col].values].values  # shape (n_rows, p)
    weights = df[value_col].values  # shape (n_rows,)

    # Factor column names come from the provided kf_df
    kf_cols = list(kf_df.columns)

    # Attach weighted factor columns
    df = df.copy()
    for i, col in enumerate(kf_cols):
        df[col] = kf_values[:, i] * weights

    # Columns defining each group
    meta_cols = ["patient", "sample", "genotype"]
    agg_cols = [c for c in index_columns if c != aa_col]
    group_cols = agg_cols + meta_cols

    # Sum weighted KFs and weights per group, then divide to get weighted mean
    sum_cols = {col: "sum" for col in kf_cols}
    sum_cols[value_col] = "sum"
    agg = df.groupby(group_cols, as_index=False).agg(sum_cols)

    for col in kf_cols:
        agg[col] = agg[col] / agg[value_col]

    return agg.drop(columns=[value_col])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate amino-acid IMGT frequencies to Kidera-factor weighted "
            "means per individual per position group."
        )
    )
    parser.add_argument("--input",  required=True,
                        help="Path to test-ready input TSV.")
    parser.add_argument("--output", required=True,
                        help="Path to output TSV.")
    parser.add_argument(
        "--kidera_file", required=True,
        help="Path to Kidera factor CSV/TSV with columns 'AA' and KF1..KF10",
    )
    parser.add_argument(
        "--index_columns",
        nargs="+",
        default=["length", "AA", "IMGT_position", "cells"],
        help=(
            "Index columns (must include --aa_col). "
            "Default: length AA IMGT_position cells."
        ),
    )
    parser.add_argument(
        "--aa_col",
        default="AA",
        help="Name of the amino-acid column inside index_columns. Default: AA.",
    )
    parser.add_argument(
        "--value_col",
        default="value",
        help="Name of the frequency / weight column. Default: value.",
    )
    args = parser.parse_args()

    print(f"Reading: {args.input}")
    df = pd.read_csv(args.input, sep="\t")

    print(f"Loading Kidera factors: {args.kidera_file}")
    kf_df = load_kidera_df(args.kidera_file)

    needed = list(set(args.index_columns + ["patient", "sample", "genotype", args.value_col]))
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing columns: {missing}")

    if args.aa_col not in args.index_columns:
        raise ValueError(
            f"--aa_col '{args.aa_col}' must be one of the --index_columns."
        )

    print(f"  {len(df)} input rows, "
          f"index columns: {args.index_columns}, "
          f"aggregating within: {[c for c in args.index_columns if c != args.aa_col]}")

    result = aggregate_kidera(
        df=df,
        index_columns=args.index_columns,
        aa_col=args.aa_col,
        value_col=args.value_col,
        kf_df=kf_df,
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Written {len(result)} rows  →  {args.output}")


if __name__ == "__main__":
    main()
