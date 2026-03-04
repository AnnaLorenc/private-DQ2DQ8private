#!/usr/bin/env python3
"""Script: imgt_wide_to_properties.py

Summary
-------
Starting from a wide-format IMGT amino-acid count/frequency table (rows =
index-column combinations, columns = samples), computes a weighted mean of
each physicochemical property (e.g. VHSE or Kidera factors) per
(non-AA index columns) × sample.

For each group defined by the non-AA index columns and a single sample column,
the weighted mean of property P is:

    P_sample = Σ_aa ( count(aa) * P(aa) ) / Σ_aa count(aa)

where the sum runs over all amino-acid rows in that group.

Input
-----
Wide TSV such as:
    results/imgt_aa_combined/comb_aa_imgt_subs_rows.tsv

Expected structure:
    - Index columns (e.g. IMGT_position, AA, aminoAcid_length) identified via
      --index_columns.  One of them must be the amino-acid column (--aa_col).
    - Remaining columns are sample columns containing numeric counts.
IMGT_position   AA      aminoAcid_length        F20790E_subs_rows       F34137N_subs_rows       F15625E_subs_rows
104     C       8       3.0     2.0     3.0     3.0     4.0     2.0     2.0     2.0 


Property file
-------------
Same format as accepted by aa_properties_aggregate.py:
    - Case 1 (AA-as-rows): first column "AA", then one column per property.
    - Case 2 (transposed): first column "prop"/"property", remaining columns
      are single-letter amino acids.
See sandbox/VHS.tsv for an example of Case 2.

Output
------
Long-format TSV with one row per (group_cols × sample):
    <non-AA index columns>  sample  <property columns>

Usage
-----
    python imgt_wide_to_properties.py \\
        --input  results/imgt_aa_combined/comb_aa_imgt_subs_rows.tsv \\
        --property_file sandbox/VHS.tsv \\
        --output sandbox/properties_wide_out.tsv \\
        [--index_columns IMGT_position AA aminoAcid_length] \\
        [--aa_col AA]
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Property file loader (shared logic with aa_properties_aggregate.py)
# ---------------------------------------------------------------------------

def load_property_df(path: str) -> pd.DataFrame:
    """Load property table from CSV/TSV.

    Accepted layouts:
    - AA-as-rows  : column 'AA' (index) + one column per property.
    - Transposed  : column 'prop' or 'property' (property names as rows),
                    remaining columns are single-letter amino acids.

    Returns a DataFrame indexed by uppercase single-letter AA codes with
    property columns.
    """
    raw = pd.read_csv(path, sep=None, engine="python")
    raw.columns = [str(c).strip() for c in raw.columns]

    # Case 1: AA column present
    if "AA" in raw.columns:
        prop_cols = [c for c in raw.columns if c != "AA"]
        kf = raw.set_index("AA")[prop_cols]
        kf.index = [str(x).strip().upper() for x in kf.index]
        kf = kf.apply(pd.to_numeric, errors="coerce")
        return kf.loc[:, kf.notna().any(axis=0)]

    # Case 2: transposed with 'prop' / 'property'
    for candidate in ("prop", "property"):
        if candidate in raw.columns:
            kf_t = raw.set_index(candidate).T
            kf_t.index = [str(x).strip().upper() for x in kf_t.index]
            kf_t.columns = [str(c).strip() for c in kf_t.columns]
            kf_t = kf_t.apply(pd.to_numeric, errors="coerce")
            kf_t = kf_t.loc[:, kf_t.notna().any(axis=0)]
            if kf_t.shape[1] == 0:
                raise ValueError(f"No numeric property columns found in {path}")
            return kf_t

    raise ValueError(
        f"{path} must have an 'AA' column or a 'prop'/'property' column (transposed)"
    )


# ---------------------------------------------------------------------------
# Core aggregation
# ---------------------------------------------------------------------------

def aggregate_properties(
    df: pd.DataFrame,
    index_columns: list[str],
    aa_col: str,
    prop_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compute weighted-mean property vectors per (group × sample).

    Parameters
    ----------
    df            : wide-format DataFrame (index cols + sample cols)
    index_columns : list of index column names including aa_col
    aa_col        : name of the AA column within index_columns
    prop_df       : DataFrame indexed by AA letter, columns = properties

    Returns
    -------
    Long-format DataFrame: one row per (non-AA group cols × sample),
    columns = non-AA index cols + 'sample' + property cols.
    """
    # Normalise AA
    df = df.copy()
    df[aa_col] = df[aa_col].astype(str).str.strip().str.upper()

    # Drop AAs not in property table
    unknown = set(df[aa_col].unique()) - set(prop_df.index)
    if unknown:
        n = df[aa_col].isin(unknown).sum()
        print(f"  Dropping {n} rows with unknown AAs: {sorted(unknown)}")
        df = df[~df[aa_col].isin(unknown)].copy()

    sample_cols = [c for c in df.columns if c not in index_columns]
    if not sample_cols:
        raise ValueError("No sample columns found (all columns are index columns)")

    prop_cols = list(prop_df.columns)
    group_cols = [c for c in index_columns if c != aa_col]  # non-AA index cols

    # Property values for each row (aligned to df order)
    prop_values = prop_df.loc[df[aa_col].values].values  # (n_rows, n_props)

    rows = []
    for sample in sample_cols:
        counts = pd.to_numeric(df[sample], errors="coerce").fillna(0).values
        # weight × property for each row
        weighted = prop_values * counts[:, np.newaxis]   # (n_rows, n_props)

        # attach to a temp frame to group
        tmp = df[group_cols].copy()
        for i, p in enumerate(prop_cols):
            tmp[p] = weighted[:, i]
        tmp["_w"] = counts

        agg = tmp.groupby(group_cols, as_index=False).sum()

        # divide by total weight to get weighted mean
        w = agg["_w"].values
        for p in prop_cols:
            agg[p] = np.where(w > 0, agg[p].values / w, np.nan)

        agg = agg.drop(columns=["_w"])
        agg.insert(len(group_cols), "sample", sample)
        rows.append(agg)

    result = pd.concat(rows, ignore_index=True)
    # reorder: group_cols | sample | prop_cols
    return result[group_cols + ["sample"] + prop_cols]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute weighted-mean physicochemical property vectors from a "
            "wide IMGT AA count table."
        )
    )
    parser.add_argument("--input", required=True,
                        help="Wide-format IMGT AA TSV (rows = AA × position, "
                             "columns = sample counts).")
    parser.add_argument("--property_file", required=True,
                        help="Property table (Kidera/VHSE etc.) in AA-rows or "
                             "transposed format.")
    parser.add_argument("--output", required=True,
                        help="Output TSV path.")
    parser.add_argument(
        "--index_columns", nargs="+",
        default=["IMGT_position", "AA", "aminoAcid_length"],
        help="Index columns in the input (must include --aa_col). "
             "Default: IMGT_position AA aminoAcid_length.",
    )
    parser.add_argument(
        "--aa_col", default="AA",
        help="Name of the amino-acid column. Default: AA.",
    )
    args = parser.parse_args()

    print(f"Reading input: {args.input}")
    df = pd.read_csv(args.input, sep="\t", dtype={"IMGT_position": str})

    missing = [c for c in args.index_columns if c not in df.columns]
    if missing:
        raise ValueError(f"Index columns not found in input: {missing}")
    if args.aa_col not in args.index_columns:
        raise ValueError(f"--aa_col '{args.aa_col}' must be in --index_columns")

    print(f"Loading property file: {args.property_file}")
    prop_df = load_property_df(args.property_file)
    print(f"  Properties: {list(prop_df.columns)}")

    n_samples = len([c for c in df.columns if c not in args.index_columns])
    print(f"  Input: {len(df)} rows × {n_samples} sample columns")

    result = aggregate_properties(
        df=df,
        index_columns=args.index_columns,
        aa_col=args.aa_col,
        prop_df=prop_df,
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Written {len(result)} rows → {args.output}")


if __name__ == "__main__":
    main()
