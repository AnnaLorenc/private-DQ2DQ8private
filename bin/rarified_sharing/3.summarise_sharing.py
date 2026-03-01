#!/usr/bin/env python3
"""
summarise_sharing.py
--------------------
Summarise rarefaction-based sharing matrices produced by rare_publicity.py.

Each input file has the structure:
    <index_cols...>   k1   k2   ...   kK
where each k-column holds an integer = number of individuals in a group that
shared a clonotype in that subsampling round.

For each input file this script computes summarise_sharing_group():
    A clonotype × N matrix where the value at (clonotype, N) is the count of
    subsampling rounds in which sharing was >= N.
    N runs from 2 to the maximum observed sharing value in that file.

The script then:
    1. Identifies "public" clonotypes in each group: those detected shared by
       >= A individuals in at least one round (i.e., the N=A column >= 1).
    2. Takes the UNION of public clonotypes across all groups.
    3. Subsets every group summary to that union.
    4. Writes one TSV per group and one merged TSV with a 'group' column.
       The merged file contains columns A..max_N (zero-filled where needed).

Usage
-----
python summarise_sharing.py \\
    --input  results/DQ2_sharing.tsv results/DQ8_sharing.tsv results/hete_sharing.tsv \\
    --A 3 \\
    --output results/summary \\
    [--index_cols aminoAcid,vFamilyName,jGeneName] \\
    [--names DQ2,DQ8,hete]          # optional; defaults to file basenames

Arguments
---------
--input        One or more sharing TSV files (output of rare_publicity.py).
--A            Minimum sharing threshold for a clonotype to be called "public".
               Clonotypes with sharing >= A in at least 1 round are retained.
--output       Output directory (created if absent).
--index_cols   Comma-separated index column names. Auto-detected if omitted
               (all columns whose names do NOT match k{integer} are treated as
               index columns).
--names        Comma-separated group labels matching the order of --input files.
               Defaults to the file basename (without _sharing.tsv / .tsv).
"""

import argparse
import os
import re
import sys

import numpy as np
import polars as pl


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def summarise_sharing_group(
    df: pl.DataFrame,
    index_cols: list[str],
) -> pl.DataFrame:
    """
    Compute, for each clonotype, how many subsampling rounds had sharing >= N
    for every N in 2..max_observed_sharing.

    Parameters
    ----------
    df : pl.DataFrame
        Input sharing file as loaded from TSV.  Must contain the index columns
        and at least one k-column (k1, k2, ...).
    index_cols : list of str
        Column names that identify a clonotype (not k-columns).

    Returns
    -------
    pl.DataFrame
        Columns: index_cols + [str(N) for N in 2..max_val].
        Each numeric column contains the count of rounds with sharing >= N.
        Returns None if max sharing value < 2 (no multi-individual sharing).
    """
    k_cols = _detect_k_cols(df, index_cols)
    if not k_cols:
        raise ValueError("No k-columns found in DataFrame.")

    counts = df.select(k_cols).to_numpy()          # shape (n_clonotypes, K)
    max_val = int(counts.max())

    if max_val < 2:
        print(
            "  WARNING: maximum sharing value is < 2; "
            "no multi-individual sharing observed.",
            file=sys.stderr,
        )
        return None

    N_range = range(2, max_val + 1)

    # Vectorised: for each N, sum boolean (counts >= N) across K axis
    summary_arrays = {
        str(N): (counts >= N).sum(axis=1).astype(np.int32)
        for N in N_range
    }

    result = df.select(index_cols)
    for col_name, arr in summary_arrays.items():
        result = result.with_columns(pl.Series(col_name, arr))

    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _detect_k_cols(df: pl.DataFrame, index_cols: list[str]) -> list[str]:
    """Return columns that look like k{integer} and are not in index_cols."""
    pattern = re.compile(r"^k\d+$")
    return [c for c in df.columns if pattern.match(c) and c not in index_cols]


def _auto_index_cols(df: pl.DataFrame) -> list[str]:
    """Return all columns that are NOT k{integer} columns."""
    pattern = re.compile(r"^k\d+$")
    return [c for c in df.columns if not pattern.match(c)]


def _load_file(path: str) -> pl.DataFrame:
    """Load a sharing TSV file (handles plain .tsv and .tsv.gz)."""
    print(f"  Reading: {path}")
    df = pl.read_csv(path, separator="\t")
    print(f"    {len(df):,} clonotypes")
    return df


def _group_name_from_path(path: str) -> str:
    """Derive a group label from a file path (strip _sharing.tsv suffix)."""
    base = os.path.basename(path)
    base = re.sub(r"_sharing\.tsv(\.gz)?$", "", base)
    base = re.sub(r"\.tsv(\.gz)?$", "", base)
    return base


def _public_clonotype_mask(
    summary: pl.DataFrame,
    index_cols: list[str],
    A: int,
) -> pl.Series:
    """
    Return a boolean Series: True where the clonotype has sharing >= A in
    at least 1 round (i.e., the column str(A) >= 1).
    """
    col_A = str(A)
    if col_A not in summary.columns:
        # No column for A means max sharing < A → no public clonotypes
        return pl.Series("public", [False] * len(summary))
    return summary[col_A] >= 1


def _clonotype_key(df: pl.DataFrame, index_cols: list[str]) -> pl.Series:
    """Concatenate index columns into a single string key for set operations."""
    sep = "\x00"  # null-byte separator unlikely to appear in sequences
    if len(index_cols) == 1:
        return df[index_cols[0]].cast(pl.String)
    key = df[index_cols[0]].cast(pl.String)
    for col in index_cols[1:]:
        key = key + sep + df[col].cast(pl.String)
    return key


def _subset_to_keys(
    df: pl.DataFrame,
    index_cols: list[str],
    keys: set[str],
) -> pl.DataFrame:
    """Filter df to rows whose clonotype key is in keys."""
    sep = "\x00"
    key_col = _clonotype_key(df, index_cols)
    mask = key_col.is_in(list(keys))
    return df.filter(mask)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarise rarefaction sharing matrices.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", nargs="+", required=True,
        help="One or more sharing TSV files (output of rare_publicity.py).",
    )
    parser.add_argument(
        "--A", type=int, required=True,
        help=(
            "Minimum sharing threshold. Clonotypes shared by >= A individuals "
            "in at least 1 round are kept."
        ),
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Output directory (created if absent).",
    )
    parser.add_argument(
        "--index_cols", default=None,
        help=(
            "Comma-separated index column names. "
            "Auto-detected if omitted (all non-k{int} columns)."
        ),
    )
    parser.add_argument(
        "--names", default=None,
        help=(
            "Comma-separated group labels matching the order of --input. "
            "Defaults to file basenames (stripping _sharing.tsv)."
        ),
    )

    args = parser.parse_args()

    # ---- Validate arguments ----
    if args.A < 2:
        sys.exit("ERROR: --A must be >= 2 (sharing by at least 2 individuals).")

    input_files = args.input

    if args.names:
        names = [n.strip() for n in args.names.split(",")]
        if len(names) != len(input_files):
            sys.exit(
                f"ERROR: --names has {len(names)} entries but "
                f"--input has {len(input_files)} files."
            )
    else:
        names = [_group_name_from_path(f) for f in input_files]

    for path in input_files:
        if not os.path.exists(path):
            sys.exit(f"ERROR: Input file not found: {path}")

    os.makedirs(args.output, exist_ok=True)

    # ---- Parse index cols (once, from the first file) ----
    if args.index_cols:
        index_cols = [c.strip() for c in args.index_cols.split(",") if c.strip()]
    else:
        # Auto-detect from first file header
        first_df = pl.read_csv(input_files[0], separator="\t", n_rows=1)
        index_cols = _auto_index_cols(first_df)
        print(f"Auto-detected index columns: {index_cols}")

    print(f"\n{'='*60}")
    print(f"summarise_sharing.py  |  A={args.A}")
    print(f"{'='*60}")
    print(f"  Input files : {len(input_files)}")
    print(f"  Groups      : {names}")
    print(f"  Index cols  : {index_cols}")
    print(f"  Output dir  : {args.output}")
    print()

    # ---- Step 1: Load + summarise each file ----
    summaries: dict[str, pl.DataFrame] = {}

    for path, name in zip(input_files, names):
        print(f"[{name}]")
        df = _load_file(path)

        # Validate index cols
        missing = [c for c in index_cols if c not in df.columns]
        if missing:
            sys.exit(
                f"ERROR: Index column(s) {missing} not found in {path}.\n"
                f"Available columns: {df.columns[:10]}..."
            )

        summary = summarise_sharing_group(df, index_cols)
        if summary is None:
            print(f"  Skipping {name}: no sharing >= 2 observed.")
            continue

        n_cols = len(summary.columns) - len(index_cols)
        print(f"  Summary columns (N=2..{2 + n_cols - 1}): {n_cols} levels")
        summaries[name] = summary

    if not summaries:
        sys.exit("ERROR: No valid summaries produced.")

    # ---- Step 2: Find public clonotypes per group ----
    print(f"\nIdentifying clonotypes with sharing >= A={args.A} in >= 1 round ...")
    public_keys_per_group: dict[str, set[str]] = {}

    for name, summary in summaries.items():
        mask = _public_clonotype_mask(summary, index_cols, args.A)
        keys = set(_clonotype_key(summary.filter(mask), index_cols).to_list())
        public_keys_per_group[name] = keys
        print(f"  [{name}]: {len(keys):,} public clonotypes")

    # ---- Step 3: Union of public clonotypes ----
    union_keys: set[str] = set()
    for keys in public_keys_per_group.values():
        union_keys |= keys
    print(f"\n  Union: {len(union_keys):,} clonotypes (shared >= {args.A} in >= 1 group)")

    if not union_keys:
        sys.exit(
            f"ERROR: No clonotypes found with sharing >= {args.A}. "
            "Try a smaller --A value."
        )

    # ---- Step 4: Subset each summary to union, write per-group files ----
    subsetted: dict[str, pl.DataFrame] = {}

    # Determine the global N range for the merged file (A .. global_max)
    global_max_N = max(
        int(c) for s in summaries.values() for c in s.columns if c.isdigit()
    )

    print(f"\nWriting per-group summaries ...")
    for name, summary in summaries.items():
        sub = _subset_to_keys(summary, index_cols, union_keys)

        # Ensure all N columns from A..global_max exist (zero-fill missing)
        existing_n_cols = {c for c in sub.columns if c.isdigit()}
        for N in range(args.A, global_max_N + 1):
            col = str(N)
            if col not in existing_n_cols:
                sub = sub.with_columns(pl.lit(0).cast(pl.Int32).alias(col))

        # Keep only columns from A onward (drop 2..A-1) and sort N columns
        n_cols_to_keep = sorted(
            [c for c in sub.columns if c.isdigit() and int(c) >= args.A],
            key=int,
        )
        sub = sub.select(index_cols + n_cols_to_keep)

        subsetted[name] = sub

        out_path = os.path.join(args.output, f"{name}_summary.tsv")
        sub.write_csv(out_path, separator="\t")
        print(f"  Written: {out_path}  ({len(sub):,} clonotypes × {len(n_cols_to_keep)} N-levels)")

    # ---- Step 5: Merged file ----
    print(f"\nBuilding merged output ...")
    merged_parts = []
    for name, sub in subsetted.items():
        part = sub.with_columns(pl.lit(name).alias("group"))
        merged_parts.append(part)

    merged = pl.concat(merged_parts, how="diagonal_relaxed")

    # Sort by index cols for readability
    merged = merged.sort(index_cols)

    merged_path = os.path.join(args.output, "merged_summary.tsv")
    merged.write_csv(merged_path, separator="\t")
    print(f"  Written: {merged_path}  ({len(merged):,} rows)")

    print("\nDone.")


if __name__ == "__main__":
    main()
