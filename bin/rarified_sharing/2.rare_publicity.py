#!/usr/bin/env python3
"""
Script: rare_publicity.py

Description
-----------
Rarefaction-based clonotype publicity scoring.

Starting from a matrix of ALL_clonotypes × ALL_samples (counts), this script:
  1. Subsets to the requested samples S.
  2. Filters to clonotypes present (count > 0) in at least one sample in S.
  3. Repeats K times:
       - For each sample, draws M clonotypes with replacement, weighted by their
         original counts (multinomial draw, same logic as subsample_merged_sequences.py).
       - Records a binary presence/absence vector per sample.
       - For each sharing group G (a named subset of S), counts how many samples
         in G detected each clonotype in this round.
  4. Outputs, per group, a matrix of shape (n_clonotypes × K), where each column k
     holds the sharing counts for that subsampling round.
     The row labels (index columns identifying each clonotype) are included in the output.

The output matrices can be used downstream to compute:
  - Mean sharing across rounds
  - Probability that sharing >= threshold (public clone calling)

Usage
-----
# Groups on the command line:
python rare_publicity.py \\
    --input  /path/to/merged.tsv.gz \\
    --K 25 \\
    --M 5000 \\
    --groups homo:F28396E,F28396N,F20790N hete:T04918N,T06306E,T13372E \\
    --index_cols aminoAcid,vFamilyName,jGeneName \\
    --output results/rare_out \\
    --seed 42

# Groups from a YAML file (--samples derived automatically as union of groups):
python rare_publicity.py \\
    --input  /path/to/merged.tsv.gz \\
    --K 25 --M 5000 \\
    --groups_file groups.yml \\
    --output results/rare_out

YAML groups file format (groups.yml)
-------------------------------------
# Each top-level key is a group name; its value is a list of sample names.
homo:
  - F28396E
  - F28396N
  - F20790N
hete:
  - T04918N
  - T06306E
  - T13372E

Arguments
---------
--input        Path to input file (TSV or TSV.gz). Columns: index_cols + sample columns.
--samples      Optional comma-separated list of sample names (set S). If omitted,
               S is set to the union of all samples listed across groups.
--K            Number of subsampling rounds (default: 25).
--M            Subsampling depth: number of clonotypes to draw per sample per round.
--groups       One or more group definitions: NAME:sample1,sample2,...
               Mutually exclusive with --groups_file.
--groups_file  Path to a YAML file defining the groups (see format above).
               Mutually exclusive with --groups.
--index_cols   Comma-separated index column names (default: aminoAcid,vFamilyName,jGeneName).
--output       Output directory. Created if it does not exist.
--seed         Optional random seed for reproducibility.
"""

import argparse
import os
import sys

import numpy as np
import polars as pl


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_groups(group_strings: list[str]) -> dict[str, list[str]]:
    """
    Parse group definitions from command-line tokens.

    Parameters
    ----------
    group_strings : list of str
        Each element has the form  "NAME:sample1,sample2,..."

    Returns
    -------
    dict mapping group_name -> list of sample names
    """
    groups: dict[str, list[str]] = {}
    for token in group_strings:
        if ":" not in token:
            raise ValueError(
                f"Invalid group definition '{token}'. "
                "Expected format: NAME:sample1,sample2,..."
            )
        name, samples_str = token.split(":", 1)
        name = name.strip()
        samples = [s.strip() for s in samples_str.split(",") if s.strip()]
        if not samples:
            raise ValueError(f"Group '{name}' has no samples.")
        if name in groups:
            raise ValueError(f"Duplicate group name: '{name}'")
        groups[name] = samples
    return groups


def parse_groups_from_yaml(path: str) -> dict[str, list[str]]:
    """
    Load group definitions from a YAML file.

    Expected format::

        homo:
          - F28396E
          - F28396N
        hete:
          - T04918N
          - T06306E

    Each top-level key is a group name; its value is a list of sample names.

    Parameters
    ----------
    path : str
        Path to the YAML file.

    Returns
    -------
    dict mapping group_name -> list of sample names
    """
    try:
        import yaml  # PyYAML
    except ImportError:
        sys.exit(
            "ERROR: PyYAML is required to use --groups_file. "
            "Install with: pip install pyyaml"
        )

    with open(path) as fh:
        data = yaml.safe_load(fh)

    if not isinstance(data, dict):
        raise ValueError(
            f"YAML groups file must contain a top-level mapping, got: {type(data).__name__}"
        )

    groups: dict[str, list[str]] = {}
    for name, members in data.items():
        if not isinstance(members, list):
            raise ValueError(
                f"Group '{name}' must be a list of sample names, "
                f"got: {type(members).__name__}"
            )
        members = [str(s).strip() for s in members if str(s).strip()]
        if not members:
            raise ValueError(f"Group '{name}' has no samples.")
        if name in groups:
            raise ValueError(f"Duplicate group name: '{name}'")
        groups[name] = members

    return groups


def load_and_subset(
    input_file: str,
    index_cols: list[str],
    samples: list[str],
) -> pl.DataFrame:
    """
    Load the input file and return a DataFrame restricted to:
      - The given index columns and sample columns
      - Rows (clonotypes) that are non-zero in at least one selected sample

    Parameters
    ----------
    input_file : str
        Path to TSV or TSV.gz file.
    index_cols : list of str
        Column names identifying each clonotype.
    samples : list of str
        Sample columns to retain.

    Returns
    -------
    Polars DataFrame with columns: index_cols + samples
    """
    print(f"  Reading: {input_file}")
    df = pl.read_csv(input_file, separator="\t")

    # Validate index columns
    missing_idx = [c for c in index_cols if c not in df.columns]
    if missing_idx:
        raise ValueError(f"Index columns not found in input: {missing_idx}")

    # Validate sample columns
    missing_s = [s for s in samples if s not in df.columns]
    if missing_s:
        raise ValueError(f"Sample columns not found in input: {missing_s}")

    # Restrict to relevant columns
    df = df.select(index_cols + samples)

    # Cast sample columns to Int64 (they may be float after some joins)
    df = df.with_columns([
        pl.col(s).cast(pl.Int64) for s in samples
    ])

    # Filter to clonotypes non-zero in at least one selected sample
    any_nonzero = pl.any_horizontal([pl.col(s) > 0 for s in samples])
    df = df.filter(any_nonzero)

    return df


def precompute_sparse(count_matrix: np.ndarray) -> list:
    """
    Pre-compute per-column non-zero indices and normalised probabilities.

    This is done once before the K-round loop so that the sparse structure
    (non-zero set per sample) is not recomputed on every round. For a
    6.6M-clonotype matrix with ~330k non-zeros per sample this reduces
    per-round time from ~12s to ~2s.

    Parameters
    ----------
    count_matrix : np.ndarray, shape (n_clonotypes, n_samples)

    Returns
    -------
    List of (nz_indices, probabilities) tuples, one per sample column.
    """
    print("  Pre-computing sparse structure (done once)...")
    cols = []
    for j in range(count_matrix.shape[1]):
        col = count_matrix[:, j]
        nz = np.flatnonzero(col)
        if len(nz) == 0:
            cols.append((nz, np.array([])))
            continue
        p = col[nz].astype(np.float64)
        p /= p.sum()
        cols.append((nz, p))
    return cols


def run_rarefaction(
    count_matrix: np.ndarray,
    sample_names: list[str],
    groups: dict[str, list[str]],
    K: int,
    M: int,
) -> dict[str, np.ndarray]:
    """
    Run K rounds of rarefaction subsampling and accumulate per-group sharing
    counts across rounds.

    Sparse optimisation: non-zero indices and normalised probabilities per
    sample column are precomputed once before the K-round loop. Only columns
    belonging to the union of all groups are subsampled — samples in S but
    not in any group are skipped entirely. The inner loop then calls
    np.random.multinomial only on the non-zero subset (~330k entries instead
    of 6.6M for the real data), reducing per-round time from ~12s to ~2s
    (≈6.8 min for K=200, M=10000, 120 samples).

    Parameters
    ----------
    count_matrix : np.ndarray, shape (n_clonotypes, n_samples)
        Integer counts; columns aligned with sample_names.
    sample_names : list of str
        Sample names corresponding to columns of count_matrix.
    groups : dict
        {group_name: [sample_name, ...]}  — all samples must be in sample_names.
    K : int
        Number of subsampling rounds.
    M : int
        Subsampling depth per sample per round.

    Returns
    -------
    dict {group_name: np.ndarray of shape (n_clonotypes, K)}
        Each column k contains the sharing count (number of samples in the
        group that detected the clonotype) for subsampling round k.
    """
    n_clonotypes, n_samples = count_matrix.shape
    sample_idx = {s: i for i, s in enumerate(sample_names)}

    # Validate that group samples exist in sample_names
    for gname, gsamples in groups.items():
        missing = [s for s in gsamples if s not in sample_idx]
        if missing:
            raise ValueError(
                f"Group '{gname}' references samples not in --samples S: {missing}"
            )

    # Only subsample columns that appear in at least one group (union of G)
    union_g_samples = []
    seen = set()
    for gsamples in groups.values():
        for s in gsamples:
            if s not in seen:
                union_g_samples.append(s)
                seen.add(s)

    union_g_idx = np.array([sample_idx[s] for s in union_g_samples], dtype=np.intp)
    union_sample_pos = {s: i for i, s in enumerate(union_g_samples)}

    # Precompute sparse structure only for the union(G) columns
    sparse_cols = precompute_sparse(count_matrix[:, union_g_idx])

    # Remap group indices to positions within the reduced (union-only) matrix
    group_indices = {
        gname: np.array([union_sample_pos[s] for s in gsamples], dtype=np.intp)
        for gname, gsamples in groups.items()
    }

    n_union = len(union_g_samples)
    if n_union < n_samples:
        print(f"  Subsampling {n_union}/{n_samples} columns (union of groups)")

    # Result accumulation matrices
    results = {
        gname: np.zeros((n_clonotypes, K), dtype=np.int32)
        for gname in groups
    }

    for k in range(K):
        if (k + 1) % 5 == 0 or k == 0:
            print(f"    Round {k + 1}/{K} ...")

        # Build presence/absence matrix for this round: (n_clonotypes, n_union)
        presence = np.zeros((n_clonotypes, n_union), dtype=np.int8)
        for j, (nz, p) in enumerate(sparse_cols):
            if len(nz) == 0:
                continue
            drawn = np.random.multinomial(M, p)
            presence[nz[drawn > 0], j] = 1

        # Accumulate sharing per group
        for gname, g_idx in group_indices.items():
            results[gname][:, k] = presence[:, g_idx].sum(axis=1)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rarefaction-based clonotype publicity scoring.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to input file (TSV or TSV.gz).",
    )
    parser.add_argument(
        "--samples", default=None,
        help=(
            "Comma-separated list of sample names to use (set S). "
            "If omitted, S is derived from the union of all group samples."
        ),
    )
    parser.add_argument(
        "--K", type=int, default=25,
        help="Number of subsampling rounds (default: 25).",
    )
    parser.add_argument(
        "--M", type=int, required=True,
        help="Subsampling depth: number of clonotypes drawn per sample per round.",
    )

    groups_src = parser.add_mutually_exclusive_group(required=True)
    groups_src.add_argument(
        "--groups", nargs="+",
        help=(
            "Group definitions: NAME:sample1,sample2,... "
            "e.g. --groups homo:F28396E,F28396N hete:T04918N,T06306E"
        ),
    )
    groups_src.add_argument(
        "--groups_file",
        help=(
            "Path to a YAML file defining groups. "
            "Format: top-level keys are group names, values are lists of sample names."
        ),
    )
    parser.add_argument(
        "--index_cols", default="aminoAcid,vFamilyName,jGeneName",
        help="Comma-separated index columns (default: aminoAcid,vFamilyName,jGeneName).",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory (created if absent).",
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility.",
    )

    args = parser.parse_args()

    # ---- Parse and validate parameters ----
    if args.seed is not None:
        np.random.seed(args.seed)
        print(f"Random seed set to {args.seed}")

    index_cols = [c.strip() for c in args.index_cols.split(",") if c.strip()]

    if args.K < 1:
        sys.exit("ERROR: --K must be >= 1.")
    if args.M < 1:
        sys.exit("ERROR: --M must be >= 1.")

    # Load groups
    if args.groups_file:
        groups = parse_groups_from_yaml(args.groups_file)
    else:
        groups = parse_groups(args.groups)

    if not groups:
        sys.exit("ERROR: No groups defined.")

    # Derive sample set S
    if args.samples is not None:
        samples = [s.strip() for s in args.samples.split(",") if s.strip()]
        if not samples:
            sys.exit("ERROR: --samples is empty.")
        # Validate that all group samples are present in S
        for gname, gsamples in groups.items():
            not_in_S = [s for s in gsamples if s not in samples]
            if not_in_S:
                sys.exit(
                    f"ERROR: Group '{gname}' contains samples not listed in --samples: "
                    f"{not_in_S}"
                )
    else:
        # S = union of all group samples, preserving first-seen order
        seen: set[str] = set()
        samples = []
        for gsamples in groups.values():
            for s in gsamples:
                if s not in seen:
                    samples.append(s)
                    seen.add(s)
        print(f"  --samples not provided; using union of groups ({len(samples)} samples).")

    print("=" * 60)
    print("rare_publicity.py")
    print("=" * 60)
    print(f"  Input         : {args.input}")
    print(f"  Samples (S)   : {samples}")
    print(f"  K (rounds)    : {args.K}")
    print(f"  M (depth)     : {args.M}")
    for gname, gsamples in groups.items():
        print(f"  Group '{gname}' ({len(gsamples)} samples): {gsamples}")
    print(f"  Index cols    : {index_cols}")
    print(f"  Output dir    : {args.output}")
    print()

    # ---- Load data ----
    df = load_and_subset(args.input, index_cols, samples)

    n_clonotypes = len(df)
    print(f"  Clonotypes non-zero in S : {n_clonotypes}")

    if n_clonotypes == 0:
        sys.exit("ERROR: No clonotypes remaining after filtering. Check --samples.")

    # ---- Extract count matrix ----
    count_matrix = df.select(samples).to_numpy().astype(np.float64)
    index_df = df.select(index_cols)

    # ---- Run rarefaction ----
    print(f"\nRunning {args.K} subsampling rounds (M={args.M} per sample)...")
    results = run_rarefaction(count_matrix, samples, groups, args.K, args.M)

    # ---- Write outputs ----
    os.makedirs(args.output, exist_ok=True)

    for gname, matrix in results.items():
        out_df = index_df.with_columns([
            pl.Series(f"k{k + 1}", matrix[:, k].tolist())
            for k in range(args.K)
        ])
        out_path = os.path.join(args.output, f"{gname}_sharing.tsv")
        out_df.write_csv(out_path, separator="\t")
        print(f"  Written: {out_path}  ({n_clonotypes} clonotypes × {args.K} rounds)")

    print("\nDone.")


if __name__ == "__main__":
    main()
