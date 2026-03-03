#!/usr/bin/env python3
"""
Script: rare_pairwise_sharing.py

Description
-----------
Rarefaction-based pairwise clonotype sharing between samples - better than bin/compute_overlaps_subs.py

For each subsampling round (K rounds, M clonotypes drawn per sample):
  - within_group entries: assess all pairwise sharing within each named group.
  - between_group entries: assess pairwise sharing for every (A-sample, B-sample)
    cross-product between two named groups.

Output: one TSV per declared comparison (within or between), with columns:
    sample_a, sample_b, k1, k2, ..., kK
where each k_i holds the number of clonotypes shared by that pair in round i.

Usage
-----
python rare_pairwise_sharing.py \\
    --input  /path/to/merged.tsv.gz \\
    --K 25 \\
    --M 5000 \\
    --groups_file sandbox/pairwise_share/pairwise_groups_E.yml \\
    --index_cols aminoAcid,vFamilyName,jGeneName \\
    --output results/pairwise_out \\
    --seed 42

YAML groups file format
-----------------------
# Named sample groups
groups:
  homo:
    - F28396E
    - F28396N
    - F20790N
  hete:
    - T04918N
    - T06306E
    - T13372E

# Pairwise comparisons to run
within_group:
  - homo
  - hete

between_group:
  - [homo, hete]
  - [homo, homo]   # cross between two differently-named copies is also supported

Arguments
---------
--input        Path to input file (TSV or TSV.gz).
--K            Number of subsampling rounds (default: 25).
--M            Subsampling depth per sample per round.
--groups_file  Path to YAML file (see format above).
--index_cols   Comma-separated index column names
               (default: aminoAcid,vFamilyName,jGeneName).
--output       Output directory (created if absent).
--seed         Optional random seed for reproducibility.
"""

import argparse
import os
import sys
from itertools import combinations

import numpy as np
import polars as pl


# ---------------------------------------------------------------------------
# YAML loading
# ---------------------------------------------------------------------------

def load_groups_file(path: str) -> tuple[dict, list, list]:
    """
    Parse the YAML groups file.

    Returns
    -------
    groups        : dict  {group_name: [sample, ...]}
    within_group  : list  [group_name, ...]
    between_group : list  [(group_a, group_b), ...]

    Supports between_group entries as either:
      - space-separated string:  "DQ2 DQ8"
      - YAML list:               [DQ2, DQ8]
    """
    try:
        import yaml
    except ImportError:
        sys.exit("ERROR: PyYAML is required. Install with: pip install pyyaml")

    with open(path) as fh:
        data = yaml.safe_load(fh)

    if not isinstance(data, dict):
        sys.exit(f"ERROR: YAML file must be a mapping, got {type(data).__name__}")

    # --- groups block ---
    raw_groups = data.get("groups", {})
    if not raw_groups:
        sys.exit("ERROR: YAML must contain a 'groups' block.")
    groups = {}
    for name, members in raw_groups.items():
        members = [str(s).strip() for s in (members or []) if str(s).strip()]
        if not members:
            sys.exit(f"ERROR: Group '{name}' has no samples.")
        groups[name] = members

    # --- within_group block ---
    within_raw = data.get("within_group", []) or []
    within_group = []
    for entry in within_raw:
        entry = str(entry).strip()
        if entry not in groups:
            sys.exit(f"ERROR: within_group entry '{entry}' not in groups block.")
        within_group.append(entry)

    # --- between_group block ---
    # Accepts both "GroupA GroupB" (space-separated string) and [GroupA, GroupB]
    between_raw = data.get("between_group", []) or []
    between_group = []
    for pair in between_raw:
        if isinstance(pair, str):
            parts = pair.split()
            if len(parts) != 2:
                sys.exit(
                    f"ERROR: between_group string entry must be 'GroupA GroupB', got: '{pair}'"
                )
            a, b = parts[0].strip(), parts[1].strip()
        elif isinstance(pair, (list, tuple)) and len(pair) == 2:
            a, b = str(pair[0]).strip(), str(pair[1]).strip()
        else:
            sys.exit(
                f"ERROR: Each between_group entry must be a 2-element list or "
                f"space-separated string, got: {pair}"
            )
        for g in (a, b):
            if g not in groups:
                sys.exit(f"ERROR: between_group references unknown group '{g}'.")
        between_group.append((a, b))

    if not within_group and not between_group:
        sys.exit("ERROR: YAML must define at least one within_group or between_group entry.")

    return groups, within_group, between_group


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_and_subset(
    input_file: str,
    index_cols: list[str],
    samples: list[str],
) -> pl.DataFrame:
    """Load input TSV, keep index_cols + sample columns, filter empty rows."""
    print(f"  Reading: {input_file}")
    df = pl.read_csv(input_file, separator="\t")

    missing_idx = [c for c in index_cols if c not in df.columns]
    if missing_idx:
        sys.exit(f"ERROR: Index columns not found: {missing_idx}")

    missing_s = [s for s in samples if s not in df.columns]
    if missing_s:
        sys.exit(f"ERROR: Sample columns not found: {missing_s}")

    df = df.select(index_cols + samples)
    df = df.with_columns([pl.col(s).cast(pl.Int64) for s in samples])

    any_nonzero = pl.any_horizontal([pl.col(s) > 0 for s in samples])
    df = df.filter(any_nonzero)
    return df


# ---------------------------------------------------------------------------
# Sparse pre-computation (reused from rare_publicity.py)
# ---------------------------------------------------------------------------

def precompute_sparse(count_matrix: np.ndarray) -> list:
    """Pre-compute non-zero indices and normalised probabilities per column."""
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


# ---------------------------------------------------------------------------
# Rarefaction + pairwise sharing
# ---------------------------------------------------------------------------

def run_pairwise_rarefaction(
    count_matrix: np.ndarray,
    sample_names: list[str],
    groups: dict[str, list[str]],
    within_group: list[str],
    between_group: list[tuple[str, str]],
    K: int,
    M: int,
    n_clonotypes: int,
) -> dict[str, dict]:
    """
    Run K rarefaction rounds and accumulate pairwise sharing counts.

    Returns
    -------
    dict keyed by comparison label (e.g. "within_homo", "between_homo_hete"),
    each value is a dict:
        {
          "pairs": [(sample_a, sample_b), ...],
          "matrix": np.ndarray shape (n_pairs, K)   # sharing counts
        }
    """
    sample_idx = {s: i for i, s in enumerate(sample_names)}

    # Build pair lists for each comparison
    # Each entry: label -> {pairs, matrix, type_of_sharing, group_label}
    comparisons: dict[str, dict] = {}

    for gname in within_group:
        gsamples = groups[gname]
        pairs = list(combinations(gsamples, 2))
        if not pairs:
            print(f"  WARNING: within_group '{gname}' has < 2 samples; skipping.")
            continue
        comparisons[f"within_{gname}"] = {
            "pairs": pairs,
            "type_of_sharing": "within",
            "group_label": gname,
        }

    for ga, gb in between_group:
        pairs = [(a, b) for a in groups[ga] for b in groups[gb] if a != b]
        label = f"between_{ga}_{gb}"
        if not pairs:
            print(f"  WARNING: between_group ({ga}, {gb}) produced no pairs; skipping.")
            continue
        comparisons[label] = {
            "pairs": pairs,
            "type_of_sharing": "between",
            "group_label": f"{ga}_{gb}",
        }

    if not comparisons:
        sys.exit("ERROR: No valid pairs found across all comparisons.")

    # Accumulate result matrices
    results = {
        label: {
            "pairs": info["pairs"],
            "type_of_sharing": info["type_of_sharing"],
            "group_label": info["group_label"],
            "matrix": np.zeros((len(info["pairs"]), K), dtype=np.int32),
        }
        for label, info in comparisons.items()
    }

    # Only subsample columns actually needed
    needed_samples: list[str] = []
    seen: set[str] = set()
    for comp_info in comparisons.values():
        for a, b in comp_info["pairs"]:
            for s in (a, b):
                if s not in seen:
                    needed_samples.append(s)
                    seen.add(s)

    needed_idx = np.array([sample_idx[s] for s in needed_samples], dtype=np.intp)
    needed_pos = {s: i for i, s in enumerate(needed_samples)}

    sparse_cols = precompute_sparse(count_matrix[:, needed_idx])

    print(f"\nRunning {K} subsampling rounds (M={M} per sample)...")
    for k in range(K):
        if (k + 1) % 5 == 0 or k == 0:
            print(f"    Round {k + 1}/{K} ...")

        # Build presence/absence for needed samples: shape (n_clonotypes, n_needed)
        presence = np.zeros((n_clonotypes, len(needed_samples)), dtype=np.int8)
        for j, (nz, p) in enumerate(sparse_cols):
            if len(nz) == 0:
                continue
            drawn = np.random.multinomial(M, p)
            presence[nz[drawn > 0], j] = 1

        # Accumulate pairwise sharing
        for label, info in results.items():
            for p_idx, (sa, sb) in enumerate(info["pairs"]):
                ja = needed_pos[sa]
                jb = needed_pos[sb]
                shared = int(np.sum(presence[:, ja] & presence[:, jb]))
                info["matrix"][p_idx, k] = shared

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rarefaction-based pairwise clonotype sharing.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True,
                        help="Path to input TSV or TSV.gz.")
    parser.add_argument("--K", type=int, default=25,
                        help="Number of subsampling rounds (default: 25).")
    parser.add_argument("--M", type=int, required=True,
                        help="Subsampling depth per sample per round.")
    parser.add_argument("--groups_file", required=True,
                        help="Path to YAML file defining groups and comparisons.")
    parser.add_argument("--index_cols",
                        default="aminoAcid,vFamilyName,jGeneName",
                        help="Comma-separated index columns.")
    parser.add_argument("--output", required=True,
                        help="Output directory (created if absent).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility.")

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        print(f"Random seed set to {args.seed}")

    if args.K < 1:
        sys.exit("ERROR: --K must be >= 1.")
    if args.M < 1:
        sys.exit("ERROR: --M must be >= 1.")

    index_cols = [c.strip() for c in args.index_cols.split(",") if c.strip()]

    # Load YAML
    groups, within_group, between_group = load_groups_file(args.groups_file)

    # Derive full sample set (union of all referenced samples)
    seen: set[str] = set()
    all_samples: list[str] = []
    for gname in list(within_group) + [s for pair in between_group for s in pair]:
        for s in groups[gname]:
            if s not in seen:
                all_samples.append(s)
                seen.add(s)

    print("=" * 60)
    print("rare_pairwise_sharing.py")
    print("=" * 60)
    print(f"  Input         : {args.input}")
    print(f"  K (rounds)    : {args.K}")
    print(f"  M (depth)     : {args.M}")
    print(f"  within_group  : {within_group}")
    print(f"  between_group : {between_group}")
    for gname, gsamples in groups.items():
        print(f"  Group '{gname}' ({len(gsamples)}): {gsamples}")
    print(f"  Index cols    : {index_cols}")
    print(f"  Output dir    : {args.output}")
    print()

    # Load data
    df = load_and_subset(args.input, index_cols, all_samples)
    n_clonotypes = len(df)
    print(f"  Clonotypes non-zero in sample set: {n_clonotypes}")

    if n_clonotypes == 0:
        sys.exit("ERROR: No clonotypes remaining after filtering.")

    count_matrix = df.select(all_samples).to_numpy().astype(np.float64)
    index_df = df.select(index_cols)

    # Run rarefaction
    results = run_pairwise_rarefaction(
        count_matrix=count_matrix,
        sample_names=all_samples,
        groups=groups,
        within_group=within_group,
        between_group=between_group,
        K=args.K,
        M=args.M,
        n_clonotypes=n_clonotypes,
    )

    # Write all results into one combined output file
    os.makedirs(args.output, exist_ok=True)

    sharing_by = "_".join(index_cols)
    all_frames: list[pl.DataFrame] = []

    for label, info in results.items():
        pairs  = info["pairs"]
        matrix = info["matrix"]   # shape (n_pairs, K)

        rows: dict = {
            "type_of_sharing": [info["type_of_sharing"]] * len(pairs),
            "groups":          [info["group_label"]]     * len(pairs),
            "sharing_by":      [sharing_by]              * len(pairs),
            "sample_a":        [p[0] for p in pairs],
            "sample_b":        [p[1] for p in pairs],
        }
        for k in range(args.K):
            rows[f"k{k + 1}"] = matrix[:, k].tolist()

        all_frames.append(pl.DataFrame(rows))

    combined_df = pl.concat(all_frames)
    out_path = os.path.join(args.output, "pairwise_sharing.tsv")
    combined_df.write_csv(out_path, separator="\t")

    total_pairs = sum(len(info["pairs"]) for info in results.values())
    print(f"  Written: {out_path}  ({total_pairs} pairs × {args.K} rounds)")
    print("\nDone.")


if __name__ == "__main__":
    main()
