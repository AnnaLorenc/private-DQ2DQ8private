#!/usr/bin/env python3
"""ad_length_distribution.py

Compare CDR3 length distributions across multiple groups using
Anderson-Darling k-sample test, with optional permutation p-value.

Input
-----
--freqs_dir   : directory containing *_summ.tsv files (one per sample)
--sample_map  : CSV with sample IDs and grouping columns
--group_col   : grouping column used as AD groups (supports >2 groups)
--sample_col  : sample ID column in sample_map
--output      : output TSV path
--freqs_col   : 'orig' (col 3) or 'mean' (col 4) from *_summ.tsv files
--cells       : optional sample suffix filter (N or E)
--split_col   : optional subgroup column; test runs independently per subgroup
--groups      : optional list of group labels to include; default: all present
--min_samples : minimum samples per group to include in test
--n_perm      : number of label permutations for empirical p-value (0 disables)
--pseudo_n    : pseudo-count used to convert frequency distributions to integer
                observations required by scipy.stats.anderson_ksamp
--seed        : RNG seed for reproducibility

Output columns
--------------
cells, [split_col], n_groups, groups, n_perm,
ad_stat, ad_pvalue_asymptotic, ad_pvalue_permutation,
n_samples_total, n_per_group
"""

import argparse
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import anderson_ksamp


def dist_to_observations(dist: dict, pseudo_n: int) -> np.ndarray:
    """Convert a length->frequency dictionary to integer observations.

    Anderson-Darling k-sample expects raw observations, not weighted frequencies.
    We approximate each distribution by pseudo_n draws using deterministic
    rounding of frequency masses.
    """
    lengths = np.array(sorted(dist.keys()), dtype=int)
    probs = np.array([dist[l] for l in lengths], dtype=float)

    total = probs.sum()
    if total <= 0:
        return np.array([], dtype=int)

    probs = probs / total
    counts = np.floor(probs * pseudo_n).astype(int)

    # Distribute residual counts by largest remainders for stable totals.
    remainder = pseudo_n - counts.sum()
    if remainder > 0:
        frac = probs * pseudo_n - counts
        top_idx = np.argsort(frac)[::-1][:remainder]
        counts[top_idx] += 1

    obs = np.repeat(lengths, counts)
    return obs.astype(int)


def _arr_to_observations(row: np.ndarray, lengths: np.ndarray, pseudo_n: int) -> np.ndarray:
    """Fast variant of dist_to_observations operating on pre-aligned numpy arrays.

    Avoids dict creation and key sorting — used in the hot permutation path.
    `lengths` must match the column order of the frequency matrix.
    """
    total = row.sum()
    if total <= 0:
        return np.array([], dtype=int)
    probs = row / total
    counts = np.floor(probs * pseudo_n).astype(int)
    remainder = pseudo_n - counts.sum()
    if remainder > 0:
        frac = probs * pseudo_n - counts
        counts[np.argsort(frac)[::-1][:remainder]] += 1
    return np.repeat(lengths, counts)


def _compute_ad_stat(obs_list: list) -> tuple[float, float]:
    """Anderson-Darling k-sample statistic from pre-built observation arrays."""
    res = anderson_ksamp(obs_list, variant="midrank")
    return float(res.statistic), float(res.pvalue)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Anderson-Darling k-sample test for CDR3 length distributions."
    )
    parser.add_argument("--freqs_dir", required=True, help="Directory with *_summ.tsv files")
    parser.add_argument("--sample_map", required=True, help="CSV with sample IDs and grouping columns")
    parser.add_argument("--group_col", default="genotype_short", help="Grouping column in sample_map")
    parser.add_argument("--sample_col", default="shortname", help="Sample ID column in sample_map")
    parser.add_argument("--output", default="ad_results.tsv", help="Output TSV path")
    parser.add_argument("--freqs_col", choices=["orig", "mean"], default="mean")
    parser.add_argument("--cells", choices=["N", "E"], default=None)
    parser.add_argument("--split_col", default=None,
                        help="Optional column to run separate tests per subgroup")
    parser.add_argument("--groups", nargs="+", default=None,
                        help="Optional list of group labels to include")
    parser.add_argument("--min_samples", type=int, default=3)
    parser.add_argument("--n_perm", type=int, default=10000,
                        help="Permutation count (0 disables permutation p-value)")
    parser.add_argument("--pseudo_n", type=int, default=10000,
                        help="Pseudo-observation count per group distribution")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--n_jobs", type=int, default=-1,
        help="Parallel workers for permutation loop (-1 = all cores, 1 = serial)",
    )
    args = parser.parse_args()

    freq_col_idx = 2 if args.freqs_col == "orig" else 3

    # Suppress expected scipy warnings once for the entire run (thread-safe).
    warnings.filterwarnings("ignore", message="p-value capped:.*")
    warnings.filterwarnings("ignore", message="p-value floored:.*")

    # ---- load sample map ----
    smap_df = pd.read_csv(args.sample_map)
    smap_df.columns = [c.strip() for c in smap_df.columns]

    for col in [args.sample_col, args.group_col]:
        if col not in smap_df.columns:
            raise ValueError(f"sample_map missing required column: '{col}'")
    if args.split_col and args.split_col not in smap_df.columns:
        raise ValueError(f"sample_map missing split column: '{args.split_col}'")

    group_map = smap_df.set_index(args.sample_col)[args.group_col].to_dict()
    split_map = (
        smap_df.set_index(args.sample_col)[args.split_col].to_dict()
        if args.split_col else None
    )

    # ---- read per-sample length distributions ----
    all_rows = {}
    for tsv in sorted(Path(args.freqs_dir).glob("*_summ.tsv")):
        sample = tsv.name.replace(".tsv_summ.tsv", "").replace("_summ.tsv", "")
        if args.cells and not sample.endswith(args.cells):
            continue

        row = {}
        with open(tsv) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                if parts[0].strip() != "cdr3_length":
                    continue
                try:
                    length = int(parts[1].strip())
                    freq = float(parts[freq_col_idx].strip())
                except ValueError:
                    continue
                row[length] = freq

        if row:
            all_rows[sample] = row

    if not all_rows:
        raise RuntimeError(f"No cdr3_length rows found in {args.freqs_dir}")

    freq_df = pd.DataFrame(all_rows).T.fillna(0.0)
    freq_df.index.name = "sample"

    # Match sample IDs from file names to sample_map IDs (base and full forms).
    def map_sample(sample: str, mapping: dict):
        base = sample[:-1] if sample[-1] in ("N", "E") else sample
        return mapping.get(base, mapping.get(sample, None))

    freq_df["group"] = [map_sample(s, group_map) for s in freq_df.index]
    freq_df = freq_df.dropna(subset=["group"])

    if split_map:
        freq_df["_split"] = [map_sample(s, split_map) for s in freq_df.index]
        freq_df = freq_df.dropna(subset=["_split"])
        subgroup_values = sorted(freq_df["_split"].unique())
    else:
        subgroup_values = [None]

    length_cols = [c for c in freq_df.columns if c not in ("group", "_split")]
    rng = np.random.default_rng(args.seed)
    results = []

    for subgroup in subgroup_values:
        sub_df = freq_df[freq_df["_split"] == subgroup] if subgroup is not None else freq_df

        groups_present = sorted(sub_df["group"].unique())
        if args.groups:
            groups = [g for g in args.groups if g in groups_present]
        else:
            groups = groups_present

        # Keep groups with enough samples.
        groups = [
            g for g in groups
            if (sub_df["group"] == g).sum() >= args.min_samples
        ]

        if len(groups) < 2:
            tag = f"{args.split_col}={subgroup}" if subgroup is not None else "all"
            print(f"Skipping {tag}: fewer than 2 eligible groups")
            continue

        feat = sub_df[length_cols]
        samples_by_group = {
            g: sub_df.index[sub_df["group"] == g].tolist()
            for g in groups
        }

        # Pre-compute per-sample observation arrays once (reused in every permutation).
        feat_np = feat.values
        sample_to_idx = {s: i for i, s in enumerate(feat.index)}
        lengths_arr = np.array([int(c) for c in length_cols], dtype=int)
        all_samples = [s for g in groups for s in samples_by_group[g]]
        obs_per_sample = {
            s: _arr_to_observations(feat_np[sample_to_idx[s]], lengths_arr, args.pseudo_n)
            for s in all_samples
        }

        def _build_obs_list(assignment: dict) -> list:
            return [
                np.concatenate([obs_per_sample[s] for s in assignment[g]])
                for g in groups
            ]

        # Observed AD statistic
        ad_stat, ad_p_asym = _compute_ad_stat(_build_obs_list(samples_by_group))

        # Permutation p-value by shuffling group labels among eligible samples.
        if args.n_perm > 0:
            group_sizes = [len(samples_by_group[g]) for g in groups]

            # Pre-generate all permuted label sequences before launching workers.
            all_perms = [rng.permutation(all_samples) for _ in range(args.n_perm)]

            def _eval_perm(perm) -> float:
                start = 0
                obs_list = []
                for size in group_sizes:
                    obs_list.append(
                        np.concatenate([obs_per_sample[s] for s in perm[start:start + size]])
                    )
                    start += size
                stat, _ = _compute_ad_stat(obs_list)
                return stat

            tag = f"{args.split_col}={subgroup}" if subgroup is not None else "all"
            max_workers = args.n_jobs if args.n_jobs > 0 else None
            print(f"Running {args.n_perm} permutations ({tag}) ...", flush=True)
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                perm_stats = np.fromiter(
                    pool.map(_eval_perm, all_perms, chunksize=100),
                    dtype=float,
                    count=args.n_perm,
                )

            ad_p_perm = float(((perm_stats >= ad_stat).sum() + 1) / (args.n_perm + 1))
        else:
            ad_p_perm = np.nan

        n_per_group = {g: len(samples_by_group[g]) for g in groups}
        row = {
            "cells": args.cells if args.cells else "all",
            "n_groups": len(groups),
            "groups": "|".join(groups),
            "n_perm": args.n_perm,
            "ad_stat": round(ad_stat, 6),
            "ad_pvalue_asymptotic": ad_p_asym,
            "ad_pvalue_permutation": ad_p_perm,
            "n_samples_total": sum(n_per_group.values()),
            "n_per_group": ";".join([f"{g}:{n_per_group[g]}" for g in groups]),
        }
        if subgroup is not None:
            row[args.split_col] = subgroup

        results.append(row)

        tag = f"{args.split_col}={subgroup}" if subgroup is not None else "all"
        print(
            f"AD ({tag}): groups={','.join(groups)} stat={ad_stat:.4f} "
            f"p_perm={ad_p_perm:.4f}"
        )

    out = pd.DataFrame(results)
    if args.split_col and args.split_col in out.columns:
        cols = [
            "cells", args.split_col, "n_groups", "groups", "n_perm",
            "ad_stat", "ad_pvalue_asymptotic", "ad_pvalue_permutation",
            "n_samples_total", "n_per_group",
        ]
        out = out[[c for c in cols if c in out.columns]]

    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nResults written -> {args.output}")


if __name__ == "__main__":
    main()
