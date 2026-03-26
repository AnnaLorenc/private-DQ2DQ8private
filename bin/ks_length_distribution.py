#!/usr/bin/env python3
"""ks_length_distribution.py

Compare CDR3 length distributions between sample groups using
a Kolmogorov-Smirnov statistic with a permutation-based null distribution.

High-level workflow:
1. Read per-sample *_summ.tsv files and extract rows where param == cdr3_length.
2. Build a sample x length frequency matrix.
3. Map each sample to a comparison group using sample_map.
4. Optionally split the analysis into independent subgroups (--split_col).
5. For each comparison, compute observed KS and permutation p-value.

Input
-----
--freqs_dir   : directory containing *_summ.tsv files (one per sample).
                Rows with param == "cdr3_length" are used.
                Expected columns (tab-separated, no header by default):
                    param  group(=length)  orig_freq  mean_subsampled_freq
--sample_map  : CSV with at least columns: shortname (sample ID) and one
                grouping column (see --group_col; default: genotype_short)
--output      : output TSV path                     (default: ks_results.tsv)
--freqs_col   : which frequency column to use:
                  "orig"  → column 3 (original frequencies)
                  "mean"  → column 4 (mean of subsamples, default)
--n_perm      : number of permutations               (default: 10000)
--groups      : two genotype_short values to compare (default: all pairs)
--cells       : cell-type suffix filter, e.g. N or E (default: both)
--split_col   : optional column in sample_map to split samples into subgroups;
                the test is run independently within each subgroup value
                (e.g. --split_col genotype_short to test within each genotype)
--min_samples : minimum samples per group to run test (default: 3)

Output
------
TSV with columns:
    group_a, group_b, cells, subgroup (if --split_col used), n_a, n_b, ks_stat, perm_p_value, n_perm
"""

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def weighted_ecdf(lengths: np.ndarray, weights: np.ndarray):
    """Return sorted lengths and cumulative weights.

    Note: kept as a utility for future extensions; current KS implementation
    computes cumulative sums directly from aligned vectors.
    """
    order = np.argsort(lengths)
    l_sorted = lengths[order]
    w_sorted = weights[order] / weights.sum()
    cdf = np.cumsum(w_sorted)
    return l_sorted, cdf


def ks_stat_weighted(dist_a: dict, dist_b: dict) -> float:
    """KS statistic between two length->frequency dicts.

    dist_a and dist_b are expected to represent group-level mean length
    distributions (averaged across samples in each group).
    """
    all_lengths = np.array(sorted(set(dist_a) | set(dist_b)), dtype=float)
    # Build mean frequency vectors over all seen lengths
    freq_a = np.array([dist_a.get(l, 0.0) for l in all_lengths])
    freq_b = np.array([dist_b.get(l, 0.0) for l in all_lengths])
    # Normalise (should already sum to ~1, but guard against rounding)
    fa = freq_a / freq_a.sum() if freq_a.sum() > 0 else freq_a
    fb = freq_b / freq_b.sum() if freq_b.sum() > 0 else freq_b
    # KS stat = max |CDF_A - CDF_B|
    cdf_a = np.cumsum(fa)
    cdf_b = np.cumsum(fb)
    return float(np.max(np.abs(cdf_a - cdf_b)))


def group_mean_dist(freq_matrix: pd.DataFrame, samples: list) -> dict:
    """Return mean length distribution across selected samples.

    The output is a dict keyed by CDR3 length with mean frequency values.
    """
    sub = freq_matrix.loc[samples]
    means = sub.mean(axis=0)
    return means.to_dict()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="KS test on CDR3 length distributions with permutation null."
    )
    parser.add_argument("--freqs_dir", required=True,
                        help="Directory with *_summ.tsv files.")
    parser.add_argument("--sample_map", required=True,
                        help="CSV with columns: shortname and a grouping column.")
    parser.add_argument("--group_col", default="genotype_short",
                        help="Column in sample_map to use as group label. "
                             "Default: genotype_short.")
    parser.add_argument("--sample_col", default="shortname",
                        help="Column in sample_map containing sample IDs. "
                             "Default: shortname.")
    parser.add_argument("--output", default="ks_results.tsv")
    parser.add_argument("--freqs_col", choices=["orig", "mean"], default="mean",
                        help="orig = column 3; mean = column 4 (default).")
    parser.add_argument("--n_perm", type=int, default=10000)
    parser.add_argument("--groups", nargs=2, default=None,
                        metavar=("GROUP_A", "GROUP_B"),
                        help="Two genotype_short values to compare. "
                             "Default: all pairs.")
    parser.add_argument("--cells", choices=["N", "E"], default=None,
                        help="Restrict to samples ending with N or E.")
    parser.add_argument("--split_col", default=None,
                        help="Column in sample_map to split into subgroups. "
                             "The KS test is run independently within each subgroup value.")
    parser.add_argument("--min_samples", type=int, default=3)
    args = parser.parse_args()


    freq_col_idx = 2 if args.freqs_col == "orig" else 3  # freq column in the input files : 0=param,1=length,2=orig,3=mean

    # ---- load sample map ----
    # sample_map can provide both:
    # - comparison label (--group_col)
    # - optional subgroup label (--split_col)
    smap_df = pd.read_csv(args.sample_map)
    smap_df.columns = [c.strip() for c in smap_df.columns]
    if args.sample_col not in smap_df.columns:
        raise ValueError(f"sample_map is missing sample column: '{args.sample_col}'")
    if args.group_col not in smap_df.columns:
        raise ValueError(f"sample_map is missing group column: '{args.group_col}' "
                         f"(available: {list(smap_df.columns)})")
    if args.split_col and args.split_col not in smap_df.columns:
        raise ValueError(f"sample_map is missing split column: '{args.split_col}' "
                         f"(available: {list(smap_df.columns)})")
    smap = smap_df.set_index(args.sample_col)[args.group_col].to_dict()
    # Mapping: sample ID -> subgroup value (only used when split_col is set)
    split_map = smap_df.set_index(args.sample_col)[args.split_col].to_dict() \
        if args.split_col else None

    # ---- load all freq files ----
    # Extract only rows that describe length distributions and select either
    # original or subsampled mean frequency column.
    freqs_dir = Path(args.freqs_dir)
    all_rows = {}  # sample → {length: freq}
    for tsv in sorted(freqs_dir.glob("*_summ.tsv")):
        sample = tsv.name.replace(".tsv_summ.tsv", "").replace("_summ.tsv", "")
        # cell-type filter
        if args.cells and not sample.endswith(args.cells):
            continue
        row = {}
        with open(tsv) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                if parts[0].strip() == "cdr3_length":
                    try:
                        length = int(parts[1].strip())
                        freq   = float(parts[freq_col_idx].strip())
                        row[length] = freq
                    except ValueError:
                        continue
        if row:
            all_rows[sample] = row

    if not all_rows:
        raise RuntimeError(f"No cdr3_length rows found in {freqs_dir}")

    # ---- build frequency matrix (samples x lengths) ----
    freq_df = pd.DataFrame(all_rows).T.fillna(0.0)
    freq_df.index.name = "sample"

    # ---- map samples to group labels ----
    # Input files often use sample IDs with trailing cell suffix (e.g. F5302N,
    # F5302E), while sample_map may store base IDs (F5302). Try both.
    def lookup_genotype(sample):
        base = sample[:-1] if sample[-1] in ("N", "E") else sample
        return smap.get(base, smap.get(sample, None))

    freq_df["genotype"] = [lookup_genotype(s) for s in freq_df.index]
    missing = freq_df["genotype"].isna().sum()
    if missing:
        print(f"Warning: {missing} samples could not be mapped to a genotype and will be skipped.")
    freq_df = freq_df.dropna(subset=["genotype"])

    # Attach subgroup values if requested and keep only mappable samples.
    if split_map:
        freq_df["_split"] = [
            split_map.get(s[:-1] if s[-1] in ("N", "E") else s,
                          split_map.get(s, None))
            for s in freq_df.index
        ]
        n_missing_split = freq_df["_split"].isna().sum()
        if n_missing_split:
            print(f"Warning: {n_missing_split} samples could not be mapped to split column and will be skipped.")
        freq_df = freq_df.dropna(subset=["_split"])
        subgroups = sorted(freq_df["_split"].unique())
    else:
        subgroups = [None]

    length_cols = [c for c in freq_df.columns if c not in ("genotype", "_split")]

    # ---- determine group pairs ----
    # Either use explicit --groups A B, or all pairwise combinations.
    if args.groups:
        pairs = [tuple(args.groups)]
    else:
        all_genotypes = sorted(freq_df["genotype"].unique())
        pairs = list(itertools.combinations(all_genotypes, 2))

    rng = np.random.default_rng(42)
    results = []
    cell_label = args.cells if args.cells else "all"

    for subgroup in subgroups:
        # Each subgroup is analyzed independently with its own permutation null.
        if subgroup is not None:
            sub_df = freq_df[freq_df["_split"] == subgroup]
        else:
            sub_df = freq_df
        feat = sub_df[length_cols]

        for gA, gB in pairs:
            samples_a = sub_df.index[sub_df["genotype"] == gA].tolist()
            samples_b = sub_df.index[sub_df["genotype"] == gB].tolist()

            tag = f"{gA} vs {gB} ({cell_label}" + (f", {args.split_col}={subgroup}" if subgroup else "") + ")"

            if len(samples_a) < args.min_samples or len(samples_b) < args.min_samples:
                print(f"Skipping {tag}: n_a={len(samples_a)}, n_b={len(samples_b)} < {args.min_samples}")
                continue

            n_a, n_b = len(samples_a), len(samples_b)
            all_samples = samples_a + samples_b

            obs_ks = ks_stat_weighted(
                group_mean_dist(feat, samples_a),
                group_mean_dist(feat, samples_b)
            )

            # Permutation test: shuffle sample labels between groups while
            # preserving original group sizes, then compare null KS to observed.
            perm_ks = np.empty(args.n_perm)
            for i in range(args.n_perm):
                perm = rng.permutation(all_samples)
                perm_a, perm_b = perm[:n_a].tolist(), perm[n_a:].tolist()
                perm_ks[i] = ks_stat_weighted(
                    group_mean_dist(feat, perm_a),
                    group_mean_dist(feat, perm_b)
                )

            # +1 correction avoids exact zero p-values with finite permutations.
            perm_p = float((perm_ks >= obs_ks).sum() + 1) / (args.n_perm + 1)

            row = {
                "group_a": gA, "group_b": gB, "cells": cell_label,
                "n_a": n_a, "n_b": n_b,
                "ks_stat": round(obs_ks, 6),
                "perm_p_value": round(perm_p, 6),
                "n_perm": args.n_perm
            }
            if subgroup is not None:
                row[args.split_col] = subgroup
            results.append(row)
            print(f"{tag}: KS={obs_ks:.4f}, p={perm_p:.4f} (n_a={n_a}, n_b={n_b})")

    out = pd.DataFrame(results)
    # Keep output columns in a stable order when split_col is used.
    if args.split_col and args.split_col in out.columns:
        cols = ["group_a", "group_b", "cells", args.split_col, "n_a", "n_b", "ks_stat", "perm_p_value", "n_perm"]
        out = out[[c for c in cols if c in out.columns]]
    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nResults written → {args.output}")


if __name__ == "__main__":
    main()
