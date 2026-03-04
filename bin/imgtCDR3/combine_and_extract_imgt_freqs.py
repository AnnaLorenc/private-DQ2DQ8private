#!/usr/bin/env python3
"""Combine multiple IMGT amino‑acid frequency tables and extract a column subset.

This script:
- reads one or more tab‑separated IMGT AA frequency files (e.g. per‑sample *_freq.tsv
    or *_freq_WL.tsv outputs),
- annotates each row with its source file name,
- concatenates all inputs into a single combined table, and
- optionally writes a second file containing only a user‑specified subset of
    columns (plus the source identifier) for downstream plotting or analysis.

Typical usage:

        python combine_and_extract_imgt_freqs.py \
                --input_files sample1_freq.tsv sample2_freq.tsv \
                --output_combined combined_freqs.tsv \
                --output_subset combined_freqs_subset.tsv \
                --subset_columns IMGT_position AA median_counts_freq median_rows_freq

"""

import pandas as pd
import argparse
import sys
from pathlib import Path
from typing import List, Dict, Optional





def collate_metrics_by_index(
    input_files: List[str],
    output_prefix: str,
    index_columns: List[str],
) -> None:
    """Collate *_orig/*_subs metrics from multiple med files into 8 wide tables.

    For each metric suffix:
      - orig_counts, orig_rows, orig_counts_freq, orig_rows_freq,
      - subs_counts, subs_rows, subs_counts_freq, subs_rows_freq

    this function builds a wide table where:
      - rows are defined by the shared index columns (e.g. aminoAcid_length, AA, IMGT_position)
      - columns are all sample-specific columns whose names end with the metric
        suffix (e.g. F5302E_orig_counts, F16018N_orig_counts, ...),
      - missing values are filled with 0.

    One TSV file is written per metric, with name:

        {output_prefix}_{metric}.tsv

    where metric is one of:
        orig_counts, orig_rows, orig_counts_freq, orig_rows_freq,
        subs_counts, subs_rows, subs_counts_freq, subs_rows_freq.
    """

    metric_suffixes = [
        "orig_counts",
        "orig_rows",
        "orig_counts_freq",
        "orig_rows_freq",
        "subs_counts",
        "subs_rows",
        "subs_counts_freq",
        "subs_rows_freq",
    ]

    # First pass: read all files and determine common index columns actually present
    dfs: List[Dict[str, object]] = []
    common_idx: Optional[set] = None

    print(f"Collating metrics from {len(input_files)} input files")

    for file_path in input_files:
        print(f"Reading: {file_path}")
        try:
            df = pd.read_csv(file_path, sep="\t")
        except Exception as e:
            print(f"  ! Error reading {file_path}: {e}")
            continue

        present_idx = [c for c in index_columns if c in df.columns]
        if not present_idx:
            print(f"  ! Skipping {file_path} (no requested index columns found)")
            continue

        if common_idx is None:
            common_idx = set(present_idx)
        else:
            common_idx &= set(present_idx)

        dfs.append({"path": file_path, "df": df})

    if not dfs:
        print("No valid files to process for metric collation!")
        sys.exit(1)

    if not common_idx:
        print("No common index columns found across input files; cannot collate.")
        sys.exit(1)

    idx_cols = [c for c in index_columns if c in common_idx]
    print(f"Using index columns: {idx_cols}")

    # Prepare accumulators per metric
    collated: Dict[str, Optional[pd.DataFrame]] = {m: None for m in metric_suffixes}

    # Second pass: build wide tables per metric
    for entry in dfs:
        file_path = entry["path"]
        df = entry["df"]
        print(f"Processing metrics from: {file_path}")

        for metric in metric_suffixes:
            suffix = f"_{metric}"
            metric_cols = [c for c in df.columns if c.endswith(suffix)]
            if not metric_cols:
                continue

            sub = df[idx_cols + metric_cols].copy()

            if collated[metric] is None:
                collated[metric] = sub
            else:
                collated[metric] = collated[metric].merge(
                    sub,
                    on=idx_cols,
                    how="outer",
                )

    # Write one file per metric, filling missing values with 0 in value columns
    for metric in metric_suffixes:
        metric_df = collated.get(metric)
        if metric_df is None:
            print(f"No columns found for metric '{metric}', writing empty table.")
            # Still write an empty file with just the index header
            empty_df = pd.DataFrame(columns=idx_cols)
            out_path = f"{output_prefix}_{metric}.tsv"
            empty_df.to_csv(out_path, sep="\t", index=False)
            continue

        value_cols = [c for c in metric_df.columns if c not in idx_cols]
        metric_df[value_cols] = metric_df[value_cols].fillna(0)

        # Sort rows by index columns for consistency
        metric_df = metric_df.sort_values(by=idx_cols)

        out_path = f"{output_prefix}_{metric}.tsv"
        print(f"Writing {metric} table to: {out_path}")
        metric_df.to_csv(out_path, sep="\t", index=False)

    print("Metric collation complete!")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Collate IMGT AA frequency med files into 8 wide tables "
            "(orig/subs counts/rows and their frequencies) joined on index columns."
        )
    )
    parser.add_argument(
        "--input_files",
        nargs="+",
        required=True,
        help="Input TSV med/freq files to process",
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help=(
            "Prefix for metric-collated outputs. "
            "Files will be named {output_prefix}_{metric}.tsv."
        ),
    )
    parser.add_argument(
        "--index_columns",
        nargs="+",
        default=["aminoAcid_length", "AA", "IMGT_position"],
        help=(
            "Index columns to join on when collating metrics; must be present "
            "in all input files. Default: aminoAcid_length AA IMGT_position"
        ),
    )

    args = parser.parse_args()

    print(f"Input files: {len(args.input_files)}")
    print(f"Collating metrics with index columns: {args.index_columns}")

    collate_metrics_by_index(
        args.input_files,
        args.output_prefix,
        args.index_columns,
    )

if __name__ == "__main__":
    main()