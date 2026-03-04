#!/usr/bin/env python
# rarefied_diversity_metrics.py
#
# Usage:
#   python rarefied_diversity_metrics.py \
#       --file_path <input.tsv.gz> \
#       --output <output.tsv> \
#       --group_by 'sample_short' \
#       --iterations ${iterations} \
#       --sequence_cols ${sequence_cols.join(',')} \
#       --count_col "${count_col}" \
#       --rarefy ${rarefy} \
#       --min_count ${min_count}

import argparse
import sys
from pathlib import Path

import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import assets.diversity.changed_diversity_AL as mod


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate diversity indices for AIRR-seq repertoires."
    )
    parser.add_argument("--file_path", required=True,
                        help="Path to input TSV (optionally gzip-compressed) file.")
    parser.add_argument("--output", required=True,
                        help="Path to output TSV file.")
    parser.add_argument("--group_by", default="repertoire_id",
                        help="Column to group repertoires by (default: repertoire_id).")
    parser.add_argument("--sequence_cols", default="junction_aa",
                        help="Comma-separated column(s) used as sequence index (default: junction_aa).")
    parser.add_argument("--count_col", default="duplicate_count",
                        help="Name of the abundance/count column (default: duplicate_count).")
    parser.add_argument("--rarefy", action="store_true", default=False,
                        help="Perform rarefaction analysis (default: False).")
    parser.add_argument("--min_count", type=int, default=None,
                        help="Minimum sequence count for rarefaction.")
    parser.add_argument("--iterations", type=int, default=100,
                        help="Number of rarefaction iterations (default: 100).")
    return parser.parse_args()


def main():
    args = parse_args()
    sequence_cols = [c.strip() for c in args.sequence_cols.split(",")]

    print(f"Reading input: {args.file_path}")
    

    data = pl.read_csv(args.file_path, separator="\t")
    print(f"Loaded Polars DataFrame: {data.shape[0]} rows x {data.shape[1]} cols")


    print(
        f"Computing diversity indices | group_by={args.group_by} | "
        f"sequence_cols={sequence_cols} | count_col={args.count_col} | "
        f"rarefy={args.rarefy} | min_count={args.min_count} | iterations={args.iterations}"
    )

    result = mod.diversity_indices(
        data,
        group_by=args.group_by,
        sequence_cols=sequence_cols,
        count_col=args.count_col,
        rarefy=args.rarefy,
        min_count=args.min_count,
        iterations=args.iterations,
    )

    print(f"Writing output: {args.output}")
    if isinstance(result, pl.DataFrame):
        result.write_csv(args.output, separator="\t")
    else:
        result.to_csv(args.output, sep="\t", index=False)

    print("Done.")


if __name__ == "__main__":
    main()

