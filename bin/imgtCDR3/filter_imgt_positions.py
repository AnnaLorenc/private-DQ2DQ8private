#!/usr/bin/env python3
"""Filter IMGT test results by IMGT_position range and p-value.

Usage
-----
python filter_imgt_positions.py \
    --input results_man/img_test/imgt_lm_homhom.tsv \
    [--output filtered.tsv]

Keeps rows where:
- IMGT_position is in 108–111 or 116–117, allowing letter extensions
  (e.g. 111A, 116B), and
- p_value < 0.01.

If --output is not provided, results are printed to stdout.
"""

import argparse
import sys
import pandas as pd
import re


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter IMGT test results to positions 108–111 and 116–117 "
            "(with extensions) where p_value < 0.01."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input TSV file, e.g. results_man/img_test/imgt_lm_homhom.tsv.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Optional output TSV file. If omitted, matching rows are "
            "written to stdout."
        ),
    )
    return parser.parse_args()


def imgt_position_in_ranges(pos: str) -> bool:
    """Return True if IMGT_position is 108–111 or 116–117 (with extensions).

    Examples that match: "108", "108A", "111C", "116", "116B", "117A".
    """

    if not isinstance(pos, str):
        pos = str(pos)

    # Capture leading numeric part and ignore any trailing letters
    m = re.match(r"^(\d+)", pos)
    if not m:
        return False
    num = int(m.group(1))

    return (108 <= num <= 111) or (116 <= num <= 117)


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.input, sep="\t")

    if "IMGT_position" not in df.columns:
        raise ValueError("Input file lacks IMGT_position column")
    if "p_value" not in df.columns:
        raise ValueError("Input file lacks p_value column")

    # Filter by position
    mask_pos = df["IMGT_position"].astype(str).apply(imgt_position_in_ranges)

    # Filter by p-value
    mask_p = df["p_value"].astype(float) < 0.01

    df_filt = df[mask_pos & mask_p]

    if args.output:
        df_filt.to_csv(args.output, sep="\t", index=False)
    else:
        df_filt.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
