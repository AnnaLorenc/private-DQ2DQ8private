#!/usr/bin/env python3
"""Filter rows from a per_factor or hotelling TSV where any raw p-value < threshold.

Detects p-value columns automatically (any column named p_value_* but NOT q_value_*
and NOT p_hotelling_*).

Usage
-----
    python filter_significant.py --input <file.tsv>
    python filter_significant.py --input <file.tsv> --pval 0.01 --output sig.tsv
"""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True, help="Input TSV file.")
    parser.add_argument("--output", default=None,
                        help="Output TSV. Default: <input_stem>_sig.tsv next to input.")
    parser.add_argument("--pval",   type=float, default=0.05,
                        help="P-value threshold (default: 0.05).")
    parser.add_argument("--digits", type=int,   default=4,
                        help="Decimal digits to round numeric columns to (default: 4).")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # Detect raw p-value columns: named p_value_* (exclude q_value and p_hotelling)
    p_cols = [c for c in df.columns
              if c.startswith("p_value_") and not c.startswith("q_")]
    if not p_cols:
        raise ValueError(f"No p_value_* columns found in {args.input}. Columns: {list(df.columns)}")

    print(f"Filtering on p-value columns: {p_cols}  (threshold < {args.pval})")

    mask = pd.Series(False, index=df.index)
    for col in p_cols:
        mask |= df[col].fillna(1.0) < args.pval

    sig = df[mask].copy()
    print(f"Significant rows: {len(sig)} / {len(df)}")

    # Round all numeric columns
    num_cols = sig.select_dtypes(include="number").columns
    sig[num_cols] = sig[num_cols].round(args.digits)

    if args.output:
        out_path = Path(args.output)
    else:
        inp = Path(args.input)
        out_path = inp.parent / f"{inp.stem}_sig{inp.suffix}"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(out_path, sep="\t", index=False)
    print(f"Written → {out_path}")


if __name__ == "__main__":
    main()
