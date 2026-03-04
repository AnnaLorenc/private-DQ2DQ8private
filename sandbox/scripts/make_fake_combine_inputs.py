#!/usr/bin/env python3
"""Generate 3 fake input TSV files for imgtCDR3/combine_and_extract_imgt_freqs.py.

Each file has:
  - Index columns: AA, IMGT_position, aminoAcid_length
  - 8 metric columns per sample: orig_counts, orig_rows, orig_counts_freq,
    orig_rows_freq, subs_counts, subs_rows, subs_counts_freq, subs_rows_freq

Output: F10001E_freq.tsv, F10002N_freq.tsv, F10003E_freq.tsv
        (written next to this script)

Usage
-----
    python make_fake_combine_inputs.py
    python make_fake_combine_inputs.py --seed 7 --outdir combine_and_extract_imgt_freqs
"""
import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

POSITIONS  = ["104", "107"]
AAS        = ["A", "C", "G", "T", "V"]
LENGTHS    = [8, 9]
SAMPLES    = ["F10001E", "F10002N", "F10003E"]


def make_file(sample: str, rng: np.random.Generator, outdir: Path) -> None:
    rows = []
    for pos, aa, length in itertools.product(POSITIONS, AAS, LENGTHS):
        orig_counts = int(rng.integers(4, 20))
        orig_rows   = int(rng.integers(2, orig_counts + 1))
        subs_counts = int(rng.integers(2, orig_counts + 1))
        subs_rows   = int(rng.integers(1, subs_counts + 1))
        total       = max(orig_counts, 1)
        rows.append({
            "AA":                       aa,
            "IMGT_position":            pos,
            "aminoAcid_length":         length,
            f"{sample}_orig_counts":    orig_counts,
            f"{sample}_orig_rows":      orig_rows,
            f"{sample}_orig_counts_freq": round(orig_counts / total, 4),
            f"{sample}_orig_rows_freq":   round(orig_rows   / total, 4),
            f"{sample}_subs_counts":    subs_counts,
            f"{sample}_subs_rows":      subs_rows,
            f"{sample}_subs_counts_freq": round(subs_counts / total, 4),
            f"{sample}_subs_rows_freq":   round(subs_rows   / total, 4),
        })
    df = pd.DataFrame(rows)
    # Keep exactly 10 rows
    df = df.sample(n=10, random_state=int(rng.integers(0, 999))).reset_index(drop=True)
    out = outdir / f"{sample}_freq.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"Written {len(df)} rows → {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed",   type=int,  default=42)
    parser.add_argument("--outdir", default=None,
                        help="Output directory. Default: combine_and_extract_imgt_freqs/ "
                             "next to this script.")
    args = parser.parse_args()

    here = Path(__file__).parent
    outdir = Path(args.outdir) if args.outdir else here / "combine_and_extract_imgt_freqs"
    outdir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    for sample in SAMPLES:
        make_file(sample, rng, outdir)


if __name__ == "__main__":
    main()
