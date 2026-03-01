#!/usr/bin/env python3
"""Generate a fake input TSV for imgt_freqs.py.

Produces a file with:
  - aminoAcid  : valid CDR3 sequences (start with C, end with F, valid AA only)
  - vFamilyName: e.g. TCRBV05
  - jGeneName  : e.g. TCRBJ02-03
  - {SAMPLE}_orig          : integer count (original repertoire)
  - {SAMPLE}_subsample_1..3: integer subsampled counts

Output: F99999E_subsampled_test.tsv  (next to this script)

Usage
-----
    python make_fake_imgt_freqs_input.py
    python make_fake_imgt_freqs_input.py --seed 7 --nrows 20 --nsubsamples 3
"""
import argparse
import random
from pathlib import Path

import pandas as pd

VALID_AA = list("ACDEFGHIKLMNPQRSTVWY")
VFAMILIES = [f"TCRBV{i:02d}" for i in range(1, 8)]
JGENES    = [f"TCRBJ02-{i:02d}" for i in range(1, 8)]
SAMPLE    = "F99999E"


def random_cdr3(rng: random.Random, length: int) -> str:
    """Valid CDR3: starts with C, ends with F, middle = random valid AAs."""
    middle = "".join(rng.choices(VALID_AA, k=length - 2))
    return "C" + middle + "F"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed",        type=int, default=42)
    parser.add_argument("--nrows",       type=int, default=20)
    parser.add_argument("--nsubsamples", type=int, default=3)
    parser.add_argument("--outdir",      default=None,
                        help="Output directory. Default: imgt_freqs/ next to this script.")
    args = parser.parse_args()

    here = Path(__file__).parent
    outdir = Path(args.outdir) if args.outdir else here / "imgt_freqs"
    outdir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)

    rows = []
    for _ in range(args.nrows):
        length  = rng.randint(8, 14)
        cdr3    = random_cdr3(rng, length)
        orig    = rng.randint(1, 10)
        subs    = {f"{SAMPLE}_subsample_{i}": rng.randint(0, orig)
                   for i in range(1, args.nsubsamples + 1)}
        row = {
            "aminoAcid":   cdr3,
            "vFamilyName": rng.choice(VFAMILIES),
            "jGeneName":   rng.choice(JGENES),
            f"{SAMPLE}_orig": orig,
            **subs,
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    out = outdir / f"{SAMPLE}_subsampled_fake.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"Written {len(df)} rows, {args.nsubsamples} subsamples → {out}")
    print(f"Columns: {list(df.columns)}")


if __name__ == "__main__":
    main()
