#!/usr/bin/env python3
"""Generate a fake aa_properties_test.py input with clearly separated groups.

hoDQ2   → VHSE1 ≈ +2,  VHSE2 ≈ -1,  VHSE3 ≈ +1,  VHSE4 ≈ -0.5
hoDQ8   → VHSE1 ≈ -2,  VHSE2 ≈ +1,  VHSE3 ≈ -1,  VHSE4 ≈ +0.5
heDQ2DQ8→ intermediate (≈ 0)

With 10 samples per group and these large effect sizes the oo Hotelling
p-value should be < 0.001 at every position × cell combination.

Usage
-----
    python make_test_input.py            # writes test_properties_input.tsv
    python make_test_input.py --seed 99  # reproducible with a different seed
"""
import argparse
import itertools

import numpy as np
import pandas as pd

POSITIONS = ["104", "107", "111"]
CELLS     = ["CD4", "CD8"]

# Property means per genotype  [VHSE1, VHSE2, VHSE3, VHSE4]
MEANS = {
    "hoDQ2":     np.array([ 2.0, -1.0,  1.0, -0.5]),
    "hoDQ8":     np.array([-2.0,  1.0, -1.0,  0.5]),
    "heDQ2DQ8":  np.array([ 0.1, -0.1,  0.1, -0.1]),
}
SD   = 0.35   # small SD → very clear separation
N    = {      # samples per genotype
    "hoDQ2":    10,
    "hoDQ8":    10,
    "heDQ2DQ8":  6,
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed",   type=int, default=42)
    parser.add_argument("--output", default="test_properties_input.tsv")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    rows = []

    patient_counter = itertools.count(1)

    for genotype, mean_vec in MEANS.items():
        n = N[genotype]
        patients = [f"P{next(patient_counter):03d}" for _ in range(n)]
        for pos, cell in itertools.product(POSITIONS, CELLS):
            for pat in patients:
                values = rng.normal(loc=mean_vec, scale=SD)
                rows.append({
                    "IMGT_position": pos,
                    "cells":         cell,
                    "patient":       pat,
                    "sample":        f"{pat}_{cell}_sample",
                    "genotype":      genotype,
                    "VHSE1":         round(values[0], 4),
                    "VHSE2":         round(values[1], 4),
                    "VHSE3":         round(values[2], 4),
                    "VHSE4":         round(values[3], 4),
                })

    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Written {len(df)} rows → {args.output}")
    print(f"  Positions : {POSITIONS}")
    print(f"  Cells     : {CELLS}")
    print(f"  Genotypes : {list(MEANS.keys())} (n={list(N.values())})")
    print(f"  Properties: VHSE1–VHSE4  (SD={SD}, large group separation)")


if __name__ == "__main__":
    main()
