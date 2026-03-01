#!/usr/bin/env python3
"""
Generates a small test TSV for rare_publicity.py.
30 clonotypes, 10 samples (5 ending in E, 5 in N).
A few "public" clonotypes are seeded with counts in most samples.
"""

import numpy as np
import pandas as pd

np.random.seed(99)

samples_E = ["S01E", "S02E", "S03E", "S04E", "S05E"]
samples_N = ["S06N", "S07N", "S08N", "S09N", "S10N"]
all_samples = samples_E + samples_N

v_families = ["TRBV12-3", "TRBV20-1", "TRBV7-2", "TRBV5-1", "TRBV6-5"]
j_genes    = ["TRBJ1-1", "TRBJ1-2", "TRBJ2-1", "TRBJ2-3", "TRBJ2-7"]

# Generate 30 CDR3 amino acid sequences (length 12-16)
def rand_cdr3(length):
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return "C" + "".join(np.random.choice(list(aa), length - 2)) + "F"

clonotypes = []
for i in range(30):
    clonotypes.append({
        "aminoAcid":   rand_cdr3(np.random.randint(12, 17)),
        "vFamilyName": np.random.choice(v_families),
        "jGeneName":   np.random.choice(j_genes),
    })

df = pd.DataFrame(clonotypes)

# Add sample count columns - mostly sparse (rare clonotypes)
for s in all_samples:
    df[s] = np.random.choice([0, 0, 0, 1, 2, 5, 10], size=30)

# Make 3 "public-E" clonotypes: present in all E samples with higher counts
for i in [0, 1, 2]:
    for s in samples_E:
        df.loc[i, s] = np.random.randint(20, 100)

# Make 2 "public-N" clonotypes: present in all N samples
for i in [3, 4]:
    for s in samples_N:
        df.loc[i, s] = np.random.randint(20, 100)

# Make 1 "fully public" clonotype: present in all samples
for s in all_samples:
    df.loc[5, s] = np.random.randint(15, 50)

df.to_csv("test_data.tsv", sep="\t", index=False)
print(f"Written test_data.tsv ({len(df)} rows, {len(all_samples)} samples)")
print(df.head(8).to_string())
