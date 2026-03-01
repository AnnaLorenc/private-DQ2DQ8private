#!/usr/bin/env python3
"""
make_test_sharing.py
--------------------
Generate toy sharing files to test summarise_sharing.py.

Produces two files:
  test_sharing/groupA_sharing.tsv  -- 30 clonotypes, 10 samples, K=15 rounds
  test_sharing/groupB_sharing.tsv  -- same clonotypes, 8 samples, K=15 rounds

Design:
  - Clonotypes 1-5:   always public (high sharing in both groups)
  - Clonotypes 6-10:  public in group A only
  - Clonotypes 11-15: public in group B only
  - Clonotypes 16-30: private (sharing 0 or 1 only)
"""

import os
import numpy as np
import polars as pl

RNG = np.random.default_rng(seed=99)

# ---- Clonotype definitions ----
clonotypes = [
    # aminoAcid          vFamilyName   jGeneName
    ("CASSLAPGATNEKLFF",  "TCRBV06",   "TCRBJ01-04"),  # 1  always public
    ("CASSLGQGDTQYF",     "TCRBV06",   "TCRBJ02-03"),  # 2  always public
    ("CASRTGQETQYF",      "TCRBV05",   "TCRBJ02-05"),  # 3  always public
    ("CASSQETQYF",        "TCRBV06",   "TCRBJ02-05"),  # 4  always public
    ("CASSWDTGELFF",      "TCRBV06",   "TCRBJ02-02"),  # 5  always public
    ("CASSLDRGSYNEQFF",   "TCRBV06",   "TCRBJ02-01"),  # 6  A-only public
    ("CASSLDRGSYEQYF",    "TCRBV06",   "TCRBJ02-07"),  # 7  A-only public
    ("CASSPGTAYEQYF",     "TCRBV06",   "TCRBJ02-07"),  # 8  A-only public
    ("CASSSGRSYNEQFF",    "TCRBV06",   "TCRBJ02-01"),  # 9  A-only public
    ("CASSTRAGNTEAFF",    "TCRBV06",   "TCRBJ01-01"),  # 10 A-only public
    ("CASSLGQGNTEAFF",    "TCRBV07",   "TCRBJ01-01"),  # 11 B-only public
    ("CASSLSGQETQYF",     "TCRBV07",   "TCRBJ02-05"),  # 12 B-only public
    ("CASSLVGGSNQPQHF",   "TCRBV07",   "TCRBJ01-05"),  # 13 B-only public
    ("CASSQGGDTQYF",      "TCRBV07",   "TCRBJ02-03"),  # 14 B-only public
    ("CASSWGQGDTQYF",     "TCRBV07",   "TCRBJ02-03"),  # 15 B-only public
    ("CASRGDTQYF",        "TCRBV09",   "TCRBJ02-03"),  # 16 private
    ("CASSADGYTF",        "TCRBV09",   "TCRBJ01-02"),  # 17 private
    ("CASSAGQNTEAFF",     "TCRBV09",   "TCRBJ01-01"),  # 18 private
    ("CASSARETQYF",       "TCRBV09",   "TCRBJ02-05"),  # 19 private
    ("CASSDRTNTEAFF",     "TCRBV09",   "TCRBJ01-01"),  # 20 private
    ("CASSGDRDTQYF",      "TCRBV10",   "TCRBJ02-03"),  # 21 private
    ("CASSGQGNTEAFF",     "TCRBV10",   "TCRBJ01-01"),  # 22 private
    ("CASSHGQETQYF",      "TCRBV10",   "TCRBJ02-05"),  # 23 private
    ("CASSLGDTQYF",       "TCRBV10",   "TCRBJ02-03"),  # 24 private
    ("CASSRGTNTEAFF",     "TCRBV10",   "TCRBJ01-01"),  # 25 private
    ("CASSPGTGELFF",      "TCRBV12",   "TCRBJ02-02"),  # 26 private
    ("CASSRLGDTQYF",      "TCRBV12",   "TCRBJ02-03"),  # 27 private
    ("CASSSDRGNTEAFF",    "TCRBV12",   "TCRBJ01-01"),  # 28 private
    ("CASSTRGQETQYF",     "TCRBV12",   "TCRBJ02-05"),  # 29 private
    ("CASSWGGNTEAFF",     "TCRBV12",   "TCRBJ01-01"),  # 30 private
]

K = 15       # rounds
nA = 10      # samples in group A
nB = 8       # samples in group B
n_clono = len(clonotypes)


def make_sharing_matrix(n_clono: int, n_samples: int, K: int, profile: str) -> np.ndarray:
    """
    Generate a (n_clono, K) integer sharing matrix.

    profile: one of 'high', 'A_only', 'B_only', 'private' — controls the
    distribution of per-draw sharing counts.
    """
    mat = np.zeros((n_clono, K), dtype=np.int32)
    for i in range(n_clono):
        p = profile[i]
        for k in range(K):
            if p == "high":
                mat[i, k] = RNG.integers(3, n_samples + 1)
            elif p == "A_only" and "A" in locals().get("_current_group", "A"):
                mat[i, k] = RNG.integers(2, n_samples + 1)
            elif p == "B_only" and "B" in locals().get("_current_group", "B"):
                mat[i, k] = RNG.integers(2, n_samples + 1)
            else:
                mat[i, k] = RNG.integers(0, 2)  # 0 or 1
    return mat


def make_group_matrix(profiles: list[str], n_samples: int) -> np.ndarray:
    """Generate (n_clono, K) sharing matrix from a per-clonotype profile list."""
    mat = np.zeros((n_clono, K), dtype=np.int32)
    for i, p in enumerate(profiles):
        for k in range(K):
            if p == "high":
                mat[i, k] = int(RNG.integers(max(2, n_samples // 2), n_samples + 1))
            elif p == "group":
                mat[i, k] = int(RNG.integers(2, n_samples + 1))
            else:  # private
                mat[i, k] = int(RNG.integers(0, 2))
    return mat


# Profiles per clonotype for each group
profiles_A = (
    ["high"] * 5      # 1-5: always public
    + ["group"] * 5   # 6-10: A-only public
    + ["private"] * 5 # 11-15: B-only (private in A)
    + ["private"] * 15 # 16-30: private
)

profiles_B = (
    ["high"] * 5      # 1-5: always public
    + ["private"] * 5 # 6-10: A-only (private in B)
    + ["group"] * 5   # 11-15: B-only public
    + ["private"] * 15 # 16-30: private
)

matrix_A = make_group_matrix(profiles_A, nA)
matrix_B = make_group_matrix(profiles_B, nB)

index_df = pl.DataFrame({
    "aminoAcid":   [c[0] for c in clonotypes],
    "vFamilyName": [c[1] for c in clonotypes],
    "jGeneName":   [c[2] for c in clonotypes],
})

k_cols_A = {f"k{k+1}": matrix_A[:, k].tolist() for k in range(K)}
k_cols_B = {f"k{k+1}": matrix_B[:, k].tolist() for k in range(K)}

df_A = index_df.with_columns([pl.Series(name, vals) for name, vals in k_cols_A.items()])
df_B = index_df.with_columns([pl.Series(name, vals) for name, vals in k_cols_B.items()])

os.makedirs("test_sharing", exist_ok=True)
df_A.write_csv("test_sharing/groupA_sharing.tsv", separator="\t")
df_B.write_csv("test_sharing/groupB_sharing.tsv", separator="\t")

print("Generated test sharing files:")
print(f"  test_sharing/groupA_sharing.tsv  ({len(df_A)} clonotypes, {K} rounds, {nA} max sharing)")
print(f"  test_sharing/groupB_sharing.tsv  ({len(df_B)} clonotypes, {K} rounds, {nB} max sharing)")
print()
print("Preview of groupA_sharing.tsv (first 8 rows, first 6 k-cols):")
print(df_A.head(8).select(["aminoAcid","vFamilyName","jGeneName","k1","k2","k3","k4","k5","k6"]))
