import polars as pl
import os
import numpy as np
from itertools import combinations
from scipy import sparse

# ------------------------------
# Step 1: Helpers
# ------------------------------

def sample_prefix(s):
    base = s.split("_subsample_")[0]
    return base[:-1]

def sample_letter(s):
    base = s.split("_subsample_")[0]
    return base[-1]


# ------------------------------
# Step 2: Global clonotype mapping
# ------------------------------

def build_clonotype_mapping(sample_list, subsample_dir, id_columns):
    all_clonotypes = set()

    for s in sample_list:
        filepath = os.path.join(subsample_dir, f"{s}_subsampled.tsv.gz")
        df = pl.read_csv(filepath, separator="\t")

        ids = (
            df.with_columns(
                pl.concat_str(
                    [pl.col(col).str.strip_chars().str.to_uppercase() for col in id_columns],
                    separator="|"
                ).alias("clonotype_id")
            )["clonotype_id"]
            .to_list()
        )

        all_clonotypes.update(ids)

    print("Total unique clonotypes in global mapping:", len(all_clonotypes))

    return {c: i for i, c in enumerate(all_clonotypes)}


# ------------------------------
# Step 3: Load sparse matrix per sample
# ------------------------------

def load_subsamples_matrix_sparse(filepath, clonotype2id, id_columns):
    df = pl.read_csv(filepath, separator="\t")

    # Create clonotype_id based on id_columns
    df = df.with_columns(
        pl.concat_str(
            [pl.col(col).str.strip_chars().str.to_uppercase() for col in id_columns],
            separator="|"
        ).alias("clonotype_id")
    )

    # Collapse rows with identical clonotype_id by summing their counts
    subsample_cols = [c for c in df.columns if "subsample_" in c]
    df = df.group_by("clonotype_id", maintain_order=True).agg(
        [pl.sum(col).alias(col) for col in subsample_cols]
    )

    n_clonotypes = len(clonotype2id)
    n_subs = len(subsample_cols)

    rows = []
    cols = []

    for i, col in enumerate(subsample_cols):

        sub_df = df.filter(pl.col(col) != 0)

 #       print(f"{os.path.basename(filepath)} | {col} nonzero clonotypes:", sub_df.height)

        indices = [
            clonotype2id[c]
            for c in sub_df["clonotype_id"].to_list()
        ]

        rows.extend([i] * len(indices))
        cols.extend(indices)

    data = np.ones(len(rows), dtype=np.uint8)

    matrix = sparse.csr_matrix(
        (data, (rows, cols)),
        shape=(n_subs, n_clonotypes),
        dtype=np.uint8
    )

 #   print("Subsample sizes (row sums):", np.array(matrix.sum(axis=1)).ravel())

    return subsample_cols, matrix


# ------------------------------
# Step 4: Main computation
# ------------------------------

def compute_all_overlaps_sparse(sample_list, subsample_dir, output_file, id_columns):

    samples_E = [s for s in sample_list if sample_letter(s) == "E"]
    samples_N = [s for s in sample_list if sample_letter(s) == "N"]

    prefix_to_E = {}
    prefix_to_N = {}

    for s in samples_E:
        prefix_to_E.setdefault(sample_prefix(s), []).append(s)

    for s in samples_N:
        prefix_to_N.setdefault(sample_prefix(s), []).append(s)

    print(f"Number of E samples: {len(samples_E)}")
    print(f"Number of N samples: {len(samples_N)}")

    clonotype2id = build_clonotype_mapping(sample_list, subsample_dir, id_columns)

    # Load matrices
    all_bitsets = {}
    for s in sample_list:
        filepath = os.path.join(subsample_dir, f"{s}_subsampled.tsv.gz")
        sub_cols, matrix = load_subsamples_matrix_sparse(filepath, clonotype2id, id_columns)
        all_bitsets[s] = (sub_cols, matrix)

    results = []

    def process_pair(s1, s2, comp_type):

        subs1, M1 = all_bitsets[s1]
        subs2, M2 = all_bitsets[s2]

        # 25x25 full overlap matrix
        inter_matrix = (M1 @ M2.T).toarray()

        unique1 = np.array(M1.sum(axis=1)).ravel()
        unique2 = np.array(M2.sum(axis=1)).ravel()

        print(f"\nProcessing {s1} vs {s2}")
        print("Max overlap in matrix:", inter_matrix.max())

        for i in range(len(subs1)):
            for j in range(len(subs2)):

                inter = int(inter_matrix[i, j])
                union = unique1[i] + unique2[j] - inter
                jaccard = inter / union if union > 0 else 0

                # Generate short_id for the pair
                short_id = "_".join([col[:4].upper() for col in id_columns])

                results.append({
                    "sample1": s1,
                    "subsample1": subs1[i],
                    "sample2": s2,
                    "subsample2": subs2[j],
                    "unique_overlap": inter,
                    "unique1": int(unique1[i]),
                    "unique2": int(unique2[j]),
                    "jaccard": float(jaccard),
                    "comparison_type": comp_type,
                    "short_id": short_id
                })

    # ---- Within E
    for s1, s2 in combinations(samples_E, 2):
        process_pair(s1, s2, "within_E")

    # ---- Within N
    for s1, s2 in combinations(samples_N, 2):
        process_pair(s1, s2, "within_N")

    # ---- E_vs_N (1-to-1 by prefix)
    common_prefixes = set(prefix_to_E.keys()) & set(prefix_to_N.keys())

    for prefix in common_prefixes:
        sE = prefix_to_E[prefix][0]
        sN = prefix_to_N[prefix][0]
        process_pair(sE, sN, "E_vs_N")

    pl.DataFrame(results).write_csv(output_file, separator="\t")
    print(f"\nSaved overlaps to {output_file}")


# ------------------------------
# Step 5: CLI
# ------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--samples_file", required=True)
    parser.add_argument("--subsample_dir", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument(
        "--id_columns",
        nargs="+",
        default=["aminoAcid", "vFamilyName", "jGeneName"],
        help="Columns to use for building clonotype IDs (default: aminoAcid, vFamilyName, jGeneName)",
    )

    args = parser.parse_args()

    sample_list = [line.strip() for line in open(args.samples_file)]

    compute_all_overlaps_sparse(
        sample_list,
        args.subsample_dir,
        args.output_file,
        args.id_columns
    )
