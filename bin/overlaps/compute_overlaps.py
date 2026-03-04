import polars as pl
import argparse
import numpy as np
from scipy import sparse

def compute_metrics(input_file, output_file, sample_start_column, sample_end_column=None, id_columns=None):

    df = pl.read_csv(input_file, separator="\t")

    if sample_end_column is None:
        sample_end_column = len(df.columns) - 1

    sample_columns = df.columns[sample_start_column:sample_end_column + 1]

    if id_columns:
        # Group by id_columns and summarize only sample columns
        df = df.group_by(id_columns).agg(
            [pl.sum(col).alias(col) for col in sample_columns]
        )

    if not id_columns:
    # Use all columns before the sample columns as id_columns
        id_columns = df.columns[:sample_start_column]

    # Convert to numpy
    X = df.select(sample_columns).to_numpy()

    # Build sparse count matrix
    X_sparse = sparse.csr_matrix(X)

    # Binary presence matrix (sparse)
    B_sparse = X_sparse.copy()
    B_sparse.data = np.ones_like(B_sparse.data)

    # ---- Sparse matrix multiplications ----

    # Unique overlaps
    unique_overlap = (B_sparse.T @ B_sparse).toarray()

    # Counts restricted to shared rows
    AB_counts = (X_sparse.T @ B_sparse).toarray()

    count_per_sample = np.array(X_sparse.sum(axis=0)).ravel()
    unique_per_sample = np.array(B_sparse.sum(axis=0)).ravel()

    n_samples = len(sample_columns)

    results = []
    for i in range(n_samples):
        for j in range(i + 1, n_samples):

            # Compute Jaccard index
            inter = unique_overlap[i, j]
            union = unique_per_sample[i] + unique_per_sample[j] - inter
            jaccard = inter / union if union > 0 else 0

            # Generate short_id by joining the first 4 letters of all id_columns
            short_id = "_".join([col[:4].upper() for col in id_columns])

            results.append({
                "A": sample_columns[i],
                "B": sample_columns[j],
                "unique_AB": int(unique_overlap[i, j]),
                "unique_A": int(unique_per_sample[i]),
                "unique_B": int(unique_per_sample[j]),
                "AB_countA": int(AB_counts[i, j]),
                "AB_countB": int(AB_counts[j, i]),
                "count_A": int(count_per_sample[i]),
                "count_B": int(count_per_sample[j]),
                "jaccard_index": jaccard,
                "short_id": short_id
            })

    pl.DataFrame(results).write_csv(output_file, separator="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--sample_start_column", type=int, required=True)
    parser.add_argument("--sample_end_column", type=int, default=None)
    parser.add_argument(
        "--id_columns",
        nargs="+",
        help="Columns to group by and summarize before processing",
    )

    args = parser.parse_args()

    compute_metrics(
        args.input_file,
        args.output_file,
        args.sample_start_column,
        args.sample_end_column,
        args.id_columns,
    )
