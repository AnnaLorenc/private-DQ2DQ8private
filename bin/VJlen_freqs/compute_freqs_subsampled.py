import polars as pl
import gzip
import os

def compute_frequencies(input_file, index_columns, cdr3_column):
    """
    Compute frequencies for each level of index columns and the CDR3 length.

    Parameters:
        input_file (str): Path to the input file (tsv.gz).
        index_columns (list of str): List of index columns to use as they are.
        cdr3_column (str): Name of the CDR3 column.

    Returns:
        pl.DataFrame: Resulting table with computed frequencies.
    """
    # Read the input file
    with gzip.open(input_file, "rt") as f:
        df = pl.read_csv(f, separator="\t")

    # Check if all required columns are present
    missing_columns = [col for col in index_columns + [cdr3_column] if col not in df.columns]
    if missing_columns:
        print(f"Warning: The following required columns are missing in the input file: {missing_columns}")
        print(f"File header: {df.columns}")
        return None

    df = df.with_columns(
        pl.Series([len(s) if s is not None else 0 for s in df[cdr3_column]]).alias("cdr3_length")
    )


    # Define the new index columns (including the CDR3 length)
    all_index_columns = index_columns + ["cdr3_length"]

    print(f"These are index columns: {all_index_columns}")

    # Extract the sample columns (all columns except the index columns and the raw CDR3 column)
    excluded = set(index_columns) | {cdr3_column, "cdr3_length"}
    sample_columns = [col for col in df.columns if col not in excluded]
    print(f"These are sample columns: {sample_columns}")

    # Prepare the output DataFrame
    results = []

    for sample in sample_columns:
        # Ensure the sample column is cast to a numeric type
        df = df.with_columns(pl.col(sample).cast(pl.Float64))

        # Filter rows where the sample column > 0
        filtered_df = df.filter(pl.col(sample) > 0)

        # Group by the index columns and count occurrences
        grouped = (
            filtered_df
            .group_by(all_index_columns)
            .agg(pl.len().alias("non_zero_count"))
        )

        # Compute the total number of rows with non-zero values for the sample column
        total_non_zero = filtered_df.shape[0]

        # Compute the frequency
        grouped = grouped.with_columns(
            (pl.col("non_zero_count") / total_non_zero).alias("value")
        )

        # Add the sample name as a column
        grouped = grouped.with_columns(pl.lit(sample).alias("sample"))

        # Keep only the index columns, sample, and frequency columns
        grouped = grouped.select(all_index_columns + ["sample", "value"])

        results.append(grouped)

    # Combine results for all sample columns
    result_df = pl.concat(results)

    # Combine index column names and values into two new columns
    result_df = result_df.with_columns(
        pl.lit("_".join(all_index_columns)).alias("param"),
        pl.concat_str([pl.col(col) for col in all_index_columns], separator="_").alias("group")
    )

    # Select only the new columns and the unique sample column
    reshaped = result_df.select(["param", "group", "sample", "value"])

    # Compute frequencies for each index column independently
    independent_results = []

    for index_col in all_index_columns:
        for sample in sample_columns:
            # Ensure the sample column is cast to a numeric type
            df = df.with_columns(pl.col(sample).cast(pl.Float64))

            # Filter rows where the sample column > 0
            filtered_df = df.filter(pl.col(sample) > 0)

            # Group by the current index column and count occurrences
            grouped = (
                filtered_df
                .group_by(index_col)
                .agg(pl.len().alias("non_zero_count"))
            )

            # Compute the total number of rows with non-zero values for the sample column
            total_non_zero = filtered_df.shape[0]

            # Compute the frequency
            grouped = grouped.with_columns(
                (pl.col("non_zero_count") / total_non_zero).alias("value")
            )

            # Add the sample name and index column name as columns
            grouped = grouped.with_columns(
                pl.lit(sample).alias("sample"),
                pl.lit(index_col).alias("param")
            )

            # Rename the index column to "group"
            grouped = grouped.rename({index_col: "group"})

            # Keep only the "param", "group", "sample", and "value" columns
            grouped = grouped.select(["param", "group", "sample", "value"])

            # Ensure consistent types for all columns before concatenation
            grouped = grouped.with_columns(
                pl.col("group").cast(pl.Utf8),
                pl.col("value").cast(pl.Float64),
                pl.col("sample").cast(pl.Utf8),
                pl.col("param").cast(pl.Utf8)
            )

            independent_results.append(grouped)

    # Combine results for all index columns and samples
    independent_df = pl.concat(independent_results)

    # Append the independent frequencies to the reshaped DataFrame
    reshaped = pl.concat([reshaped, independent_df])


    #just v/J
    vj_results = []

    for sample in sample_columns:
        # Ensure the sample column is cast to a numeric type
        df = df.with_columns(pl.col(sample).cast(pl.Float64))

        # Filter rows where the sample column > 0
        filtered_df = df.filter(pl.col(sample) > 0)

        # Group by the index columns and count occurrences
        grouped = (
            filtered_df
            .group_by(index_columns)
            .agg(pl.len().alias("non_zero_count"))
        )

        # Compute the total number of rows with non-zero values for the sample column
        total_non_zero = filtered_df.shape[0]

        # Compute the frequency
        grouped = grouped.with_columns(
            (pl.col("non_zero_count") / total_non_zero).alias("value")
        )

        # Add the sample name as a column
        grouped = grouped.with_columns(pl.lit(sample).alias("sample"))

        # Keep only the index columns, sample, and frequency columns
        grouped = grouped.select(index_columns + ["sample", "value"])

        vj_results.append(grouped)

    # Combine results for all sample columns
    result_df = pl.concat(vj_results)

    # Combine index column names and values into two new columns
    result_df = result_df.with_columns(
        pl.lit("_".join(index_columns)).alias("param"),
        pl.concat_str([pl.col(col) for col in index_columns], separator="_").alias("group")
    )

    # Select only the new columns and the unique sample column
    reshaped_vj = result_df.select(["param", "group", "sample", "value"])

    reshaped = pl.concat([reshaped, reshaped_vj])


    # Pivot the reshaped DataFrame so each sample has its own column
    reshaped = reshaped.pivot(
        values="value",
        index=["param", "group"],
        on="sample"
    )

    return reshaped




def summarize_frequencies(reshaped, output_file_summ):
    """
    Summarize frequencies for each level of param and group.

    Parameters:
        reshaped (pl.DataFrame): The reshaped DataFrame with frequencies.
        output_file_summ (str): Path to save the summarized output file.
    """
    orig_columns = [col for col in reshaped.columns if col.endswith("_orig")]
    subsample_columns = [col for col in reshaped.columns if "_subsample_" in col]

    # Row-wise mean of subsample columns
    if subsample_columns:
        mean_expr = pl.mean_horizontal([pl.col(c) for c in subsample_columns]).fill_null(0.0).alias("mean_subsamples")
    else:
        mean_expr = pl.lit(0.0).alias("mean_subsamples")

    summary = reshaped.select([
        pl.col("param"),
        pl.col("group"),
        *[pl.col(c).fill_null(0.0).alias(c) for c in orig_columns],
        mean_expr
    ])
        # Write the summary to the output file
    summary.write_csv(output_file_summ, separator="\t")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Compute frequencies for each level of index columns and CDR3 length.")
    parser.add_argument("--input_file", required=True, help="Path to the input file (tsv.gz).")
    parser.add_argument("--output_file", required=True, help="Path to save the output file.")
    parser.add_argument("--index_columns", nargs='+', required=True, help="List of index columns to use as they are.")
    parser.add_argument("--cdr3_column", required=True, help="Name of the CDR3 column.")

    args = parser.parse_args()

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    result = compute_frequencies(args.input_file, args.index_columns, args.cdr3_column)

    # Save the reshaped result to the output file
    result.write_csv(args.output_file, separator="\t", compression="gzip")

    # Generate the summary output file name by appending '_summ' to the base name of the output file
    base_name, ext = os.path.splitext(args.output_file)
    output_file_summ = f"{base_name}_summ.tsv"

    # Summarize the reshaped result and save to the summary file
    summarize_frequencies(result, output_file_summ)


if __name__ == "__main__":
    main()