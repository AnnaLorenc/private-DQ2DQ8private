import polars as pl
import os
import argparse
import gzip

def merge_trb_data(input_files, output_file, merge_columns, count_column):
    """
    Merge TRB data from multiple files into a single file with specified merge and count columns. Created in the context of merging TRB data for the DQ2DQ8 project,
     from all samples into one big db-like file.
    It first collapses rows with the same merge_columns values by summing their count_column values.
The output is a wide matrix with:

Rows: one per unique combination of merge_columns values found across all input files
Columns: len(merge_columns) + len(input_files)
the merge columns come first (e.g. aminoAcid, vFamilyName)
then one column per input sample, named after the first _-delimited part of each filename
So for example, with --merge_columns aminoAcid vFamilyName and 60 input files, the output is N × 62, where N is the number of distinct clonotype definitions seen across all samples.

    Parameters:
        input_files (list): List of file paths to ingest.
        output_file (str): Path to the output file (gzipped).
        merge_columns (list): List of column names to merge on (e.g., ["aminoacid", "VFamily"]).
        count_column (str): Name of the count column in the input files.
    """
 
    total_lines_ingested = 0
    numeric_columns = []

    dfs = []

    for file_index, file_path in enumerate(input_files):
        # Extract the first part of the file name for the sample name
        sample_name = os.path.basename(file_path).split("_")[0]
        numeric_columns.append(sample_name)

        with gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path, "r") as f:
            df = pl.read_csv(
                f,
                separator="\t",
                columns=merge_columns + [count_column]
            )

        # Collapse rows with the same merge_columns values by summing their count_column values
        df = (
            df.group_by(merge_columns)
              .agg(pl.col(count_column).sum())
              .rename({count_column: sample_name})
        )
        
        dfs.append(df)
        n_lines = df.height
        print("in this file ingested:", n_lines)
        total_lines_ingested += n_lines

    # Stack all files vertically
    combined = pl.concat(dfs, how="diagonal")

    # Group by the merge columns and sum the numeric columns - should not be necessary
    aggregated_data = (
        combined
        .group_by(merge_columns)
        .agg(pl.col(pl.NUMERIC_DTYPES).sum())
        .sort(merge_columns)
        .fill_null(0)
)



    # Write the output to a gzipped file
   
    aggregated_data.write_csv(
        output_file,
        separator="\t",
        compression="gzip"
)

    # Print statistics
    print(f"Total lines ingested: {total_lines_ingested}")
    print(f"Total lines in final output: {aggregated_data.shape[0]}")
  #  print(f"Rows with non-zero counts in more than one sample: {non_zero_rows}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Merge TRB data from multiple files into a single file.")
    parser.add_argument("--input_files", nargs="+", help="List of input tab-separated or gzipped files.")
    parser.add_argument("--input_dir", help="Directory containing input tab-separated or gzipped files.")
    parser.add_argument("--output_file", required=True, help="Path to the output tab-separated or gzipped file.")
    parser.add_argument("--merge_columns", required=True, help="Names of the columns to merge on (e.g., aminoacid VFamily).", nargs="+")
    parser.add_argument("--count_column", required=True, help="Name of the count column in the input files.")
    args = parser.parse_args()

    if args.input_files:
        input_files = args.input_files
    elif args.input_dir:
        input_files = [
            os.path.join(args.input_dir, f)
            for f in os.listdir(args.input_dir)
            if f.endswith(".tsv") or f.endswith(".tsv.gz")
        ]
    else:
        raise ValueError("Either --input_files or --input_dir must be provided.")

    output_file = args.output_file
    merge_columns = args.merge_columns
    count_column = args.count_column

    # Merge the data
    merge_trb_data(input_files, output_file, merge_columns, count_column)
    print(f"Merged data written to {output_file}")

if __name__ == "__main__":
    main()