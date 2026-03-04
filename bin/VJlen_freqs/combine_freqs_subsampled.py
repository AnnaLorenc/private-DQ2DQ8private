import polars as pl
import os
import argparse

def combine_freqs_subsampled(input_files, output_file):
    """
    Combine frequency files by joining them on 'param' and 'group'.

    Parameters:
        input_files (list of str): List of input file paths.
        output_file (str): Path to save the combined output file.
    """
    
    dfs = []
    for file in input_files:
        # Extract the sample name from the '_orig' column
        sample_name = None
        with open(file, 'r') as f:
            header = f.readline().strip().split("\t")
            sample_name = [col.replace("_orig", "") for col in header if col.endswith("_orig")][0]

        # Read the file
        df = pl.read_csv(file, separator="\t")

        # Rename the 'mean_subsamples' column to '<SAMPLE>_subs'
        df = df.rename({"mean_subsamples": f"{sample_name}_subs"})

        # Merge the dataframes

        dfs.append(df)

    combined_df = pl.concat(dfs, how="align")
    

    # Reorder columns: param, group, then _orig columns, then _subs columns
    orig_columns = [col for col in combined_df.columns if col.endswith("_orig")]
    subs_columns = [col for col in combined_df.columns if col.endswith("_subs")]
    combined_df = combined_df.select(["param", "group"] + orig_columns + subs_columns)

    

    # Write the combined dataframe to the output file
    combined_df.write_csv(output_file, separator="\t")


def main():
    parser = argparse.ArgumentParser(description="Combine frequency files by joining them on 'param' and 'group'.")
    parser.add_argument("--input_files", nargs='+', required=True, help="List of input file paths.")
    parser.add_argument("--output_file", required=True, help="Path to save the combined output file.")

    args = parser.parse_args()

    combine_freqs_subsampled(args.input_files, args.output_file)


if __name__ == "__main__":
    main()