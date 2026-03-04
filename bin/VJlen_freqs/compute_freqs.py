#!/usr/bin/env python
"""
This script processes an ImmunoSEQ file to calculate V, J, CDR3len  and save the results as a CSV file.

Inputs:
1. --file_path: Path to the ImmunoSEQ file.
2. --sample: Sample name.
3. --output_loc: Directory to save the diversity metrics CSV file.

Outputs:
- A CSV file containing the diversity metrics for the given sample.
"""

import argparse
import polars as pl
import pandas as pd
import lymphoseq as ls
from pathlib import Path

# Add error handling for required columns
def validate_columns(df, required_columns):
    """Ensure the DataFrame contains the required columns."""

    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

def immunoseq_to_airr(df: pl.DataFrame) -> pl.DataFrame:


    # Define the column mapping
    col_map = {
        "nucleotide": "sequence",
        "aminoAcid": "junction_aa",
        'cdr3Length':'junction_aa_length',
        "vGeneName": "v_call",
        "jGeneName": "j_call",
        "sequenceStatus": "frame_type",
        "count (templates/reads)": "duplicate_count",
        "vFamilyName": "v_family",
        "jFamilyName": "j_family",
        "vGeneAllele": "v_allele",
        "jGeneAllele": "j_allele",
        "dGeneName": "d_gene"
    }

    # Only keep mappings where the column exists in df
    existing_map = {old: new for old, new in col_map.items() if old in df.columns}
    print(f"Existing columns to rename: {list(existing_map.keys())}")

    target_cols = list(existing_map.values())
    cols_to_drop = [col for col in target_cols if col in df.columns]
    cols_to_drop += ["sequence_aa","junction","junction_length","d_call","c_call","productive","vj_in_frame","stop_codon","locus","duplicate_frequency","consensus_count","d_family","cdr1","cdr1_aa","cdr2","cdr2_aa","cdr3","cdr3_aa","clone_id","cell_id"]
    
    if cols_to_drop:
        print(f"Dropping existing columns to avoid conflicts: {cols_to_drop}")
        df = df.drop(cols_to_drop)

    # Perform renaming
    df = df.rename(existing_map)

   # Cast columns with 'length' or 'count' in the name to integers
    int_cols = [col for col in df.columns if "length" in col.lower() or "count" in col.lower()]
    if int_cols:
        print(f"Casting columns to integers: {int_cols}")
        df = df.with_columns([pl.col(col).cast(pl.Int64, strict=False) for col in int_cols])

    # Add a productive column based on frame_type
    if "frame_type" in df.columns:
        df = df.with_columns(
            (pl.col("frame_type") == "In").alias("productive")
        )

    return df

def compute_freqs(df: pd.DataFrame) ->  pd.DataFrame:
    summary_params = ['vFamilyName', 'jFamilyName', 'cdr3Length', 'jGeneName']
    all_summaries = []

    for param in summary_params:
        # Row count and fraction
        param_summary = df[param].value_counts().reset_index()
        param_summary.columns = ['group', 'row_count']
        param_summary['row_fraction'] = param_summary['row_count'] / param_summary['row_count'].sum()

        # Cell count and fraction based on 'count (templates/reads)'
        cell_counts = df.groupby(param)['count (templates/reads)'].sum().reset_index()
        cell_counts.columns = [param, 'cell_count']
        cell_counts['cell_fraction'] = cell_counts['cell_count'] / cell_counts['cell_count'].sum()

        # Merge summaries
        param_summary = param_summary.merge(cell_counts, left_on='group', right_on=param).drop(columns=[param])

        # Add parameter column
        param_summary['param'] = param

        # Append to list
        all_summaries.append(param_summary)

    # Concatenate all summaries
    final_summary = pd.concat(all_summaries, ignore_index=True)

    # Compute summary for grouping by ['vFamilyName', 'jFamilyName', 'cdr3Length']
    grouped_summary = df.groupby(['vFamilyName', 'jFamilyName', 'cdr3Length']).agg(
        row_count=('vFamilyName', 'size'),
        cell_count=('count (templates/reads)', 'sum')
    ).reset_index()

    # Add row_fraction and cell_fraction
    grouped_summary['row_fraction'] = grouped_summary['row_count'] / grouped_summary['row_count'].sum()
    grouped_summary['cell_fraction'] = grouped_summary['cell_count'] / grouped_summary['cell_count'].sum()

    # Add parameter column
    grouped_summary['param'] = 'vFamilyName-jFamilyName-cdr3Length'
    grouped_summary['group'] = (
        grouped_summary['vFamilyName'].astype(str) + "-" +
        grouped_summary['jFamilyName'].astype(str) + "-" +
        grouped_summary['cdr3Length'].astype(str)
    )
    grouped_summary.drop(columns=['vFamilyName', 'jFamilyName', 'cdr3Length'], inplace=True)

    # Append grouped summary to final summary
    final_summary = pd.concat([final_summary, grouped_summary], ignore_index=True)

    return final_summary



def main():
    import argparse
    import polars as pl
    import lymphoseq as ls
    from pathlib import Path

    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Extract diversity metrics from an ImmunoSEQ file.")
    parser.add_argument("--file_path", type=str, required=True, help="Path to the ImmunoSEQ file.")
    parser.add_argument("--sample", type=str, required=True, help="Sample name.")
    parser.add_argument("--output_loc", type=str, required=True, help="Directory to save the diversity metrics CSV file.")
    args = parser.parse_args()

    # Assign arguments to variables
    file_path = args.file_path
    sample = args.sample
    output_loc = args.output_loc

    # Read the ImmunoSEQ file
    print(f"Reading ImmunoSEQ file: {file_path}")
    data = ls.read_immunoseq(file_path)
    pd_df = data if isinstance(data, pd.DataFrame) else data.to_pandas()



    # Update main to validate columns
    required_columns = [
        "nucleotide", "aminoAcid", "cdr3Length", "vGeneName", "jGeneName",
        "sequenceStatus", "count (templates/reads)", "vFamilyName", "jFamilyName"
    ]
    validate_columns(pd_df, required_columns)
    

    # Ensure columns "count (templates/reads)" and "cdr3Length" are integers
    if not pd.api.types.is_integer_dtype(pd_df['count (templates/reads)']):
        pd_df['count (templates/reads)'] = pd_df['count (templates/reads)'].astype(int)

    if not pd.api.types.is_integer_dtype(pd_df['cdr3Length']):
        pd_df['cdr3Length'] = pd_df['cdr3Length'].astype(int)

    # Add detailed logging for each step
    print("Validating required columns...")
    print("Columns validated successfully.")
    print("Processing data...")
    print("Data processing completed.")

    # Calculate diversity metrics
    print("Calculating frequencies...")
    freqs = compute_freqs(pd_df)

    # Save the diversity metrics as a CSV file
    Path(output_loc).mkdir(parents=True, exist_ok=True)
    output_file = Path(output_loc) / f"{sample}_freqs.csv"
    freqs.to_csv(output_file)

    print(f"Frequencies saved to: {output_file}")

if __name__ == "__main__":
    main()