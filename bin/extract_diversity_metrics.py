#!/usr/bin/env python
"""
This script processes an ImmunoSEQ file to calculate diversity metrics and save the results as a CSV file.

Inputs:
1. --file_path: Path to the ImmunoSEQ file.
2. --sample: Sample name.
3. --output_loc: Directory to save the diversity metrics CSV file.

Outputs:
- A CSV file containing the diversity metrics for the given sample.
"""

import argparse
import polars as pl
import lymphoseq as ls
from pathlib import Path
import sys
sys.path.append('/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/assets')
import changed_diversity

# changed_diversity is my version of changed_diversity from lymphoseq, which is a bit more flexible and can handle polars DataFrames. It also includes the shannon_diversity function, which I implemented separately to ensure it works with polars DataFrames. The compute_diversity_metrics function combines the results from changed_diversity.diversity_metrics with the Shannon diversity calculation, returning a single DataFrame with all the metrics.


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

def compute_diversity_metrics(df: pl.DataFrame) -> pl.DataFrame:
    def shannon_diversity(df: pl.DataFrame, count_col="duplicate_count") -> float:
        if count_col in df.columns:
            total_count = df[count_col].sum()
            shannon_div = (
                -df
                .with_columns((pl.col(count_col) / total_count).alias("proportion"))
                .with_columns((pl.col("proportion") * pl.col("proportion").log()).alias("shannon_term"))
                ["shannon_term"]
                .sum()
            )
            return shannon_div
        return 0.0

    # Calculate Shannon diversity as a single number
    shannon_div = shannon_diversity(df)

    # Calculate other diversity metrics
    diversity = changed_diversity.diversity_metrics(df)
    isinstance(diversity, pl.DataFrame)

    # Add shannon_div as a new column to diversity (handles polars DataFrame, pandas DataFrame/Series, or scalars)
    diversity = diversity.with_columns(pl.lit(shannon_div).alias("shannon_diversity"))
    return diversity

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
    airr_df = immunoseq_to_airr(data)

    # Calculate diversity metrics
    print("Calculating diversity metrics...")
    diversity = compute_diversity_metrics(airr_df)

    # Save the diversity metrics as a CSV file
    Path(output_loc).mkdir(parents=True, exist_ok=True)
    output_file = Path(output_loc) / f"{sample}_diversity.csv"
    diversity.write_csv(output_file)

    print(f"Diversity metrics saved to: {output_file}")

if __name__ == "__main__":
    main()