#!/usr/bin/env python
"""
This script processes immunoseq data by removing specified sequences, splitting the data into productive and non-productive sequences,
and generating a cleanup summary. The script takes the following command-line arguments:

1. --seqs_to_remove: Path to the file containing sequences to remove.
2. --output_loc: Directory to save the cleaned results.
3. --sample: Sample name.
4. --file_path: Path to the sample file.

The outputs include:
- Productive sequences file
- Non-productive sequences file
- Cleanup summary file
"""

def main():
    import argparse
    import pandas as pd
    from pathlib import Path

    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Process immunoseq data.")
    parser.add_argument("--seqs_to_remove", type=str, required=True, help="Path to the file containing sequences to remove.")
    parser.add_argument("--output_loc", type=str, required=True, help="Directory to save the cleaned results.")
 #   parser.add_argument("--samplesheet", type=str, required=True, help="Path to the samplesheet CSV file.")
    parser.add_argument("--sample", type=str, required=True, help="sample name")
    parser.add_argument("--file_path", type=str, required=True, help="Path to the sample file")
    args = parser.parse_args()

    # Assign arguments to variables
    seqs_to_remove = args.seqs_to_remove
    output_loc = args.output_loc
 #   samplesheet_path = args.samplesheet
    filepath = args.file_path
    sample_short = args.sample

    # Read the samplesheet
   # samplesheet = pd.read_csv(samplesheet_path)
   # row0 = samplesheet.iloc[0]
  #  filepath = Path(row0['LOC']) / row0['SAMPLE']
  #  sample_short = row0['sample_short']

    # Read the sample data
    sample_df = pd.read_csv(filepath, sep='\t', compression='gzip', low_memory=False)
    sample_df['sample_short'] = sample_short

    # Read in the file seqs_to_remove, which contains a list of CDR3 sequences to remove from the data
    seqs_to_remove_df = pd.read_csv(seqs_to_remove, sep='\t', names=['nucleotide', 'aminoAcid'], header=None)

    # Remove sequences by matching their nucleotide and aminoAcid
    merged_df = sample_df.merge(seqs_to_remove_df, on=['nucleotide', 'aminoAcid'], how='left', indicator=True)
    removed_sequences = merged_df[merged_df['_merge'] == 'both']
    remaining_sequences = merged_df[merged_df['_merge'] == 'left_only']

    # Write info about how many sequences were found and removed
    num_removed_rows = len(removed_sequences)
    removed_counts_sum = removed_sequences['count (templates/reads)'].sum()
    print(f"Number of sequences removed: {num_removed_rows}")
    print(f"Sum of counts removed: {removed_counts_sum}")

    # Take only rows where jGeneName and vFamilyName are not NA
    filtered_df = remaining_sequences[~remaining_sequences[['jGeneName', 'vFamilyName']].isna().any(axis=1)]

    # Split into productive and non-productive sequences
    productive_df = filtered_df[filtered_df['sequenceStatus'] == 'In']
    nonproductive_df = filtered_df[filtered_df['sequenceStatus'] != 'In']

    # Save the outputs
    Path(output_loc).mkdir(parents=True, exist_ok=True)
    productive_output_path = f"{output_loc}/{sample_short}_productive.tsv.gz"
    nonproductive_output_path = f"{output_loc}/{sample_short}_nonproductive.tsv.gz"
    productive_df.to_csv(productive_output_path, sep='\t', index=False, compression='gzip')
    nonproductive_df.to_csv(nonproductive_output_path, sep='\t', index=False, compression='gzip')

    print(f"Productive sequences saved to: {productive_output_path}")
    print(f"Non-productive sequences saved to: {nonproductive_output_path}")

    # Create a cleanup summary DataFrame
    cleanup_summary = pd.DataFrame({
        'Category': ['Removed Sequences', 'Productive Sequences', 'Non-Productive Sequences'],
        'Num Sequences': [num_removed_rows, len(productive_df), len(nonproductive_df)],
        'Sum of Counts': [removed_counts_sum, productive_df['count (templates/reads)'].sum(), nonproductive_df['count (templates/reads)'].sum()]
    })

    # Save the cleanup summary
    cleanup_summary_path = f"{output_loc}/{sample_short}_cleanup_summary.tsv"
    cleanup_summary.to_csv(cleanup_summary_path, sep='\t', index=False)
    print(f"Cleanup summary saved to: {cleanup_summary_path}")

if __name__ == "__main__":
    main()