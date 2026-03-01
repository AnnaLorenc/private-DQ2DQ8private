"""
Script: imgt_freqs.py

Description:
This script processes CDR3 sequences to compute summary and frequency tables based on IMGT numbering. It filters valid sequences based on length, starting and ending amino acids, and valid amino acid patterns. The script outputs valid and invalid sequences, as well as summary and frequency tables.

Expected Inputs:
- A tab-separated input file containing columns for CDR3 sequences, VFamily, Jgene, and sample-specific counts (in the current usage SAMPLE_orig, SAMPLE_subsample1 SAMPLE_subsample_2...).
VFamily and Jgene are NOT used for grouping

Outputs, in a directory XXXX_subsampled_aa_freqs:
SAMPLE_subsampled_valid.tsv
- A file containing valid CDR3 sequences, exploaded to every aa within every CDR3.

SAMPLE_subsampled_invalid.tsv
- A file containing invalid CDR3 sequences. XXXX_subsampled_invalid.tsv - just rows from the input file. with invalid sequences (mostly length now)

NOT ANYMORE
SAMPLE_subsampled_summ.tsv
- A summary table with aggregated counts and rows grouped by amino acid, IMGT position, and CDR3 length.
# AA      IMGT_position   aminoAcid_length  F15625E_orig_counts     F15625E_subsample_1_counts      F15625E_subsample_2_counts...   F15625E_orig_rows       F15625E_subsample_1_rows ....q


SAMPLE_subsampled_freq.tsv
- A frequency table with normalized counts and rows grouped by amino acid and IMGT position and CDR3 length.
# AA      IMGT_position   aminoAcid_length  F15625E_orig_counts     F15625E_subsample_1_counts      F15625E_subsample_2_counts...   F15625E_orig_rows       F15625E_subsample_1_rows ....q
# + normalised

SAMPLE_subsampled_freq_withoutLength.tsv
- A frequency table with normalized counts and rows grouped by amino acid and IMGT position (filtering for length still happens!).
# AA      IMGT_position  F15625E_orig_counts     F15625E_subsample_1_counts      F15625E_subsample_2_counts...   F15625E_orig_rows       F15625E_subsample_1_rows ....q
# + normalised

Usage:
python imgt_freqs.py --input <input_file> --summary_output <summary_file> --freq_output <freq_file> \
                     --valid_output <valid_file> --invalid_output <invalid_file> [--min_len <min_length>] [--max_len <max_length>]

Arguments:
--input: Path to the input file (required).
--summary_output: Path to the summary output file (required).
--freq_output: Path to the frequency output file (required).
--valid_output: Path to save valid CDR3 sequences (required).
--invalid_output: Path to save invalid CDR3 sequences (required).
--min_len: Minimum length of valid CDR3 sequences (default: 8).
--max_len: Maximum length of valid CDR3 sequences (default: 25).
"""

#!/usr/bin/env python3

import argparse
import polars as pl
import re
import os


VALID_AA_REGEX = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


# --------------------------------------------------
# IMGT numbering (vector-safe helper)
# --------------------------------------------------
def imgt_number_trb(seq: str):

    positions = []
    positions.append(("104", seq[0]))

    core = seq[1:-1]
    core_len = len(core)

    canonical = 13

    left = ["105","106","107","108","109","110","111"]
    right = ["112","113","114","115","116","117"]

    if core_len == canonical:
        core_positions = left + right

    elif core_len > canonical:
        extra = core_len - canonical
        right_extra = extra // 2
        left_extra = extra - right_extra
        left_insertions = [f"111{chr(65+i)}" for i in range(left_extra)]
        right_insertions = [f"112{chr(65+right_extra-1-i)}" for i in range(right_extra)]
        core_positions = left + left_insertions + right_insertions + right

    else:
        deficit = canonical - core_len
        print(f"deficit: {deficit}")
        remove_right = deficit // 2
        
        remove_left = deficit - remove_right
        # print(f"remove_left: {remove_left}")
        # print(f"remove_right: {remove_right}")
        trimmed_left = left[:len(left)-remove_left]  # Remove last remove_left elements
        # print(f"trimmed_left: {trimmed_left}")
        trimmed_right = right[remove_right:]
        # print(f"trimmed_right: {trimmed_right}")
        core_positions = trimmed_left + trimmed_right
        # print(f"core_positions: {core_positions}")    

    for pos, aa in zip(core_positions, core):
        positions.append((pos, aa))

    positions.append(("118", seq[-1]))

    return positions


# --------------------------------------------------
# Main
# --------------------------------------------------
def process(input_file, output_base, min_len, max_len):

    # Derive base path and output directory
    base_root = os.path.splitext(output_base)[0]
    output_dir = f"{base_root}_aa_freqs"
    os.makedirs(output_dir, exist_ok=True)

    file_stem = os.path.basename(base_root)

    # summary_out = os.path.join(output_dir, f"{file_stem}_summ.tsv")
    freq_out = os.path.join(output_dir, f"{file_stem}_freq.tsv")
    freq_out_Vfam = os.path.join(output_dir, f"{file_stem}_freq_Vfam.tsv")
    freq_out_withoutLength = os.path.join(output_dir, f"{file_stem}_freq_WL.tsv")
    freq_out_withoutLength_Vfam = os.path.join(output_dir, f"{file_stem}_freq_WL_Vfam.tsv")
    valid_out = os.path.join(output_dir, f"{file_stem}_valid.tsv")
    invalid_out = os.path.join(output_dir, f"{file_stem}_invalid.tsv")

    df = pl.read_csv(input_file, separator="\t")

    sample_cols = [c for c in df.columns if c not in ["aminoAcid", "vFamilyName", "jGeneName"]]

    total = df.height

    df = (
        df
        .with_columns(
            pl.col("aminoAcid")
            .str.strip_chars()
            .str.to_uppercase()
            .alias("aminoAcid")
        )
        .with_columns(
            (
                (pl.col("aminoAcid").str.len_chars() >= min_len) &
                (pl.col("aminoAcid").str.len_chars() <= max_len) &
                (pl.col("aminoAcid").str.starts_with("C")) &
                (
                    pl.col("aminoAcid").str.ends_with("F") |
                    pl.col("aminoAcid").str.ends_with("W") |
                    pl.col("aminoAcid").str.ends_with("YV")
                ) &
                (pl.col("aminoAcid").str.contains(VALID_AA_REGEX))
            ).alias("valid")
        )
    )

    invalid_df = df.filter(~pl.col("valid"))
    invalid_df.write_csv(invalid_out, separator="\t")

    df = df.filter(pl.col("valid"))

    print("Total sequences:", total)
    print("Filtered out:", invalid_df.height)
    print("Retained:", df.height)

    # --------------------------------------------------
    # Explode into IMGT positions
    # --------------------------------------------------
    exploded = []

    for row in df.iter_rows(named=True):

        seq = row["aminoAcid"]
        seq_len = len(seq)
        positions = imgt_number_trb(seq)

        for pos, aa in positions:
            record = {
                "aminoAcid": seq,
                "AA": aa,
                "IMGT_position": pos,
                "aminoAcid_length": seq_len,
                "vFamilyName": row["vFamilyName"],
            }

            for s in sample_cols:
                record[s] = row[s]

            exploded.append(record)

    exploded_df = pl.DataFrame(exploded)
    exploded_df.write_csv(valid_out, separator="\t")


    # --------------------------------------------------
    # SUMMARY TABLE
    # --------------------------------------------------
    # summary = (
    #     exploded_df
    #     .group_by(["AA", "IMGT_position", "aminoAcid_length"])
    #     .agg(
    #         [
    #             pl.col(s).sum().alias(f"{s}_counts")
    #             for s in sample_cols
    #         ] +
    #         [
    #             (pl.col(s) > 0)
    #             .cast(pl.Int32)
    #             .sum()
    #             .alias(f"{s}_rows")
    #             for s in sample_cols
    #         ]
    #     )
    # )

    # summary.write_csv(summary_out, separator="\t")

    # --------------------------------------------------
    # FREQUENCY TABLE (with length)
    # --------------------------------------------------
    freq_counts = (
        exploded_df
        .group_by(["AA", "IMGT_position", "aminoAcid_length"])
        .agg(
            [
                pl.col(s).sum().alias(f"{s}_counts")
                for s in sample_cols
            ] +
            [
                (pl.col(s) > 0)
                .cast(pl.Int32)
                .sum()
                .alias(f"{s}_rows")
                for s in sample_cols
            ]
        )
    )

    totals = (
        freq_counts
        .group_by(["IMGT_position", "aminoAcid_length"])
        .agg(
            [
                pl.col(f"{s}_counts").sum().alias(f"{s}_counts_total")
                for s in sample_cols
            ] +
            [
                pl.col(f"{s}_rows").sum().alias(f"{s}_rows_total")
                for s in sample_cols
            ]
        )
    )

    freq = freq_counts.join(totals, on=["IMGT_position", "aminoAcid_length"])

    # compute normalized frequencies
    for s in sample_cols:
        freq = freq.with_columns(
            pl.when(pl.col(f"{s}_counts_total") == 0)
            .then(0)
            .otherwise(pl.col(f"{s}_counts") / pl.col(f"{s}_counts_total"))
            .alias(f"{s}_counts_freq")
        ).with_columns(
            pl.when(pl.col(f"{s}_rows_total") == 0)
            .then(0)
            .otherwise(pl.col(f"{s}_rows") / pl.col(f"{s}_rows_total"))
            .alias(f"{s}_rows_freq")
        )
    # remove total columns
    drop_cols = [c for c in freq.columns if c.endswith("_total")]
    freq = freq.drop(drop_cols)

    freq.write_csv(freq_out, separator="\t")

    # --------------------------------------------------
    # FREQUENCY TABLE (with length, including vFamilyName)
    # --------------------------------------------------
    freq_counts_Vfam = (
        exploded_df
        .group_by(["AA", "IMGT_position", "aminoAcid_length", "vFamilyName"])
        .agg(
            [
                pl.col(s).sum().alias(f"{s}_counts")
                for s in sample_cols
            ] +
            [
                (pl.col(s) > 0)
                .cast(pl.Int32)
                .sum()
                .alias(f"{s}_rows")
                for s in sample_cols
            ]
        )
    )

    totals_Vfam = (
        freq_counts_Vfam
        .group_by(["IMGT_position", "aminoAcid_length", "vFamilyName"])
        .agg(
            [
                pl.col(f"{s}_counts").sum().alias(f"{s}_counts_total")
                for s in sample_cols
            ] +
            [
                pl.col(f"{s}_rows").sum().alias(f"{s}_rows_total")
                for s in sample_cols
            ]
        )
    )

    freq_Vfam = freq_counts_Vfam.join(totals_Vfam, on=["IMGT_position", "aminoAcid_length", "vFamilyName"])

    for s in sample_cols:
        freq_Vfam = freq_Vfam.with_columns(
            pl.when(pl.col(f"{s}_counts_total") == 0)
            .then(0)
            .otherwise(pl.col(f"{s}_counts") / pl.col(f"{s}_counts_total"))
            .alias(f"{s}_counts_freq")
        ).with_columns(
            pl.when(pl.col(f"{s}_rows_total") == 0)
            .then(0)
            .otherwise(pl.col(f"{s}_rows") / pl.col(f"{s}_rows_total"))
            .alias(f"{s}_rows_freq")
        )
    drop_cols = [c for c in freq_Vfam.columns if c.endswith("_total")]
    freq_Vfam = freq_Vfam.drop(drop_cols)

    freq_Vfam.write_csv(freq_out_Vfam, separator="\t")

    # --------------------------------------------------
    # FREQUENCY TABLE (without length)
    # --------------------------------------------------
    freq_counts_wl = (
        exploded_df
        .group_by(["AA", "IMGT_position"])
        .agg(
            [
                pl.col(s).sum().alias(f"{s}_counts")
                for s in sample_cols
            ] +
            [
                (pl.col(s) > 0)
                .cast(pl.Int32)
                .sum()
                .alias(f"{s}_rows")
                for s in sample_cols
            ]
        )
    )

    totals_wl = (
        freq_counts_wl
        .group_by("IMGT_position")
        .agg(
            [
                pl.col(f"{s}_counts").sum().alias(f"{s}_counts_total")
                for s in sample_cols
            ] +
            [
                pl.col(f"{s}_rows").sum().alias(f"{s}_rows_total")
                for s in sample_cols
            ]
        )
    )

    freq_withoutLength = freq_counts_wl.join(totals_wl, on="IMGT_position")

    for s in sample_cols:
        freq_withoutLength = freq_withoutLength.with_columns(
            (
                pl.col(f"{s}_counts") /
                pl.col(f"{s}_counts_total")
            ).alias(f"{s}_counts_freq")
        ).with_columns(
            (
                pl.col(f"{s}_rows") /
                pl.col(f"{s}_rows_total")
            ).alias(f"{s}_rows_freq")
        )

    drop_cols = [c for c in freq_withoutLength.columns if c.endswith("_total")]
    freq_withoutLength = freq_withoutLength.drop(drop_cols)

    freq_withoutLength.write_csv(freq_out_withoutLength, separator="\t")

    # --------------------------------------------------
    # FREQUENCY TABLE (without length, including vFamilyName)
    # --------------------------------------------------
    freq_counts_wl_Vfam = (
        exploded_df
        .group_by(["AA", "IMGT_position", "vFamilyName"])
        .agg(
            [
                pl.col(s).sum().alias(f"{s}_counts")
                for s in sample_cols
            ] +
            [
                (pl.col(s) > 0)
                .cast(pl.Int32)
                .sum()
                .alias(f"{s}_rows")
                for s in sample_cols
            ]
        )
    )

    totals_wl_Vfam = (
        freq_counts_wl_Vfam
        .group_by(["IMGT_position", "vFamilyName"])
        .agg(
            [
                pl.col(f"{s}_counts").sum().alias(f"{s}_counts_total")
                for s in sample_cols
            ] +
            [
                pl.col(f"{s}_rows").sum().alias(f"{s}_rows_total")
                for s in sample_cols
            ]
        )
    )

    freq_withoutLength_Vfam = freq_counts_wl_Vfam.join(totals_wl_Vfam, on=["IMGT_position", "vFamilyName"])

    for s in sample_cols:
        freq_withoutLength_Vfam = freq_withoutLength_Vfam.with_columns(
            pl.when(pl.col(f"{s}_counts_total") == 0)
            .then(0)
            .otherwise(pl.col(f"{s}_counts") / pl.col(f"{s}_counts_total"))
            .alias(f"{s}_counts_freq")
        ).with_columns(
            pl.when(pl.col(f"{s}_rows_total") == 0)
            .then(0)
            .otherwise(pl.col(f"{s}_rows") / pl.col(f"{s}_rows_total"))
            .alias(f"{s}_rows_freq")
        )

    drop_cols = [c for c in freq_withoutLength_Vfam.columns if c.endswith("_total")]
    freq_withoutLength_Vfam = freq_withoutLength_Vfam.drop(drop_cols)

    freq_withoutLength_Vfam.write_csv(freq_out_withoutLength_Vfam, separator="\t")

# --------------------------------------------------
# CLI
# --------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", dest="output")
    parser.add_argument("--min_len", type=int, default=8)
    parser.add_argument("--max_len", type=int, default=25)

    args = parser.parse_args()

    process(
        args.input,
        args.output,
        args.min_len,
        args.max_len
    )
