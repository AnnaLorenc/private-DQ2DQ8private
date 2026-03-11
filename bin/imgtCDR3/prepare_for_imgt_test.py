#!/usr/bin/env python3
#prepare_for_imgt_test.py
"""Reshape IMGT AA median frequency tables to long format for testing.

This script:
- reads a wide IMGT AA table such as
    results/imgt_aa_combined/comb_aa_imgt_full_subs_rows_freq.tsv
  where rows are defined by index columns (e.g. IMGT_position, AA,
  aminoAcid_length) and columns after the index are sample-specific
  frequency columns (e.g. F16018E_subs_rows_freq, F20790N_subs_rows_freq, ...),
- optionally filters rows using a matching "_subs_rows" table, keeping only
  index combinations where at least N non-index columns are above a minimum
  value,
- converts the table to a long format with one row per
    (sample, length, AA, IMGT_position, patient, cells),
- joins genotype and sample annotation from a collated_info.csv-style file,
  and
- writes a tab-separated output suitable as input for bin/test_imgt.py.

Example output columns:
    sample  length  AA  IMGT_position   patient cells   value   genotype
    F16018E_subs_rows_freq  8   C   104 F16018  E   1   hoDQ2

Filtering logic (if --filter_file is provided):
- read filter_file (e.g. comb_aa_imgt_full_subs_rows.tsv), which shares the
  same index columns as input_file but has counts instead of frequencies,
- for each row, count how many non-index columns have value > --min_value,
- keep rows where that count >= --min_non_index,
- restrict input_file to those index combinations before reshaping.

The script prints:
- number of rows in the original input_file,
- if filtering is used, number of rows after filtering,
- a message indicating that join with the annotation file is performed.
"""

import argparse
from pathlib import Path
from typing import List

import polars as pl


DEFAULT_INDEX_COLUMNS: List[str] = ["IMGT_position", "AA", "aminoAcid_length"]


def recode_genotype(genotype_short: pl.Series) -> pl.Series:
    """Map genotype_short values to compact labels used for testing.

    Example mappings (adjust as needed for new data):
      - "homoDQ2"      -> "hoDQ2"
      - "homoDQ8"      -> "hoDQ8"
      - "heteroDQ2DQ8" -> "heDQ2DQ8"
    Other values are passed through unchanged.
    """

    return (
        pl.when(genotype_short == "homoDQ2").then("hoDQ2")
        .when(genotype_short == "homoDQ8").then("hoDQ8")
        .when(genotype_short == "heteroDQ2DQ8").then("heDQ2DQ8")
        .otherwise(genotype_short)
    )


def apply_filter(
    df_input: pl.DataFrame,
    filter_path: Path,
    index_columns: List[str],
    min_non_index: int,
    min_value: float,
) -> pl.DataFrame:
    """Filter input rows using a separate filter file.

    Returns a new DataFrame with only rows whose index combination appears in
    filter_file and satisfies the N / min criteria.
    """

    print(f"Reading filter file: {filter_path}")
    df_filter = pl.read_csv(
        filter_path,
        separator="\t",
        schema_overrides={"IMGT_position": pl.Utf8},
    )

    # Ensure all requested index columns are present in the filter file
    missing_idx = [c for c in index_columns if c not in df_filter.columns]
    if missing_idx:
        raise ValueError(
            f"Index columns {missing_idx} not found in filter_file {filter_path}"
        )

    non_index_cols = [c for c in df_filter.columns if c not in index_columns]
    if not non_index_cols:
        print("No non-index columns found in filter_file; skipping filtering.")
        return df_input

    # Count how many non-index columns exceed min_value for each row
    # Build a row-wise counter of how many non-index columns exceed min_value
    cond_exprs = [
        (pl.col(c) > min_value).cast(pl.Int64) for c in non_index_cols
    ]
    df_filter = df_filter.with_columns(
        sum_above_min=pl.sum_horizontal(cond_exprs),
    )

    df_kept = df_filter.filter(pl.col("sum_above_min") >= min_non_index)
    print(
        "Filtering: kept "
        f"{df_kept.height} of {df_filter.height} rows in filter_file"
    )

    # Restrict input dataframe to index combinations present in df_kept
    keys = df_kept.select(index_columns).unique()
    df_joined = df_input.join(keys, on=index_columns, how="inner")

    return df_joined


def prepare_for_imgt_test(
    input_file: Path,
    annotation_file: Path,
    output_file: Path,
    index_columns: List[str],
    filter_file: Path | None = None,
    min_non_index: int = 1,
    min_value: float = 0.0,
) -> None:
    print(f"Reading input file: {input_file}")
    # IMGT_position can contain non-numeric labels (e.g. "111A"), so
    # force it to be read as a string to avoid parse errors.
    df = pl.read_csv(
        input_file,
        separator="\t",
        schema_overrides={"IMGT_position": pl.Utf8},
    )

    missing_idx = [c for c in index_columns if c not in df.columns]
    if missing_idx:
        raise ValueError(
            f"Index columns {missing_idx} not found in input_file {input_file}"
        )

    n_input_rows = df.height
    print(f"Input rows: {n_input_rows}")

    if filter_file is not None:
        df = apply_filter(
            df_input=df,
            filter_path=filter_file,
            index_columns=index_columns,
            min_non_index=min_non_index,
            min_value=min_value,
        )
        print(f"Rows after filtering: {df.height}")

    # Melt wide table to long format: one row per sample column
    value_cols = [c for c in df.columns if c not in index_columns]
    if not value_cols:
        raise ValueError("No non-index (sample) columns found in input_file")

    print(f"Reshaping to long format with {len(value_cols)} sample columns")

    # Use unpivot (replacement for deprecated melt)
    df_long = df.unpivot(
        index=index_columns,
        on=value_cols,
        variable_name="sample",
        value_name="value",
    )

    # Derive sample_short from the sample column (strip known suffix if present)
    df_long = df_long.with_columns(
        sample_short=pl.col("sample").str.replace("_subs_rows_freq$", ""),
    )

    # Read annotation / collated_info
    print(f"Joining with annotation file: {annotation_file}")
    annot = pl.read_csv(annotation_file, separator=",")

    required_annot_cols = ["sample_short", "shortname", "cells", "genotype_short"]
    missing_annot = [c for c in required_annot_cols if c not in annot.columns]
    if missing_annot:
        raise ValueError(
            "Annotation file is missing required columns: " f"{missing_annot}"
        )

    annot_small = annot.select(required_annot_cols)

    df_joined = df_long.join(annot_small, on="sample_short", how="left")

    # Recode genotype_short into compact labels used for testing
    df_joined = df_joined.with_columns(
        genotype=(
            pl.when(pl.col("genotype_short") == "homoDQ2").then(pl.lit("hoDQ2"))
            .when(pl.col("genotype_short") == "homoDQ8").then(pl.lit("hoDQ8"))
            .when(pl.col("genotype_short") == "heteroDQ2DQ8").then(pl.lit("heDQ2DQ8"))
            .otherwise(pl.col("genotype_short"))
        )
    )

    # Build final output columns in the requested order
    out_cols = []
    out_cols.append(pl.col("sample"))

    # Map aminoAcid_length -> length if present
    if "aminoAcid_length" in df_joined.columns:
        out_cols.append(pl.col("aminoAcid_length").alias("length"))

    if "vFamilyName" in df_joined.columns:
        out_cols.append(pl.col("vFamilyName").alias("vFamilyName"))    

    out_cols.extend(
        [
            pl.col("AA"),
            pl.col("IMGT_position"),
            pl.col("shortname").alias("patient"),
            pl.col("cells"),
            pl.col("value"),
            pl.col("genotype"),
        ]
    )

    df_out = df_joined.select(out_cols)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing long-format output to: {output_file}")
    df_out.write_csv(output_file, separator="\t")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reshape IMGT AA median frequency table to long format and "
            "join with sample annotations for IMGT testing."
        )
    )
    parser.add_argument(
        "--input_file",
        required=True,
        help=(
            "Path to wide IMGT AA table (e.g. comb_aa_imgt_full_subs_rows_freq.tsv) "
            "with index columns followed by sample-specific columns."
        ),
    )
    parser.add_argument(
        "--annotation_file",
        required=True,
        help=(
            "Annotation CSV file like data/collated_info.csv containing "
            "at least columns: genotype_short,sample_short,shortname,cells."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output TSV with long-format data for testing.",
    )
    parser.add_argument(
        "--filter_file",
        default=None,
        help=(
            "Optional filter file (e.g. comb_aa_imgt_full_subs_rows.tsv). "
            "If provided, rows are kept only when at least N non-index "
            "columns exceed min_value."
        ),
    )
    parser.add_argument(
        "--index_columns",
        nargs="+",
        default=DEFAULT_INDEX_COLUMNS,
        help=(
            "Index columns defining unique IMGT rows. "
            "Default: IMGT_position AA aminoAcid_length."
        ),
    )
    parser.add_argument(
        "--min_non_index",
        type=int,
        default=1,
        help=(
            "Minimum number of non-index columns that must exceed min_value "
            "to keep a row when --filter_file is used. Default: 1."
        ),
    )
    parser.add_argument(
        "--min_value",
        type=float,
        default=0.0,
        help=(
            "Minimum value threshold used with --filter_file. "
            "A row is kept if at least min_non_index non-index columns are "
            "> min_value. Default: 0.0."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_path = Path(args.input_file)
    annot_path = Path(args.annotation_file)
    output_path = Path(args.output)
    filter_path = Path(args.filter_file) if args.filter_file else None


    prepare_for_imgt_test(
        input_file=input_path,
        annotation_file=annot_path,
        output_file=output_path,
        index_columns=args.index_columns,
        filter_file=filter_path,
        min_non_index=args.min_non_index,
        min_value=args.min_value,
    )


if __name__ == "__main__":
    main()

