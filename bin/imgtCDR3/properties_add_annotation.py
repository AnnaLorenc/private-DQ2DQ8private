#!/usr/bin/env python3
"""Script: properties_add_annotation.py

Summary
-------
Adds patient, cells, and genotype columns to the output of
imgtCDR3/imgt_wide_to_properties.py, producing a file with the same shape as
sandbox/properties_aggregated.tsv.

Input
-----
1. Long-format TSV from imgtCDR3/imgt_wide_to_properties.py:
       <group cols>  sample  <property cols>
   where 'sample' contains values like 'F15625E_subs_rows'.
# Script: properties_add_annotation.py




2. Annotation CSV (collated_info-style) with at least:
       sample_short, shortname, cells, genotype_short
   where sample_short matches the stripped sample name.

Output
------
TSV with columns:
    <non-AA group cols (e.g. aminoAcid_length → length, IMGT_position)>
    cells  patient  sample  genotype  <property cols>

matching the format of sandbox/properties_aggregated.tsv.


    # Recodes genotype to match properties_aggregated.tsv logic
    python properties_add_annotation.py \\
        --input  sandbox/properties_wide_out.tsv \\
        --annotation data/collated_info.csv \\
        --output sandbox/properties_annotated.tsv \\
        [--sample_col sample] \\
        [--sample_suffix _subs_rows] \\
        [--length_col aminoAcid_length]
"""

import argparse
from pathlib import Path

import pandas as pd


GENOTYPE_MAP = {
    "homoDQ2":      "hoDQ2",
    "homoDQ8":      "hoDQ8",
    "heteroDQ2DQ8": "heDQ2DQ8",
}


def recode_genotype(series: pd.Series) -> pd.Series:
    return series.map(lambda x: GENOTYPE_MAP.get(x, x))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Add patient, cells and genotype annotation to imgtCDR3/imgt_wide_to_properties.py "
            "output, producing properties_aggregated.tsv-style output."
        )
    )
    parser.add_argument("--input", required=True,
                        help="Output TSV from imgtCDR3/imgt_wide_to_properties.py.")
    parser.add_argument("--annotation", required=True,
                        help="Annotation CSV with columns: sample_short, "
                             "shortname, cells, genotype_short.")
    parser.add_argument("--output", required=True,
                        help="Output TSV path.")
    parser.add_argument(
        "--sample_col", default="sample",
        help="Name of the sample column in --input. Default: sample.",
    )
    parser.add_argument(
        "--sample_suffix", default=None,
        help="Suffix to strip from sample names to get sample_short "
             "(e.g. '_subs_rows'). Auto-detected if not provided.",
    )
    parser.add_argument(
        "--length_col", default="aminoAcid_length",
        help="Name of the length column in --input to rename to 'length'. "
             "Default: aminoAcid_length. Set to '' to skip renaming.",
    )
    parser.add_argument(
        "--index_col", default=None,
        help="Optional index column name (e.g. a row ID or group identifier). "
             "Will be preserved but excluded from property columns.",
    )
    parser.add_argument(
        "--exclude_cols",
        default="sample_short,genotype_short,cells,patient,genotype,IMGT_position,AA,aminoAcid_length,length,vFamilyName",
        help="Comma-separated list of column names to exclude from property columns. "
             "Default includes annotation and grouping columns.",
    )
    args = parser.parse_args()

    print(f"Reading: {args.input}")
    df = pd.read_csv(args.input, sep="\t", dtype={"IMGT_position": str})

    if args.sample_col not in df.columns:
        raise ValueError(f"Column '{args.sample_col}' not found in input")

    # --- Strip suffix to get sample_short ---
    if args.sample_suffix:
        suffix = args.sample_suffix
        df["sample_short"] = df[args.sample_col].str.replace(
            suffix + "$", "", regex=True
        )
    else:
        # Auto-detect: try to extract FxxxxN/E pattern
        # e.g. F15625E_subs_rows -> F15625E
        extracted = df[args.sample_col].str.extract(r'^([A-Z]\d+[NE])', expand=False)
        if extracted.notna().all():
            df["sample_short"] = extracted
        else:
            # fallback: strip everything after last _ run that is all lowercase
            df["sample_short"] = df[args.sample_col].str.replace(
                r'[_a-z].*$', '', regex=True
            )
        print(f"  Auto-detected sample_short examples: "
              f"{df['sample_short'].unique()[:5].tolist()}")

    # --- Load annotation ---
    print(f"Reading annotation: {args.annotation}")
    # try tab then comma
    try:
        annot = pd.read_csv(args.annotation, sep="\t")
        if annot.shape[1] < 2:
            raise ValueError
    except Exception:
        annot = pd.read_csv(args.annotation, sep=",")

    required = ["sample_short", "shortname", "cells", "genotype_short"]
    missing = [c for c in required if c not in annot.columns]
    if missing:
        raise ValueError(f"Annotation file missing columns: {missing}")

    annot_small = annot[required].drop_duplicates()

    # --- Join ---
    before = len(df)
    df = df.merge(annot_small, on="sample_short", how="left")
    n_unmatched = df["shortname"].isna().sum()
    if n_unmatched:
        print(f"  Warning: {n_unmatched} rows could not be matched to annotation")
    print(f"  Joined {before} rows; matched {before - n_unmatched}")

    # Recode genotype
    df["genotype"] = recode_genotype(df["genotype_short"])
    df = df.rename(columns={"shortname": "patient"})

    # Rename length column if requested
    if args.length_col and args.length_col in df.columns:
        df = df.rename(columns={args.length_col: "length"})
        length_out = "length"
    else:
        length_out = args.length_col if args.length_col in df.columns else None

    # Identify property columns (all remaining numeric cols that are not
    # in the exclude list)
    exclude_list = set(args.exclude_cols.split(",")) if args.exclude_cols else set()
    exclude_list.add(args.sample_col)
    if args.index_col:
        exclude_list.add(args.index_col)
    
    prop_cols = [c for c in df.columns
                 if c not in exclude_list
                 and pd.api.types.is_numeric_dtype(df[c])]

    # Build output column order:
    # [index_col]  [length]  [all exclude_cols that exist]  <prop_cols>
    ordered = []
    if args.index_col and args.index_col in df.columns:
        ordered.append(args.index_col)
    
    if length_out and length_out in df.columns:
        ordered.append(length_out)
    
    # Add all exclude_cols that exist in the dataframe
    exclude_cols_list = [c.strip() for c in args.exclude_cols.split(",")] if args.exclude_cols else []
    for c in exclude_cols_list:
        if c in df.columns and c not in ordered and c != length_out:
            ordered.append(c)
    
    # Add property columns
    for c in prop_cols:
        if c not in ordered:
            ordered.append(c)

    df_out = df[ordered]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.output, sep="\t", index=False)
    print(f"Written {len(df_out)} rows → {args.output}")


if __name__ == "__main__":
    main()
