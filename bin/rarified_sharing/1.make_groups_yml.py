#!/usr/bin/env python3
"""
make_groups_yml.py
------------------
Generate a YAML groups file for rare_publicity.py from a delimited table.

Each unique value of the --group_col column becomes one YAML group, and the
corresponding values in --sample_col are the group members.

The table can come from:
  • a file (TSV, CSV, or any single-character delimiter)
  • stdin — so you can pre-filter with awk, grep, or csvkit before piping in

Examples
--------
# Basic usage (reads a file):
python make_groups_yml.py \\
    --input data/collated_info.csv \\
    --group_col genotype_short \\
    --sample_col sample_short \\
    --output groups.yml

# Auto-detect delimiter, write to stdout:
python make_groups_yml.py \\
    --input ../../data/collated_info.csv \\
    --group_col genotype_short \\
    --sample_col sample_short

# Filter rows first, then pipe to this script (stdin):
awk -F',' 'NR==1 || $3=="E"' data/collated_info.csv | \
    python ./sandbox/rare/make_groups_yml.py \
        --group_col genotype_short \
        --sample_col sample_short \
        --output groups_E.yml \
        --rename "heteroDQ2DQ8=hete,homoDQ8=DQ8,homoDQ2=DQ2"

# Using csvkit for more complex filtering:
csvgrep -c cells -m E data/collated_info.csv | \\
    python make_groups_yml.py \\
        --group_col genotype_short \\
        --sample_col sample_short \\
        --output groups_E_only.yml

# Rename groups with --rename (applied after grouping):
python make_groups_yml.py \\
    --input data/collated_info.csv \\
    --group_col genotype_short \\
    --sample_col sample_short \\
    --rename "heteroDQ2DQ8=hete,homoDQ8=homo" \\
    --output groups.yml
"""

import argparse
import csv
import io
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sniff_delimiter(sample: str) -> str:
    """Sniff the delimiter from a small sample of text."""
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        return dialect.delimiter
    except csv.Error:
        return ","  # default


def read_table(
    source,
    delimiter: str | None,
    group_col: str,
    sample_col: str,
) -> dict[str, list[str]]:
    """
    Read *source* (file-like object) and return an ordered dict:
        { group_value: [sample1, sample2, ...], ... }
    Insertion order of groups and samples is preserved.
    """
    raw = source.read()

    # Auto-detect delimiter if not provided
    sep = delimiter if delimiter else sniff_delimiter(raw[:4096])

    reader = csv.DictReader(io.StringIO(raw), delimiter=sep)

    if reader.fieldnames is None:
        sys.exit("ERROR: Could not parse table header.")

    missing = [c for c in (group_col, sample_col) if c not in reader.fieldnames]
    if missing:
        sys.exit(
            f"ERROR: Column(s) not found in table: {missing}\n"
            f"Available columns: {list(reader.fieldnames)}"
        )

    groups: dict[str, list[str]] = {}
    for row in reader:
        grp = row[group_col].strip()
        smp = row[sample_col].strip()
        if not grp or not smp:
            continue
        groups.setdefault(grp, [])
        if smp not in groups[grp]:          # deduplicate within group
            groups[grp].append(smp)

    return groups


def apply_renames(
    groups: dict[str, list[str]],
    rename_str: str,
) -> dict[str, list[str]]:
    """
    Apply user-supplied renames of the form  'old1=new1,old2=new2'.
    Group names not present in the mapping are kept as-is.
    """
    renames: dict[str, str] = {}
    for token in rename_str.split(","):
        token = token.strip()
        if "=" not in token:
            sys.exit(
                f"ERROR: Invalid --rename token '{token}'. "
                "Expected format: old_name=new_name"
            )
        old, new = token.split("=", 1)
        renames[old.strip()] = new.strip()

    renamed: dict[str, list[str]] = {}
    for name, members in groups.items():
        renamed[renames.get(name, name)] = members
    return renamed


def groups_to_yaml(groups: dict[str, list[str]]) -> str:
    """
    Serialise the groups dict to a YAML string without requiring PyYAML
    (plain Python formatting is sufficient for this simple structure).
    """
    lines: list[str] = []
    for name, members in groups.items():
        lines.append(f"{name}:")
        for m in members:
            lines.append(f"  - {m}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a YAML groups file for rare_publicity.py from a delimited table.\n"
            "Input can be a file or piped through stdin."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", default="-",
        help=(
            "Path to a CSV/TSV file, or '-' to read from stdin (default: stdin). "
            "Delimiter is auto-detected unless --delimiter is set."
        ),
    )
    parser.add_argument(
        "--group_col", "-g", required=True,
        help="Column whose unique values become group names.",
    )
    parser.add_argument(
        "--sample_col", "-s", required=True,
        help="Column whose values are sample names within each group.",
    )
    parser.add_argument(
        "--delimiter", "-d", default=None,
        help="Field delimiter (auto-detected if omitted).",
    )
    parser.add_argument(
        "--rename", "-r", default=None,
        help=(
            "Rename groups after grouping. "
            "Comma-separated list of old=new pairs, "
            "e.g. 'heteroDQ2DQ8=hete,homoDQ8=homo'."
        ),
    )
    parser.add_argument(
        "--output", "-o", default="-",
        help="Output YAML file path, or '-' to write to stdout (default: stdout).",
    )

    args = parser.parse_args()

    # ---- Read input ----
    if args.input == "-":
        if sys.stdin.isatty():
            sys.exit(
                "ERROR: No --input file specified and nothing piped to stdin.\n"
                "       Provide --input FILE or pipe data: cat file.csv | python make_groups_yml.py ..."
            )
        source = sys.stdin
    else:
        try:
            source = open(args.input, newline="")
        except FileNotFoundError:
            sys.exit(f"ERROR: File not found: {args.input}")

    groups = read_table(source, args.delimiter, args.group_col, args.sample_col)

    if args.input != "-":
        source.close()

    if not groups:
        sys.exit("ERROR: No groups found. Check --group_col and --sample_col values.")

    # ---- Optional renames ----
    if args.rename:
        groups = apply_renames(groups, args.rename)

    # ---- Serialise ----
    yaml_str = groups_to_yaml(groups)

    # ---- Write output ----
    if args.output == "-":
        sys.stdout.write(yaml_str)
    else:
        with open(args.output, "w") as fh:
            fh.write(yaml_str)
        # Print summary to stderr so it doesn't pollute stdout when piping
        print(
            f"Written: {args.output}  "
            f"({len(groups)} groups, "
            f"{sum(len(v) for v in groups.values())} samples total)",
            file=sys.stderr,
        )
        for name, members in groups.items():
            print(f"  {name!r}: {members}", file=sys.stderr)


if __name__ == "__main__":
    main()
