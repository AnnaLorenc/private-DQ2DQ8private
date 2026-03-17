#!/usr/bin/env python3
"""Extract productive/non-productive sequence counts from cleanup summary files."""

import csv
import re
from pathlib import Path

cleanup_dir = Path(__file__).parent.parent / "results" / "cleanup"
out_file = cleanup_dir / "cleanup_summary_all.csv"

rows = []
for tsv in sorted(cleanup_dir.rglob("*_cleanup_summary.tsv")):
    sample = tsv.parent.name
    data = {"sample": sample, "num_seq_prod": "", "num_seq_nonprod": "", "counts_prod": "", "counts_nonprod": ""}
    with open(tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            cat = parts[0].strip()
            if re.search(r"Productive Sequences", cat, re.IGNORECASE) and not re.search(r"Non", cat, re.IGNORECASE):
                data["num_seq_prod"] = parts[1].strip()
                data["counts_prod"] = parts[2].strip()
            elif re.search(r"Non.Productive Sequences", cat, re.IGNORECASE):
                data["num_seq_nonprod"] = parts[1].strip()
                data["counts_nonprod"] = parts[2].strip()
    rows.append(data)

with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["sample", "num_seq_prod", "num_seq_nonprod", "counts_prod", "counts_nonprod"])
    writer.writeheader()
    writer.writerows(rows)

print(f"Written {len(rows)} samples → {out_file}")
