#!/usr/bin/env python

#Helper script for preparing a useful sample sheet with better sample names for downstream analyses. It reads in the original sample sheet, applies the specified rules to generate new sample names, and saves the updated information to a new CSV file.

import csv

def process_samples(input_file, output_file):
    """Process the input CSV file and add  new columns with better sample names based on the specified rules."""
    with open(input_file, mode='r') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        # Read the header and add the new column names
        # header = next(reader)
        # header.extend(["newname", "shortname"])
        writer.writerow(["newname", "shortname", "genotype_short"])

        # Process each row
        header = next(reader, None)
        for row in reader:
            newrow = []
            first_col = row[0]
            second_col = row[2] if len(row) > 1 else ""
            last_col = row[-1] if len(row) > 0 else ""
            short_genotype_col = row[-2].replace(" ", "") if len(row) > 0 else ""

            # Determine the newname value
            if first_col == "Finland":
                if "IN ERROR SEND AS 6828" in last_col:
                    newname = "FIN_6826"
                else:
                    newname = f"FIN_{second_col}"
            elif first_col == "TrialNet":
                newname = second_col
            else:
                newname = ""

            # Determine the shortname value
            if newname.startswith("FIN_"):
                shortname = newname.replace("IN_", "")
            elif newname.startswith("TN"):
                shortname = f"T{newname[4:-3]}"
            else:
                shortname = ""

            # Append the new column values to the row
            newrow.extend([newname, shortname, short_genotype_col])
            writer.writerow(newrow)

def main():
    input_file = "/Users/ania/Documents/DQ2DQ8/data_link/D20210208A_1/metadata/aim_1_samples.csv"
    output_file = "./data/aim_1_samples_with_newname.csv"
    process_samples(input_file, output_file)
    print(f"Processed file saved to {output_file}")

if __name__ == "__main__":
    main()