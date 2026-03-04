#!/usr/bin/env python


import pandas as pd

# File paths
individual_info_file = "data/aim_1_samples_with_newname.csv"
samplesheet_file = "validated_samplesheet.tsv"
output_file = "data/collated_info.csv"

def collate_info():
    # Read the individual_info.csv file
    individual_info = pd.read_csv(individual_info_file, delimiter=",")

    # Read the samplesheet.csv file
    samplesheet = pd.read_csv(samplesheet_file, delimiter="\t")

    # Extract the base newname from the SAMPLE field in samplesheet
    samplesheet['base_newname'] = samplesheet['SAMPLE'].str.extract(r'^(.*)_(?:NAIVE|AgXP)_TCRB\.tsv\.gz$')[0]
     # Handle exceptions
    samplesheet['base_newname'] = samplesheet['base_newname'].replace({
        'FIN_6828': 'FIN_6826',
        'TN0013371F01': 'TN0013372F01'
    })
    samplesheet['cells'] = samplesheet['SAMPLE'].str.extract(r'^(.*)_(NAIVE|AgXP)_(TCRB\.tsv\.gz)$')[1]



    # Merge the two dataframes on base_newname (samplesheet) and newname (individual_info)
    collated = samplesheet.merge(individual_info, left_on="base_newname", right_on="newname", how="left")

    # add short_samplename as shortname + cells (joined with underscore)
    # map cells values: AgXP -> E, NAIVE -> N
    collated['cells'] = collated['cells'].replace({'AgXP': 'E', 'NAIVE': 'N'})
    collated['sample_short'] = collated['shortname'].astype(str) +  collated['cells'].astype(str)

    # Drop the temporary base_newname column
    collated.drop(columns=['base_newname'], inplace=True)

    # Save the collated dataframe to a new CSV file
    collated.to_csv(output_file, index=False)
    print(f"Collated file saved to {output_file}")

if __name__ == "__main__":
    collate_info()