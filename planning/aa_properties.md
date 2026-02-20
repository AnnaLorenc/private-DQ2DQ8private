You are a seasoned data analyst, specialising in understanding how aminoacids in CDR3 TRB are influenced by their physychomechemical properties.
You will write scripts continuing the previous analysis and some new approaches.
Write scripts in python, making them executable from CLI. add header explaining usage, argparse, run with main(). Explain functions.
You can read all the files in /Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8, but write only in /Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/sandbox - you can create there all the directories, files you need.

1. First approach:
Summary: Aminoacids proerties instead of specific aminoacids
Translate aminoacids into their physichochemical properties, according to several schemes, described below. Then tests for specific positions being enriched for some aminoacid properties. Below I am describing this on the example of Kidera factors - makee it general for any set of properties given as a matrix in a tsv/csv file.

script A (invent a good name)
The input file is a file like results/imgt_aa_fortest/comb_aa_imgt_full_subs_rows_freq_5_1_test_ready.tsv. 
Its index columns are length,	AA,	IMGT_position, cells - but could be different, make these default but nor hardcoded.
Columns to aggregate within are index columns (minus AA, as we are summarising AAs).
Convert each aminoacid to 10 Kidera factors. For each combination of index columns (minus AA, as we are summarising AAs), aggregate the Kidera factors per individual as mean of their Kidera factors, weighted by the aminoacid frequencies (field "value"). Take care to keep the sample and genotype fields. Output an intermediate file. 

script B (invent a good name)
Perform statistical testing on the output of script A. First, use Multivariate Hotteling's T. Then compare for each Kidera factor independently, for a posthoc analysis.
Comparisons should be done within a cell type (E or N)
Groups to compare: hoDQ8 versus hoDQ2 (oo), heDQ2DQ8 versus (hoDQ8 and hoDQ2), (eo).
Output files: file 1 with Hotteling's T results (one row per a combination of columns to aggregate within), side by side oo and eo tests.
file 2: A row per Kidera  factor per per a combination of columns to aggregate within, side by side oo and eo tests, also add pvalue from Hotelling test. 

In all statistical test, provide also effect size whenever possible. Report raw pvalues and FDR corrected (q values). Describe how you onbtained these in the file aa_properties_reasoning.md in the sandbox. 

Kidera factors:
   KF        A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V
 1 KF1   -1.56  0.22  1.14  0.58  0.12 -0.47 -1.45  1.46 -0.41 -0.73 -1.04 -0.34 -1.4  -0.21  2.06  0.81  0.26  0.3   1.38 -0.74
 2 KF2   -1.67  1.27 -0.07 -0.22 -0.89  0.24  0.19 -1.96  0.52 -0.16  0     0.82  0.18  0.98 -0.33 -1.08 -0.7   2.1   1.48 -0.71
 3 KF3   -0.97  1.37 -0.12 -1.58  0.45  0.07 -1.61 -0.23 -0.28  1.79 -0.24 -0.23 -0.42 -0.36 -1.15  0.16  1.21 -0.72  0.8   2.04
 4 KF4   -0.27  1.87  0.81  0.81 -1.05  1.1   1.17 -0.16  0.28 -0.77 -1.1   1.7  -0.73 -1.43 -0.75  0.42  0.63 -1.57 -0.56 -0.4 
 5 KF5   -0.93 -1.7   0.18 -0.92 -0.71  1.1  -1.31  0.1   1.61 -0.54 -0.55  1.54  2     0.22  0.88 -0.21 -0.1  -1.16  0     0.5 
 6 KF6   -0.78  0.46  0.37  0.15  2.41  0.59  0.4  -0.11  1.01  0.03 -2.05 -1.62  1.52 -0.81 -0.45 -0.43  0.21  0.57 -0.68 -0.81
 7 KF7   -0.2   0.92 -0.09 -1.52  1.52  0.84  0.04  1.32 -1.85 -0.83  0.96  1.15  0.26  0.67  0.3  -1.89  0.24 -0.48 -0.31 -1.07
 8 KF8   -0.08 -0.39  1.23  0.47 -0.69 -0.71  0.38  2.36  0.47  0.51 -0.76 -0.08  0.11  1.1  -2.3  -1.15 -1.15 -0.4   1.03  0.06
 9 KF9    0.21  0.23  1.1   0.76  1.13 -0.03 -0.35 -1.66  1.13  0.66  0.45 -0.48 -1.27  1.71  0.74 -0.97 -0.56 -2.3  -0.05 -0.46
10 KF10  -0.48  0.93 -1.73  0.7   1.1  -2.33 -0.12  0.46  1.63 -1.78  0.93  0.6   0.27 -0.44 -0.28 -0.23  0.19 -0.6   0.53  0.65

Describe your decisions and planning in the file aa_properties_reasoning.md in the sandbox.