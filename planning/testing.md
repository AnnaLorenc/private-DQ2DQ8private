You are a seasoned data analyst. Write the following scripts, test them by producing output in the directory results_man/imgt_test. 
The only files you are allowed to write to outside of results_man/imgt_test are:


1. Write a python script bin/imgtCDR3/prepare_for_imgt_test.py to restructure the input_file of a shape like ./results/imgt_aa_combined/comb_aa_imgt_full_subs_rows_freq.tsv:
IMGT_position is a string.
Add as an option filtering before this processing, by a filter_file like  results/imgt_aa_combined/comb_aa_imgt_full_subs_rows.tsv.
Take from filter_file index columns for the rows where at least N non-index columns  (Index columns: by default IMGT_position   AA      aminoAcid_length , should be possible to indicate) have value > min. 

input_file
IMGT_position   AA      aminoAcid_length        F16018E_subs_rows_freq  F20790N_subs_rows_freq  F18072N_subs_rows_freq ...
104     C       8       1.0     0.8     0.3    ...
Into a long format like:
sample length AA IMGT_position patient cells      value
F16018E_subs_rows_freq 8 C 104 F16018 E 1
F20790N_subs_rows_freq 8 C 104 F20790 N 0.8
F18072N_subs_rows_freq 8 C 104 F18072 N 0.3
...


The scipt then joins this this table two last columns of an annotation file  like data/collated_info.csv (genotype_short,sample_short) via sample==sample_short.
Use current version of polars or other fast implementations when possible.
Implement similarily to other scripts in bin/*.py, so it is runnable with CLI, has a header, gets main() call.

arguments: input_file, filter_file (optional), index_columns ( default IMGT_position   AA      aminoAcid_length), annotation_file, output
Print out the number of input_file rows, if filtering - how many rows left after filtering, info that joining with annotation happens

2. Write a python script bin/imgtCDR3/test_imgt.py. 
Use current version of polars or other fast implementations when possible.
Implement similarily to other scripts in bin/*.py, so it is runnable with CLI, has a header, gets main() call.

Input: file of the shape:
sample length AA IMGT_position patient cells   value genotype
F16018E_subs_rows_freq 8 C 104 F16018 E 1  homoDQ2
F20790N_subs_rows_freq 8 C 104 F20790 N 0.8 homoDQ8
F18072N_subs_rows_freq 8 C 104 F18072 N 0.3  heteroDQ2DQ8

The script should test with a linear model, fitted in several ways.
It should test separately for each combination of cells, length, IMGT_position, AA.
It should test groups: 
    grouping1, "homhom" :homoDQ2 versus homoDQ8;
    grouping2, "hethom": heteroDQ2DQ8 versus (any of homoDQ2, homoDQ8)


If this would be in R, I would use with geno the column indicating group, ref_geno indicating :
input%>%
   group_by(cells, length, IMGT_position, AA)%>%
  group_modify(~ broom::tidy(lm(value~geno, data=.x)))

and to get Cohen's D
input%>%
   group_by(cells, length, IMGT_position, AA)%>%
  group_modify(~ cohens_d(formula = value~geno, data=.x, ref.group=ref_geno))


Add another testing to perform. This should be a paired test between naive (N) and experienced (E) cells in each individual. seperate test for each  set of index columns (by default IMGT position, AA, length) Perform it for all 61 individuals together and then also by splitting into 3 genotype groups. Output results of this test in a similar format as the previous one, in a separate file with _cells in the name.