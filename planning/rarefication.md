I want to write a fast, polars and matrix algebra based, based python script for judging publicity of clonotypes.

 It will start from a matrix ALL_clonotypes x  ALL_samples, the values are counts in these samples. This matrix has to be build from an input file; clonotypes are identified in this file by a combination of index_columns. 
This input should be subset to the samples list clonotypes_non0_in_samples_S x S 
Then, for each sample M clonotypes should be sampled with replacement in a similar way as in the script bin/subsample_merged_sequences.py - using the original counts of a given clonotype in the starting matrix as weight for sampling
	Per one subsampling:Then, for each of the sharing group G, script counts in how many sub-samples given clonotype is shared. So the output of the first round should be length(G) 1-column matrices, each one clonotypes_non0_in_samples_S x how many individuals shared.
	This step is repeated K times, with sharing values added in each round.
	The final output are length(G) matrices, each one of size clonotypes_non0_in_samples_S x K.
	
input:file with ALL_clonotypes x  samples, with index columns specifying clonotypes
/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/merged/merged_aminoAcid_vFamilyName_jGeneName_productive.tsv.gz
parameters: 
	S samples to use, (str names)
	K how many subsamples to perform (25, int)
	M how deep to subsample
	G sharing groups (list of sample names, subsets of S, with names, for example --groups homo:F28396E,F28396N,F20790N hete:T04918N,T06306E,T13372E would be 2 groups, group names homo and hete, samples in homo:F28396E,F28396N,F20790, in  hete:T04918N,T06306E,T13372E )
	index_columns index columns to identify clonotypes, like aminoAcid       vFamilyName     jGeneName 
	  
output: directory with matrices. Matrices names from group names.
Include a header, explanations for functions, checks of the input and parameters.
 Produce also a small play file, with 20-40 clonotypes and 10 samples to test the script.