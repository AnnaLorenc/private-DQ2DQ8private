
Write a script which summarizes sharing: it will read a tsv file. Input file has a structure: index columns (identifyling a clonotype) and then columns with number of individuals sharing a clonotype in a random subsample.

Input:
file like hete_sharing.tsv - one or more
parameter: A

Obtain a summary matrix, clonotypes x shared by at least N individuals, where values are in how many draws.N should be in the range from 2 to the max value in the input matrix. For example,if there were 4 draws and a clonotype was shared by 3,5,7,0 individuals in each, the matrix would be:
aminoAcid       vFamilyName     jGeneName 	2 3 4 5 6 7 8 9 10
CAAGGAGGTGELFF  TCRBV02 TCRBJ02-02	3	3	2	1	0	0	0	0	0 

Make it a function summarise_sharing_group.


Make it possible for the script to read several input files, compute summarise_sharing_group for each of them.
For each output, identify clonotypes shared by A or more individuals at least 1.
subset each of the  summarise_sharing_group outputs to the union of those clonotypes.
Merge these outputs, adding sharing group name, for example for A=3:

aminoAcid       vFamilyName     jGeneName 	 3 4 5 6 7 8 9 10	group
CAAGGAGGTGELFF  TCRBV02 TCRBJ02-02	3	2	1	0	0	0	0	0 het
CAAGGAGGTGELFF  TCRBV02 TCRBJ02-02	3	2	2	1	0	0	0	0 hom

Output should be a set of files, one for each output of summarise_sharing_group and a common file with merged outputs.

Include a header, explanations for functions, checks of the input and parameters.
 Produce also a small play file, with 20-40 clonotypes and 10 samples to test the script.