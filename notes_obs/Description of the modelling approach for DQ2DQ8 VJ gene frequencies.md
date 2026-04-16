I used 'standard' modelling for classification of repertoires based on /J gene frequency. The rationale was that it is easily visible on the PCA, so linear discriminant analysis should do the job.
Modelling was done in R (tidymodels etc)

1. Considered input data
	1. Data
		1. V genes detected in at least 41 samples (removes "TCRBV26" "TCRBV08" "TCRBV17" "TCRBV22" "TCRBVA" , last one more prevalent, detected in 14% of samples )
		2. Averaged frequencies from subsampling, all cells, all stages (122 samples)
		3. Subsamples 25 per individual- for cross-vaildation split is by individual, so the model will see a new individual in test/validation data
		4. Info about celltype, stage and genotype
	2. Data preprocessing
		1. for all predictors:  removing predictors without variance and with tiny variance; normalise numerical variables (V/J gene frequencies) (centered and scaled)
	3. Predictors
		1. V genes (family)
		2. V (family) and J genes
		3. V (family) a, pca-preprocessed
		4. V (family) and J genes, pca-preprocessed
		5. celltype, stage and V (family) and J genes
		6. celltype, stage and V (family) and J genes, genes pca-preprocessed
			
2. Considered modelling approaches
	1. LDA - originally I wanted to try more, but as results were really good and this model is easier in interpretation than others, I stayed with it
3. Modelling
	1. Training data: 80% of the data, selected by individual
	2. 10x crossvalidation (data is randomly split in 10 folds; training is done separately 10 times, with 9 folds each time and the final fold used as validation)
	 3. Final crossvalidation check with the test data 


## Results


### Data: 2.
	![[geno_from_vj_modelchoice_122B.png]]
	![[geno_from_vj_modelchoice_122 1.png]]
All versions have overlapping ranges; therefore I decided to test with subsamples. Also, looked on inference on the test data with the model with predictor set 2 (V and J gene)
##### Final crossvalidation 
![[geno_from_vj_model_122 1.png]]

### Data: 3.
![[geno_from_vj_modelchoice_3050A 1.png]]
V and J gene data (just normalised), enriched by stage and celltype info are the best model - this is seen in almost all folds of the 10-fold crossvalidation.
![[geno_from_vj_modelchoice_3050B.png]]
##### Final crossvalidation
With held out 24 samples
![[geno_from_vj_model.png]]

However, there is an issue. The LDs are very corrrelated: (Pearson correlaton R=0.99; SVD 86 and 34)
![[LD_from_vj_model.png]]

When we look at top 5 predictors, DQ2DQ8 has values in between of homozygotes (TCRBV28, TCRBV21, TCRBV15) or outside of the range of homozygotes (TCRBV05, TCRBV10)
![[mean_contrib_of_top5preds.png]]
But the catual values are not so convincing...
![[actual_values_of_top5preds.png]]

But I should have used relative top:
![[actual_values_of_top5RELATIVEpreds.png]]