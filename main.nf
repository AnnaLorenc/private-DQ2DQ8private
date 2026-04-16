#!/usr/bin/env nextflow


nextflow.enable.dsl=2

// Parameters
params.input = './data/collated_info.csv'
params.outdir = './results'
params.seqs_to_remove = './data/D20210208D_1-overrepresented-sequences.txt'
params.N = 10000 
params.M = 25
params.index_columns = ["aminoAcid", "vFamilyName", "jGeneName"] //to join samples in merging
params.min_len = 8
params.max_len = 25	// CDR3 length filters for IMGT freqs
params.annotation_file = './data/collated_info.csv'
params.min_non_index = 5
params.min_value = 1
params.prop_VHS = 'assets/VHSE.tsv'
params.prop_kidera = 'assets/kidera_factors.tsv'


include {INITIAL_CLEANUP_SPLIT } from './modules/read_sheet.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_PROD } from './modules/extract_diversity_metrics.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_NONPROD } from './modules/extract_diversity_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_PROD} from './modules/compute_freqs_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_NONPROD} from './modules/compute_freqs_metrics.nf'
include{ COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_PROD } from './modules/combine_diversity_metrics.nf'
include{ COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_NONPROD } from './modules/combine_diversity_metrics.nf'
include{ COMBINE_FREQS as COMBINE_FREQS_PROD } from './modules/combine_freqs.nf'
include{ COMBINE_FREQS as COMBINE_FREQS_NONPROD } from './modules/combine_freqs.nf'
include { MERGE_SAMPLES as MERGE_AMINO_VFAM } from './modules/merge_samples.nf'
include { SUBSAMPLE_FROM_MERGED } from './modules/subsample_from_merged.nf'
include {COMPUTE_FREQS_SUBSAMPLED} from './modules/compute_freqs_subsampled.nf'
include { COMBINE_FREQS_SUBSAMPLED } from './modules/combine_freqs_subsampled.nf'
include {COMPUTE_OVERLAPS} from './modules/compute_overlaps.nf'
include {COMPUTE_OVERLAPS_SUBSAMPLED} from './modules/compute_overlaps_subsampled.nf'
include {COMPUTE_IMGT_AA_FREQS_SUBS} from './modules/imgt_aa_freqs_subsampled.nf'
include { COMBINE_IMGT_AA_FREQS_MED } from './modules/combine_IMGT_AA_freqs.nf'
include { COMBINE_IMGT_AA_FREQS_MED as COMBINE_IMGT_AA_FREQS_MED_WL } from './modules/combine_IMGT_AA_freqs.nf'
include { COMBINE_IMGT_AA_FREQS_MED as COMBINE_IMGT_AA_FREQS_VFAM_MED_WL } from './modules/combine_IMGT_AA_freqs.nf'
include { COMBINE_IMGT_AA_FREQS_MED as COMBINE_IMGT_AA_FREQS_VFAM_MED } from './modules/combine_IMGT_AA_freqs.nf'
include {PREP_COMB_IMGT_TEST } from './modules/prep_comb_IMGT_test.nf'
include {PREP_COMB_IMGT_TEST as PREP_COMB_IMGT_TEST_WL } from './modules/prep_comb_IMGT_test.nf'
include {PREP_COMB_IMGT_TEST as PREP_COMB_IMGT_TEST_VFAM } from './modules/prep_comb_IMGT_test.nf'
include {PREP_COMB_IMGT_TEST as PREP_COMB_IMGT_TEST_WL_VFAM	 } from './modules/prep_comb_IMGT_test.nf'
include {RUN_IMGT_TEST } from './modules/run_imgt_test.nf'
include {RUN_IMGT_TEST as RUN_IMGT_TEST_WL } from './modules/run_imgt_test.nf'
include {RUN_IMGT_TEST as RUN_IMGT_TEST_VFAM } from './modules/run_imgt_test.nf'
include {RUN_IMGT_TEST as RUN_IMGT_TEST_WL_VFAM } from './modules/run_imgt_test.nf'

include {AA_PROPERTIES_TEST} from './modules/aa_properties.nf'
include {AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL } from './modules/aa_properties.nf'
include {AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_VFAM } from './modules/aa_properties.nf'
include {AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL_VFAM } from './modules/aa_properties.nf'
include {AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_kidera } from './modules/aa_properties.nf'
include {AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL_kidera } from './modules/aa_properties.nf'
include {AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_VFAM_kidera } from './modules/aa_properties.nf'
include {AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL_VFAM_kidera } from './modules/aa_properties.nf'

include {RAREFIED_SHARING} from './modules/rarefied_sharing.nf'
include {RAREFIED_SHARING as RAREFIED_SHARING_VFAM } from './modules/rarefied_sharing.nf'
include {RAREFIED_SHARING as RAREFIED_SHARING_AA } from './modules/rarefied_sharing.nf'

include {RAREFIED_DIVERSITY_METRICS}	from './modules/rarefied_diversity_metrics.nf'
include {COMBINE_RAREFIED_DIVERSITY_METRICS} from './modules/combine_rarefied_diversity_metrics.nf'

include { MERGE_SAMPLES as MERGE_AMINO_VFAM_NONPROD } from './modules/merge_samples.nf'
include {SUBSAMPLE_FROM_MERGED as SUBSAMPLE_FROM_MERGED_NONPROD } from './modules/subsample_from_merged.nf'
include {COMPUTE_FREQS_SUBSAMPLED as COMPUTE_FREQS_SUBSAMPLED_NONPROD} from './modules/compute_freqs_subsampled.nf'
include {COMBINE_FREQS_SUBSAMPLED as COMBINE_FREQS_SUBSAMPLED_NONPROD } from './modules/combine_freqs_subsampled.nf'
include {AD_LENGTH_DISTRIBUTION } from './modules/ad_length_distribution.nf'

// Print a header for your pipeline 
log.info """\

=======================================================================================
DQ2DQ8 
=======================================================================================

Created by Ania


=======================================================================================
Workflow run parameters 
=======================================================================================
input       : ${params.input}
results     : ${params.outdir}
workDir     : ${workflow.workDir}
=======================================================================================

"""


workflow {

// Show help message if --help is run or (||) a required parameter (input) is not provided
		
		log.info "seqs_to_remove: ${params.seqs_to_remove}"
		log.info "DEBUG: params.input = '${params.input}'"  // Add this debug line
    	log.info "DEBUG: params.input type = ${params.input?.getClass()}"  // Check type
		// Channel 1: sample_short + genotype_short
		channel.fromPath(params.input)
			.splitCsv(header: true, sep: ',')
			.map { row -> tuple(row.sample_short, row.genotype_short) }
			.set { sample_genotype_ch }

		// Channel 2: sample_short + filepath
		channel.fromPath(params.input)
			.splitCsv(header: true, sep: ',')
			.map { row -> tuple(row.sample_short, "${row.LOC}/${row.SAMPLE}", params.seqs_to_remove, params.outdir) }
			.set { sample_filepath_ch}

		// Use channels
		
		INITIAL_CLEANUP_SPLIT (sample_filepath_ch)

		// INITIAL_CLEANUP_SPLIT.out.productive.view { sample, file -> "Sample ${sample}: ${file}"    }
		EXTRACT_DIVERSITY_METRICS_PROD(INITIAL_CLEANUP_SPLIT.out.productive)
		EXTRACT_DIVERSITY_METRICS_NONPROD(INITIAL_CLEANUP_SPLIT.out.nonproductive)
		COMPUTE_FREQS_PROD(INITIAL_CLEANUP_SPLIT.out.productive)
		COMPUTE_FREQS_NONPROD(INITIAL_CLEANUP_SPLIT.out.nonproductive)

		// Collect productive CSV files
		productive_data = EXTRACT_DIVERSITY_METRICS_PROD.out.diversity_metrics
			.unique()  // Remove duplicates
			.toList()
			.map { items -> 
				def samples = items.collect { it[0] }
				def files = items.collect { it[1] }
				tuple('productive', samples, files)
			}
		
		// Collect nonproductive with deduplication
		nonproductive_data = EXTRACT_DIVERSITY_METRICS_NONPROD.out.diversity_metrics
			.unique()  // Remove duplicates
			.toList()
			.map { items -> 
				def samples = items.collect { it[0] }
				def files = items.collect { it[1] }
				tuple('nonproductive', samples, files)
			}
		
		// Combine each type
		COMBINE_DIVERSITY_METRICS_PROD(productive_data)
		COMBINE_DIVERSITY_METRICS_NONPROD(nonproductive_data)


		// Collect productive CSV files for freqs
		productive_freqs = COMPUTE_FREQS_PROD.out.freqs
		.unique()  // Remove duplicates
		.toList()
		.map { items -> 
			def samples = items.collect { it[0] }
			def files = items.collect { it[1] }
			tuple('productive', samples, files)
		}
	
	// Collect nonproductive 
		nonproductive_freqs = COMPUTE_FREQS_NONPROD.out.freqs
			.unique()  // Remove duplicates
			.toList()
			.map { items -> 
				def samples = items.collect { it[0] }
				def files = items.collect { it[1] }
				tuple('nonproductive', samples, files)
			}

	COMBINE_FREQS_PROD(productive_freqs)
	COMBINE_FREQS_NONPROD(nonproductive_freqs)

	INITIAL_CLEANUP_SPLIT.out.productive
        .map { sample_name, file_path -> file_path }  // Extract just the files
        .collect()  // Collect all files into a single list
        .map { files -> 
            tuple(files, params.index_columns, "count (templates/reads)","productive") 
        }
        .set { prod_merge_ch }



	// prepare one database of productive sequences
    MERGE_AMINO_VFAM(prod_merge_ch)

	//subsample Vfam, Jgene, CDR3
	MERGE_AMINO_VFAM.out.merged
		.map { merged_file -> tuple(merged_file, params.N, params.M, 3, "productive") } // Example values for N, M, num_index_columns .set { subsample_input_ch }
		.set { subsample_input_ch }
	//sample M times N sequences from each sample
	SUBSAMPLE_FROM_MERGED( subsample_input_ch)



	SUBSAMPLE_FROM_MERGED.out.subsampled_files
		.flatten()
    	.map { subsampled_file -> 
        def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
        tuple(subsampled_file, sample_name, "productive", "aminoAcid", [ "vFamilyName","jGeneName"])
    } 
    .set { compute_freqs_subsampled_input_ch }
	
	COMPUTE_FREQS_SUBSAMPLED(compute_freqs_subsampled_input_ch)

	COMPUTE_FREQS_SUBSAMPLED.out.subsampled
		.map { big_file, summ_file -> summ_file}
		.collect()
		.map{ files ->
		tuple(files,  "productive_subs_freqs.tsv") }
		.set { combine_freqs_subsampled_input_ch  }	

	COMBINE_FREQS_SUBSAMPLED(combine_freqs_subsampled_input_ch)

	def first_column_with_sample_in_merged = params.index_columns.size()

	MERGE_AMINO_VFAM.out.merged
		.map { merged_file -> tuple(merged_file,  "productive", first_column_with_sample_in_merged) }
		.set { compute_overlaps_input_ch }
		
	COMPUTE_OVERLAPS(compute_overlaps_input_ch) 

	// Coompute overlaps on subsampled
	SUBSAMPLE_FROM_MERGED.out.sample_names
		.collect()
		.combine(SUBSAMPLE_FROM_MERGED.out.subsampled_dir.collect())
		.map { sample_names, subsampled_dirs -> 
			tuple(sample_names, subsampled_dirs,  "productive")  // Assuming single files
		}
	.set { compute_overlaps_subsampled_input_ch }
	COMPUTE_OVERLAPS_SUBSAMPLED(compute_overlaps_subsampled_input_ch  )

	//----aminoacid usage in CDR3
	SUBSAMPLE_FROM_MERGED.out.subsampled_files
		.flatten() 
		.map { subsampled_file -> 
		// println "Processing file: ${subsampled_file}"
        // println "File name: ${subsampled_file.name}"
        // println "Base name: ${subsampled_file.baseName}"
		def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
		// println "Extracted sample name: ${sample_name}"
		tuple(subsampled_file, sample_name, "productive", params.min_len, params.max_len)
	} 
	.set { compute_imgt_freqs_subsampled_input_ch }	

	COMPUTE_IMGT_AA_FREQS_SUBS(compute_imgt_freqs_subsampled_input_ch)

	//with length, no Vfam
	// Collect all aa_imgt_freq_full_med files
	COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_full_med
		.map { sample_name, file -> file }  // Extract just the files
		.collect()  // Collect all files into a single list
		.set { freq_med_files_ch }

	def index_cols = params.imgt_index_columns ?: ["IMGT_position", "AA", "aminoAcid_length"]
	COMBINE_IMGT_AA_FREQS_MED(freq_med_files_ch, index_cols,"_full")
	
	//without length, no Vfam
	// Collect all aa_imgt_freq_full_med files
	COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_WL_med
	.map { sample_name, file -> file }  // Extract just the files
	.collect()  // Collect all files into a single list
	.set { freq_med_WL_files_ch }

	def index_cols_WL = params.imgt_index_columns_WL ?: ["IMGT_position", "AA"]
	COMBINE_IMGT_AA_FREQS_MED_WL(freq_med_WL_files_ch, index_cols_WL, "_WL")


	//with length, Vfam
	// Collect all aa_imgt_freq_full_med files
	COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_full_Vfam_med	
		.map { sample_name, file -> file }  // Extract just the files
		.collect()  // Collect all files into a single list
		.set { freq_med_vfam_files_ch }

	def index_cols_vfam = params.imgt_index_columns ?: ["IMGT_position", "AA", "aminoAcid_length","vFamilyName"]
	COMBINE_IMGT_AA_FREQS_VFAM_MED(freq_med_vfam_files_ch, index_cols_vfam,"_full_vfam")
	
	//without length, Vfam
	// Collect all aa_imgt_freq_full_med files
	COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_WL_Vfam_med
	.map { sample_name, file -> file }  // Extract just the files
	.collect()  // Collect all files into a single list
	.set { freq_med_WL_vfam_files_ch }

	def index_cols_WL_vfam = params.imgt_index_columns_WL ?: ["IMGT_position", "AA","vFamilyName"]
	COMBINE_IMGT_AA_FREQS_VFAM_MED_WL(freq_med_WL_vfam_files_ch, index_cols_WL_vfam, "_WL_vfam")




    def combined = COMBINE_IMGT_AA_FREQS_MED.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_MED.out.combined_subs_counts)

	//prep for testing - with length, no Vfam &testing
	PREP_COMB_IMGT_TEST(combined, index_cols, params.annotation_file,  params.min_non_index, params.min_value)
	
	def index_cols1 = params.imgt_index_columns ?: ["IMGT_position", "AA", "length"]
	RUN_IMGT_TEST(PREP_COMB_IMGT_TEST.out.imgt_test_input, index_cols1)



	def combined_WL = COMBINE_IMGT_AA_FREQS_MED_WL.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_MED_WL.out.combined_subs_counts)

	//prep for testing - with length, no Vfam &testing
	PREP_COMB_IMGT_TEST_WL(combined_WL, index_cols_WL, params.annotation_file,  params.min_non_index, params.min_value)
	
	RUN_IMGT_TEST_WL(PREP_COMB_IMGT_TEST_WL.out.imgt_test_input, index_cols_WL)



    def combined_vfam = COMBINE_IMGT_AA_FREQS_VFAM_MED.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_VFAM_MED.out.combined_subs_counts)

	//prep for testing - with length,Vfam &testing
	PREP_COMB_IMGT_TEST_VFAM(combined_vfam, index_cols_vfam, params.annotation_file,  params.min_non_index, params.min_value)
	
	def index_cols1_vfam = params.imgt_index_columns ?: ["IMGT_position", "AA", "length", "vFamilyName"]
	RUN_IMGT_TEST_VFAM(PREP_COMB_IMGT_TEST_VFAM.out.imgt_test_input, index_cols1_vfam)



	def combined_WL_vfam	 = COMBINE_IMGT_AA_FREQS_VFAM_MED_WL.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_VFAM_MED_WL.out.combined_subs_counts)

	//prep for testing - with length, no Vfam &testing
	PREP_COMB_IMGT_TEST_WL_VFAM(combined_WL_vfam, index_cols_WL_vfam, params.annotation_file,  params.min_non_index, params.min_value)
	
	RUN_IMGT_TEST_WL_VFAM(PREP_COMB_IMGT_TEST_WL_VFAM.out.imgt_test_input, index_cols_WL_vfam)




	//testing for properties
	

	AA_PROPERTIES_TEST(COMBINE_IMGT_AA_FREQS_MED.out.combined_subs_rows,
	 index_cols,
	  params.prop_VHS,
	  params.annotation_file,
	 "VHS",
	 ["IMGT_position", "cells","length"])

	AA_PROPERTIES_TEST_WL(COMBINE_IMGT_AA_FREQS_MED_WL.out.combined_subs_rows,
	 index_cols_WL,
	  params.prop_VHS,
	  params.annotation_file,
	 "VHS",
	 ["IMGT_position", "cells"])

	AA_PROPERTIES_TEST_VFAM(COMBINE_IMGT_AA_FREQS_VFAM_MED.out.combined_subs_rows,
	 index_cols_vfam,
	  params.prop_VHS,
	  params.annotation_file,
	 "VHS",
	 ["IMGT_position", "cells","length","vFamilyName"])

	AA_PROPERTIES_TEST_WL_VFAM(COMBINE_IMGT_AA_FREQS_VFAM_MED_WL.out.combined_subs_rows,
	 index_cols_WL_vfam,
	  params.prop_VHS,
	  params.annotation_file,
	 "VHS",
	 ["IMGT_position", "cells","vFamilyName"])

	AA_PROPERTIES_TEST_kidera(COMBINE_IMGT_AA_FREQS_MED.out.combined_subs_rows,
	 index_cols,
	  params.prop_kidera,
	  params.annotation_file,
	 "KF",
	 ["IMGT_position", "cells","length"])

	AA_PROPERTIES_TEST_WL_kidera(COMBINE_IMGT_AA_FREQS_MED_WL.out.combined_subs_rows,
	 index_cols_WL,
	  params.prop_kidera,
	  params.annotation_file,
	 "KF",
	 ["IMGT_position", "cells"])

	AA_PROPERTIES_TEST_VFAM_kidera(COMBINE_IMGT_AA_FREQS_VFAM_MED.out.combined_subs_rows,
	 index_cols_vfam,
	  params.prop_kidera,
	  params.annotation_file,
	 "KF",
	 ["IMGT_position", "cells","length","vFamilyName"])

	AA_PROPERTIES_TEST_WL_VFAM_kidera(COMBINE_IMGT_AA_FREQS_VFAM_MED_WL.out.combined_subs_rows,
	 index_cols_WL_vfam,
	  params.prop_kidera,
	  params.annotation_file,
	 "KF",
	 ["IMGT_position", "cells","vFamilyName"])



    def K = params.K ?: 200
    def M = params.N ?: 10000
    def cols_clonotype = ["aminoAcid", "vFamilyName","jGeneName"]
	

	MERGE_AMINO_VFAM.out.merged
    .filter { file -> file.name =~ /.*productive.*/ }
    .set { merged_file }

	RAREFIED_SHARING (cols_clonotype, merged_file, params.annotation_file, K, M)
	RAREFIED_SHARING_VFAM (["aminoAcid", "vFamilyName"], merged_file, params.annotation_file, K, M)
	RAREFIED_SHARING_AA (["aminoAcid"], merged_file, params.annotation_file, K, M)

	def rarefy=true
	def sequence_cols= params.index_columns
	def count_col = 'count (templates/reads)'
	def min_count = params.N
	def iterations = params.M

	RAREFIED_DIVERSITY_METRICS(INITIAL_CLEANUP_SPLIT.out.productive,
        iterations,
        sequence_cols,
        count_col,
        rarefy,
        min_count)

	    // Collect just the summary files (drop sample_id)
    all_summaries = RAREFIED_DIVERSITY_METRICS.out.diversity_metrics.map { sample_id, file -> file }.collect()    
    COMBINE_RAREFIED_DIVERSITY_METRICS(all_summaries)


	
		INITIAL_CLEANUP_SPLIT.out.nonproductive
        .map { sample_name, file_path -> file_path }  // Extract just the files
        .collect()  // Collect all files into a single list
        .map { files -> 
            tuple(files, params.index_columns, "count (templates/reads)","nonproductive") 
        }
        .set { nonprod_merge_ch }

		    MERGE_AMINO_VFAM_NONPROD(nonprod_merge_ch)

		//subsample Vfam, Jgene, CDR3
		MERGE_AMINO_VFAM_NONPROD.out.merged
			.map { merged_file -> tuple(merged_file, 5000, params.M, 3, "nonproductive") } // Example values for N, M, num_index_columns .set { subsample_input_ch }
			.set { nonprod_subsample_input_ch }
		//sample M times N sequences from each sample
		SUBSAMPLE_FROM_MERGED_NONPROD( nonprod_subsample_input_ch)

		SUBSAMPLE_FROM_MERGED_NONPROD.out.subsampled_files
		.flatten()
    	.map { subsampled_file -> 
        def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
        tuple(subsampled_file, sample_name, "nonproductive", "aminoAcid", [ "vFamilyName","jGeneName"])
    } 
    .set { nonprod_compute_freqs_subsampled_input_ch }
	
	COMPUTE_FREQS_SUBSAMPLED_NONPROD(nonprod_compute_freqs_subsampled_input_ch)

	COMPUTE_FREQS_SUBSAMPLED_NONPROD.out.subsampled
		.map { big_file, summ_file -> summ_file}
		.collect()
		.map{ files ->
		tuple(files,  "nonproductive_subs_freqs",'nonproductive') }
		.set { nonprod_combine_freqs_subsampled_input_ch  }	

	COMBINE_FREQS_SUBSAMPLED_NONPROD(nonprod_combine_freqs_subsampled_input_ch)




	//Anderson tests for length differences

	    // Define common parameters
    def sample_col_val = params.sample_col ?: 'shortname'
    def n_perm_val = params.N ?: 10000
	def type = "prod"
    params.annotation_file = "${projectDir}/data/collated_info.csv"
    // Create a channel with your freqs_dir (this would come from an earlier process)
    freqs_dir_ch = channel.fromPath("${params.outdir}/freqs_subsampled/productive", type: 'dir')
    
    // Define your different analysis configurations
    // Format: [group_col, cells, split_col, groups]
    def analysis_configs = channel.of(
        ['genotype_short', 'N', null, null],                    // ad_length_results_N_geno.tsv
        ['genotype_short', 'N', 'source', null],                // ad_length_results_N_geno_byTF.tsv
        ['source', 'E', 'genotype_short', ['T', 'F']],          // ad_length_results_E_TF_bygeno.tsv
        ['source', 'N', null, null],                             // ad_length_results_N_TF.tsv
        ['genotype_short', 'E', null, null],                    // ad_length_results_E_geno.tsv
        ['genotype_short', 'E', 'source', null]                 // ad_length_results_E_geno_byTF.tsv
    )
    
    // Run the process for each configuration
    AD_LENGTH_DISTRIBUTION(
        freqs_dir_ch,
        analysis_configs.map { cfg -> cfg[0] },  // group_col
        channel.of(sample_col_val),
        analysis_configs.map { cfg -> cfg[1] },  // cells
        channel.of(n_perm_val),
        analysis_configs.map { cfg -> [cfg[2], cfg[3]] },  // [split_col, groups]
		channel.of(type),
		params.annotation_file
    )
    
    // Collect all results
    AD_LENGTH_DISTRIBUTION.out.results.collect().view { files ->
        "Generated AD length distribution results: ${files}"
    }



	}






// Print workflow execution summary 
workflow.onComplete {
summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
results     : ${params.outdir}

=======================================================================================
  """
println summary

}
