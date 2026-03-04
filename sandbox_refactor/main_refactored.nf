#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.input = './data/collated_info.csv'
params.outdir = './results'
params.seqs_to_remove = './data/D20210208D_1-overrepresented-sequences.txt'
params.N = 10000 
params.M = 25
params.index_columns = ["aminoAcid", "vFamilyName", "jGeneName"] // to join samples in merging
params.min_len = 8 // CDR3 length filters for IMGT freqs
params.max_len = 25  // CDR3 length filters for IMGT freqs
params.annotation_file = './data/collated_info.csv' // Annotation file 
params.min_non_index = 5
params.min_value = 1
params.prop_VHS = 'assets/VHSE.tsv'
params.prop_kidera = 'assets/kidera_factors.tsv'
params.K = 200
params.imgt_index_columns = ["IMGT_position", "AA", "aminoAcid_length"]
params.imgt_index_columns_WL = ["IMGT_position", "AA"]

// Include initial cleanup module
include { INITIAL_CLEANUP_SPLIT } from './modules/read_sheet.nf'

// Include sub-workflows
include { DIVERSITY_METRICS_WORKFLOW } from './subworkflows/diversity_metrics.nf'
include { FREQUENCY_ANALYSIS_WORKFLOW } from './subworkflows/frequency_analysis.nf'
include { MERGING_SUBSAMPLING_WORKFLOW } from './subworkflows/merging_subsampling.nf'
include { IMGT_AA_FREQUENCY_WORKFLOW } from './subworkflows/imgt_aa_frequency.nf'
include { AA_PROPERTIES_TESTING_WORKFLOW } from './subworkflows/aa_properties_testing.nf'
include { RAREFIED_SHARING_WORKFLOW } from './subworkflows/rarefied_sharing.nf'

// Print a header 
log.info """
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
    // Log parameters
    log.info "seqs_to_remove: ${params.seqs_to_remove}"
    log.info "DEBUG: params.input = '${params.input}'"
    log.info "DEBUG: params.input type = ${params.input?.getClass()}"

    // Channel 1: sample_short + genotype_short
    channel.fromPath(params.input)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_short, row.genotype_short) }
        .set { sample_genotype_ch }

    // Channel 2: sample_short + filepath
    channel.fromPath(params.input)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row.sample_short, "${row.LOC}/${row.SAMPLE}", params.seqs_to_remove, params.outdir) }
        .set { sample_filepath_ch }

    // Initial cleanup and split
    INITIAL_CLEANUP_SPLIT(sample_filepath_ch)

    // Sub-workflow 1: Diversity Metrics
    DIVERSITY_METRICS_WORKFLOW(
        INITIAL_CLEANUP_SPLIT.out.productive,
        INITIAL_CLEANUP_SPLIT.out.nonproductive
    )

    // Sub-workflow 2: Frequency Analysis
    FREQUENCY_ANALYSIS_WORKFLOW(
        INITIAL_CLEANUP_SPLIT.out.productive,
        INITIAL_CLEANUP_SPLIT.out.nonproductive
    )

    // Sub-workflow 3: Merging and Subsampling
    MERGING_SUBSAMPLING_WORKFLOW(
        INITIAL_CLEANUP_SPLIT.out.productive,
        params.index_columns,
        "count (templates/reads)",
        params.N,
        params.M
    )

    // Sub-workflow 4: IMGT AA Frequency Analysis
    IMGT_AA_FREQUENCY_WORKFLOW(
        MERGING_SUBSAMPLING_WORKFLOW.out.subsampled_files,
        params.min_len,
        params.max_len,
        params.annotation_file,
        params.min_non_index,
        params.min_value
    )

    // Sub-workflow 5: AA Properties Testing
    AA_PROPERTIES_TESTING_WORKFLOW(
        IMGT_AA_FREQUENCY_WORKFLOW.out.combined_subs_rows_full,
        IMGT_AA_FREQUENCY_WORKFLOW.out.combined_subs_rows_WL,
        IMGT_AA_FREQUENCY_WORKFLOW.out.combined_subs_rows_vfam,
        IMGT_AA_FREQUENCY_WORKFLOW.out.combined_subs_rows_WL_vfam,
        params.prop_VHS,
        params.prop_kidera,
        params.annotation_file
    )

    // Sub-workflow 6: Rarefied Sharing
    MERGING_SUBSAMPLING_WORKFLOW.out.merged_file
        .filter { file -> file.name =~ /.*productive.*/ }
        .set { merged_file }

    RAREFIED_SHARING_WORKFLOW(
        merged_file,
        params.annotation_file,
        params.K,
        params.N
    )
}

// Print workflow execution summary 
workflow.onComplete {
    def summary = """
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
