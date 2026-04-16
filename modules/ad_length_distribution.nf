#!/usr/bin/env nextflow

nextflow.preview.strict = true

/*
 * Process for computing Anderson test for length distributions 
 * Handles multiple parameter combinations for different analysis types
 */


process AD_LENGTH_DISTRIBUTION {
    tag "Anderson test"
    publishDir "${params.outdir}/ad_length_results/${type}", mode: 'copy'
    
    input:
    path freqs_dir
    val group_col
    val sample_col
    val cells
    val n_perm
    tuple val(split_col), val(groups)
    val type
    val annotation_file
    
    output:
    path "${output_name}", emit: results
    
    script:
    // Build output filename based on parameters
    def split_suffix = split_col ? "_by${split_col}" : ""
    def groups_suffix = groups ? "_${groups.join('')}" : ""
    output_name = "ad_length_results_${cells}_${group_col}${split_suffix}${groups_suffix}.tsv"
    
    // Build optional arguments
    def split_arg = split_col ? "--split_col ${split_col}" : ""
    def groups_arg = groups ? "--groups ${groups.join(' ')}" : ""
    
    """
    python ${projectDir}/bin/ad_length_distribution.py \\
        --freqs_dir ${freqs_dir} \\
        --sample_map ${annotation_file} \\
        --group_col ${group_col} \\
        --sample_col ${sample_col} \\
        --cells ${cells} \\
        --n_perm ${n_perm} \\
        ${split_arg} \\
        ${groups_arg} \\
        --output ${output_name}
    """
}