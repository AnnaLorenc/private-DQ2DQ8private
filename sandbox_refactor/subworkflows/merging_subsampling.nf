#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MERGE_SAMPLES as MERGE_AMINO_VFAM } from '../modules/merge_samples.nf'
include { SUBSAMPLE_FROM_MERGED } from '../modules/subsample_from_merged.nf'
include { COMPUTE_FREQS_SUBSAMPLED } from '../modules/compute_freqs_subsampled.nf'
include { COMBINE_FREQS_SUBSAMPLED } from '../modules/combine_freqs_subsampled.nf'
include { COMPUTE_OVERLAPS } from '../modules/compute_overlaps.nf'
include { COMPUTE_OVERLAPS_SUBSAMPLED } from '../modules/compute_overlaps_subsampled.nf'

workflow MERGING_SUBSAMPLING_WORKFLOW {
    take:
    productive_files_ch // channel: files
    index_columns       // list
    count_column        // string
    N                   // integer
    M                   // integer

    main:
    // Prepare merge input
    def prod_merge_ch = productive_files_ch
        .map { sample_name, file_path -> file_path }  // Extract just the files
        .collect()  // Collect all files into a single list
        .map { files -> 
            tuple(files, index_columns, count_column, "productive") 
        }

    // Merge samples
    MERGE_AMINO_VFAM(prod_merge_ch)

    // Subsample
    def first_column_with_sample_in_merged = index_columns.size()
    
    MERGE_AMINO_VFAM.out.merged
        .map { merged_file -> tuple(merged_file, N, M, 3, "productive") }
        .set { subsample_input_ch }

    SUBSAMPLE_FROM_MERGED(subsample_input_ch)

    // Compute frequencies on subsampled data
    SUBSAMPLE_FROM_MERGED.out.subsampled_files
        .flatten()
        .map { subsampled_file -> 
            def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
            tuple(subsampled_file, sample_name, "productive", "aminoAcid", ["vFamilyName", "jGeneName"])
        }
        .set { compute_freqs_subsampled_input_ch }
    
    COMPUTE_FREQS_SUBSAMPLED(compute_freqs_subsampled_input_ch)

    // Combine subsampled frequencies
    COMPUTE_FREQS_SUBSAMPLED.out.subsampled
        .map { big_file, summ_file -> summ_file }
        .collect()
        .map { files ->
            tuple(files, "productive_subs_freqs.tsv")
        }
        .set { combine_freqs_subsampled_input_ch }

    COMBINE_FREQS_SUBSAMPLED(combine_freqs_subsampled_input_ch)

    // Compute overlaps on merged data
    MERGE_AMINO_VFAM.out.merged
        .map { merged_file -> tuple(merged_file, "productive", first_column_with_sample_in_merged) }
        .set { compute_overlaps_input_ch }
    
    COMPUTE_OVERLAPS(compute_overlaps_input_ch)

    // Compute overlaps on subsampled data
    SUBSAMPLE_FROM_MERGED.out.sample_names
        .collect()
        .combine(SUBSAMPLE_FROM_MERGED.out.subsampled_dir.collect())
        .map { sample_names, subsampled_dirs -> 
            tuple(sample_names, subsampled_dirs, "productive")
        }
        .set { compute_overlaps_subsampled_input_ch }
    
    COMPUTE_OVERLAPS_SUBSAMPLED(compute_overlaps_subsampled_input_ch)

    emit:
    merged_file = MERGE_AMINO_VFAM.out.merged
    subsampled_files = SUBSAMPLE_FROM_MERGED.out.subsampled_files
    combined_freqs = COMBINE_FREQS_SUBSAMPLED.out
    overlaps = COMPUTE_OVERLAPS.out
    overlaps_subsampled = COMPUTE_OVERLAPS_SUBSAMPLED.out
}
