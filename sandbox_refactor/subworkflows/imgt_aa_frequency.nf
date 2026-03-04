#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COMPUTE_IMGT_AA_FREQS_SUBS } from '../modules/imgt_aa_freqs_subsampled.nf'
include { COMBINE_IMGT_AA_FREQS_MED } from '../modules/combine_IMGT_AA_freqs.nf'
include { COMBINE_IMGT_AA_FREQS_MED as COMBINE_IMGT_AA_FREQS_MED_WL } from '../modules/combine_IMGT_AA_freqs.nf'
include { COMBINE_IMGT_AA_FREQS_MED as COMBINE_IMGT_AA_FREQS_VFAM_MED_WL } from '../modules/combine_IMGT_AA_freqs.nf'
include { COMBINE_IMGT_AA_FREQS_MED as COMBINE_IMGT_AA_FREQS_VFAM_MED } from '../modules/combine_IMGT_AA_freqs.nf'
include { PREP_COMB_IMGT_TEST } from '../modules/prep_comb_IMGT_test.nf'
include { PREP_COMB_IMGT_TEST as PREP_COMB_IMGT_TEST_WL } from '../modules/prep_comb_IMGT_test.nf'
include { PREP_COMB_IMGT_TEST as PREP_COMB_IMGT_TEST_VFAM } from '../modules/prep_comb_IMGT_test.nf'
include { PREP_COMB_IMGT_TEST as PREP_COMB_IMGT_TEST_WL_VFAM } from '../modules/prep_comb_IMGT_test.nf'
include { RUN_IMGT_TEST } from '../modules/run_imgt_test.nf'
include { RUN_IMGT_TEST as RUN_IMGT_TEST_WL } from '../modules/run_imgt_test.nf'
include { RUN_IMGT_TEST as RUN_IMGT_TEST_VFAM } from '../modules/run_imgt_test.nf'
include { RUN_IMGT_TEST as RUN_IMGT_TEST_WL_VFAM } from '../modules/run_imgt_test.nf'

workflow IMGT_AA_FREQUENCY_WORKFLOW {
    take:
    subsampled_files_ch // channel: files
    min_len             // integer
    max_len             // integer
    annotation_file     // path
    min_non_index       // integer
    min_value           // integer

    main:
    // Compute IMGT AA frequencies on subsampled data
    subsampled_files_ch
        .flatten()
        .map { subsampled_file -> 
            def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
            tuple(subsampled_file, sample_name, "productive", min_len, max_len)
        }
        .set { compute_imgt_freqs_subsampled_input_ch }

    COMPUTE_IMGT_AA_FREQS_SUBS(compute_imgt_freqs_subsampled_input_ch)

    // Define index columns
    def index_cols = ["IMGT_position", "AA", "aminoAcid_length"]
    def index_cols_WL = ["IMGT_position", "AA"]
    def index_cols_vfam = ["IMGT_position", "AA", "aminoAcid_length", "vFamilyName"]
    def index_cols_WL_vfam = ["IMGT_position", "AA", "vFamilyName"]

    // With length, no Vfam
    def freq_med_files_ch = COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_full_med
        .map { sample_name, file -> file }
        .collect()

    COMBINE_IMGT_AA_FREQS_MED(freq_med_files_ch, index_cols, "_full")

    def combined = COMBINE_IMGT_AA_FREQS_MED.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_MED.out.combined_subs_counts)

    PREP_COMB_IMGT_TEST(combined, index_cols, annotation_file, min_non_index, min_value)
    
    def index_cols1 = ["IMGT_position", "AA", "length"]
    RUN_IMGT_TEST(PREP_COMB_IMGT_TEST.out.imgt_test_input, index_cols1)

    // Without length, no Vfam
    def freq_med_WL_files_ch = COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_WL_med
        .map { sample_name, file -> file }
        .collect()

    COMBINE_IMGT_AA_FREQS_MED_WL(freq_med_WL_files_ch, index_cols_WL, "_WL")

    def combined_WL = COMBINE_IMGT_AA_FREQS_MED_WL.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_MED_WL.out.combined_subs_counts)

    PREP_COMB_IMGT_TEST_WL(combined_WL, index_cols_WL, annotation_file, min_non_index, min_value)
    RUN_IMGT_TEST_WL(PREP_COMB_IMGT_TEST_WL.out.imgt_test_input, index_cols_WL)

    // With length, Vfam
    def freq_med_vfam_files_ch = COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_full_Vfam_med
        .map { sample_name, file -> file }
        .collect()

    COMBINE_IMGT_AA_FREQS_VFAM_MED(freq_med_vfam_files_ch, index_cols_vfam, "_full_vfam")

    def combined_vfam = COMBINE_IMGT_AA_FREQS_VFAM_MED.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_VFAM_MED.out.combined_subs_counts)

    PREP_COMB_IMGT_TEST_VFAM(combined_vfam, index_cols_vfam, annotation_file, min_non_index, min_value)
    
    def index_cols1_vfam = ["IMGT_position", "AA", "length", "vFamilyName"]
    RUN_IMGT_TEST_VFAM(PREP_COMB_IMGT_TEST_VFAM.out.imgt_test_input, index_cols1_vfam)

    // Without length, Vfam
    def freq_med_WL_vfam_files_ch = COMPUTE_IMGT_AA_FREQS_SUBS.out.aa_imgt_freq_WL_Vfam_med
        .map { sample_name, file -> file }
        .collect()

    COMBINE_IMGT_AA_FREQS_VFAM_MED_WL(freq_med_WL_vfam_files_ch, index_cols_WL_vfam, "_WL_vfam")

    def combined_WL_vfam = COMBINE_IMGT_AA_FREQS_VFAM_MED_WL.out.combined_subs_rows_freq
        .combine(COMBINE_IMGT_AA_FREQS_VFAM_MED_WL.out.combined_subs_counts)

    PREP_COMB_IMGT_TEST_WL_VFAM(combined_WL_vfam, index_cols_WL_vfam, annotation_file, min_non_index, min_value)
    RUN_IMGT_TEST_WL_VFAM(PREP_COMB_IMGT_TEST_WL_VFAM.out.imgt_test_input, index_cols_WL_vfam)

    emit:
    combined_subs_rows_full = COMBINE_IMGT_AA_FREQS_MED.out.combined_subs_rows
    combined_subs_rows_WL = COMBINE_IMGT_AA_FREQS_MED_WL.out.combined_subs_rows
    combined_subs_rows_vfam = COMBINE_IMGT_AA_FREQS_VFAM_MED.out.combined_subs_rows
    combined_subs_rows_WL_vfam = COMBINE_IMGT_AA_FREQS_VFAM_MED_WL.out.combined_subs_rows
    test_results_full = RUN_IMGT_TEST.out
    test_results_WL = RUN_IMGT_TEST_WL.out
    test_results_vfam = RUN_IMGT_TEST_VFAM.out
    test_results_WL_vfam = RUN_IMGT_TEST_WL_VFAM.out
}
