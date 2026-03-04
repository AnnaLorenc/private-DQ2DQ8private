#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { RAREFIED_SHARING } from '../modules/rarefied_sharing.nf'
include { RAREFIED_SHARING as RAREFIED_SHARING_VFAM } from '../modules/rarefied_sharing.nf'
include { RAREFIED_SHARING as RAREFIED_SHARING_AA } from '../modules/rarefied_sharing.nf'

workflow RAREFIED_SHARING_WORKFLOW {
    take:
    merged_file_ch   // channel: file (filtered productive merged file)
    annotation_file  // path
    K                // integer
    M                // integer

    main:
    // Define clonotype columns
    def cols_clonotype = ["aminoAcid", "vFamilyName", "jGeneName"]
    def cols_vfam = ["aminoAcid", "vFamilyName"]
    def cols_aa = ["aminoAcid"]

    // Run rarefied sharing with different column configurations
    RAREFIED_SHARING(cols_clonotype, merged_file_ch, annotation_file, K, M)
    RAREFIED_SHARING_VFAM(cols_vfam, merged_file_ch, annotation_file, K, M)
    RAREFIED_SHARING_AA(cols_aa, merged_file_ch, annotation_file, K, M)

    emit:
    sharing_full = RAREFIED_SHARING.out
    sharing_vfam = RAREFIED_SHARING_VFAM.out
    sharing_aa = RAREFIED_SHARING_AA.out
}
