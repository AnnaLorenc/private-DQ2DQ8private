#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AA_PROPERTIES_TEST } from '../modules/aa_properties.nf'
include { AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL } from '../modules/aa_properties.nf'
include { AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_VFAM } from '../modules/aa_properties.nf'
include { AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL_VFAM } from '../modules/aa_properties.nf'
include { AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_kidera } from '../modules/aa_properties.nf'
include { AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL_kidera } from '../modules/aa_properties.nf'
include { AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_VFAM_kidera } from '../modules/aa_properties.nf'
include { AA_PROPERTIES_TEST as AA_PROPERTIES_TEST_WL_VFAM_kidera } from '../modules/aa_properties.nf'

workflow AA_PROPERTIES_TESTING_WORKFLOW {
    take:
    combined_subs_rows_full      // channel
    combined_subs_rows_WL        // channel
    combined_subs_rows_vfam      // channel
    combined_subs_rows_WL_vfam   // channel
    prop_VHS                     // path
    prop_kidera                  // path
    annotation_file              // path

    main:
    // Define index columns
    def index_cols = ["IMGT_position", "AA", "aminoAcid_length"]
    def index_cols_WL = ["IMGT_position", "AA"]
    def index_cols_vfam = ["IMGT_position", "AA", "aminoAcid_length", "vFamilyName"]
    def index_cols_WL_vfam = ["IMGT_position", "AA", "vFamilyName"]

    // VHS properties testing - with length, no Vfam
    AA_PROPERTIES_TEST(
        combined_subs_rows_full,
        index_cols,
        prop_VHS,
        annotation_file,
        "VHS",
        ["IMGT_position", "cells", "length"]
    )

    // VHS properties testing - without length, no Vfam
    AA_PROPERTIES_TEST_WL(
        combined_subs_rows_WL,
        index_cols_WL,
        prop_VHS,
        annotation_file,
        "VHS",
        ["IMGT_position", "cells"]
    )

    // VHS properties testing - with length, Vfam
    AA_PROPERTIES_TEST_VFAM(
        combined_subs_rows_vfam,
        index_cols_vfam,
        prop_VHS,
        annotation_file,
        "VHS",
        ["IMGT_position", "cells", "length", "vFamilyName"]
    )

    // VHS properties testing - without length, Vfam
    AA_PROPERTIES_TEST_WL_VFAM(
        combined_subs_rows_WL_vfam,
        index_cols_WL_vfam,
        prop_VHS,
        annotation_file,
        "VHS",
        ["IMGT_position", "cells", "vFamilyName"]
    )

    // Kidera factors testing - with length, no Vfam
    AA_PROPERTIES_TEST_kidera(
        combined_subs_rows_full,
        index_cols,
        prop_kidera,
        annotation_file,
        "KF",
        ["IMGT_position", "cells", "length"]
    )

    // Kidera factors testing - without length, no Vfam
    AA_PROPERTIES_TEST_WL_kidera(
        combined_subs_rows_WL,
        index_cols_WL,
        prop_kidera,
        annotation_file,
        "KF",
        ["IMGT_position", "cells"]
    )

    // Kidera factors testing - with length, Vfam
    AA_PROPERTIES_TEST_VFAM_kidera(
        combined_subs_rows_vfam,
        index_cols_vfam,
        prop_kidera,
        annotation_file,
        "KF",
        ["IMGT_position", "cells", "length", "vFamilyName"]
    )

    // Kidera factors testing - without length, Vfam
    AA_PROPERTIES_TEST_WL_VFAM_kidera(
        combined_subs_rows_WL_vfam,
        index_cols_WL_vfam,
        prop_kidera,
        annotation_file,
        "KF",
        ["IMGT_position", "cells", "vFamilyName"]
    )

    emit:
    vhs_full = AA_PROPERTIES_TEST.out
    vhs_WL = AA_PROPERTIES_TEST_WL.out
    vhs_vfam = AA_PROPERTIES_TEST_VFAM.out
    vhs_WL_vfam = AA_PROPERTIES_TEST_WL_VFAM.out
    kidera_full = AA_PROPERTIES_TEST_kidera.out
    kidera_WL = AA_PROPERTIES_TEST_WL_kidera.out
    kidera_vfam = AA_PROPERTIES_TEST_VFAM_kidera.out
    kidera_WL_vfam = AA_PROPERTIES_TEST_WL_VFAM_kidera.out
}
