#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COMPUTE_FREQS as COMPUTE_FREQS_PROD } from '../modules/compute_freqs_metrics.nf'
include { COMPUTE_FREQS as COMPUTE_FREQS_NONPROD } from '../modules/compute_freqs_metrics.nf'
include { COMBINE_FREQS as COMBINE_FREQS_PROD } from '../modules/combine_freqs.nf'
include { COMBINE_FREQS as COMBINE_FREQS_NONPROD } from '../modules/combine_freqs.nf'

workflow FREQUENCY_ANALYSIS_WORKFLOW {
    take:
    productive_ch    // channel: tuple(sample, file)
    nonproductive_ch // channel: tuple(sample, file)

    main:
    // Compute frequencies
    COMPUTE_FREQS_PROD(productive_ch)
    COMPUTE_FREQS_NONPROD(nonproductive_ch)

    // Collect productive CSV files for freqs
    productive_freqs = COMPUTE_FREQS_PROD.out.freqs
        .unique()  // Remove duplicates
        .toList()
        .map { items -> 
            def samples = items.collect { it[0] }
            def files = items.collect { it[1] }
            tuple('productive', samples, files)
        }
    
    // Collect nonproductive with deduplication
    nonproductive_freqs = COMPUTE_FREQS_NONPROD.out.freqs
        .unique()  // Remove duplicates
        .toList()
        .map { items -> 
            def samples = items.collect { it[0] }
            def files = items.collect { it[1] }
            tuple('nonproductive', samples, files)
        }

    // Combine frequencies
    COMBINE_FREQS_PROD(productive_freqs)
    COMBINE_FREQS_NONPROD(nonproductive_freqs)

    emit:
    productive_combined = COMBINE_FREQS_PROD.out
    nonproductive_combined = COMBINE_FREQS_NONPROD.out
}
