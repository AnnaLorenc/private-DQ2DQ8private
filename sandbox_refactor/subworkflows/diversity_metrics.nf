#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_PROD } from '../modules/extract_diversity_metrics.nf'
include { EXTRACT_DIVERSITY_METRICS as EXTRACT_DIVERSITY_METRICS_NONPROD } from '../modules/extract_diversity_metrics.nf'
include { COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_PROD } from '../modules/combine_diversity_metrics.nf'
include { COMBINE_DIVERSITY_METRICS as COMBINE_DIVERSITY_METRICS_NONPROD } from '../modules/combine_diversity_metrics.nf'

workflow DIVERSITY_METRICS_WORKFLOW {
    take:
    productive_ch    // channel: tuple(sample, file)
    nonproductive_ch // channel: tuple(sample, file)

    main:
    // Extract diversity metrics
    EXTRACT_DIVERSITY_METRICS_PROD(productive_ch)
    EXTRACT_DIVERSITY_METRICS_NONPROD(nonproductive_ch)

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

    emit:
    productive_combined = COMBINE_DIVERSITY_METRICS_PROD.out
    nonproductive_combined = COMBINE_DIVERSITY_METRICS_NONPROD.out
}
