process COMBINE_RAREFIED_DIVERSITY_METRICS {
   tag "Collecting rarefied diversity metrics}"
   publishDir "${params.outdir}/diversity_metrics_rarefied/", mode: 'copy'
    
    input:
    path summaries, stageAs: 'input_*.tsv'
    
    output:
    path 'diversities_summary.txt'
    
    script:
    """


       # Get header from first file
    head -n 1 ${summaries[0]} > diversities_summary.txt
    
    # Append all files (skipping their headers)
    for file in ${summaries}; do
        tail -n +2 \$file >> diversities_summary.txt
    done
    """
}
//   // cat ${summaries} > diversities_summary.txt