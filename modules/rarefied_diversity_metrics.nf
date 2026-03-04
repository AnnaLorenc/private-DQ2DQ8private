process RAREFIED_DIVERSITY_METRICS {
    tag "Computing rarefied diversity metrics for sample ${sample_name}"
    publishDir "${params.outdir}/diversity_metrics_rarefied/per_sample/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(sample_file)
    val(iterations)
    val(sequence_cols)
    val(count_col)
    val(rarefy)
    val(min_count)

    output:
    tuple val(sample_name), path("${sample_name}_diversity.tsv"), emit: diversity_metrics

    script:
    def rarefy_flag = rarefy ? '--rarefy' : ''
    """
# Extract the type (productive or nonproductive) from the file name
if echo "${sample_file}" | grep -q "productive"; then
    if echo "${sample_file}" | grep -q "nonproductive"; then
        sample_type="nonproductive"
    else
        sample_type="productive"
    fi
else
    sample_type="unknown"
fi

# Run the Python script with the modified output name
python ${projectDir}/bin/diversity/rarefied_diversity_metrics.py \
    --file_path ${sample_file} \
    --group_by 'sample_short' \
    --iterations ${iterations} \
    --sequence_cols ${sequence_cols.join(',')} \
    --count_col "${count_col}" \
    ${rarefy_flag} \
    --min_count ${min_count} \
    --output ${sample_name}_diversity.tsv
    """
}




// workflow {
//     samples_ch = channel.fromPath('samples/*.txt')
//         .map { file -> tuple(file.baseName, file) }
    
//     summaries = SUMMARIZE_SAMPLE(samples_ch)
    
//     // Collect just the summary files (drop sample_id)
//     all_summaries = summaries.map { sample_id, file -> file }.collect()
    
//     COMBINE_SUMMARIES(all_summaries)
// }