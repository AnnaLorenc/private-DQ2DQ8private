process EXTRACT_DIVERSITY_METRICS {
    tag "Extracting diversity metrics for sample ${sample_name}"
    publishDir "${params.outdir}/diversity_metrics/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(sample_file)

    output:
    tuple val(sample_name), path("${sample_name}_*_diversity.csv"), emit: diversity_metrics

    script:
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
python ${projectDir}/bin/diversity/extract_diversity_metrics.py \
    --file_path ${sample_file} \
    --sample ${sample_name}_\${sample_type} \
    --output_loc .
    """
}