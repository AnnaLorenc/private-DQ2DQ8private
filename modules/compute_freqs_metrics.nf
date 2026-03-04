process COMPUTE_FREQS {
    tag "Computing V, J, CDR3 len ferquencies from ${sample_name}"
    publishDir "${params.outdir}/freqs/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), path(sample_file)

    output:
    tuple val(sample_name), path("${sample_name}_*_freqs.csv"), emit: freqs

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
python ${projectDir}/bin/VJlen_freqs/compute_freqs.py \
    --file_path ${sample_file} \
    --sample ${sample_name}_\${sample_type} \
    --output_loc .
    """
}