process PREP_COMB_IMGT_TEST {
    tag "preparing for IMGT test"
    publishDir "${params.outdir}/imgt_aa_fortest/", mode: 'copy'

    input:
    tuple path(freqs), path(counts)
    val(index_columns)
    val(annotation_file)
    val(min_non_index)
    val(min_value)



    output:
    path "*_test_ready.tsv", emit: imgt_test_input

    script:
    def index_cols_arg = index_columns ? "--index_columns ${index_columns.join(' ')}" : ""

        // Get base name of the freqs file (without extension)
    def base_name = freqs.baseName

    // Construct an output filename using that base name
    def out_name = "${base_name}_${min_non_index}_${min_value}_test_ready.tsv"

    """
    # Filter to have at least min_non_index samples with at least min_value counts of given position (established based on --filter_file with counts)
    # Then reshape into long format: Index columns, celltype,
    #combine with the annotation file 
    python ${projectDir}/bin/imgtCDR3/prepare_for_imgt_test.py \
        --input_file ${freqs} \
        --annotation_file ${projectDir}/${annotation_file} \
        --output ${out_name} \
        --filter_file ${counts} \
        ${index_cols_arg} \
        --min_non_index ${min_non_index} \
        --min_value ${min_value}
    """
}