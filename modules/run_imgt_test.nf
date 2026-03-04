process RUN_IMGT_TEST{
    tag "running IMGT test"
    publishDir "${params.outdir}/imgt_aa_test/", mode: 'copy'

    input:
    path(imgt_test_input_file)
    val(index_columns)



    output:
    path "*_res/*", emit: imgt_tests

    script:

    def base_name = imgt_test_input_file.baseName.replaceAll('_test_ready$', '')
    def output_dir = "${base_name}_res"
    def index_cols_arg = index_columns ? "--index_columns ${index_columns.join(' ')}" : ""

    """
    # Filter to have at least min_non_index samples with at least min_value counts of given position (established based on --filter_file with counts)
    # Then reshape into long format: Index columns, celltype,
    #combine with the annotation file 
    python ${projectDir}/bin/imgtCDR3/test_imgt.py \
        --input ${imgt_test_input_file} \
        --output_dir ${output_dir} \
        ${index_cols_arg}

    """
}