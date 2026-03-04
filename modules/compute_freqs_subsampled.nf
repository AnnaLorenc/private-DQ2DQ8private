process COMPUTE_FREQS_SUBSAMPLED {
    tag "frequencies of subsampled"
    publishDir "${params.outdir}/freqs_subsampled/", mode: 'copy'



    input:
    tuple path(input_file), val(sample_name), val(type), val(cdr3_column), val(index_columns)

    output:
    tuple path ("${type}/*.tsv.gz"),path ("${type}/*_summ.tsv") , emit: subsampled
   

    script:
    def index_columns_arg = index_columns.join(' ')
    """
    python ${projectDir}/bin/VJlen_freqs/compute_freqs_subsampled.py \
        --input_file ${input_file} \
        --output_file ${type}/${sample_name}.tsv.gz \
        --index_columns ${index_columns_arg } \
        --cdr3_column ${cdr3_column}
    """
}

    // parser.add_argument("--input_file", required=True, help="Path to the input file (tsv.gz).")
    // parser.add_argument("--output_file", required=True, help="Path to save the output file.")
    // parser.add_argument("--index_columns", nargs='+', required=True, help="List of index columns to use as they are.")
    // parser.add_argument("--cdr3_column", required=True, help="Name of the CDR3 column.")