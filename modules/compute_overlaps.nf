process COMPUTE_OVERLAPS {
    tag "overlaps"
    publishDir "${params.outdir}/overlaps_raw/", mode: 'copy'



    input:
    tuple path(input_file),  val(type), val(sample_start_column)

    output:
    path ("overlaps_*.tsv"), emit: overlaps
   


"""
python ${projectDir}/bin/overlaps/compute_overlaps.py \
    --input_file ${input_file} \
    --output_file overlaps_${type}.tsv \
    --sample_start_column ${sample_start_column} 
"""
}

    // parser.add_argument("--input_file", required=True)
    // parser.add_argument("--output_file", required=True)
    // parser.add_argument("--sample_start_column", type=int, required=True)
    // parser.add_argument("--sample_end_column", type=int, default=None)
    // parser.add_argument(
    //     "--id_columns",
    //     nargs="+",
    //     help="Columns to group by and summarize before processing",
    // )