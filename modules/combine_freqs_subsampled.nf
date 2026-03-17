process COMBINE_FREQS_SUBSAMPLED {
    tag "joined frequencies of subsampled"
    publishDir "${params.outdir}/combined_freqs/${type}", mode: 'copy'



    input:
    tuple val(input_files), val(output_file), val(type)

    output:
    path ("*.tsv"), emit: all_ssamples_freqs
   

    script:
    def input_files_all = input_files.join(' ')
    """
    python ${projectDir}/bin/VJlen_freqs/combine_freqs_subsampled.py \
        --input_files ${input_files_all} \
        --output_file ${output_file}.tsv
    """
}