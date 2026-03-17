process SUBSAMPLE_FROM_MERGED {
    tag "subsampling ${type}"
    publishDir "${params.outdir}/subsampled/${type}", mode: 'copy'

    input:
    tuple val(input_file), val(N), val(M), val(num_index_columns), val(type)

    output:
    path "${type}/*_subsampled.tsv.gz", emit: subsampled_files
    path "sample_list.txt", emit: sample_names
    path "${type}/", emit: subsampled_dir
   
    script:
    """
    python ${projectDir}/bin/subsample_merged_sequences.py \
        --input_file ${input_file} \
        --output_dir ${type} \
        --N  ${N} \
        --M ${M} \
        --num_index_columns ${num_index_columns}
    
    # Extract sample names and create a list
    for file in ${type}/*_subsampled.tsv.gz; do
        basename "\$file" _subsampled.tsv.gz
    done > sample_list.txt
    """
}

        // input_file (str): Path to the input merged file (tsv.gz).
        // output_dir (str): Directory to save the output files.
        // N (int): Number of indexes to sample per column.
        // M (int): Number of repetitions for subsampling.
        // num_index_columns (int): Number of columns to treat as indexes.