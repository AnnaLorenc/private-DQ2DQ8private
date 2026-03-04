process INITIAL_CLEANUP_SPLIT {
    tag "${sample}"  // Show sample name in logs
    publishDir "${params.outdir}/cleanup/${sample}/", mode: 'copy'


    input:
    tuple val(sample), val(filepath), val(seqs_to_remove), val(outdir)

    output:
    tuple val(sample), path("${sample}_productive.tsv.gz"), emit: productive
    tuple val(sample), path("${sample}_nonproductive.tsv.gz"), emit: nonproductive
    tuple val(sample), path("${sample}_cleanup_summary.tsv"), emit: summary

    script:
    """
    python  ${projectDir}/bin/preprocess/clean_split_immunoseq.py \
        --seqs_to_remove ${projectDir}/${seqs_to_remove} \
        --output_loc . \
        --sample ${sample} \
        --file_path ${filepath}
    """
}