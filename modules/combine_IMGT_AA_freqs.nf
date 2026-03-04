process COMBINE_IMGT_AA_FREQS_MED {
    tag "combine_imgt_aa_freqs_med"
    publishDir "${params.outdir}/imgt_aa_combined/", mode: 'copy'

    input:
    path(freq_med_files)
    val(index_columns)
    val(suffix)

    output:
    path "comb_aa_imgt${suffix}_orig_counts.tsv",  emit: combined_orig_counts
    path "comb_aa_imgt${suffix}_orig_counts_freq.tsv", emit: combined_orig_counts_freq
    path "comb_aa_imgt${suffix}_orig_rows.tsv", emit: combined_orig_rows
    path "comb_aa_imgt${suffix}_orig_rows_freq.tsv", emit: combined_orig_rows_freq
    path "comb_aa_imgt${suffix}_subs_counts.tsv", emit: combined_subs_counts
    path "comb_aa_imgt${suffix}_subs_counts_freq.tsv", emit: combined_subs_counts_freq
    path "comb_aa_imgt${suffix}_subs_rows.tsv", emit: combined_subs_rows
    path "comb_aa_imgt${suffix}_subs_rows_freq.tsv", emit: combined_subs_rows_freq


    script:
    def index_cols_arg = index_columns ? "--index_columns ${index_columns.join(' ')}" : ""
    """
    echo "Processing ${freq_med_files.size()} frequency median files:"
    ls -la *.tsv
    
    # Combine and extract subset of columns
    python ${projectDir}/bin/imgtCDR3/combine_and_extract_imgt_freqs.py \
        --input_files ${freq_med_files.join(' ')} \
        --output_prefix comb_aa_imgt${suffix} \
        ${index_cols_arg}
    """
}