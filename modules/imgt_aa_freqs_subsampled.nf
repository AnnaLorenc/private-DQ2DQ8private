process COMPUTE_IMGT_AA_FREQS_SUBS {
    tag "imgt_aa_freqs_${sample_name}"
    publishDir "${params.outdir}/imgt_aa/", mode: 'copy'

    input:
    tuple path(subsampled_file), val(sample_name), val(type), val(min_len), val(max_len)

    output:
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_invalid.tsv"), emit: invalid_CDR3
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq.tsv"), emit: aa_imgt_freq_full  
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_valid.tsv"), emit: exploded_df
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq_WL.tsv"), emit: aa_imgt_freq_WL
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq_med.tsv"), emit: aa_imgt_freq_full_med
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq_WL_med.tsv"), emit: aa_imgt_freq_WL_med 
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq_Vfam_med.tsv"), emit: aa_imgt_freq_full_Vfam_med
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq_WL_Vfam_med.tsv"), emit: aa_imgt_freq_WL_Vfam_med  
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq_Vfam.tsv"), emit: aa_imgt_freq_full_Vfam
    tuple val(sample_name), path("${sample_name}_aa_freqs/${sample_name}_freq_WL_Vfam.tsv"), emit: aa_imgt_freq_WL_Vfam  

    script:
    """
    echo "Processing file: ${subsampled_file}"
    echo "Sample name: ${sample_name}"
    echo "File contents preview:"
    gunzip -c ${subsampled_file} | head -5
    
    # Run the first script
    python ${projectDir}/bin/imgtCDR3/imgt_freqs.py \
        --input ${subsampled_file} \
        --min_len ${min_len} \
        --max_len ${max_len} \
        --output ${sample_name}

    echo "Files created by imgtCDR3/imgt_freqs.py:"
    ls -la ${sample_name}_aa_freqs/
    
    echo "Content of freq_WL file:"
    head -5 ${sample_name}_aa_freqs/${sample_name}_freq_WL.tsv
    
    # Run median computation scripts
    python ${projectDir}/bin/imgtCDR3/compute_medians.py \
        --input ${sample_name}_aa_freqs/${sample_name}_freq_WL.tsv \
        --sample_name ${sample_name} \
        --output ${sample_name}_aa_freqs/${sample_name}_freq_WL_med.tsv \
        --index_columns IMGT_position AA    
    
    python ${projectDir}/bin/imgtCDR3/compute_medians.py \
        --input ${sample_name}_aa_freqs/${sample_name}_freq.tsv \
        --output ${sample_name}_aa_freqs/${sample_name}_freq_med.tsv  \
        --sample_name ${sample_name} \

    python ${projectDir}/bin/imgtCDR3/compute_medians.py \
    --input ${sample_name}_aa_freqs/${sample_name}_freq_WL_Vfam.tsv \
    --sample_name ${sample_name} \
    --output ${sample_name}_aa_freqs/${sample_name}_freq_WL_Vfam_med.tsv \
    --index_columns IMGT_position AA vFamilyName   
    
    python ${projectDir}/bin/imgtCDR3/compute_medians.py \
        --input ${sample_name}_aa_freqs/${sample_name}_freq_Vfam.tsv \
        --output ${sample_name}_aa_freqs/${sample_name}_freq_Vfam_med.tsv  \
        --sample_name ${sample_name} \
        --index_columns IMGT_position aminoAcid_length AA vFamilyName  

    """
}

    // parser.add_argument("--input", required=True)
    // parser.add_argument("--freq_output", required=False)
    // parser.add_argument("--valid_output", required=False, help="Path to save valid aminoAcid sequences")
    // parser.add_argument("--invalid_output", required=False, help="Path to save invalid sequences")
    // parser.add_argument("--min_len", type=int, default=8)
    // parser.add_argument("--max_len", type=int, default=25)
// SUBSAMPLE_FROM_MERGED.out.subsampled_files
// 		.map { subsampled_file -> 
// 		def sample_name = subsampled_file.baseName.replaceAll(/_subsampled\.tsv$/, '')
// 		tuple(subsampled_file, sample_name, "productive", params.min_len, params.max_len)