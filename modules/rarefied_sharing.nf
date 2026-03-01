process RAREFIED_SHARING {

    publishDir "${params.outdir}/rarified_sharing/", mode: 'copy'

    input:
    val index_cols
    path merged_file
    val annotation_file
    val K
    val M
    
    output: 
    path "groups_*.yml", emit: yml_files
    path "rare_groupE_K${K}_M${M}/", emit: rare_groupE_sharing_dir
    path "rare_groupN_K${K}_M${M}/", emit: rare_groupN_sharing_dir
    path "sharing_summary_E_A3/", emit: sharing_summary_E
    path "sharing_summary_N_A3/", emit: sharing_summary_N

    script:
    """
    echo "preparing yaml.."
    awk -F',' 'NR==1 || \$3=="E"' ${projectDir}/${annotation_file} | \
        python ${projectDir}/bin/rarified_sharing/1.make_groups_yml.py \
            --group_col genotype_short \
            --sample_col sample_short \
            --output groups_E.yml \
            --rename "heteroDQ2DQ8=hete,homoDQ8=DQ8,homoDQ2=DQ2"

    awk -F',' 'NR==1 || \$3=="N"' ${projectDir}/${annotation_file} | \
        python ${projectDir}/bin/rarified_sharing/1.make_groups_yml.py \
            --group_col genotype_short \
            --sample_col sample_short \
            --output groups_N.yml \
            --rename "heteroDQ2DQ8=hete,homoDQ8=DQ8,homoDQ2=DQ2"

    echo "publicity E..."
    python ${projectDir}/bin/rarified_sharing/2.rare_publicity.py \
        --input ${merged_file} \
        --K ${K} --M ${M} \
        --groups_file groups_E.yml \
        --output rare_groupE_K${K}_M${M} \
        --seed 42         

    echo "publicity N..."
    python ${projectDir}/bin/rarified_sharing/2.rare_publicity.py \
        --input ${merged_file} \
        --K ${K} --M ${M} \
        --groups_file groups_N.yml \
        --output rare_groupN_K${K}_M${M} \
        --seed 42           

    echo "summarising sharing..."
    python ${projectDir}/bin/rarified_sharing/3.summarise_sharing.py \
        --input rare_groupE_K${K}_M${M}/DQ2_sharing.tsv \
                rare_groupE_K${K}_M${M}/DQ8_sharing.tsv \
                rare_groupE_K${K}_M${M}/hete_sharing.tsv \
        --A 3 \
        --output sharing_summary_E_A3

    python ${projectDir}/bin/rarified_sharing/3.summarise_sharing.py \
        --input rare_groupN_K${K}_M${M}/DQ2_sharing.tsv \
                rare_groupN_K${K}_M${M}/DQ8_sharing.tsv \
                rare_groupN_K${K}_M${M}/hete_sharing.tsv \
        --A 3 \
        --output sharing_summary_N_A3


    """
}

   

// Output files:
// - `results/sharing_summary_A3/DQ2_summary.tsv`   