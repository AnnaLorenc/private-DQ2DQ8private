process AA_PROPERTIES_TEST{

  tag "translate aminoacids to AA properties, derive summary stats per site and individual, test for differences between groups"
    publishDir "${params.outdir}/aa_properties/${propname}", mode: 'copy'

    input:
    path(file_with_counts)
    val(index_columns)
    val(property_file)
    val(annotation_file)
    val(propname)
    val(agg_columns)


    output:
    tuple path("*_anno.tsv"), path("*_properties_wide_out.tsv") ,path("*_properties_results/*"), emit: prop_tests

    script:
    def index_cols_arg = index_columns ? "--index_columns ${index_columns.join(' ')}" : ""
    def agg_columns_arg = agg_columns ? "--agg_columns ${agg_columns.join(' ')}" : ""
    def exclude_cols_arg = index_columns ? "--exclude_cols ${agg_columns.join(',')}" : ""
    def baseName = file_with_counts.baseName
    """
        python ${projectDir}/bin/imgtCDR3/imgt_wide_to_properties.py \
        --input ${file_with_counts} \
        --property_file  ${projectDir}/${property_file} \
        --output ${baseName}_properties_wide_out.tsv \
        ${index_cols_arg}  \
        --aa_col "AA" 

        python ${projectDir}/bin/imgtCDR3/properties_add_annotation.py \
        --input ${baseName}_properties_wide_out.tsv \
        --annotation  ${projectDir}/${annotation_file} \
        --output ${baseName}_${propname}_anno.tsv \
    #    ${exclude_cols_arg}

        python ${projectDir}/bin/imgtCDR3/aa_properties_test.py \
        --input ${baseName}_${propname}_anno.tsv \
        ${agg_columns_arg} \
        --propname ${propname} 
    """
}

