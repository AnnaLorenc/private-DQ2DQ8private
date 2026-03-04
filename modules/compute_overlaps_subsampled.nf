process COMPUTE_OVERLAPS_SUBSAMPLED {
    tag "overlaps_subs"
    publishDir "${params.outdir}/overlaps_subs/", mode: 'copy'



    input:
    tuple path(sample_names_file), val(subsampled_dirs), val(type)

    output:
    path ("overlaps_subs_*.tsv"), emit: overlaps
   


    """
    python ${projectDir}/bin/overlaps/compute_overlaps_subs.py \
        --samples_file ${sample_names_file} \
        --subsample_dir ${subsampled_dirs} \
        --output_file overlaps_subs_${type}.tsv 
    """
}

//   parser = argparse.ArgumentParser()
//     parser.add_argument("--samples_file", required=True)
//     parser.add_argument("--subsample_dir", required=True)
//     parser.add_argument("--output_file", required=True)
//     parser.add_argument(
//         "--id_columns",
//         nargs="+",
//         default=["aminoAcid", "vFamilyName", "jGeneName"],
//         help="Columns to use for building clonotype IDs (default: aminoAcid, vFamilyName, jGeneName)",