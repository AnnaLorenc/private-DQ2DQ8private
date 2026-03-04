Efforts to move from by aminoacid testing to testing by properties (Kidera, VHS, crucianii.... - available in asstes)

1. Getting properties from the frequencies, version_via_freqs/aa_kidera_aggregate.py. output:  version_via_freqs/kid_trans.tsv
Then testing with aa_properties_test.py . output:  version_via_freqs/kidera_results

2. Then I realised that the filtering, done when preparing the final frequency files, leads to errors (single AA might be rare,so not testable, but their property could be used)
Solution: to prepare the input file for testing from imgtCDR3/imgt_wide_to_properties.py, from results/imgt_aa_combined/comb_aa_imgt_subs_rows.tsv.
Then the annotation is added with imgtCDR3/properties_add_annotation.py.
Testing again with aa_properties_test.py

