#! bin/bash

for cells in E N
do
    echo "Processing ${cells} cells..."
    python bin/ad_length_distribution.py \
  --freqs_dir results/freqs_subsampled/productive \
  --sample_map data/collated_info_new.csv \
  --group_col genotype_short \
  --sample_col shortname \
  --cells ${cells} \
  --n_perm 10000 \
  --output results_man/ad_length_results_${cells}_geno.tsv

python bin/ad_length_distribution.py \
  --freqs_dir results/freqs_subsampled/productive \
  --sample_map data/collated_info_new.csv \
  --group_col genotype_short \
  --split_col source \
  --sample_col shortname \
  --cells ${cells} \
  --n_perm 10000 \
  --output results_man/ad_length_results_${cells}_geno_byTF.tsv



python bin/ad_length_distribution.py \
  --freqs_dir results/freqs_subsampled/productive \
  --sample_map data/collated_info_new.csv \
  --group_col source \
  --sample_col shortname \
  --split_col genotype_short \
  --groups T F \
  --cells ${cells} \
  --n_perm 10000 \
  --output results_man/ad_length_results_${cells}_TF_bygeno.tsv


python bin/ad_length_distribution.py \
  --freqs_dir results/freqs_subsampled/productive \
  --sample_map data/collated_info_new.csv \
  --group_col source \
  --sample_col shortname \
  --cells ${cells} \
  --n_perm 10000 \
  --output results_man/ad_length_results_${cells}_TF.tsv
done

awk 'FNR==1 && NR!=1 {next} 1' results_man/ad_length_results_*by*  > results_man/combined_ad_length_results_all_by.tsv

awk -F'\t' '
FNR==1 {
  keep = (NF==9)
  if (!keep) nextfile
  if (seen_header) next
  seen_header=1
}
keep
' results_man/ad_length_results_* > results_man/combined_ad_length_results_all_simple.tsv