library(tidyverse)

input <- read_tsv("sandbox/test_inputs/imgt_freqs/F5302E_subsampled_test.tsv" )

# changing in R to imgt positions

cdr3_to_imgt <- function(cdr3) {
  n <- nchar(cdr3)
  aa <- strsplit(cdr3, "", fixed = TRUE)[[1]]
  # IMGT: 104 (C), 105-111 (core), 111A/B... (insertions), 112-117 (core), 118 (F, if present)
  canonical_core <- 13
  pos <- c("104")
  core_len <- n - 2
  left <- as.character(105:111)
  right <- as.character(112:117)
  if (core_len <= canonical_core) {
    # No insertions
    n_right <- core_len - 7
    if (n_right > 0) {
      pos <- c(pos, left, right[1:n_right])
    } else {
      pos <- c(pos, left[1:core_len])
    }
    pos <- c(pos, "117")
  } else {
    # Insertions after 111: 111A, 111B, ...
    n_insert <- core_len - canonical_core
    insertions <- paste0("111", LETTERS[1:n_insert])
    pos <- c(pos, left, insertions, right)
    pos <- pos[1:(n-1)]
    pos <- c(pos, ifelse(n >= 15, "118", "117"))
  }
  tibble::tibble(aa = aa, pos = pos)
}

cdr3_vec <- c("CAAAIGGAQDTQYF", "CASSLGTDTQYF", "CASSLGGGSYNEQFF", "CASSLGGGGSYNEQFF")

input_imgt <- input %>%
  mutate(rowid = row_number()) %>%
  left_join(
    map_dfr(input$aminoAcid, ~cdr3_to_imgt(.x) %>% mutate(cdr3 = .x), .id = "rowid") %>%
      mutate(rowid = as.integer(rowid)),
    by = c( "rowid", "aminoAcid"="cdr3")
  )%>%
  mutate(aa_len=nchar(aminoAcid))

input_imgt %>%
  pivot_longer(cols=starts_with("F5302"), names_to = "sample", values_to = "value")%>%
  group_by( aa,    pos,   sample )%>%
  summarise(counts_uniq_TRBs =n_distinct(aminoAcid, vFamilyName) )%>%
  ungroup()%>%
  pivot_wider(names_from=sample, values_from=counts_uniq_TRBs)%>%
  group_by( aa, pos)%>%
  summarise(across(starts_with("F5302E"), ~.x/sum(.x)))

 to 
####=== I played with these to know whether it was doing the correct thing
raw <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/subsampled/productive/F5302E_subsampled.tsv.gz")
F5302E_freq_med.tsv <-read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq_med.tsv")
F5302E_freq_WL_med.tsv <-read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq_WL_med.tsv")
F5302E_freq.tsv<- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq.tsv")
F5302E_valid <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_valid.tsv")
F5302E_freq_Vfam<- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq_Vfam.tsv")



test <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/rarified_sharing/sharing_summary_E_A3/DQ8_summary.tsv")

test%>%mutate(across(`3`:`15`, ~if_else(.x>=100,.x,0)))%>%filter(`6`>0)

### old plotting functions with new input
pval_cutoff <-0.05
all_together  <- read_csv( file="~/sangermac/sr_projects/DR3DR4/analysis_rusted/processed_data/aa_pos_freq/all_together.csv.gz")
kidera_properties_for_plotting <- sapply(names(AAdata$kideraFactors),function(i){
  a=rep(c("blue","green","orange","red"), each=5)
  names(a) <-  AAdata$kideraFactors[[i]]%>%sort()%>%names()
  return(a)
}, USE.NAMES=TRUE, simplify=FALSE)
plot_enrichment_depletion2(all_together%>%filter(cells=="AgXP",pval_DQ2_vs_DQ8<pval_cutoff)%>%
                             mutate(imgt=as.character(imgt), length=paste0("l", length)),
                           length_col="length", estim_col="estim_DQ2_vs_DQ8", estim_expr_for_plotting="estim_DQ2_vs_DQ8>0",
                           inset_colour="red", aminoacid_mapping=kidera_properties_for_plotting[["KF1"]])+
  ggtitle("AgXp cells, enriched in DQ2", subtitle="")

change_imgt_format <- function(f){
  changef <- function(x){
    prefix <- substr(x, 1, nchar(x)-1)
    suffix <- substr(x, nchar(x), nchar(x))%>%
      match(., LETTERS)
    paste(prefix, suffix, sep=".")
  }
  g <-ifelse(grepl(pat="[ABCDEF]$",f),
             changef(f),
           f)
  g
  }

nextflow_output <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa_test/comb_aa_imgt_full_subs_rows_freq_5_1_res/imgt_lm_combined.tsv")
#Kidera7,8 - hightat the start, other and at the end
nextflow_output%>%filter(p_value_oo<0.001, !(IMGT_position%in%as.character(c(105:109, 115))))%>%
  mutate(imgt=change_imgt_format(IMGT_position),
                                                   aa=AA,
                                                   length=paste0("l", length))%>%
  plot_enrichment_depletion2(., 
                                                
                                                   length_col="length",
                                                   estim_col="cohens_d_oo",
                                                   estim_expr_for_plotting="cohens_d_oo>0",
                                                   inset_colour="red", aminoacid_mapping=kidera_properties_for_plotting[["KF1"]])



####paper figures
colors=c(DQ2="red",
         DQ2DQ8="purple",
         DQ8="blue")


####fig 1

diversity_next <- read_csv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_metrics/combined_productive_diversity_metrics.csv")

next_for_fig1 <- inner_join(anno_df%>%
                              select(file, newname, patient, cells,source, geno,anonym_sample_id),  
                       DR3DR4_samples_info%>%
                         select(`Specimen Name`,total_templates,total_rearrangements), by=c("file"="Specimen Name"))%>%
  select(-c(file))%>%
  inner_join(., 
             diversity_next%>%select(sample_name, unique_seqs=unique_sequences, total_count, clonality, gini_coefficient),
             by=c("newname"="sample_name"))%>%
  pivot_longer(cols = c(unique_seqs:gini_coefficient),
               names_to = "measure")%>%
  mutate(geno=gsub("^.* ",rep="",geno))

p<-list()
for(j in c("NAIVE","MEMORY")){
  for(i in c("unique_seqs","total_count","clonality","gini_coefficient")){
      p[[paste(i,j, sep="_")]] <- plot_paper_plot (i=i, cells=j, next_for_fig1)+ggtitle(paste(j, "cells, ",i))
  }
}


############ PCA
freqs <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_freqs/productive_subs_freqs.tsv.tsv")

pca_as_previous_one <- freqs%>%filter(param%in%c("jGeneName","vFamilyName"))%>%
  mutate(group=paste(param, group, sep="_"))%>%
  select(group, ends_with("subs"))
  
#filtering: eliinating levels where most of entries are 0
combinations_to_exclude <-
pca_as_previous_one$group[which(pca_as_previous_one%>%rowwise()%>%summarise(across(ends_with("subs"), ~sum(is.na(.x)|.x<0.001))) %>%rowSums()>10)]

df <- pca_as_previous_one%>%
  rename_with(., ~gsub(pat="_subs", rep="",.x))%>%
  filter(!group%in%combinations_to_exclude)%>%
  data.frame(row.names = "group")%>%
  
  as.matrix()%>%
  t()

library(ggplot2)

pca_res   <- prcomp(df, scale. = TRUE)
scores    <- as.data.frame(pca_res$x)

rownames(scores) <- gsub(pat="F5302", rep="F05302",rownames(scores) )%>%gsub(pat="F6826", rep="F06826" )
df_meta <- inner_join(data.frame(sample=rownames(scores)), anno_df%>%select(cells,newname,geno), by=c("sample"="newname"))

scores$Group <- df_meta$geno

ggplot(scores, aes(PC1, PC2, color = Group)) +
  geom_point(alpha = 0.7) +
  stat_ellipse(type = "t", linetype = 2) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
  )

###should be fine

