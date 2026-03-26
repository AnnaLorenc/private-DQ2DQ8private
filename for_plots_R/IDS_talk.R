
library(tidyverse)
library(kableExtra)
library(scales)
library(ggpubr)
library(rstatix)


###missing plot - raw numbers



source("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/plot_cdr3_functions.R") 


annotation_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/data/collated_info.csv"

diversities_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/diversity_metrics_rarefied/cleaner_diversities_summary.txt"
freqs_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_freqs/productive_subs_freqs.tsv.tsv"

diversities_raw_loc <- '/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_metrics/combined_productive_diversity_metrics.csv'
  
  plot_dir <- "/Users/ania/Documents/DQ2DQ8/pipeline/260304"

# colors=c(DQ2="red",
#          DQ2DQ8="purple",
#          DQ8="blue")

colors =c(DQ2="#FF4D4D",
          DQ2DQ8="#B455E0",
          DQ8="#2E86DE")

colors_two <- c("#E69F00", "#009E73")


anno <- read_csv(annotation_loc)%>%
  mutate(source=substr(newname,1,1))%>%
  rename(patient=shortname)%>%
  unique

#adding anonymous IDs
anonym_patient_ids <- 
  anno %>%
  select(patient, source)%>%
  unique()%>%
  group_by(source)%>%
  mutate(anonym_patient_id=1:n()%>%str_pad(.,side = "left",width=2,"0")%>%paste0(substr(source,1,1),.))


anno <- anno %>%
  inner_join(., anonym_patient_ids%>%ungroup%>%select(patient, anonym_patient_id), by=c("patient"))%>%
  mutate(anonym_sample_id=paste0(anonym_patient_id, cells))%>%
  mutate(stage=ifelse(source=="F","stage 3","stage 2"),
         cells_long=ifelse(cells=="N","naive","memory"),
         geno=gsub(pat=".*o", rep="", genotype_short))

imgt_tests_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa_test"
imgt_tests <- c("WL",  "full",  "WL_vfam",     "full_vfam")
imgt_prefixes <-c("comb_aa_imgt_", "_subs_rows_freq_5_1_res")
file <- "imgt_lm_combined.tsv"

properties_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/aa_properties"
properties <- c("KF", "VHS")
properties_tests <- c("WL",  "full",  "WL_vfam",     "full_vfam")
properties_res <- c("hotelling", "per_factor")
properties_prefixes <-c("comb_aa_imgt_", "_subs_rows_","_anno_")
#"comb_aa_imgt_" (WL) "_subs_rows_", (VHS) "_anno_" (VHS) (hotelling).tsv


##############################
############### RAW diversities
##############################


diversity_dir <- file.path(plot_dir, "diversity_raw" )
dir.create(diversity_dir, recursive = TRUE)

diversities <- read_csv(diversities_raw_loc)

diversities_anno <- diversities%>%
  select(-c(repertoire_id))%>%
  left_join(., anno%>%select(c(cells,
                               patient,
                               genotype_short,
                               sample_short,
                               source,
                               anonym_patient_id,
                               anonym_sample_id)), by=c("sample_name"="sample_short"))%>%
  pivot_longer(cols =total_sequences:shannon_diversity, names_to = "measure" )%>%
  unique()



for( i in c("unique_sequences", "total_count","gini_coefficient")){
  for(j in c("MEMORY","NAIVE")){
    plot_paper_plot (i=i, cells=j, diversities_anno)+ggtitle(paste(i, "in ",j, "cells raw data"))
    
    ggsave(file.path(diversity_dir, paste0("otherDiv_",substr(i,1,4),"_",substr(j,1,1) ,".png")), width=5, height = 7.5)
    
  }
}


##############################
############### LENGTHS
##############################

##### #### 
#### by stages
length_db_raref <- read_tsv(freqs_loc )%>%
  filter(param=="cdr3_length")%>%
  mutate(group=as.numeric(group))%>%
  filter(group>7, group<24)%>%
  select(aa_len=group, ends_with("subs")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))%>%
  pivot_longer(cols=ends_with("subs"),names_to = "sample_short")%>%
  mutate(sample_short=gsub(pat="_subs",rep="",sample_short))%>%
  left_join(anno%>%select(c(cells,
                            patient,
                            genotype_short,
                            sample_short,source,
                            anonym_patient_id,
                            anonym_sample_id,
                            cells_long,
                            geno,
                            stage))%>%unique(), by=c("sample_short"))



length_db_raref %>%
  arrange(geno, stage,sample_short, aa_len) %>%   # ensure proper ordering
  filter(aa_len>=11,aa_len<19)%>%
  group_by(stage, sample_short) %>%
  mutate(value=value/sum(value),
         cells=ifelse(cells=='N','NAIVE','MEMORY'))%>%
  ungroup() %>%
  ggplot(aes(x = aa_len, y = value, col = stage, fill=stage, group = interaction(stage, aa_len))) +
  geom_boxplot(alpha=0.5, outliers=FALSE) +
  theme_classic() +
  xlab("CDR3 length (aa)") +
  ggtitle("CDR3 length (mean frequencies)  across genotypes and stages") +
  facet_wrap(~cells )+
  # Significance stars: compares stage within each aa_len
  # stat_compare_means(
  #   aes(group = stage),
  #   method = "wilcox.test",   # or t.test if appropriate
  #   label = "p.signif",
  #   hide.ns = TRUE
  # )+
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)

ggsave(file.path(length_dir, "cons_length_cumsum_bystage.png"), width = 7.5, height = 6)




#### by genotypes, 
length_db_raref %>%
  arrange(geno, stage,sample_short, aa_len) %>%   # ensure proper ordering
  filter(aa_len>=11,aa_len<19)%>%
  group_by(geno, sample_short) %>%
  mutate(value=value/sum(value),
         cells=ifelse(cells=='N','NAIVE','MEMORY'))%>%
  ungroup() %>%
  ggplot(aes(x = aa_len, y = value, col = geno, fill=geno, group = interaction(geno, aa_len))) +
  geom_boxplot(alpha=0.5, outliers=FALSE) +
  theme_classic() +
  xlab("CDR3 length (aa)") +
  ggtitle("CDR3 length (mean frequencies) across genotypes and stages") +
  facet_wrap(~cells )+
  # Significance stars: compares stage within each aa_len
  # stat_compare_means(
  #   aes(group = geno),
  #   method = "wilcox.test",   # or t.test if appropriate
  #   label = "p.signif",
  #   hide.ns = TRUE
  # )+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)

ggsave(file.path(length_dir, "cons_length_cumsum_bygeno.png"), width = 7.5, height = 6)

#### within genotypes, by stages
length_db_raref %>%
  arrange(geno, stage,sample_short, aa_len) %>%   # ensure proper ordering
  filter(aa_len>=11,aa_len<19)%>%
  group_by(geno, stage,sample_short) %>%
  mutate(value=value/sum(value),
         cells=ifelse(cells=='N','NAIVE','MEMORY'))%>%
  ungroup() %>%
  ggplot(aes(x = aa_len, y = value, col = stage,fill=stage, group = interaction(stage, aa_len))) +
  geom_boxplot(alpha=0.5, outliers=FALSE) +
  theme_classic() +
  xlab("CDR3 length (aa)") +
  ggtitle("CDR3 length (mean frequencies) across genotypes and stages") +
  facet_wrap(cells ~ geno)+
  # Significance stars: compares stage within each aa_len
  # stat_compare_means(
  #   aes(group = stage),
  #   method = "wilcox.test",   # or t.test if appropriate
  #   label = "p.signif",
  #   hide.ns = TRUE
  # )+
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)

ggsave(file.path(length_dir, "cons_length_cumsum_bystage_within_geno.png"), width = 7.5, height = 9)



#### within  stages,by genotypes,
length_db_raref %>%
  arrange(geno, stage,sample_short, aa_len) %>%   # ensure proper ordering
  filter(aa_len>=11, aa_len<19)%>%
  group_by(geno, stage, sample_short) %>%
  mutate(value=value/sum(value),
         cells=ifelse(cells=='N','NAIVE','MEMORY'))%>%
  ungroup() %>%
  ggplot(aes(x = aa_len, y = value, col = geno,fill=geno, group = interaction(geno, aa_len))) +
  geom_boxplot(alpha=0.5, outliers=FALSE) +
  theme_classic() +
  xlab("CDR3 length (aa)") +
  ggtitle("CDR3 length (mean frequencies) across genotypes and stages") +
  facet_wrap(cells ~ stage)+
  # Significance stars: compares stage within each aa_len
  # stat_compare_means(
  #   aes(group = geno),
  #   method = "wilcox.test",   # or t.test if appropriate
  #   label = "p.signif",
  #   hide.ns = TRUE
  # )+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)

ggsave(file.path(length_dir, "cons_length_cumsum_bygeno_withinstage.png"), width = 7.5, height = 9)



####################
###### PCA
####################

pca_dir <- file.path(plot_dir, "pca_raref" )


V_J <- read_tsv(freqs_loc )%>%
  filter(param%in%c("vFamilyName","jGeneName"))%>%
  select(vj=group, ends_with("subs")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))


###### PCA naive (V and J separatedly)
######

pca_mat <- V_J%>%select(c(vj, ends_with("N_subs")))%>%
  filter(!grepl(pat="TCRBJ", vj))%>%
  rename_with(., ~gsub(pat="_subs", rep="",.x))%>%
  data.frame(row.names = "vj")%>%
  as.matrix()%>%
  t()

pca_res   <- prcomp(pca_mat, scale. = TRUE)
scores    <- as.data.frame(pca_res$x)

df_meta <- inner_join(data.frame(sample_short=rownames(scores)),
                      anno%>%unique()%>%select(cells,sample_short, genotype_short), by=c("sample_short"))%>%
  mutate(geno=case_when(genotype_short =="homoDQ2" ~ "DQ2",
                        genotype_short =="heteroDQ2DQ8" ~"DQ2DQ8",
                        genotype_short =="homoDQ8" ~ "DQ8"))
scores$Group <- df_meta$geno

ggplot(scores, aes(PC1, PC2, color = Group, shape=Group)) +
  geom_point(size=3) +
  stat_ellipse(type = "t", linetype = 1) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
  ) +coord_fixed() +scale_color_manual(values=colors)+
  scale_shape_manual(values=c(19,8,19))+
  ggtitle("Naive cells, V usage")

ggsave(file.path(pca_dir, "pca_naive_V.png"), width=7.5, height = 7)


###### PCA memory (V and J separatedly)
############

pca_mat <- V_J%>%select(c(vj, ends_with("E_subs")))%>%
  filter(!grepl(pat="TCRBJ", vj))%>%
  rename_with(., ~gsub(pat="_subs", rep="",.x))%>%
  data.frame(row.names = "vj")%>%
  as.matrix()%>%
  t()

pca_res   <- prcomp(pca_mat, scale. = TRUE)
scores    <- as.data.frame(pca_res$x)

df_meta <- inner_join(data.frame(sample_short=rownames(scores)),
                      anno%>%unique()%>%select(cells,sample_short, genotype_short), by=c("sample_short"))%>%
  mutate(geno=case_when(genotype_short =="homoDQ2" ~ "DQ2",
                        genotype_short =="heteroDQ2DQ8" ~"DQ2DQ8",
                        genotype_short =="homoDQ8" ~ "DQ8"))
scores$Group <- df_meta$geno

ggplot(scores, aes(PC1, PC2, color = Group, shape=Group)) +
  geom_point(size=3) +
  stat_ellipse(type = "t", linetype = 1) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
  ) +coord_fixed() +scale_color_manual(values=colors)+
  scale_shape_manual(values=c(19,8,19))+
  ggtitle("Memory cells, V usage")

ggsave(file.path(pca_dir, "pca_memory_V.png"), width=7.5, height = 7)




##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### Better Seq logo for heterozygotes
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

results_imgt_WL  <- read_tsv(file.path(imgt_tests_loc, file.path(paste0(imgt_prefixes[1], imgt_tests[1],imgt_prefixes[2]),file))) 
cols <-list(imgt = "imgt",
            aa = "aa",
            mean_a = "mean_group_a_eo",
            mean_b = "mean_group_b_eo", 
            cohens_d = "cohens_d_eo",
            p_value = "p_value_eo",
            p_value_do="p_value_do",cohens_d_do="cohens_d_do",
            label_a = "label_group_a_eo",
            label_b = "label_group_b_eo")

test_version <- "imgt_het_hom"
min_freq = 0.01
p_cutoff = 0.05
d_cutoff = 0.8
q_cutoff = 0.1

cell_fullnames <- c(N = "naive", E = "memory")

for(cells_sel in c("N","E")){
  
  p <- plot_aa_composition_logo( results_imgt_WL%>%filter(cells==cells_sel)%>%
                                  
                                   mutate(imgt =change_imgt_format(IMGT_position))%>%
                                   filter(imgt%in%c(as.character(105:111),'111.1',"111.2","112.1",as.character(112:118)))%>%
                                   rename("aa"=AA)%>%select(imgt, aa,ends_with("eo"),ends_with("do")),
                                 min_freq =min_freq, use_seqlogo = TRUE,
                                 cols=cols, p_cutoff = p_cutoff, d_cutoff = d_cutoff,
                                 add_to_title = paste(cell_fullnames[cells_sel],"cells"))
  ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL", cells_sel, sep="_"), ".png"))
  
  
  phet_log = plot_aa_composition_logo( results_imgt_WL%>%filter(cells==cells_sel)%>%
                                         
                                         mutate(imgt =change_imgt_format(IMGT_position))%>%
                                         filter(imgt%in%c(as.character(105:111),'111.1',"111.2","112.1",as.character(112:118)))%>%
                                         rename("aa"=AA)%>%select(imgt, aa,ends_with("eo"),ends_with("do")),
                                       min_freq =min_freq, use_seqlogo = TRUE,
                                       cols=cols, p_cutoff = p_cutoff, d_cutoff = d_cutoff,
                                       add_to_title = paste(cell_fullnames[cells_sel],"cells"),scale_log = TRUE)
  ggsave(plot= phet_log, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL_log", cells_sel, sep="_"), ".png"))
  
  
  phet = plot_aa_composition_logo( results_imgt_WL%>%filter(cells==cells_sel)%>%
                                         
                                         mutate(imgt =change_imgt_format(IMGT_position))%>%
                                         filter(imgt%in%c(as.character(105:111),'111.1',"111.2","112.1",as.character(112:118)))%>%
                                         rename("aa"=AA)%>%select(imgt, aa,ends_with("eo"),ends_with("do")),
                                       min_freq =min_freq, use_seqlogo = TRUE,
                                       cols=cols, p_cutoff = p_cutoff, d_cutoff = d_cutoff,
                                       add_to_title = paste(cell_fullnames[cells_sel],"cells"),scale_log = FALSE)
  
  ggsave(plot= phet, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL_nolog", cells_sel, sep="_"), ".png"))
  
} 


##now plogo like

imgt_range <- c(as.character(106:111),'111.1',"111.2","112.1",as.character(112:117))

bonf <- results_imgt_WL%>%filter(cells==cells_sel)%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%rowwise()%>%
  filter(imgt%in%c(imgt_range), max(mean_group_a_oo, mean_group_b_oo, mean_group_a_eo)>0.01)%>%
  rename("aa"=AA)%>%select(imgt, aa,ends_with("oo"))%>%group_by(imgt)%>%summarise(a=n_distinct(aa))%>%summarise(sum(a))%>%unlist()

plot_plogo_like(df= results_imgt_WL%>%filter(cells==cells_sel)%>%
                   mutate(imgt =change_imgt_format(IMGT_position))%>%
                   filter(imgt%in%imgt_range)%>%
                   rename("aa"=AA)%>%select(imgt, aa,ends_with("eo"),ends_with("do"))%>%
                  mutate(label_group_a_eo=gsub(pat="he", rep="",label_group_a_eo), label_group_b_eo=gsub(pat="hom", rep="DQ2, DQ8",label_group_b_eo) ),
                height_by="cohens_d",
                 cols=cols, p_cutoff = p_cutoff/(bonf), d_cutoff = d_cutoff,
                 add_to_title = paste(cell_fullnames[cells_sel],"cells, Bonferroni:", bonf))




test_version <- "prop_het_hom"

cols = list(
  imgt     = "imgt",
  prop       = "prop",
  cells = "cells",
  mean_a   = "mean_a_eo",
  mean_b   = "mean_b_eo",
  cohens_d = "cohens_d_eo",
  p_value  = "p_value_eo",
  q_value = "q_value_eo",
  label_a  = "label_group_a_eo",
  label_b  = "label_group_b_eo"
)



underlying_data_KF_WL <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/aa_properties/KF/comb_aa_imgt_WL_subs_rows_KF_anno.tsv")

main_pval_cutoff <- 0.05

prepared_KF_WL <- results_KF_WL_by_fact%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%select(cells, imgt, prop,ends_with("eo"),ends_with("do"),ends_with("po"))%>%select(-c(contains('hotelling')))%>%
  mutate(label_group_a_eo="DQ2DQ8",label_group_b_eo="homo",label_group_b_po="DQ2",label_group_b_do="DQ8")%>%
  filter(imgt%in%c(as.character(105:111),'111.1',"111.2","112.1",as.character(112:118)) )%>%
  select(-c(starts_with("mean_diff"), starts_with("t_stat"), starts_with("df")))%>%
  arrange(cells, imgt, prop)

signif_KF_WL <- prepared_KF_WL%>%
  filter(imgt%in%positions,
         p_value_do<main_pval_cutoff,
         p_value_po<main_pval_cutoff,
         p_value_eo<main_pval_cutoff )


#here I am actually plottingfiltered by pvals, not qvals!
positions = c("109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113")
q_cutoff = main_pval_cutoff
d_cutoff = 0.5

cols <- list(
  imgt     = "imgt",
  prop       = "prop",
  cells = "cells",
  mean_a   = "mean_b_po",
  mean_b   = "mean_b_do",
  mean_c = "mean_a_eo",
  cohens_d = "cohens_d_eo",
  q_value = "p_value_eo",       
  cohens_d = "cohens_d_eo",
  cohens_d_do = "cohens_d_do",
  cohens_d_po = "cohens_d_po",
  q_value_do = "p_value_do",     
  q_value_po = "p_value_po",  
  label_a = "label_group_b_po",
  label_b = "label_group_b_do",
  label_c = "label_group_a_eo"        # Third group label
)



for(cells_sel in c("N","E")){
  
  p <-  plot_property_logo(
    prepared_KF_WL,
    positions = positions,
    q_cutoff   = p_cutoff,
    d_cutoff   = d_cutoff,
    cols = cols,
    cell_types = cells_sel,
    color_a    = "#e41a1c",   
    color_b    = "#377eb8" ,  
    color_c = "purple",
    add_to_title = paste(cell_fullnames[cells_sel],"cells"))
  
  ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL", cells_sel, sep="_"), ".png"))
  


###now boxplots



p <- underlying_data_KF_WL %>%
  pivot_longer(cols = paste0('KF',1:10),names_to = 'prop')%>%
  mutate(imgt=change_imgt_format(IMGT_position),
         genotype=gsub(pat="h.D", rep="D",genotype))%>%
  inner_join(., signif_KF_WL, by=c('cells','imgt', "prop"))%>%
  filter(cells==cells_sel)%>%
  ggplot(aes(x=genotype, y=value, group=genotype, col=genotype) )+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width=0.3, height =0)+
  facet_wrap(imgt~prop, scales='free_y') +theme_bw() +scale_color_manual(values=colors)+guides(color="none")+ 
  stat_compare_means(
    aes(group = genotype),
    method = "wilcox.test",  comparisons = list(
      c("DQ2", "DQ2DQ8"),
      c("DQ2DQ8", "DQ8"),
      c("DQ2", "DQ8")
    ) ,# or t.test if appropriate
    # label = "p.signif",
    # hide.ns = TRUE
  )
p <- p+ ggtitle(paste("KF properties,", cell_fullnames[cells_sel], "cells, all lengths"))



ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"_box_WL", cells_sel, sep="_"), ".png"))
}
###the same BUT as single plots





  for(len_sel in 13:16){
    p <-  
      plot_property_logo(
        results_KF_by_fact%>%
          filter(length==len_sel)%>%
          mutate(imgt =change_imgt_format(IMGT_position))%>%select(cells, imgt, prop,ends_with("eo"),ends_with("do"),ends_with("po"))%>%select(-c(contains('hotelling')))%>%
          mutate(label_group_a_eo="DQ2DQ8",label_group_b_eo="homo",label_group_b_po="DQ2",label_group_b_do="DQ8"),
        positions = c("109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"),
        
        q_cutoff   = q_cutoff,
        d_cutoff   = d_cutoff,
        cols = cols,
        cell_types = cells_sel,
        color_b    = "#e41a1c",  
        color_a    = "#377eb8"  ,  
        color_c = "purple",
        add_to_title =  paste(cell_fullnames[cells_sel],"cells, length=", len_sel))
    ggsave(plot = p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version, paste0("l",len_sel), cells_sel, sep="_"),".png"))
  }
}




