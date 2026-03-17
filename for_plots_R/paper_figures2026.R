

library(tidyverse)
library(kableExtra)
library(scales)
library(ggpubr)
library(rstatix)

#knitr::opts_chunk$set( message=FALSE, fig.width = 8, fig.height = 12, warning = FALSE)


source("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/plot_cdr3.R") 


annotation_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/data/collated_info.csv"

diversities_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/diversity_metrics_rarefied/cleaner_diversities_summary.txt"
freqs_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_freqs/productive_subs_freqs.tsv.tsv"


plot_dir <- "/Users/ania/Documents/DQ2DQ8/pipeline/260304"

colors=c(DQ2="red",
         DQ2DQ8="purple",
         DQ8="blue")

colors_two <- c("darkgrey","darkgreen")

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


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Diversity
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

diversity_dir <- file.path(plot_dir, "diversity_raref" )
dir.create(diversity_dir, recursive = TRUE)

diversities <- read_tsv(diversities_loc)

diversities_anno <- diversities%>%
  select(-c(total_sequences, rarefaction_iterations, sequence_index))%>%
left_join(., anno%>%select(c(cells,
                           patient,
                           genotype_short,
                           sample_short,source,
                           anonym_patient_id,
                           anonym_sample_id)), by=c("sample_short"))%>%
  pivot_longer(cols =unique_sequences:shannon_index, names_to = "measure" )%>%
  unique()


plot_paper_plot <- function(i, cells, for_fig1) {
  
  ## ---- Required columns check ----
  required_cols <- c("measure", "genotype_short", "cells", "value", "source")
  missing_cols <- setdiff(required_cols, colnames(for_fig1))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf("Missing required column(s) in for_fig1: %s",
              paste(missing_cols, collapse = ", "))
    )
  }
  
  ## ---- Check content of cells column ----
  allowed_raw_cells <- c("E", "N")
  observed_cells <- unique(for_fig1$cells)
  
  if (!all(observed_cells %in% allowed_raw_cells)) {
    stop(
      sprintf(
        "Unexpected values in column 'cells'. Found: %s. Allowed: %s",
        paste(observed_cells, collapse = ", "),
        paste(allowed_raw_cells, collapse = ", ")
      )
    )
  }
  
  ## ---- Validate requested cells argument ----
  allowed_final_cells <- c("MEMORY", "NAIVE")
  if (!cells %in% allowed_final_cells) {
    stop(
      sprintf(
        "Argument 'cells' must be one of: %s",
        paste(allowed_final_cells, collapse = ", ")
      )
    )
  }
  
  ## ---- Check external dependencies ----
  if (!exists("colors")) {
    stop("Object 'colors' not found in environment.")
  }
  
  if (!exists("new_param_names")) {
    stop("Object 'new_param_names' not found in environment.")
  }
  
  ## ---- Data transformation ----
  for_fig <- for_fig1 %>%
    dplyr::filter(measure == i) %>%
    dplyr::mutate(
      cells = ifelse(cells == "E", "MEMORY", "NAIVE"),
      genotype_short = gsub(".*o", "", genotype_short)
    )
  
  ## ---- Plot ----
  p <- for_fig %>%
    dplyr::filter(cells == cells) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = genotype_short,
      y = value,
      color = genotype_short,
      fill = genotype_short,
      shape = source,
      group = genotype_short
    )) +
    ggplot2::geom_boxplot(outlier.size = 0, notch = FALSE, alpha = 0.3) +
    ggplot2::geom_jitter(height = 0, width = 0.2, size = 4, alpha = 0.6) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::ggtitle(new_param_names[i]) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12)) +
    ggplot2::scale_shape_manual(
      values = c(1, 19),
      labels = c("stage 2", "stage 3"),
      name = "stage"
    ) +
    ggplot2::guides(color = "none", fill = "none")
  
  return(p)
}




plot_paper_plot (i="inverse_simpson", cells="NAIVE", diversities_anno)+ ggtitle("Inverse Simpson index, naive cells (rarefaction to 10K cells)")
ggsave(file.path(diversity_dir, "fig2A2.png"))

plot_paper_plot (i="shannon_index", cells="NAIVE", diversities_anno) + ggtitle("Shannon diversity, naive cells")
ggsave(file.path(diversity_dir, "fig2A1.png"))



plot_paper_plot (i="inverse_simpson", cells="MEMORY", diversities_anno) + ggtitle("Inversed Simpson index, memory cells (rarefaction to 10K cells)")
ggsave(file.path(diversity_dir, "fig3A2.png"))

plot_paper_plot (i="shannon_index", cells="MEMORY", diversities_anno) +ggtitle("Shannon diversity, memory cells")
ggsave(file.path(diversity_dir, "fig3A1.png"))


for( i in c("unique_sequences", "clonality","gini_coefficient")){
  for(j in c("MEMORY","NAIVE")){
    plot_paper_plot (i=i, cells=j, diversities_anno)+ggtitle(paste(i, "in ",j, "cells (rarefaction to 10K cells"))
    
    ggsave(file.path(diversity_dir, paste0("otherDiv_",substr(i,1,4),"_",substr(j,1,1) ,".png")))
    
  }
}

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### New fig 2B/3B, comparing lengths of CDR3s
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

length_dir <- file.path(plot_dir, "length_raref" )
dir.create(length_dir, recursive = TRUE)

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

##### #### Boxplots per stage
##### #### 

ggplot(
  length_db_raref,
  aes(
    x = factor(aa_len),
    y = value,
    colour = stage    # box fill by stage
  )
)  +
  # Overlay individual points
  geom_jitter(
    aes(color = stage),  # optional, same as fill
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    alpha = 0.5,
    size = 1
  ) +
  # Boxplots per aa_len × stage
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA
  )+
  # Significance stars: compares stage within each aa_len
  stat_compare_means(
    aes(group = stage),
    method = "wilcox.test",   # or t.test if appropriate
    label = "p.signif",
    hide.ns = TRUE
  ) +
  # Facet by cells (rows) × geno (columns)
  facet_grid(cells ~ geno) +
  theme_classic() +
  labs(
    x = "AA length",
    y = "Value",
    fill = "Stage"
  ) +
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)
ggsave(file.path(length_dir, "length_db.png"))

##### #### Boxplots per stage cumulative
##### #### 


length_db_raref %>%
arrange(geno, stage,sample_short, aa_len) %>%   # ensure proper ordering
  group_by(geno, stage,sample_short) %>%
  mutate(cum_value = cumsum(value)) %>%
  ungroup() %>%
  ggplot(aes(x = aa_len, y = cum_value, col = stage, group = interaction(stage, aa_len))) +
  geom_boxplot() +
  theme_classic2() +
  scale_shape_manual(values = c(1, 19)) +
  xlab("CDR3 length (aa)") +
  scale_x_discrete(
    breaks = paste(sep = ".", "homo DQ8", c(8, 12, 16, 20, 24)),
    labels = c(8, 12, 16, 20, 24)
  ) +
  ggtitle("CDR3 length (cumulative mean) across genotypes and stages\n naive cells") +
  facet_wrap(cells ~ geno)+
  # Significance stars: compares stage within each aa_len
  stat_compare_means(
    aes(group = stage),
    method = "wilcox.test",   # or t.test if appropriate
    label = "p.signif",
    hide.ns = TRUE
  )+
  scale_fill_manual(values=colors_two)+
  scale_color_manual(values=colors_two)

ggsave(file.path(length_dir, "length_cumsum.png"))

####testing
enough_TCRs <-  length_db_raref %>%group_by(cells, aa_len) %>%
  slice_min(value)%>%
  filter(value>=0.01)%>%
  select(cells, aa_len)

###I use lengths 11-18, as they comprise at least 1% of the sample
###lengths 11-18
length_db_raref %>%
  filter(aa_len>10, aa_len<19)%>%
  group_by(cells, aa_len)%>%
  filter(all(table(geno) >= 2)) %>% 
  t_test(value ~ geno,var.equal=FALSE,p.adjust.method="bonf",detailed=T)%>%
  ungroup()%>%
  rowwise()%>%
  mutate(p_fully_adjusted= min(p*8*3,1))%>% #now adjust: for 2 additional comparisons, for 8 lengths
  group_by(cells)%>%
  slice_min(order_by = p, n = 5)%>%select(-starts_with("p.adj"))
  
###ns

###stage within genotype
length_db_raref %>%
  filter(aa_len>10, aa_len<19)%>%
  group_by(cells, aa_len, geno)%>%
  filter(all(table(stage) >= 2)) %>% 
  t_test(value ~ stage,var.equal=FALSE,p.adjust.method="fdr",detailed=T)%>%
  ungroup()%>%
  rowwise()%>%
  mutate(p_fully_adjusted= min(p*8*3,1))%>% #now adjust: for 3 within genotype comparisons, for 8 lengths
  group_by(cells)%>%
  slice_min(order_by = p, n = 5)%>%select(-starts_with("p.adj"))


# aa_len cells geno   estimate estimate1 estimate2 .y.   group1  group2     n1    n2 statistic        p    df  conf.low conf.high method alternative p_fully_adjusted
# <dbl> <chr> <chr>     <dbl>     <dbl>     <dbl> <chr> <chr>   <chr>   <int> <int>     <dbl>    <dbl> <dbl>     <dbl>     <dbl> <chr>  <chr>                  <dbl>
#   1     16 E     DQ2     0.00682    0.143     0.136  value stage 2 stage 3    10    10      2.75 0.0165    13.1  0.00146   0.0122   T-test two.sided             0.396 
# 2     12 E     DQ2DQ8 -0.00580    0.0739    0.0797 value stage 2 stage 3    10    10     -2.44 0.0271    15.5 -0.0108   -0.000749 T-test two.sided             0.650 
# 3     16 E     DQ2DQ8  0.00682    0.139     0.132  value stage 2 stage 3    10    10      2.34 0.0316    17.2  0.000674  0.0130   T-test two.sided             0.758 
# 4     12 E     DQ2    -0.00409    0.0726    0.0767 value stage 2 stage 3    10    10     -1.98 0.0652    16.0 -0.00847   0.000289 T-test two.sided             1     
# 5     17 E     DQ2DQ8  0.00487    0.0664    0.0615 value stage 2 stage 3    10    10      1.81 0.0878    18.0 -0.000798  0.0105   T-test two.sided             1     
# 6     13 N     DQ2    -0.0119     0.164     0.176  value stage 2 stage 3    10    10     -4.38 0.00043   16.5 -0.0177   -0.00617  T-test two.sided             0.0103
# 7     16 N     DQ2     0.00791    0.146     0.138  value stage 2 stage 3    10    10      4.11 0.000726  17.0  0.00385   0.0120   T-test two.sided             0.0174
# 8     17 N     DQ2     0.00761    0.0744    0.0668 value stage 2 stage 3    10    10      3.45 0.00367   14.8  0.00290   0.0123   T-test two.sided             0.0881
# 9     12 N     DQ2    -0.00577    0.0706    0.0764 value stage 2 stage 3    10    10     -3.10 0.00738   14.8 -0.00974  -0.00180  T-test two.sided             0.177 
# 10     18 N     DQ2     0.00454    0.0348    0.0303 value stage 2 stage 3    10    10      2.99 0.0107    12.8  0.00125   0.00783  T-test two.sided             0.257 


###stage without genotype
length_db_raref %>%
  filter(aa_len>10, aa_len<19)%>%
  group_by(cells, aa_len)%>%
  filter(all(table(stage) >= 2)) %>% 
  t_test(value ~ stage,var.equal=FALSE,p.adjust.method="fdr",detailed=T)%>%arrange(p)%>%
  ungroup()%>%
  rowwise()%>%
  mutate(p_fully_adjusted= min(p*8,1))%>% #now adjust: for 2 additional comparisons, for 8 lengths
  group_by(cells)%>%
  slice_min(order_by = p, n = 5)%>%select(-starts_with("p.adj"))





##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### PCA
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

pca_dir <- file.path(plot_dir, "pca_raref" )
dir.create(pca_dir, recursive = TRUE)

V_J <- read_tsv(freqs_loc )%>%
  filter(param%in%c("vFamilyName","jGeneName"))%>%
  select(vj=group, ends_with("subs")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))


VJ <- read_tsv(freqs_loc )%>%
  filter(param=="vFamilyName_jGeneName")%>%
  select(vj=group, ends_with("subs")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

###### selected TRBV across genotypes
######

V_J %>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
  #  filter(value > 0) %>%
  mutate(sample_short=gsub(pat="_subs", rep="",sample_short))%>%
  left_join(anno%>%select(c(cells,
                            patient,
                            genotype_short,
                            sample_short,source,
                            anonym_patient_id,
                            anonym_sample_id,
                            geno,
                            cells_long,
                            stage))%>%unique(), by=c("sample_short")) %>%
  filter(((vj %in% c("TCRBV15", "TCRBV05"))&(cells == "N")) |
           ((vj %in% c("TCRBV18", "TCRBV19", "TCRBV20")) & cells == "E"))%>%
  mutate(
    genotype = factor(geno, levels = c("DQ2", "DQ2DQ8",  "DQ8")),
    cells = ifelse(cells == "E", "MEMORY", "NAIVE"),
  ) %>%
  ggplot(aes(
    group = genotype,
    x = genotype,
    y = value,
    col = genotype,
    shape = stage
  )) + geom_boxplot(notch = F, outliers = F) +
  geom_jitter(aes( shape = stage), height = 0,
              alpha = 0.5,
              size = 2) + facet_wrap( cells ~ vj, scales = "free_y") +
  theme_bw()  +
  stat_compare_means(method = "wilcox.test", comparisons = list(
    c("DQ2", "DQ2DQ8"),
    c("DQ2DQ8", "DQ8"),
    c("DQ2", "DQ8")
  )) +
  scale_color_manual(values = colors) +
  guides(color = "none", shape = "none") +
  scale_shape_manual(values = c(1, 19)) +
  ggtitle("Examples of TRBV usage difference")

ggsave(filename =file.path(pca_dir, "SuppFig4.png"))


#statistical testing
signif_vj_oo <- V_J %>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
  #  filter(value > 0) %>%
  mutate(sample_short=gsub(pat="_subs", rep="",sample_short))%>%
  left_join(anno%>%select(c(cells,
                            patient,
                            genotype_short,
                            sample_short,source,
                            anonym_patient_id,
                            anonym_sample_id,
                            geno,
                            cells_long,
                            stage))%>%unique(), by=c("sample_short")) %>%
  group_by(cells, vj)%>%
  t_test(value ~ geno,var.equal=FALSE,p.adjust.method="fdr",detailed=T)%>%filter(p.adj<0.05, group1%in%c("DQ2","DQ8"),  group2%in%c("DQ2","DQ8")) #%>%filter(p*20*2<0.05) 
# vj      cells estimate estimate1 estimate2 .y.   group1 group2    n1    n2 statistic        p    df conf.low conf.high method alternative    p.adj p.adj.signif
# <chr>   <chr>    <dbl>     <dbl>     <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl> <dbl>    <dbl>     <dbl> <chr>  <chr>          <dbl> <chr>       
#   1 TCRBV05 E      0.0215     0.141    0.120   value DQ2    DQ2DQ8    20    20      6.20 3.32e- 7  37.2  0.0145    0.0286  T-test two.sided   5.32e- 7 ****        
#   2 TCRBV05 E      0.0221     0.141    0.119   value DQ2    DQ8       20    21      6.13 3.55e- 7  38.7  0.0148    0.0294  T-test two.sided   5.32e- 7 ****        
#   3 TCRBV15 E      0.00544    0.0133   0.00789 value DQ2    DQ8       20    21      8.97 5.23e-11  38.9  0.00421   0.00666 T-test two.sided   1.57e-10 ****        
#   4 TCRBV15 E      0.00325    0.0111   0.00789 value DQ2DQ8 DQ8       20    21      4.75 2.91e- 5  37.9  0.00186   0.00463 T-test two.sided   4.36e- 5 ****        
#   5 TCRBV05 N      0.0226     0.138    0.115   value DQ2    DQ2DQ8    20    20      5.66 2.76e- 6  32.6  0.0144    0.0307  T-test two.sided   8.28e- 6 ****        
#   6 TCRBV05 N      0.0222     0.138    0.116   value DQ2    DQ8       20    21      5.28 6.35e- 6  36.0  0.0137    0.0308  T-test two.sided   9.53e- 6 ****        
#   7 TCRBV15 N      0.00264    0.0148   0.0121  value DQ2    DQ2DQ8    20    20      4.48 6.71e- 5  37.6  0.00145   0.00384 T-test two.sided   6.71e- 5 ****        
#   8 TCRBV15 N      0.00558    0.0148   0.00919 value DQ2    DQ8       20    21      9.09 3.63e-11  38.9  0.00434   0.00683 T-test two.sided   1.09e-10 ****        
#   9 TCRBV15 N      0.00294    0.0121   0.00919 value DQ2DQ8 DQ8       20    21      5.02 1.20e- 5  38.9  0.00175   0.00412 T-test two.sided   1.80e- 5 ****        
#   10 TCRBV18 N     -0.00556    0.0329   0.0384  value DQ2    DQ8       20    21     -3.73 6.89e- 4  34.6 -0.00859  -0.00253 T-test two.sided   2   e- 3 **          
#   11 TCRBV19 N     -0.00463    0.0456   0.0502  value DQ2    DQ2DQ8    20    20     -3.64 8.35e- 4  36.5 -0.00720  -0.00205 T-test two.sided   3   e- 3 **          
#   12 TCRBV20 N      0.0257     0.111    0.0855  value DQ2    DQ8       20    21      6.28 2.87e- 7  36.1  0.0174    0.0339  T-test two.sided   8.61e- 7 ****        
#   13 TCRBV30 N     -0.0108     0.0196   0.0304  value DQ2    DQ8       20    21     -3.67 7.21e- 4  38.9 -0.0168   -0.00486 T-test two.sided   2   e- 3 ** 




###### selected TRBV , naive versus memory
######

V_J %>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
  #  filter(value > 0) %>%
  mutate(sample_short=gsub(pat="_subs", rep="",sample_short))%>%
  left_join(anno%>%select(c(cells,
                            patient,
                            genotype_short,
                            sample_short,source,
                            anonym_patient_id,
                            anonym_sample_id,
                            geno,
                            cells_long,
                            stage))%>%unique(), by=c("sample_short")) %>%
  filter(vj %in% c("TCRBV16", "TCRBV02"))%>%
  mutate(
    genotype = factor(geno, levels = c("DQ2", "DQ2DQ8",  "DQ8")),
    cells = ifelse(cells == "E", "MEMORY", "NAIVE"),
  ) %>%
  ggplot(aes(
    group = cells,
    x = cells,
    y = value,
    col = genotype,
    shape = stage
  )) + geom_boxplot(notch = F, outliers = F) +
  geom_jitter(
    height = 0.00001,
    alpha = 0.5,
    size = 3,
    width = 0.2
  ) + facet_wrap( ~ vj, scales = "free_y") +
  theme_bw() + scale_color_manual(values = c("red", "purple", "blue")) +
  guides(color = "none") +
  scale_shape_manual(values = c(1, 19)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("MEMORY", "NAIVE")))


ggsave(filename = file.path(pca_dir, "SuppFig3B.png"))



###### PCA naive (V and J separatedly)
######

pca_mat <- V_J%>%select(c(vj, ends_with("N_subs")))%>%
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

ggplot(scores, aes(PC1, PC2, color = Group)) +
geom_point(alpha = 0.7, size=2) +
stat_ellipse(type = "t", linetype = 1) +
theme_minimal() +
labs(
x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
) +coord_fixed() +scale_color_manual(values=colors)+
  ggtitle("Naive cells, VJ usage")

ggsave(file.path(pca_dir, "fig2C.png"))


###### PCA memory (V and J separatedly)
############
pca_mat <- V_J%>%select(c(vj, ends_with("E_subs")))%>%
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

ggplot(scores, aes(PC1, PC2, color = Group)) +
  geom_point(alpha = 0.7, size=2) +
  stat_ellipse(type = "t", linetype = 1) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
  ) +coord_fixed() +scale_color_manual(values=colors)+
  ggtitle("Memory cells, VJ usage")

ggsave(file.path(pca_dir, "fig3C.png"))


###otehr parts - migrate to another file

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Pairwise sharing
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Public sequences
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Chocolate plots + chocolate bar plots (for length agnostic tests)
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

# Extract underlying values for plots
# OOF repertoire have similar V/J usage across genotypes, showing that recombination events are not linked to the DQ status TO BE DONE - OOF/ priority:2
# Public sequences:Sharing of TCRB sequences is highest between individuals from the same genotype and lowest for the individuals of two opposite homozygous phenotypes
# { pending results, as the public sequences were analysed for the memory cells only: public sequences sharing} TO BE DONE
# { Aminoacid usage: genotype-unique TRBs – where do they fall in this pattern } TO BE DONE


