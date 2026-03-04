

library(tidyverse)
library(kableExtra)
library(scales)
library(ggpubr)

#knitr::opts_chunk$set( message=FALSE, fig.width = 8, fig.height = 12, warning = FALSE)


source("/Users/ania/sangermac/sr_projects/DR3DR4/analysis_rusted/scripts/plot_cdr3.R") 


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
  rename(patient=shortname)

#adding anonymous IDs
anonym_patient_ids <- 
anno %>%
  select(patient, source)%>%
  unique()%>%
  group_by(source)%>%
  mutate(anonym_patient_id=1:n()%>%str_pad(.,side = "left",width=2,"0")%>%paste0(substr(source,1,1),.))


anno <- anno %>%
  inner_join(., anonym_patient_ids%>%ungroup%>%select(patient, anonym_patient_id), by=c("patient"))%>%
  mutate(anonym_sample_id=paste0(anonym_patient_id, cells))


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
                            anonym_sample_id))%>%unique(), by=c("sample_short"))
  
length_db_raref <- length_db_raref%>%
  mutate(stage=ifelse(source=="F","stage 3","stage 2"),
         cells=ifelse(cells=="N","naive","memory"),
         geno=gsub(pat=".*o", rep="", genotype_short))


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


V_J %>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
  #  filter(value > 0) %>%
  mutate(sample_short=gsub(pat="_subs", rep="",sample_short))%>%
  left_join(anno%>%select(c(cells,
                            patient,
                            genotype_short,
                            sample_short,source,
                            anonym_patient_id,
                            anonym_sample_id))%>%unique(), by=c("sample_short")) %>%
  filter(((vj %in% c("TCRBV15", "TCRBV05"))&(cells == "N")) |
           ((vj %in% c("TCRBV18", "TCRBV19", "TCRBV20")) & cells == "E"))%>%
  mutate(
    genotype = case_when(genotype_short =="homoDQ2" ~ "DQ2",
                         genotype_short =="heteroDQ2DQ8" ~"DQ2DQ8",
                         genotype_short =="homoDQ8" ~ "DQ8"),
    genotype = factor(genotype, levels = c("DQ2", "DQ2DQ8",  "DQ8")),
    cells = ifelse(cells == "E", "MEMORY", "NAIVE"),
    stage = c("stage 2", "stage 3")[grepl(pat="F", sample_short)+1]
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


###### frequency of VJ usage, naive versus memory
V_J %>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
  #  filter(value > 0) %>%
  mutate(sample_short=gsub(pat="_subs", rep="",sample_short))%>%
  left_join(anno%>%select(c(cells,
                            patient,
                            genotype_short,
                            sample_short,source,
                            anonym_patient_id,
                            anonym_sample_id))%>%unique(), by=c("sample_short")) %>%
  filter(vj %in% c("TCRBV16", "TCRBV02"))%>%
  mutate(
    genotype = case_when(genotype_short =="homoDQ2" ~ "DQ2",
                         genotype_short =="heteroDQ2DQ8" ~"DQ2DQ8",
                         genotype_short =="homoDQ8" ~ "DQ8"),
    genotype = factor(genotype, levels = c("DQ2", "DQ2DQ8",  "DQ8")),
    cells = ifelse(cells == "E", "MEMORY", "NAIVE"),
    stage = c("stage 2", "stage 3")[grepl(pat="F", sample_short)+1]
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

######




###### PCA naive
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


###### PCA memory
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


