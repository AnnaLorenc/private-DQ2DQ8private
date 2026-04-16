#OOF


library(tidyverse)
library(kableExtra)
library(scales)
library(ggpubr)
library(gt)
library(dplyr)
library(rstatix)

#knitr::opts_chunk$set( message=FALSE, fig.width = 8, fig.height = 12, warning = FALSE)


source("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/plot_cdr3_functions.R") 

annotation_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/data/collated_info.csv"
plot_dir <- "/Users/ania/Documents/DQ2DQ8/pipeline/260304"

colors =c(DQ2="#FF4D4D",
          DQ2DQ8="#B455E0",
          DQ8="#2E86DE")

colors_two <- c("#E69F00", "#009E73")


cell_fullnames <- c(N = "Naive", E = "Memory")
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

new_param_names <- c(total_templates="TRB_templates",
                     total_rearrangements="TRB_uniques",
                     unique_seqs="TRB_clean_uniques",
                     exp_Shannon="Shannon_diversity", 
                     inv_Simpson="inv_Simpson_index" )



plot_dir <- file.path(plot_dir, "oof" )
dir.create(plot_dir, recursive = TRUE)


diversity_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_metrics/combined_nonproductive_diversity_metrics.csv"

vj_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_freqs/combined_nonproductive_freqs.csv"

divers_oof <-read_csv(diversity_loc)

diversities_anno <- divers_oof%>%
  left_join(., anno%>%select(c(cells,
                               patient,
                               genotype_short,
                               sample_short,source,
                               anonym_patient_id,
                               anonym_sample_id)), by=c("sample_name"="sample_short"))%>%
  pivot_longer(cols =unique_sequences:shannon_diversity, names_to = "measure" )%>%
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

plots_div <-list()

for(meas in  diversities_anno$measure%>%unique()){
  for( cells_s in c("NAIVE", "MEMORY")){
    plots_div[[paste(cells_s, meas, sep="_")]] <-plot_paper_plot (i=meas, cells="NAIVE", diversities_anno)+ ggtitle(paste(meas, ", raw nonproductive, cells ", cells_s))
    ggsave(plot= plots_div[[paste(cells_s, meas, sep="_")]], path = plot_dir, filename = paste0(meas, "_raw_nonpr_", cells_s,".png"), height = 5, width = 7)
}
}

#just counts/unique sequences
prod_nonprod <-read_csv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/cleanup/cleanup_summary_all.csv")%>%
  left_join(., anno%>%select(c(cells,
                               patient,
                               genotype_short,
                               sample_short,source,
                               anonym_patient_id,
                               anonym_sample_id)), by=c("sample"="sample_short"))
  
prod_nonprod%>%group_by(cells, genotype_short)%>%slice_min(counts_nonprod, n=2)

prod_nonprod%>%filter(num_seq_nonprod>2000)%>%group_by(cells, genotype_short)%>%summarise(n=n())

###About fraction of productive vs nonproductive 
prod_nonprod%>%
  mutate(geno=gsub(genotype_short, pat=".*mo|.*ro", rep=""), cells=cell_fullnames[cells], 
  stage=c('stage 2', 'stage 3')[as.integer(source=='F')+1])%>%
  ggplot(aes(x=geno, y=num_seq_nonprod/(num_seq_nonprod+num_seq_prod), col=geno))+geom_boxplot(outliers = FALSE)+
  geom_jitter(aes(shape=stage, col=geno), width=0.1, size=3)+
  facet_wrap(~cells)+ theme_bw()+ scale_color_manual(values=colors)+ylab('Nonproductive in unique sequences')+scale_shape_manual(values=c(1,19))
ggsave(path = plot_dir, filename="nonprod_in_all_unique.png", width=7.5, height =5)


prod_nonprod%>%
  mutate(geno=gsub(genotype_short, pat=".*mo|.*ro", rep=""), cells=cell_fullnames[cells], 
         stage=c('stage 2', 'stage 3')[as.integer(source=='F')+1])%>%
  ggplot(aes(x=geno, y=counts_nonprod/(counts_nonprod+counts_prod), col=geno))+geom_boxplot(outliers = FALSE)+
  geom_jitter(aes(shape=stage, col=geno), width=0.1, size=3)+
  facet_wrap(~cells)+ theme_bw()+ scale_color_manual(values=colors)+ylab('Nonproductive in all TCRs')+scale_shape_manual(values=c(1,19))
ggsave(path = plot_dir, filename="nonprod_in_all_counts.png", width=7.5, height=5)

########## VJ
vj_oof <- read_csv(vj_loc)
vj_oof_anno <- vj_oof %>%
  left_join(., anno%>%select(c(cells,
                               patient,
                               genotype_short,
                               sample_short,source,
                               anonym_patient_id,
                               anonym_sample_id)), by=c("sample_name"="sample_short"))%>%
  unique()
#comparison with productive?
which_often_enough <- vj_oof_anno%>%filter(param=="vFamilyName")%>%group_by(group,cells)%>%
  slice_min(row_fraction)%>%
  filter(row_fraction>0.01)%>%
  select(cells, group)




vj_oof_anno%>%filter(param=="vFamilyName")%>%
group_by(cells, group)%>%
  filter(all(table(genotype_short) >= 2)) %>% 
  t_test(row_fraction ~ genotype_short,var.equal=FALSE,p.adjust.method="fdr",detailed=T)%>%
  filter(p*20*2<0.01)# Nothing stays significant after correcting for multiple testing; also for J


###subsampled results
vj_subs_loc <- '/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/combined_freqs/nonproductive/nonproductive_subs_freqs.tsv'
#I use subsampled values, select the genes present as at least 1% of repertoire in any sample and perform robust ttest
vj_subs_oof <- read_tsv(vj_subs_loc)%>%
  filter(param%in%c("vFamilyName","jGeneName"))%>%
  select(vj=group, ends_with("subs")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

signif_vj_oof <- vj_subs_oof  %>%
  rowwise%>%filter(any(across(ends_with("subs") )>0.01))%>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
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
  t_test(value ~ geno,var.equal=FALSE,p.adjust.method="fdr",detailed=T)


eff_size <- vj_subs_oof  %>%
  rowwise%>%filter(any(across(ends_with("subs") )>0.01))%>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
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
  cohens_d(value ~ geno)

signif_vj_oof <- signif_vj_oof %>%
  add_column(cohens_d=eff_size$effsize)

vj_subs_oof%>%
  filter(vj%in%signif_vj_oof$vj)%>%
  rowwise%>%filter(any(across(ends_with("subs") )>0.01))%>%
  pivot_longer(cols = c(ends_with("subs")), names_to = "sample_short") %>%
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
  ggplot(aes(x=genotype_short, y=value, col=genotype_short)) +geom_boxplot()+
  facet_wrap(~cells, scales="free_y")


signif_vj_oof %>%
  group_by(cells)%>%
  # filter(group1%in%c("DQ2","DQ8"),  group2%in%c("DQ2","DQ8"))%>%
  ungroup()%>%
  rowwise()%>%
  mutate(p_fully_adjusted= min(p*38*2,1))%>% 
  group_by(cells)%>%
  slice_min(order_by = p, n = 5)%>%select(-starts_with("p.adj"))
# vj      cells estimate estimate1 estimate2 .y.   group1 group2    n1    n2 statistic     p    df  conf.low  conf.high method alternative p_fully_adjusted
# <chr>   <chr>    <dbl>     <dbl>     <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl>     <dbl>      <dbl> <chr>  <chr>                  <dbl>
#   1 TCRBV03 E      0.00661    0.0451    0.0385 value DQ2DQ8 DQ8       20    21      3.20 0.003  37.4  0.00242   0.0108    T-test two.sided              0.228
# 2 TCRBV20 E      0.00545    0.0553    0.0498 value DQ2    DQ2DQ8    20    20      2.61 0.013  34.9  0.00121   0.00969   T-test two.sided              0.988
# 3 TCRBV20 E     -0.00416    0.0498    0.0540 value DQ2DQ8 DQ8       20    21     -2.61 0.013  37.1 -0.00740  -0.000926  T-test two.sided              0.988
# 4 TCRBV13 E     -0.00243    0.0120    0.0145 value DQ2DQ8 DQ8       20    21     -2.53 0.016  33.0 -0.00439  -0.000479  T-test two.sided              1    
# 5 TCRBV01 E      0.00249    0.0200    0.0175 value DQ2    DQ2DQ8    20    20      2.49 0.017  36.1  0.000466  0.00451   T-test two.sided              1    
# 6 TCRBV09 N     -0.00197    0.0325    0.0344 value DQ2    DQ2DQ8    20    20     -2.18 0.036  33.7 -0.00380  -0.000136  T-test two.sided              1    
# 7 TCRBV06 N      0.00557    0.0954    0.0898 value DQ2    DQ8       20    21      2.12 0.04   38.1  0.000263  0.0109    T-test two.sided              1    
# 8 TCRBV07 N     -0.00369    0.0798    0.0835 value DQ2    DQ2DQ8    20    20     -2.12 0.041  37.8 -0.00723  -0.000159  T-test two.sided              1    
# 9 TCRBV21 N     -0.00256    0.0535    0.0561 value DQ2DQ8 DQ8       20    21     -1.99 0.054  38.0 -0.00517   0.0000481 T-test two.sided              1    
# 10 TCRBV20 N     -0.00284    0.0490    0.0518 value DQ2DQ8 DQ8       20    21     -1.97 0.056  37.5 -0.00575   0.0000787 T-test two.sided              1    

signif_vj_oof %>%
  group_by(cells)%>%
 filter(group1%in%c("DQ2","DQ8"),  group2%in%c("DQ2","DQ8"))%>%
  ungroup()%>%
  rowwise()%>%
  mutate(p_fully_adjusted= min(p*38,1))%>% 
  group_by(cells)%>%
  slice_min(order_by = p, n = 5)%>%select(-starts_with("p.adj"))

# vj         cells estimate estimate1 estimate2 .y.   group1 group2    n1    n2 statistic     p    df  conf.low conf.high method alternative p_fully_adjusted
# <chr>      <chr>    <dbl>     <dbl>     <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl>     <dbl>     <dbl> <chr>  <chr>                  <dbl>
#   1 TCRBJ02-02 E      0.00252    0.0718    0.0692 value DQ2    DQ8       20    21      1.74 0.09   38.6 -0.000411  0.00545  T-test two.sided                  1
# 2 TCRBV01    E      0.00149    0.0200    0.0185 value DQ2    DQ8       20    21      1.57 0.127  34.2 -0.000444  0.00343  T-test two.sided                  1
# 3 TCRBJ02-06 E      0.00133    0.0373    0.0359 value DQ2    DQ8       20    21      1.47 0.148  38.8 -0.000494  0.00315  T-test two.sided                  1
# 4 TCRBV28    E     -0.00400    0.0332    0.0372 value DQ2    DQ8       20    21     -1.32 0.196  39.0 -0.0101    0.00215  T-test two.sided                  1
# 5 TCRBJ02-04 E      0.00157    0.0456    0.0440 value DQ2    DQ8       20    21      1.21 0.236  36.5 -0.00107   0.00421  T-test two.sided                  1
# 6 TCRBV06    N      0.00557    0.0954    0.0898 value DQ2    DQ8       20    21      2.12 0.04   38.1  0.000263  0.0109   T-test two.sided                  1
# 7 TCRBV30    N     -0.00168    0.0272    0.0289 value DQ2    DQ8       20    21     -1.76 0.086  38.9 -0.00362   0.000252 T-test two.sided                  1
# 8 TCRBV28    N     -0.00558    0.0324    0.0379 value DQ2    DQ8       20    21     -1.75 0.088  38.3 -0.0120    0.000871 T-test two.sided                  1
# 9 TCRBV18    N     -0.00172    0.0254    0.0271 value DQ2    DQ8       20    21     -1.72 0.093  39.0 -0.00373   0.000299 T-test two.sided                  1
# 10 TCRBJ02-06 N      0.00174    0.0375    0.0358 value DQ2    DQ8       20    21      1.72 0.097  29.5 -0.000332  0.00381  T-test two.sided                  1
pval=0.01
groups <- signif_vj_oof%>%.$vj%>%unique()%>%length()

signif_vj_oof%>%
  filter(cells=="N")%>%
  ggplot(aes(x=p)) +geom_histogram(fill="white",col='darkgrey',breaks=c(seq(from=0, to=1,by=pval)) ) +
  geom_vline(xintercept=pval, lty=2)+geom_hline(yintercept=round(pval*groups*3), lty=2) +theme_bw()+facet_wrap(~cells, scale="free_y") +
  ggtitle("VJ usage in OOF naives: genotypes comparison",subtitle = paste("Dashed lines:vertical sign pvalue",pval,", horizontal:expected hits"))
ggsave(path = plot_dir, filename="nonprod_VJ_pvals_naive.png", width=7.5, height =3)

vj_oof_forpca <- vj_oof%>%
  filter(param%in%c("vFamilyName","jGeneName"))%>%
  select(vj=group, row_fraction,sample_name) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))%>%left_join(anno%>%select(c(cells,
                                                                                     patient,
                                                                                     genotype_short,
                                                                                     sample_name=sample_short,source,
                                                                                     anonym_patient_id,
                                                                                     anonym_sample_id,
                                                                                     geno,
                                                                                     cells_long,
                                                                                     stage))%>%unique(), by=c("sample_name")) %>%  mutate(geno=case_when(genotype_short =="homoDQ2" ~ "DQ2",
                                                                                                                                                         genotype_short =="heteroDQ2DQ8" ~"DQ2DQ8",
                                                                                                                                                         genotype_short =="homoDQ8" ~ "DQ8"))


vj_oof_forpca  <- vj_oof_forpca%>%
  select(vj, row_fraction,sample_name, geno, cells, stage)%>%
  filter(cells=="N")%>%
  pivot_wider(names_from = vj,values_from = row_fraction,values_fill = 0, id_cols = c(sample_name, geno, cells, stage))


pca_mat <- vj_oof_forpca%>%select(c(sample_name, starts_with("TCRBV")))%>%
  data.frame(row.names = "sample_name")%>%
  as.matrix()

pca_res   <- prcomp(pca_mat, scale. = TRUE)
scores    <- as.data.frame(pca_res$x)
scores$Group <- vj_oof_forpca$geno
scores$stage <- vj_oof_forpca$stage

ggplot(scores, aes(PC1, PC2, color = Group, shape=Group)) +
  geom_point(alpha = 0.7, size=3) +
  stat_ellipse(type = "t", linetype = 1) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
  ) +coord_fixed() +scale_color_manual(values=colors)+
  ggtitle("NAIVE OOF cells, VJ usage")+
  scale_shape_manual(values=c(19,8,19))
ggsave(path = plot_dir, filename="nonprod_V_pca_naive_geno.png",  width=7.5, height = 7)

ggplot(scores, aes(PC1, PC2, color = stage)) +
  geom_point(alpha = 0.7, size=2) +
  stat_ellipse(type = "t", linetype = 1) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)")
  ) +coord_fixed() +scale_color_manual(values=colors_two)+
  ggtitle("NAIVE OOF cells, VJ usage")
ggsave(path = plot_dir, filename="nonprod_V_pca_naive_stage.png", width=7.5, height =6.5)


####NAIVE lengths - for oof lengths are in nucleotides
vj_oof_anno%>%
  filter(param=="cdr3Length")%>%
  mutate(len=as.integer(group),
         geno=gsub(pat="h.*o",rep="",genotype_short),
         type=c("out","stop")[as.integer(len%%3==0)+1],
         stage=c("stage 1", "stage 3")[(source=="T")+1])%>%
filter(cells=="N",len>=30, len<=66, type=="out")%>%
  ggplot(aes(x=len, y=row_fraction, col=geno, group=interaction(geno, len))) +geom_boxplot(outliers = FALSE)+
scale_color_manual(values=colors)+ theme_classic()+
  ggtitle("Length (nucleotides) of non-productive TRBs in naive cells")

ggsave(path = plot_dir, filename="nonprod_length_naive_geno.png",  width = 7.5, height = 6)


vj_oof_anno%>%
  filter(param=="cdr3Length")%>%
  mutate(len=as.integer(group),
         geno=gsub(pat="h.*o",rep="",genotype_short),
         type=c("out","stop")[as.integer(len%%3==0)+1],
         stage=c("stage 1", "stage 3")[(source=="T")+1])%>%
  filter(cells=="N")%>%
  ggplot(aes(x=len, y=row_fraction, col=stage, group=interaction(stage, len))) +geom_boxplot(outliers = FALSE)+
  scale_color_manual(values=colors_two)+ theme_bw()+facet_wrap(~type, scales="free",nrow=2)+
  ggtitle("Length (nucleotides) of non-productive TRBs in naive cells")

ggsave(path = plot_dir, filename="nonprod_length_naive_stage.png", width=12, height = 8)


