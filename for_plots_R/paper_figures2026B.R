
library(tidyverse)
library(kableExtra)
library(scales)
library(ggpubr)

#knitr::opts_chunk$set( message=FALSE, fig.width = 8, fig.height = 12, warning = FALSE)


source("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/plot_cdr3.R") 

annotation_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/data/collated_info.csv"

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


##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### CDR3 over/underrepresentation
##### ##### ##### ##### ##### ##### ##### ##### #####

results_imgt  <- read_tsv(file.path(imgt_tests_loc, file.path(paste0(imgt_prefixes[1], imgt_tests[2],imgt_prefixes[2]),file))) 

#  significant between homozygotes:
##summary
results_imgt%>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))),
         p_value_oo <0.005,
         )

res_imgt <- results_imgt%>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  rename("aa"=AA)

#KF4, KF10; one celltype=608 tests 
plot_enrichment_depletion3 (res_imgt%>%filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8),  length_col="length", estim_col="cohens_d_oo",
                            estim_expr_for_plotting="cohens_d_oo>0", 
                            aminoacid_mapping=color_scale_properties_property(AAdata$kideraFactors,"KF4") , inset_colour="red")
plot_enrichment_depletion3 (res_imgt%>%filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8),  length_col="length", estim_col="cohens_d_oo",
                            estim_expr_for_plotting="cohens_d_oo<0", 
                            aminoacid_mapping=color_scale_properties_property(AAdata$kideraFactors,"KF4") , inset_colour="blue")


res_imgt%>%filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8)# 53 results (6 exp)


res_imgt%>% 
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8, cohens_d_oo>0)
###the same for E cells


############ WL

results_imgt_WL  <- read_tsv(file.path(imgt_tests_loc, file.path(paste0(imgt_prefixes[1], imgt_tests[1],imgt_prefixes[2]),file))) 

results_imgt_WL%>% #one celltype=158 tests, 19 significant (2 exp)
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  rename("aa"=AA)%>%filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8)
####### PLOT

############ VFAM

#VFAM which were significantly different between DQ2, DQ8
results_imgt_VFAM  <- read_tsv(file.path(imgt_tests_loc, file.path(paste0(imgt_prefixes[1], imgt_tests[4],imgt_prefixes[2]),file))) 

#used more by DQ2
results_imgt_VFAM %>% 
filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  rename("aa"=AA)%>%filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8, cohens_d_oo>0,
                           !(vFamilyName%in%c( signif_vj_oo%>%filter(cells=="N",statistic>0)%>%.$vj)))


###compare results with, without VFAM
naive=
 results_imgt %>% 
   filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
    mutate(imgt =change_imgt_format(IMGT_position))%>%
 rename("aa"=AA)%>%filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8)%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group")))

VFAM=
  results_imgt_VFAM %>% 
 filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
    mutate(imgt =change_imgt_format(IMGT_position))%>%
   rename("aa"=AA)%>%filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8)%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group")))
VFAM_all=
  results_imgt_VFAM %>% 
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  rename("aa"=AA)%>%filter(cells=="N")%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group")))

#Attempt to show that these aminoacids are not driven by TRBV gene - for significant aminoacids, when testing is repeated AFTER splitting into TRBV gene groups, some groups remain significant (eventhough much lower power etc)
left_join(naive,
          VFAM_all,  by=c("cells","IMGT_position","imgt","aa","length", "label_group_a_oo" ,"label_group_b_oo"))%>%
  group_by( IMGT_position, aa ,   length)%>%
  summarise(TRBV_tested=length(vFamilyName),
            TRBV_sign=sum(p_value_oo.y<0.05, na.rm=T),
            TRBV_sign_and_vjdiff=sum((vFamilyName%in%(signif_vj_oo%>%filter(cells=="N")%>%.$vj))&p_value_oo.y<0.05))%>%arrange(length, IMGT_position, aa )

#######what about simple correlation Vgene versus AA (across lengths/imgt positions)?

#####################
##### PROPERITES
###################

properties_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/aa_properties"
properties <- c("KF", "VHS")
properties_tests <- c("WL",  "full",  "WL_vfam",     "full_vfam")
properties_res <- c("hotelling", "per_factor")
properties_prefixes <-c("comb_aa_imgt_", "_subs_rows_","_anno_")
#"comb_aa_imgt_" (WL) "_subs_rows_", (VHS) "_anno_" (VHS) (hotelling).tsv

results_KF <- read_tsv(file.path(properties_loc,
                                 file.path(properties[1],
                                           file.path(paste0(properties[1],"_properties_results"),
                                                     paste0(properties_prefixes[1],properties_tests[2],properties_prefixes[2], properties[1] ,properties_prefixes[3],properties[1] ,"_",properties_res[1], ".tsv")))))

results_KF %>%filter(!is.na(q_value_oo))%>%filter(!(IMGT_position%in%as.character(c(104:110, 113:118))),q_value_oo<0.05)



results_KF_by_fact <- read_tsv(file.path(properties_loc,
                                 file.path(properties[1],
                                           file.path(paste0(properties[1],"_properties_results"),
                                                     paste0(properties_prefixes[1],properties_tests[2],properties_prefixes[2], properties[1] ,properties_prefixes[3],properties[1] ,"_",properties_res[2], ".tsv")))))


results_KF_by_fact%>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8,q_hotelling_oo<0.05)%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group")))%>%arrange(cells, IMGT_position, prop, length)%>%data.frame()


results_KF_WL_VFAM_by_fact <- read_tsv(file.path(properties_loc,
                                         file.path(properties[1],
                                                   file.path(paste0(properties[1],"_properties_results"),
                                                             paste0(properties_prefixes[1],properties_tests[3],properties_prefixes[2], properties[1] ,properties_prefixes[3],properties[1] ,"_",properties_res[2], ".tsv")))))


results_KF_WL_VFAM_by_fact %>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8,q_hotelling_oo<0.05)%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group")))%>%arrange(cells, IMGT_position, prop, vFamilyName)%>%data.frame()
