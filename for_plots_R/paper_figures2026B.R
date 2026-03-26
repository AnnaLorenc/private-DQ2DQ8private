
library(tidyverse)
library(kableExtra)
library(scales)
library(ggpubr)
library(gt)
library(dplyr)

#knitr::opts_chunk$set( message=FALSE, fig.width = 8, fig.height = 12, warning = FALSE)


source("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/plot_cdr3_functions.R") 
  
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

  
#####
#####OLD plots
#####
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

#####


############ WL

results_imgt_WL  <- read_tsv(file.path(imgt_tests_loc, file.path(paste0(imgt_prefixes[1], imgt_tests[1],imgt_prefixes[2]),file))) 

######################## 
############ VFAM
########################


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



#####################
##### PROPERITES
###################

properties_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/aa_properties"
properties <- c("KF", "VHS")
properties_tests <- c("WL",  "full",  "WL_vfam",     "full_vfam")
properties_res <- c("hotelling", "per_factor")
properties_prefixes <-c("comb_aa_imgt_", "_subs_rows_","_anno_")
#"comb_aa_imgt_" (WL) "_subs_rows_", (VHS) "_anno_" (VHS) (hotelling).tsv



# most aggregated version:
results_KF_WL_by_fact <- read_tsv(file.path(properties_loc,
                                            file.path(properties[1],
                                                      file.path(paste0(properties[1],"_properties_results"),
                                                                paste0(properties_prefixes[1],properties_tests[1],properties_prefixes[2], properties[1] ,properties_prefixes[3],properties[1] ,"_",properties_res[2], ".tsv")))))



results_KF_by_fact <- read_tsv(file.path(properties_loc,
                                 file.path(properties[1],
                                           file.path(paste0(properties[1],"_properties_results"),
                                                     paste0(properties_prefixes[1],properties_tests[2],properties_prefixes[2], properties[1] ,properties_prefixes[3],properties[1] ,"_",properties_res[2], ".tsv")))))


results_KF_WL_VFAM_by_fact <- read_tsv(file.path(properties_loc,
                                         file.path(properties[1],
                                                   file.path(paste0(properties[1],"_properties_results"),
                                                             paste0(properties_prefixes[1],properties_tests[3],properties_prefixes[2], properties[1] ,properties_prefixes[3],properties[1] ,"_",properties_res[2], ".tsv")))))


results_KF_WL_VFAM_by_fact %>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8,q_hotelling_oo<0.05)%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group")))%>%arrange(cells, IMGT_position, prop, vFamilyName)%>%data.frame()


#240, 15 significant at qval 1%
results_KF_WL_by_fact%>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8,q_hotelling_oo<0.05)%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group")))%>%arrange(cells, IMGT_position, prop)%>%data.frame()


####Is it dependent on the Vgene? Examples when it is still significant withn specific V genes
#significant:
df <- left_join(
results_KF_WL_by_fact %>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_oo <0.05, cells=="N",abs(cohens_d_oo)>0.8,q_hotelling_oo<0.05)%>%
  select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group"),IMGT_position, ends_with("hotelling_oo"))),

#significant, stratified by V:
results_KF_WL_VFAM_by_fact %>%
  filter(!(IMGT_position%in%as.character(c(104:110, 113:118))))%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_oo <0.01, cells=="N",abs(cohens_d_oo)>0.8,q_hotelling_oo<0.05)%>%select(-c(ends_with("eo"), starts_with("t_stat"), starts_with("n_group"),IMGT_position, ends_with("hotelling_oo"))), by=c("imgt","cells", "prop"), suffix = c("_f", "v")                                                                                                                                                                                    
)%>%
  left_join(., 
            
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
              group_by(cells, vj)%>%
              t_test(value ~ geno,var.equal=FALSE,p.adjust.method="fdr",detailed=T)%>%
              filter(group1=="DQ2",group2=="DQ8")%>%
              select(vj,cells, freq_DQ2=estimate1,to_DQ8=estimate)%>%unique(), by=c("cells","vFamilyName"="vj"))



#######################################
################### Table printing for significant results WITH/WITHOUT VJ
#######################################

# 
df %>%
  gt(groupname_col = "cells") %>% # Groups by 'E' and 'N'
  
  # 1. Format column labels to be more readable
  cols_label(
    prop = "Property",
    mean_a_oo_f = "Mean A",
    mean_b_oo_f = "Mean B",
    cohens_d_oo_f = "Cohen's D",
    p_value_oo_f = "p-value",
    q_value_oo_f = "q-value",
    imgt = "IMGT",
    vFamilyName = "V Family",
    mean_a_oov = "Mean A (V)",
    mean_b_oov = "Mean B (V)",
    freqV_DQ2 = "Freq DQ2",
    freqV_diffto_DQ8 = "Diff DQ8"
  ) %>%
  
  # 2. Handle NA values (makes the table much cleaner)
  fmt_missing(columns = everything(), missing_text = "-") %>%
  
  # 3. Format Decimals
  fmt_number(
    columns = c(contains("mean"), contains("cohens"), contains("freq")),
    decimals = 3
  ) %>%
  
  # 4. Format P-values and Q-values with scientific notation
  fmt_scientific(
    columns = c(contains("p_value"), contains("q_value")),
    decimals = 2
  ) %>%
  
  # 5. Visual Styling
  tab_header(
    title = "Analysis of Cell Properties and V-Family Frequencies",
    subtitle = "Comparing DQ2 and DQ8 groups"
  ) %>%
  tab_options(
    table.font.size = px(12),
    column_labels.font.weight = "bold",
    row_group.background.color = "#f9f9f9"
  ) %>%
  
  # 6. Highlight significant results (optional)
  tab_style(
    style = cell_text(color = "red", weight = "bold"),
    locations = cells_body(
      columns = q_value_oo_f,
      rows = q_value_oo_f < 0.05
    )
  )


# ##examples: N, KF10, 111: across and in TRBV03,06,19
# cells prop  n_a_oo_f n_b_oo_f mean_a_oo_f mean_b_oo_f mean_diff_oo_f cohens_d_oo_f df_oo_f  p_value_oo_f q_value_oo_f imgt  vFamilyName n_a_oov n_b_oov mean_a_oov mean_b_oov mean_diff_oov cohens_d_oov df_oov p_value_oov  q_value_oov freq_DQ2
# <chr> <chr>    <dbl>    <dbl>       <dbl>       <dbl>          <dbl>         <dbl>   <dbl>         <dbl>        <dbl> <chr> <chr>         <dbl>   <dbl>      <dbl>      <dbl>         <dbl>        <dbl>  <dbl>       <dbl>        <dbl>    <dbl>
#   1 N     KF10        20       21     0.140        0.124          0.0166         2.31     38.9 0.00000000583  0.000000406 111   TCRBV03          20      21      0.170      0.141        0.0289        1.14    38.8    7.25e- 4  0.0258        0.0289
# 2 N     KF10        20       21     0.140        0.124          0.0166         2.31     38.9 0.00000000583  0.000000406 111   TCRBV06          20      21      0.139      0.100        0.0389        2.78    36.6    1.25e-10  0.000000223   0.122 
# 3 N     KF10        20       21     0.140        0.124          0.0166         2.31     38.9 0.00000000583  0.000000406 111   TCRBV19          20      21      0.146      0.122        0.0245        0.946   36.6    4.69e- 3  0.0920        0.0458
# 4 N     KF5         20       21    -0.222       -0.206         -0.0166        -1.08     36.3 0.00146        0.0135      111   TCRBV03          20      21     -0.258     -0.224       -0.0338       -1.24    36.0    3.13e- 4  0.0155        0.0289
# 5 N     KF5         20       21    -0.222       -0.206         -0.0166        -1.08     36.3 0.00146        0.0135      111   TCRBV04          20      21     -0.270     -0.243       -0.0279       -1.42    38.3    5.35e- 5  0.00576       0.0526
# 6 N     KF5         20       21    -0.222       -0.206         -0.0166        -1.08     36.3 0.00146        0.0135      111   TCRBV06          20      21     -0.208     -0.187       -0.0213       -0.938   34.3    5.31e- 3  0.102         0.122 
# 7 N     KF5         20       21    -0.222       -0.206         -0.0166        -1.08     36.3 0.00146        0.0135      111   TCRBV07          20      21     -0.221     -0.200       -0.0211       -1.06    35.5    1.76e- 3  0.0468        0.116 
# 8 N     KF5         20       21    -0.222       -0.206         -0.0166        -1.08     36.3 0.00146        0.0135      111   TCRBV19          20      21     -0.235     -0.209       -0.0265       -0.886   36.1    7.81e- 3  0.135         0.0458
# 9 N     KF3         20       21     0.00247      0.0255        -0.0230        -1.05     38.1 0.00180        0.0144      111.2 NA               NA      NA     NA         NA           NA            NA       NA     NA        NA            NA     
# 10 N     KF2         20       21    -0.671       -0.693          0.0222         1.06     38.9 0.00153        0.0135      112   TCRBV04          20      21     -0.639     -0.677        0.0377        0.870   38.1    8.45e- 3  0.142         0.0526
# 11 N     KF2         20       21    -0.671       -0.693          0.0222         1.06     38.9 0.00153        0.0135      112   TCRBV06          20      21     -0.656     -0.694        0.0385        1.32    38.8    1.31e- 4  0.01000       0.122 
# 12 N     KF2         20       21    -0.671       -0.693          0.0222         1.06     38.9 0.00153        0.0135      112   TCRBV20          20      21     -0.656     -0.686        0.0296        0.965   38.9    3.69e- 3  0.0780        0.111 
# 13 N     KF4         20       21     0.101        0.112         -0.0108        -0.845    39.0 0.0101         0.0565      112   TCRBV18          20      21      0.112      0.140       -0.0274       -1.06    39.0    1.59e- 3  0.0440        0.0328
# 14 N     KF9         20       21    -0.462       -0.475          0.0132         1.12     37.9 0.00101        0.0112      112   TCRBV06          20      21     -0.466     -0.488        0.0214        1.06    36.9    1.52e- 3  0.0423        0.122 
# 15 N     KF9         20       21    -0.462       -0.475          0.0132         1.12     37.9 0.00101        0.0112      112   TCRBV20          20      21     -0.441     -0.470        0.0286        1.32    38.7    1.41e- 4  0.0106        0.111 
# 16 N     KF3         20       21     0.0230       0.0478        -0.0248        -1.13     26.8 0.00138        0.0135      112.1 NA               NA      NA     NA         NA           NA            NA       NA     NA        NA            NA     
# 17 N     KF4         20       21     0.0897       0.119         -0.0297        -1.49     38.5 0.0000266      0.000520    112.1 TCRBV11          20      21      0.114      0.164       -0.0502       -0.890   38.1    7.17e- 3  0.126         0.0388
# 18 N     KF8         20       21     0.244        0.207          0.0372         1.16     38.9 0.000615       0.00790     112.1 NA               NA      NA     NA         NA           NA            NA       NA     NA        NA            NA     
# 19 N     KF3         20       21     0.0191       0.0536        -0.0345        -0.948    28.9 0.00560        0.0359      112.2 NA               NA      NA     NA         NA           NA            NA       NA     NA        NA            NA     

#######################################
################### PLOTS #############
#######################################

source("~/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/plot_aa_logo.R")



min_freq = 0.01
p_cutoff = 0.05
d_cutoff = 0.8
q_cutoff = 0.1

cell_fullnames <- c(N = "naive", E = "memory")


###############
#### DQ2 vs DQ8


cols <-list(imgt = "imgt",
            aa = "aa",
            mean_a = "mean_group_a_oo",
            mean_b = "mean_group_b_oo", 
            cohens_d = "cohens_d_oo",
            p_value = "p_value_oo",
            label_a = "label_group_a_oo",
            label_b = "label_group_b_oo")

test_version <- "imgt"
dir.create(file.path(plot_dir, test_version ))




for(cells_sel in c("N","E")){
  
  p <- plot_aa_composition_logo( results_imgt_WL%>%filter(cells==cells_sel)%>%
                                   mutate(imgt =change_imgt_format(IMGT_position))%>%rename("aa"=AA),
                                 min_freq = min_freq, use_seqlogo = TRUE,
                                 p_cutoff = p_cutoff, d_cutoff = d_cutoff, cols=cols,
                                 add_to_title = paste(cell_fullnames[cells_sel],"cells"))
  ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL", cells_sel, sep="_"), ".png"))
  
  ##only significant ones
  
  #for Bonferroni correction
  bonf <- results_imgt_WL%>%filter(cells==cells_sel)%>%
    mutate(imgt =change_imgt_format(IMGT_position))%>%rowwise()%>%
    filter(imgt%in%c(as.character(106:111),'111.1',"111.2","112.1",as.character(112:117)), max(mean_group_a_oo, mean_group_b_oo)>0.01)%>%
    rename("aa"=AA)%>%select(imgt, aa,ends_with("oo"))%>%group_by(imgt)%>%summarise(a=n_distinct(aa))%>%summarise(sum(a))%>%unlist()
  
  q <- plot_plogo_like(df= results_imgt_WL%>%filter(cells==cells_sel)%>%
                    mutate(imgt =change_imgt_format(IMGT_position))%>%
                    filter(imgt%in%c(as.character(106:111),'111.1',"111.2","112.1",as.character(112:116)))%>%
                    rename("aa"=AA)%>%select(imgt, aa,ends_with("oo"))%>%
                    mutate(label_group_a_oo=gsub(pat="ho", rep="",label_group_a_oo), label_group_b_oo=gsub(pat="ho", rep="",label_group_b_oo) ),
                  height_by="cohens_d",
                  cols=cols, p_cutoff = p_cutoff/bonf, d_cutoff = d_cutoff,
                  add_to_title = paste(cell_fullnames[cells_sel],"cells, Bonferroni:",bonf))
  
  ggsave(plot=q, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL_sign", cells_sel, sep="_"), ".png"))
  
  
  ###one length plots
  for(len_sel in 13:16){
    p <- plot_aa_composition_logo( results_imgt%>%filter(cells==cells_sel, length==len_sel)%>%
                                     mutate(imgt =change_imgt_format(IMGT_position))%>%rename("aa"=AA),
                                   min_freq = min_freq, use_seqlogo = TRUE,
                                   p_cutoff = p_cutoff, d_cutoff = d_cutoff, cols=cols,
                                   add_to_title =  paste(cell_fullnames[cells_sel],"cells, length=", len_sel))
    ggsave(plot = p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version, paste0("l",len_sel), cells_sel, sep="_"),".png"))
  }
} 

test_version <- "prop"
dir.create(file.path(plot_dir, test_version ))
cols = list(
  imgt     = "imgt",
  prop       = "prop",
  cells = "cells",
  mean_a   = "mean_a_oo",
  mean_b   = "mean_b_oo",
  cohens_d = "cohens_d_oo",
  p_value  = "p_value_oo",
  q_value = "q_value_oo",
  label_a  = "label_group_a_oo",
  label_b  = "label_group_b_oo"
)

for(cells_sel in c("N","E")){
  
  p <-  plot_property_logo(
    results_KF_WL_by_fact%>%
      mutate(imgt =change_imgt_format(IMGT_position))%>%select(cells, imgt, prop,ends_with("oo"))%>%select(-c(contains('hotelling')))%>%mutate(label_group_a_oo="DQ2",label_group_b_oo="DQ8"),
    q_cutoff   = q_cutoff,
    d_cutoff   = d_cutoff,
    cols = cols,
    cell_types = cells_sel,
    color_a    = "#e41a1c",   # red = group A
    color_b    = "#377eb8"    # blue = group B
    ,positions = c("108","109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"),
    add_to_title = paste(cell_fullnames[cells_sel],"cells"))
  ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL", cells_sel, sep="_"), ".png"))
  
  
  for(len_sel in 13:16){
    p <-  
      plot_property_logo(
        results_KF_by_fact%>%
          mutate(imgt =change_imgt_format(IMGT_position))%>%select(cells,imgt, prop,ends_with("oo"), length)%>%select(-c(contains('hotelling')))%>%mutate(label_group_a_oo="DQ2",label_group_b_oo="DQ8")%>%filter(length==15),
        q_cutoff   = q_cutoff,
        d_cutoff   = d_cutoff,
        cols = cols,
        cell_types = cells_sel,
        color_a    = "#e41a1c",   # red = group A
        color_b    = "#377eb8"    # blue = group B
        ,positions = c("108","109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"),
        add_to_title =  paste(cell_fullnames[cells_sel],"cells, length=", len_sel))
    ggsave(plot = p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version, paste0("l",len_sel), cells_sel, sep="_"),".png"))
  }
}


################## 
######### Now compare hetero vs hom
################## 

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
dir.create(file.path(plot_dir, test_version ))

for(cells_sel in c("N","E")){
  
  p <- plot_aa_composition_logo( results_imgt_WL%>%filter(cells==cells_sel)%>%
                                   mutate(imgt =change_imgt_format(IMGT_position))%>%rename("aa"=AA)%>%select(imgt, aa,ends_with("eo"),ends_with("do")),
                                 min_freq =min_freq, use_seqlogo = TRUE,
                                  cols=cols, p_cutoff = p_cutoff, d_cutoff = d_cutoff,
                                 add_to_title = paste(cell_fullnames[cells_sel],"cells"))
  ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL", cells_sel, sep="_"), ".png"))
  
  
  ###one length plots
  for(len_sel in 13:16){
    p <- plot_aa_composition_logo( results_imgt%>%filter(cells==cells_sel, length==len_sel)%>%
                                     mutate(imgt =change_imgt_format(IMGT_position))%>%rename("aa"=AA)%>%select(imgt, aa,ends_with("eo"),ends_with("do")),
                                   min_freq =min_freq, use_seqlogo = TRUE,
                                   p_cutoff = p_cutoff, d_cutoff = d_cutoff, cols=cols,
                                   add_to_title =  paste(cell_fullnames[cells_sel],"cells, length=", len_sel))
    ggsave(plot = p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version, paste0("l",len_sel), cells_sel, sep="_"),".png"))
  }
}


test_version <- "prop_het_hom"
dir.create(file.path(plot_dir, test_version ))
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



for(cells_sel in c("N","E")){
  
  p <-  plot_property_logo(
    results_KF_WL_by_fact%>%
      mutate(imgt =change_imgt_format(IMGT_position))%>%select(cells, imgt, prop,ends_with("eo"))%>%select(-c(contains('hotelling')))%>%mutate(label_group_a_eo="DQ2DQ8",label_group_b_eo="homo"),
    q_cutoff   = q_cutoff,
    d_cutoff   = d_cutoff,
    cols = cols,
    cell_types = cells_sel,
    color_a    = "#e41a1c",   # red = group A
    color_b    = "#377eb8"    # blue = group B
    ,positions = c("108","109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"),
    add_to_title = paste(cell_fullnames[cells_sel],"cells"))
  ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL", cells_sel, sep="_"), ".png"))
  
  
  for(len_sel in 13:16){
    p <-  
      plot_property_logo(
        results_KF_by_fact%>%
          mutate(imgt =change_imgt_format(IMGT_position))%>%select(cells,imgt, prop,ends_with("eo"), length)%>%select(-c(contains('hotelling')))%>%mutate(label_group_a_eo="DQ2DQ8",label_group_b_eo="homo")%>%filter(length==15),
        q_cutoff   = q_cutoff,
        d_cutoff   = d_cutoff,
        cols = cols,
        cell_types = cells_sel,
        color_a    = "#e41a1c",   # red = group A
        color_b    = "#377eb8"    # blue = group B
        ,positions = c("108","109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"),
        add_to_title =  paste(cell_fullnames[cells_sel],"cells, length=", len_sel))
    ggsave(plot = p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version, paste0("l",len_sel), cells_sel, sep="_"),".png"))
  }
}


cols <- list(
  imgt     = "imgt",
  prop       = "prop",
  cells = "cells",
  mean_a   = "mean_b_po",
  mean_b   = "mean_b_do",
  mean_c = "mean_a_do",
  cohens_d = "cohens_d_eo",
  q_value = "q_value_eo",       
  cohens_d = "cohens_d_eo",
  cohens_d_do = "cohens_d_do",
  cohens_d_po = "cohens_d_po",
  q_value_do = "q_value_do",     
  q_value_po = "q_value_po",  
  label_a = "label_group_b_do",
  label_b = "label_group_b_po",
  label_c = "label_group_a_eo"        # Third group label
)


for(cells_sel in c("N","E")){
  
  p <-  plot_property_logo(
  results_KF_WL_by_fact%>%
    mutate(imgt =change_imgt_format(IMGT_position))%>%select(cells, imgt, prop,ends_with("eo"),ends_with("do"),ends_with("po"))%>%select(-c(contains('hotelling')))%>%
    mutate(label_group_a_eo="DQ2DQ8",label_group_b_eo="homo",label_group_b_po="DQ2",label_group_b_do="DQ8"),
  positions = c("109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"),
  q_cutoff   = q_cutoff,
  d_cutoff   = d_cutoff,
  cols = cols,
  cell_types = cells_sel,
  color_b    = "#e41a1c",   
  color_a    = "#377eb8" ,  
  color_c = "purple",
  add_to_title = paste(cell_fullnames[cells_sel],"cells"))
  ggsave(plot=p, path = file.path(plot_dir, file.path("cdr3", test_version) ), device = "png",width = 7.25, height = 6.1, filename = paste0(paste(test_version,"WL", cells_sel, sep="_"), ".png"))
  
  
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









################## 
######### Now check whether hetero is just DQ8 or DQ2...
################## 
cutoff=0.01

results_imgt_WL%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%rename("aa"=AA)%>%
  filter(p_value_eo<cutoff, p_value_do<cutoff, p_value_po< cutoff)%>%
  filter(imgt%in%c("109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"))%>%
  select(c(cells,imgt, aa,starts_with("p_value"),starts_with("cohens_d")), mean_dq2=mean_group_a_oo, mean_dq8= mean_group_b_oo, mean_het=mean_group_a_po)%>%
  kable(digits=c(rep(0,3),rep(4,4),rep(1,4),rep(2,3)))

results_imgt%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%rename("aa"=AA)%>%
  filter(p_value_eo<cutoff, p_value_do<cutoff, p_value_po< cutoff)%>%
  filter(imgt%in%c("109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"))%>%
  select(c(cells,length, imgt, aa,starts_with("p_value"),starts_with("cohens_d")), mean_dq2=mean_group_a_oo, mean_dq8= mean_group_b_oo, mean_het=mean_group_a_po)%>%
  arrange(cells, length, imgt)%>%
  rename_with(~ .x |> gsub("_oo$", " dq2 vs dq8", x = _) |>gsub("_eo", " het vs hom", x = _) |> gsub("_po$", " het vs dq2", x = _) |>gsub("_do", " het vs dq8", x = _)
  )%>%
  kable(digits=c(rep(0,4),2,rep(4,3),rep(1,4),rep(3,3)))


results_KF_WL_by_fact%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_eo<cutoff, p_value_do<cutoff, p_value_po< cutoff)%>%
  filter(imgt%in%c("109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"))%>%
  select(c(cells,imgt, prop, starts_with("p_value"),starts_with("cohens_d")), mean_dq2=mean_a_oo, mean_dq8= mean_b_oo, mean_het=mean_a_po)%>%

  kable(digits=c(0,0,0,2,rep(4,3),rep(1,4),rep(2,3)))


results_KF_by_fact%>%
  mutate(imgt =change_imgt_format(IMGT_position))%>%
  filter(p_value_eo<cutoff, p_value_do<cutoff, p_value_po< cutoff)%>%
  filter(imgt%in%c("109","110","111", "111.1","111.2", "112", "112.1", "112.2", "113"))%>%
  select(c(cells,length, imgt, prop, starts_with("p_value"),starts_with("cohens_d"), starts_with("q_val")), mean_dq2=mean_a_oo, mean_dq8= mean_b_oo, mean_het=mean_a_po)%>%
  rename_with(~ .x |> gsub("_oo$", " dq2 vs dq8", x = _) |>gsub("_eo", " het vs hom", x = _) |> gsub("_po$", " het vs dq2", x = _) |>gsub("_do", " het vs dq8", x = _)
  )%>%gt() %>%
  fmt_number(
    columns = starts_with("mean"),
    decimals = 3
  ) %>%
  fmt_number(
    columns = starts_with("cohens_d"),
    decimals = 1
  ) %>%
  fmt_number(
    columns = starts_with("p_value"),
    decimals = 4
  ) %>%
  fmt_number(
    columns = starts_with("q_value"),
    decimals = 2
  )

