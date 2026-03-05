
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

