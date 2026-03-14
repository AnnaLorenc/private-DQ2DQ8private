library(tidyverse)
library(kableExtra)
library(scales)
library(ggpubr)


###pairwise sharing

annotation_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/data/collated_info.csv"

pairwise_loc <- "/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/overlaps_subs/overlaps_subs_productive.tsv"

colors=c(DQ2="red",
         DQ2DQ8="purple",
         DQ8="blue")


pairwise <- read_tsv(pairwise_loc)
anno <- read_csv(annotation_loc)%>%
  mutate(source=substr(newname,1,1), genotype_short=gsub(pat="heteroD|homoD", rep="D", genotype_short))%>%
  rename(patient=shortname)%>%
  unique

df <- pairwise%>%group_by(short_id, comparison_type,sample1,sample2)%>%
  summarise(jaccard=mean(jaccard))

###within genotype

plot_data <- df%>%
  filter(comparison_type=="within_N")%>%
  ungroup%>%
  left_join(., anno%>%select(sample_short, source1=source, geno_s1=genotype_short), by=c("sample1"="sample_short"))%>%
  left_join(., anno%>%select(sample_short, source2=source, geno_s2=genotype_short), by=c("sample2"="sample_short"))%>%
  mutate(geno_comp = paste0(pmin(geno_s1, geno_s2), "-", pmax(geno_s1, geno_s2))) %>%
  # Define the specific order as a factor
  mutate(geno_comp = factor(geno_comp, levels = c("DQ2-DQ2", "DQ2DQ8-DQ2DQ8", "DQ8-DQ8"))) %>%
  filter(!is.na(geno_comp))

my_comparisons <- list( 
  c("DQ2-DQ2", "DQ2DQ8-DQ2DQ8"), 
  c("DQ2DQ8-DQ2DQ8", "DQ8-DQ8"), 
  c("DQ2-DQ2", "DQ8-DQ8") 
)
  
p <- ggplot(plot_data, aes(x = geno_comp, y = jaccard, fill = geno_comp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(
    title = "Jaccard Index, naive cells",
    subtitle = "Order: DQ2, DQ2DQ8, DQ8",
    x = "Genotype Comparison",
    y = "Jaccard Index",
    fill = "Comparison"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_fill_manual(values=unname(colors))

p+ stat_compare_means(
  comparisons = my_comparisons,
  method = "wilcox.test",  # Non-parametric (recommended for Jaccard/overlaps)
  label = "p.format"      # Use "p.format" for actual numbers, "p.signif" for ****
) +
  # Add an overall ANOVA/Kruskal-Wallis test at the top
  stat_compare_means(label.y = max(plot_data$jaccard) * 1.5) 
  

###between genotypes
plot_data <- df%>%
  filter(comparison_type=="within_N")%>%
  ungroup%>%
  left_join(., anno%>%select(sample_short, source1=source, geno_s1=genotype_short), by=c("sample1"="sample_short"))%>%
  left_join(., anno%>%select(sample_short, source2=source, geno_s2=genotype_short), by=c("sample2"="sample_short"))%>%
  mutate(geno_comp = paste0(pmin(geno_s1, geno_s2), "-", pmax(geno_s1, geno_s2))) %>%
  # Define the specific order as a factor
  mutate(geno_comp = factor(geno_comp, levels = c("DQ2-DQ8", "DQ2-DQ2DQ8", "DQ2DQ8-DQ8"))) %>%
  # Filter to only keep the three combinations you requested
  filter(!is.na(geno_comp))

my_comparisons <- list( 
  c("DQ2-DQ8", "DQ2-DQ2DQ8"), 
  c("DQ2-DQ8", "DQ2DQ8-DQ8"), 
  c("DQ2DQ8-DQ8", "DQ2-DQ2DQ8") 
)

p <- ggplot(plot_data, aes(x = geno_comp, y = jaccard, fill = geno_comp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(
    title = "Jaccard Index, naive cells",
    x = "Genotype Comparison",
    y = "Jaccard Index",
    fill = "Comparison"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_fill_brewer(palette = "Set2")

p + stat_compare_means(
  comparisons = my_comparisons,
  method = "wilcox.test",  # Non-parametric (recommended for Jaccard/overlaps)
  label = "p.format"      # Use "p.format" for actual numbers, "p.signif" for ****
) +
  # Add an overall ANOVA/Kruskal-Wallis test at the top
  stat_compare_means(label.y = max(plot_data$jaccard) * 1.5) 



########Memory cells
###within genotype

plot_data <- df%>%
  filter(comparison_type=="within_E")%>%
  ungroup%>%
  left_join(., anno%>%select(sample_short, source1=source, geno_s1=genotype_short), by=c("sample1"="sample_short"))%>%
  left_join(., anno%>%select(sample_short, source2=source, geno_s2=genotype_short), by=c("sample2"="sample_short"))%>%
  mutate(geno_comp = paste0(pmin(geno_s1, geno_s2), "-", pmax(geno_s1, geno_s2))) %>%
  # Define the specific order as a factor
  mutate(geno_comp = factor(geno_comp, levels = c("DQ2-DQ2", "DQ2DQ8-DQ2DQ8", "DQ8-DQ8"))) %>%
  filter(!is.na(geno_comp))

my_comparisons <- list( 
  c("DQ2-DQ2", "DQ2DQ8-DQ2DQ8"), 
  c("DQ2DQ8-DQ2DQ8", "DQ8-DQ8"), 
  c("DQ2-DQ2", "DQ8-DQ8") 
)

p <- ggplot(plot_data, aes(x = geno_comp, y = jaccard, fill = geno_comp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(
    title = "Jaccard Index, memory cells",
    subtitle = "Order: DQ2, DQ2DQ8, DQ8",
    x = "Genotype Comparison",
    y = "Jaccard Index",
    fill = "Comparison"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_fill_manual(values=unname(colors))

p+ stat_compare_means(
  comparisons = my_comparisons,
  method = "wilcox.test",  # Non-parametric (recommended for Jaccard/overlaps)
  label = "p.format"      # Use "p.format" for actual numbers, "p.signif" for ****
) +
  # Add an overall ANOVA/Kruskal-Wallis test at the top
  stat_compare_means(label.y = max(plot_data$jaccard) * 1.5) 

###between genotypes
plot_data <- df%>%
  filter(comparison_type=="within_E")%>%
  ungroup%>%
  left_join(., anno%>%select(sample_short, source1=source, geno_s1=genotype_short), by=c("sample1"="sample_short"))%>%
  left_join(., anno%>%select(sample_short, source2=source, geno_s2=genotype_short), by=c("sample2"="sample_short"))%>%
  mutate(geno_comp = paste0(pmin(geno_s1, geno_s2), "-", pmax(geno_s1, geno_s2))) %>%
  # Define the specific order as a factor
  mutate(geno_comp = factor(geno_comp, levels = c("DQ2-DQ8", "DQ2-DQ2DQ8", "DQ2DQ8-DQ8"))) %>%
  # Filter to only keep the three combinations you requested
  filter(!is.na(geno_comp))


my_comparisons <- list( 
  c("DQ2-DQ8", "DQ2-DQ2DQ8"), 
  c("DQ2-DQ8", "DQ2DQ8-DQ8"), 
  c("DQ2DQ8-DQ8", "DQ2-DQ2DQ8") 
)

p <- ggplot(plot_data, aes(x = geno_comp, y = jaccard, fill = geno_comp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(
    title = "Jaccard Index, memory cells",
    x = "Genotype Comparison",
    y = "Jaccard Index",
    fill = "Comparison"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  scale_fill_brewer(palette = "Set2")

p+ stat_compare_means(
  comparisons = my_comparisons,
  method = "wilcox.test",  # Non-parametric (recommended for Jaccard/overlaps)
  label = "p.format"      # Use "p.format" for actual numbers, "p.signif" for ****
) +
  # Add an overall ANOVA/Kruskal-Wallis test at the top
  stat_compare_means(label.y = max(plot_data$jaccard) * 1.5) 




################
##
## PUBLIC SEQUENCES
##
################

##Improve the main title part
##add N of clonotypes in the 'leading' column
##overpotting issue???


plot_sharing <- function(sharing_summary, cfs, min_frac, clonotype_cols=c(  "aminoAcid" ,  "vFamilyName" ,"jGeneName"), for_title=""){
  max_sharing <- colnames(sharing_summary) %>%
    grep(pat="^X", val=T) %>%
    gsub(pat='X', rep='') %>%
    as.integer(.)%>%
    max()
  max_shared_within_geno <- sapply(c("DQ8", "DQ2", "hete"), function(grp)
    sharing_summary%>%
      filter(.data[[paste0("X", cfs)]] > 0, group == grp)%>%
      select(-c(paste0('X',c(3:max_sharing)[3:max_sharing <cfs])))%>%
      pivot_longer(cols=starts_with('X'), names_to="shared_by")%>%
      filter(value>= min_frac)%>%
      mutate(shared=as.integer(gsub(pat='X', rep='',shared_by)))%>%
      unite("clonotype", all_of(clonotype_cols), sep = "_", remove = FALSE)%>%
      group_by(clonotype)%>%
      slice_max(shared), USE.NAMES = TRUE, simplify=FALSE
  )
  
  top_hits <- sapply(max_shared_within_geno, function(mswg)
    left_join(mswg%>%
                select(clonotype),
              sharing_summary%>%pivot_longer(cols=starts_with('X'), names_to="shared_by")%>%
                filter(value>= min_frac)%>%
                mutate(shared=as.integer(gsub(pat='X', rep='',shared_by)))%>%
                select(-c(shared_by, value))%>%
                unite("clonotype", all_of(clonotype_cols), sep = "_", remove = FALSE)%>%
                group_by(across(all_of(c("group", "clonotype"))))%>%
                slice_max(shared),
              by ="clonotype")%>%
      pivot_wider( names_from = group,values_from = shared, values_fn=max,values_fill = 0),
    USE.NAMES = TRUE, simplify=FALSE)
  
  #now visualise top hits
  plot_one_panel <- function(panel_df, name,clonotype_cols){
    p <- panel_df%>%
      pivot_longer(c(DQ2,hete,DQ8), names_to = "geno")%>%
      mutate(geno=factor(geno, levels=c("DQ2", "hete", "DQ8")))%>%
      ggplot(aes(x=geno, y=value, group=clonotype, col=clonotype)) +geom_line() +theme_bw()+
      scale_y_discrete()+ylim(0, max_sharing)+
      guides(color="none")+
      ggtitle(name)
    p
  }
  library(patchwork)
  
  all3plots <- lapply(c("DQ2", "hete", "DQ8"), function(i)plot_one_panel (top_hits[[i]], i,clonotype_cols=clonotype_cols))
  
  p <- wrap_plots(all3plots, ncol=3)
  p+ plot_annotation(title = paste("clonotype by", paste0(sapply(clonotype_cols, substr,start=1, stop=4),collapse="_"),"min.frac=",min_frac, " cfs=",cfs," ", for_title))
}



sharing_summary <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/rarified_sharing/aminoAcid_vFamilyName_jGeneName/sharing_summary_N_A3/merged_summary.tsv")
names(sharing_summary) <- make.names(names(sharing_summary))
plot_sharing (sharing_summary, cfs=7, min_frac=20, clonotype_cols=c(  "aminoAcid" ,  "vFamilyName" ,"jGeneName"), for_title=" naive cells")

sharing_summary <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/rarified_sharing/aminoAcid_vFamilyName/sharing_summary_N_A3/merged_summary.tsv")
names(sharing_summary) <- make.names(names(sharing_summary))

plot_sharing (sharing_summary, cfs=7, min_frac=20, clonotype_cols=c(  "aminoAcid" ,  "vFamilyName"), for_title=" naive cells")


##something wrong when just aminoacid column
sharing_summary <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/rarified_sharing/aminoAcid/sharing_summary_N_A3/merged_summary.tsv")
names(sharing_summary) <- make.names(names(sharing_summary))

plot_sharing (sharing_summary, cfs=7, min_frac=20, clonotype_cols=c(  "aminoAcid" ), for_title=" naive cells")



sharing_summary <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/rarified_sharing/aminoAcid_vFamilyName_jGeneName/sharing_summary_E_A3/merged_summary.tsv")
names(sharing_summary) <- make.names(names(sharing_summary))
plot_sharing (sharing_summary, cfs=7, min_frac=20, clonotype_cols=c(  "aminoAcid" ,  "vFamilyName" ,"jGeneName"), for_title=" memory cells")

sharing_summary <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/rarified_sharing/aminoAcid_vFamilyName/sharing_summary_E_A3/merged_summary.tsv")
names(sharing_summary) <- make.names(names(sharing_summary))

plot_sharing (sharing_summary, cfs=7, min_frac=20, clonotype_cols=c(  "aminoAcid" ,  "vFamilyName"), for_title=" memory cells")


##something wrong when just aminoacid column
sharing_summary <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/rarified_sharing/aminoAcid/sharing_summary_E_A3/merged_summary.tsv")
names(sharing_summary) <- make.names(names(sharing_summary))

plot_sharing (sharing_summary, cfs=7, min_frac=10, clonotype_cols=c(  "aminoAcid" ), for_title=" memory cells")
