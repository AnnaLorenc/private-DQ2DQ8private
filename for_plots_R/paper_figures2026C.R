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
