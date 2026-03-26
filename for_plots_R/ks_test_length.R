library(dplyr)
library(readr)
library(ggplot2)
library(purrr)

#-------------------------------
# 1. Load data
#-------------------------------
counts <- read_tsv("cdr3_lengths.tsv")
meta   <- read_tsv("sample_groups.tsv")

df <- counts %>%
  left_join(meta, by = "sample")

#-------------------------------
# 2. Normalize per sample
# (convert counts → probabilities)
#-------------------------------
df_norm <- df %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

#-------------------------------
# 3. Define common AA length grid
#-------------------------------
aa_grid <- sort(unique(df$aa_len))

#-------------------------------
# 4. Compute ECDF per sample (fast)
#-------------------------------
compute_ecdf <- function(subdf) {
  vec <- setNames(subdf$freq, subdf$aa_len)
  vec_full <- vec[as.character(aa_grid)]
  vec_full[is.na(vec_full)] <- 0
  cumsum(vec_full)
}

sample_ecdfs <- df_norm %>%
  group_by(sample, group) %>%
  summarise(ecdf = list(compute_ecdf(cur_data())), .groups = "drop")

#-------------------------------
# 5. Function to compute KS between groups
#-------------------------------
compute_group_ks <- function(sample_ecdfs) {
  groups <- unique(sample_ecdfs$group)
  
  ecdf1 <- sample_ecdfs %>%
    filter(group == groups[1]) %>%
    pull(ecdf) %>%
    simplify2array() %>%
    rowMeans()
  
  ecdf2 <- sample_ecdfs %>%
    filter(group == groups[2]) %>%
    pull(ecdf) %>%
    simplify2array() %>%
    rowMeans()
  
  max(abs(ecdf1 - ecdf2))
}

# observed KS
D_obs <- compute_group_ks(sample_ecdfs)

#-------------------------------
# 6. Fast permutation test (sample-level)
#-------------------------------
ks_permutation <- function(sample_ecdfs, B = 2000) {
  groups <- sample_ecdfs$group
  D_perm <- numeric(B)
  
  for(i in seq_len(B)) {
    perm_groups <- sample(groups)
    
    perm_df <- sample_ecdfs
    perm_df$group <- perm_groups
    
    D_perm[i] <- compute_group_ks(perm_df)
  }
  
  p_value <- mean(D_perm >= D_obs)
  
  list(D_obs = D_obs, p_value = p_value, D_perm = D_perm)
}

res <- ks_permutation(sample_ecdfs, B = 2000)

cat("Observed KS:", res$D_obs, "\n")
cat("Permutation p-value:", res$p_value, "\n")

#-------------------------------
# 7. Prepare ECDF plot data
#-------------------------------
plot_df <- sample_ecdfs %>%
  mutate(ecdf_mat = map(ecdf, ~ .x)) %>%
  unnest_wider(ecdf_mat) %>%
  pivot_longer(
    cols = starts_with("..."),
    names_to = "idx",
    values_to = "ecdf"
  ) %>%
  mutate(
    aa_len = aa_grid[as.numeric(gsub("\\D", "", idx))]
  )

# group average ECDF
group_ecdf <- plot_df %>%
  group_by(group, aa_len) %>%
  summarise(ecdf = mean(ecdf), .groups = "drop")

#-------------------------------
# 8. Plot ECDFs + KS distance
#-------------------------------
g <- ggplot(group_ecdf, aes(x = aa_len, y = ecdf, color = group)) +
  geom_step(size = 1.2) +
  theme_classic() +
  labs(
    x = "CDR3 length",
    y = "ECDF",
    title = paste0("KS = ", round(res$D_obs, 3),
                   ", p = ", signif(res$p_value, 3))
  )

# mark KS max difference
diff_df <- group_ecdf %>%
  pivot_wider(names_from = group, values_from = ecdf)

ks_idx <- which.max(abs(diff_df[[2]] - diff_df[[3]]))

g <- g +
  geom_segment(
    data = diff_df[ks_idx, ],
    aes(
      x = aa_len,
      xend = aa_len,
      y = diff_df[[2]][ks_idx],
      yend = diff_df[[3]][ks_idx]
    ),
    inherit.aes = FALSE,
    linetype = "dashed"
  )

print(g)