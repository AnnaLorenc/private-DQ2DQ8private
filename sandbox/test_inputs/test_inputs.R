library(tidyverse)

input <- read_tsv("sandbox/test_inputs/imgt_freqs/F5302E_subsampled_test.tsv" )

#

cdr3_to_imgt <- function(cdr3) {
  n <- nchar(cdr3)
  aa <- strsplit(cdr3, "", fixed = TRUE)[[1]]
  # IMGT: 104 (C), 105-111 (core), 111A/B... (insertions), 112-117 (core), 118 (F, if present)
  canonical_core <- 13
  pos <- c("104")
  core_len <- n - 2
  left <- as.character(105:111)
  right <- as.character(112:117)
  if (core_len <= canonical_core) {
    # No insertions
    n_right <- core_len - 7
    if (n_right > 0) {
      pos <- c(pos, left, right[1:n_right])
    } else {
      pos <- c(pos, left[1:core_len])
    }
    pos <- c(pos, "117")
  } else {
    # Insertions after 111: 111A, 111B, ...
    n_insert <- core_len - canonical_core
    insertions <- paste0("111", LETTERS[1:n_insert])
    pos <- c(pos, left, insertions, right)
    pos <- pos[1:(n-1)]
    pos <- c(pos, ifelse(n >= 15, "118", "117"))
  }
  tibble::tibble(aa = aa, pos = pos)
}

cdr3_vec <- c("CAAAIGGAQDTQYF", "CASSLGTDTQYF", "CASSLGGGSYNEQFF", "CASSLGGGGSYNEQFF")

input_imgt <- input %>%
  mutate(rowid = row_number()) %>%
  left_join(
    map_dfr(input$aminoAcid, ~cdr3_to_imgt(.x) %>% mutate(cdr3 = .x), .id = "rowid") %>%
      mutate(rowid = as.integer(rowid)),
    by = c( "rowid", "aminoAcid"="cdr3")
  )%>%
  mutate(aa_len=nchar(aminoAcid))

input_imgt %>%
  pivot_longer(cols=starts_with("F5302"), names_to = "sample", values_to = "value")%>%
  group_by( aa,    pos,   sample )%>%
  summarise(counts_uniq_TRBs =n_distinct(aminoAcid, vFamilyName) )%>%
  ungroup()%>%
  pivot_wider(names_from=sample, values_from=counts_uniq_TRBs)%>%
  group_by( aa, pos)%>%
  summarise(across(starts_with("F5302E"), ~.x/sum(.x)))


####===
raw <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/subsampled/productive/F5302E_subsampled.tsv.gz")
F5302E_freq_med.tsv <-read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq_med.tsv")
F5302E_freq_WL_med.tsv <-read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq_WL_med.tsv")
F5302E_freq.tsv<- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq.tsv")
F5302E_valid <- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_valid.tsv")
F5302E_freq_Vfam<- read_tsv("/Users/ania/Documents/DQ2DQ8/pipeline/DQ2DQ8/results/imgt_aa/F5302E_aa_freqs/F5302E_freq_Vfam.tsv")

