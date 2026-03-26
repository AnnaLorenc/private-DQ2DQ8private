library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

format_pval <- function(p, digits = 2) {
  sapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x == 0) return(paste0("<", formatC(.Machine$double.xmin, format = "e", digits = digits)))
    
    # use significant digits
    out <- signif(x, digits)
    
    # format nicely (avoid unnecessary trailing zeros)
    if (out < 0.001) {
      formatC(out, format = "e", digits = digits)
    } else {
      sub("0+$", "", sub("\\.$", "", format(out, scientific = FALSE)))
    }
  })
}
# ---------------------------------------------------------------------------
# plot_aa_preferences_logo() 
#
# Sequence-logo-style visualisation of amino acid preferences across IMGT
# positions, based on Cohen's D from pairwise group comparisons.
#
#   Letters ABOVE x-axis  →  positive Cohen's D (enriched in group A)
#   Letters BELOW x-axis  →  negative Cohen's D (enriched in group B)
#   Letter height         ∝  |Cohen's D|  (stacked within each position)
#   Letter colour         =  aa_colors mapping
#
# Default column names in df (override any via the cols argument):
#   imgt             (cols$imgt)     - IMGT position label
#   aa               (cols$aa)       - single-letter amino acid code
#   cohens_d_oo      (cols$cohens_d) - Cohen's D for the comparison
#   p_value_oo       (cols$p_value)  - p-value for the comparison
#   label_group_a_oo (cols$label_a)  - group A name (axis annotation)
#   label_group_b_oo (cols$label_b)  - group B name (axis annotation)
#
# Args:
#   df         : data frame as described above
#   p_cutoff   : maximum p-value to include              (default 0.05)
#   d_cutoff   : minimum |Cohen's D| to include          (default 0.5)
#   positions  : character vector of imgt values to show; NULL = all present
#   aa_colors  : named character vector mapping AA letters to hex colours;
#                NULL uses a built-in chemical-property colour scheme
#   cols       : named list mapping canonical internal names to actual column
#                names in df; partial lists are merged with defaults
#
# Returns: a ggplot object
# ---------------------------------------------------------------------------

# canonical full IMGT CDR3 position order (104–118)
.imgt_order <- c(
  as.character(104:111),
  paste0("111.", 1:5),   # 111A, 111B, 111C, 111D  (insertions towards centre)
  paste0("112.", 5:1),   # 112E, 112D, 112C, 112B, 112A  (insertions, reverse)
  as.character(112:118)
)

# default chemical-property colour scheme (20 standard AAs)
.default_aa_colors <- c(
  # aliphatic / nonpolar
  A = "#80b1d3", V = "#80b1d3", I = "#80b1d3", L = "#80b1d3", M = "#80b1d3",
  # aromatic
  F = "#bc80bd", W = "#bc80bd", Y = "#fdb462",
  # small / special
  G = "#b3de69", P = "#b3de69",
  # polar uncharged
  S = "#fb8072", T = "#fb8072", C = "#fdb462", N = "#bebada", Q = "#bebada",
  # charged positive
  K = "#1f78b4", R = "#1f78b4", H = "#8dd3c7",
  # charged negative
  D = "#ff7f00", E = "#ff7f00"
)

# ---------------------------------------------------------------------------
# Internal helper – rename data-frame columns from user-supplied names to
# the canonical names expected internally by each plot function.
#
#   col_map  named list; keys   = canonical (internal) column name,
#                         values = actual column name in the user's df.
#            Key-value pairs where key == value are silently skipped.
# ---------------------------------------------------------------------------
.rename_cols <- function(df, col_map) {
  # col_map: key = internal code name (target), value = user column name (source)
  
  # 1. Filter the map to only include source columns that actually exist in df
  #    and where the source name is different from the target name.
  valid_map <- Filter(
    function(key) {
      src <- col_map[[key]]
      !is.null(src) && length(src) == 1L && !is.na(src) &&
        src != key && src %in% names(df)
    },
    names(col_map)
  )
  
  if (length(valid_map) == 0L) return(df)

  # 2. Prepare the final mapping.
  #    We must handle the case where the user specifies a source column 
  #    that is ALREADY the target name for a different key, 
  #    or if the target name exists but is a different column.
  
  # To avoid "Names must be unique" in dplyr::rename, we first drop 
  # ANY existing columns in df that match our TARGET names, 
  # UNLESS that target name is actually the source for one of our renames.
  
  all_targets <- names(valid_map)
  all_sources <- unname(vapply(valid_map, function(k) col_map[[k]], character(1L)))
  
  # Columns to drop: any target name that exists in df but is NOT being used as a source
  targets_to_drop <- all_targets[all_targets %in% names(df) & !(all_targets %in% all_sources)]
  
  if (length(targets_to_drop) > 0) {
    df <- df |> dplyr::select(-dplyr::any_of(targets_to_drop))
  }

  # 3. Perform renaming. 
  #    We use a named vector for rename: c(new_name = old_name)
  rename_vec <- vapply(valid_map, function(k) col_map[[k]], character(1L))
  
  # Final safety check: ensure the rename_vec doesn't have duplicate sources or targets
  # (valid_map keys are unique by definition of named lists)
  dplyr::rename(df, !!rename_vec)
}




# ---------------------------------------------------------------------------
# plot_aa_composition_logo()
#
# Frequency-based sequence-logo plot showing ALL amino acids at each IMGT
# position. Amino acids passing significance thresholds are shown in full
# colour; all others are shown in light grey.
#
#   Letters ABOVE x-axis  →  group A frequencies  (mean_group_a_oo)
#   Letters BELOW x-axis  →  group B frequencies  (mean_group_b_oo)
#   Letter height         ∝  frequency (stacked, significant closest to axis)
#   Letter colour         =  aa_colors if significant, grey otherwise
#   x-axis                =  IMGT position labels
#
# Significance is determined per (imgt, aa) pair using p_cutoff and d_cutoff.
#
# use_seqlogo = TRUE requires:
#   install.packages(c("ggseqlogo", "patchwork"))
#   In this mode, two separate panels are shown (group A above, group B
#   below) with upright letters in both panels (no flipping).
#   Note: ggseqlogo colour schemes are global per-AA; significance colouring
#   is therefore applied at the AA level (an AA coloured if it is significant
#   at *any* shown position).
#
# Default column names in df (override any via the cols argument):
#   imgt             (cols$imgt)     - IMGT position
#   aa               (cols$aa)       - amino acid
#   mean_group_a_oo  (cols$mean_a)   - mean frequency in group A
#   mean_group_b_oo  (cols$mean_b)   - mean frequency in group B
#   cohens_d_oo      (cols$cohens_d) - Cohen's D
#   p_value_oo       (cols$p_value)  - p-value
#   label_group_a_oo (cols$label_a)  - group A label
#   label_group_b_oo (cols$label_b)  - group B label
#
# Args:
#   df           : data frame as described above
#   p_cutoff     : max p-value to consider significant     (default 0.05)
#   d_cutoff     : min |Cohen's D| to consider significant (default 0.5)
#   positions    : character vector of imgt values to show; NULL = all present
#   aa_colors    : named character vector mapping AA letters to hex colours;
#                  NULL uses the built-in chemical-property scheme
#   min_freq     : suppress letters below this frequency   (default 0.005)
#   use_seqlogo  : if TRUE use ggseqlogo + patchwork for true letter shapes
#   cols         : named list mapping canonical internal names to actual
#                  column names in df; partial lists merged with defaults
#
# Returns: a ggplot / patchwork object
# ---------------------------------------------------------------------------

plot_aa_composition_logo <- function(
    df,
    p_cutoff    = 0.05,
    d_cutoff    = 0.5,
    positions   = NULL,
    aa_colors   = NULL,
    min_freq    = 0.005,
    use_seqlogo = FALSE,
    cols = list(),
    add_to_title = "",
    scale_log = FALSE
) {
    library(ggplot2)
    library(dplyr)
    library(tidyr)

    if (is.null(aa_colors)) aa_colors <- .default_aa_colors
    if(add_to_title!=""){add_to_title <- paste0("; ",add_to_title)}

    # ---- resolve column names (require all from cols) ----
    required_cols <- c("imgt", "aa", "mean_a", "mean_b", "cohens_d", "p_value", "label_a", "label_b")
    missing <- setdiff(required_cols, names(cols))
    if (length(missing) > 0) {
      stop(paste("The following required mapping(s) are missing from 'cols':", 
                 paste(missing, collapse = ", ")))
    }

    # Identify additional p_value and cohens_d columns
    extra_pval_cols <- setdiff(names(cols)[grepl("^p_value", names(cols))], "p_value")
    extra_d_cols    <- setdiff(names(cols)[grepl("^cohens_d", names(cols))], "cohens_d")

    # Build rename mapping for all columns
    rename_map <- list(
      imgt             = cols$imgt,
      aa               = cols$aa,
      mean_group_a_oo  = cols$mean_a,
      mean_group_b_oo  = cols$mean_b,
      cohens_d_oo      = cols$cohens_d,
      p_value_oo       = cols$p_value,
      label_group_a_oo = cols$label_a,
      label_group_b_oo = cols$label_b
    )
    # Add extra columns to rename mapping
    for (col in extra_pval_cols) {
      rename_map[[paste0(col, "_oo")]] <- cols[[col]]
    }
    for (col in extra_d_cols) {
      rename_map[[paste0(col, "_oo")]] <- cols[[col]]
    }

    df <- .rename_cols(df, rename_map)

  # ---- filter positions ----
  if (!is.null(positions)) df <- df |> filter(imgt %in% positions)
  if (nrow(df) == 0) {
    message("No data for the specified positions.")
    return(ggplot() + theme_void())
  }

  # ---- canonical position order ----
  present_pos <- unique(df$imgt)
  ordered_pos <- .imgt_order[.imgt_order %in% present_pos]
  extra_pos   <- setdiff(present_pos, ordered_pos)
  ordered_pos <- c(ordered_pos, extra_pos)

  # ---- significant (imgt, aa) pairs ----
  # Build list of all p_value/cohens_d column pairs
  pval_cols <- c("p_value_oo", paste0(extra_pval_cols, "_oo"))
  d_cols    <- c("cohens_d_oo", paste0(extra_d_cols, "_oo"))

  # Require ALL filtering criteria to be fulfilled simultaneously
  filter_expr <- rep(TRUE, nrow(df))
  for (i in seq_along(pval_cols)) {
    pcol <- pval_cols[i]
    dcol <- d_cols[i]
    if (pcol %in% colnames(df) && dcol %in% colnames(df)) {
      filter_expr <- filter_expr & (df[[pcol]] <= p_cutoff & abs(df[[dcol]]) >= d_cutoff)
    }
  }
  sig_pairs <- df[filter_expr, c("imgt", "aa")] |>
    distinct() |>
    mutate(is_sig = TRUE)

  # ---- aggregate across cells if multiple cell types present ----
  df_agg <- df |>
    group_by(imgt, aa, label_group_a_oo, label_group_b_oo) |>
    summarise(
      freq_a = mean(mean_group_a_oo, na.rm = TRUE),
      freq_b = mean(mean_group_b_oo, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(sig_pairs, by = c("imgt", "aa")) |>
    replace_na(list(is_sig = FALSE)) |>
    mutate(imgt = factor(imgt, levels = ordered_pos))

  group_a <- unique(df$label_group_a_oo)[1]
  group_b <- unique(df$label_group_b_oo)[1]

  # ==========================================================================
  # Branch A: ggseqlogo + patchwork with selective position-specific coloring
  # Requires: ggseqlogo, patchwork, geom_logo_selective.R
  # Significance colouring: Position-specific (coloured only where significant)
  # Two upright panels: group A (top) and group B (bottom), no letter flipping
  # ==========================================================================
  if (use_seqlogo) {
    if (!requireNamespace("ggseqlogo", quietly = TRUE))
      stop("Package 'ggseqlogo' is required. Install with: install.packages('ggseqlogo')")
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
    
    # Source the selective logo function if not already loaded
    if (!exists("ggseqlogo_selective")) {
      source("geom_logo_selective.R")
    }
    
    library(ggseqlogo)
    library(patchwork)

    aa_order <- c("A","C","D","E","F","G","H","I","K","L",
                  "M","N","P","Q","R","S","T","V","W","Y")

    # Build position frequency matrices
    build_pfm <- function(freq_col) {
      mat <- matrix(0, nrow = length(aa_order), ncol = length(ordered_pos),
                    dimnames = list(aa_order, as.character(ordered_pos)))
      for (i in seq_len(nrow(df_agg))) {
        aa_i  <- df_agg$aa[i]
        pos_i <- as.character(df_agg$imgt[i])
        val   <- df_agg[[freq_col]][i]
        if (aa_i %in% aa_order && pos_i %in% colnames(mat) && val >= min_freq)
          mat[aa_i, pos_i] <- val
      }
      cs2 <- colSums(mat)
      mat[, cs2 > 0] <- sweep(mat[, cs2 > 0, drop = FALSE], 2, cs2[cs2 > 0], "/")
      mat
    }

    pfm_a <- build_pfm("freq_a")
    pfm_b <- build_pfm("freq_b")
    
    # Create color dataframes for selective coloring
    # For group A: color significant AAs where group A is dominant
    color_df_a <- sig_pairs %>%
      left_join(df_agg %>% select(imgt, aa, freq_a, freq_b), by = c("imgt", "aa")) %>%
      filter(freq_a >= freq_b) %>%  # Only where group A is dominant
      select(position = imgt, aa) %>%
      mutate(position = as.character(position))
    
    # For group B: color significant AAs where group B is dominant
    color_df_b <- sig_pairs %>%
      left_join(df_agg %>% select(imgt, aa, freq_a, freq_b), by = c("imgt", "aa")) %>%
      filter(freq_b > freq_a) %>%  # Only where group B is dominant
      select(position = imgt, aa) %>%
      mutate(position = as.character(position))
    
    # Handle case where no significant results exist
    if (nrow(color_df_a) == 0) {
      color_df_a <- data.frame(position = character(0), aa = character(0), stringsAsFactors = FALSE)
    }
    if (nrow(color_df_b) == 0) {
      color_df_b <- data.frame(position = character(0), aa = character(0), stringsAsFactors = FALSE)
    }

    # Shared x-axis: IMGT labels replace ggseqlogo's default integers
    imgt_x_scale <- scale_x_continuous(
      breaks = seq_along(ordered_pos),
      labels = ordered_pos
    )

    theme_top <- theme_classic(base_size = 11) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title   = element_text(face = "bold", size = 11)
      )
    theme_bot <- theme_classic(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_text(size = 10),
        plot.title   = element_text(face = "bold", size = 11)
      )

    # Create selective logo plots with position-specific coloring using custom aa_colors
    p_above <- ggseqlogo_selective(pfm_a, color_df = color_df_a, 
                                  position_labels = ordered_pos,
                                  method = "prob", namespace = aa_order,
                                  col_scheme = aa_colors, grey_col = "#cccccc") +
      imgt_x_scale +
      labs(y = "Frequency", title = group_a) +
      theme_top

    p_below <- ggseqlogo_selective(pfm_b, color_df = color_df_b,
                                  position_labels = ordered_pos, 
                                  method = "prob", namespace = aa_order,
                                  col_scheme = aa_colors, grey_col = "#cccccc") +
      imgt_x_scale +
      labs(x = "IMGT position", y = "Frequency",
           title =  group_b) +
      theme_bot
    
    if(scale_log){
      p_above <- p_above +scale_y_log10()
      p_below <- p_below +scale_y_log10()
    }

    return(
      (p_above / p_below) +
        patchwork::plot_annotation(
          title    = paste0("AA composition: ", group_a, " vs ", group_b, add_to_title),
          subtitle = paste0(
            "Significant (p \u2264 ", p_cutoff, ", |D| \u2265 ", d_cutoff,
            "): coloured at specific positions  |  non-significant: grey  |  min_freq \u2265 ", min_freq
          ),
          theme = theme(plot.title = element_text(face = "bold"))
        )
    )
  }

  # ==========================================================================
  # Branch B: geom_text stacking with position-specific coloring
  # Significant AAs: full colour only at positions where they are significant
  # Non-significant AAs: light grey (#cccccc), stacked further from axis
  # x-axis: IMGT position labels from the imgt factor
  # ==========================================================================

  d_long <- bind_rows(
    df_agg |> mutate(direction = "pos", freq = freq_a),
    df_agg |> mutate(direction = "neg", freq = freq_b)
  ) |>
    filter(freq >= min_freq) |>
    mutate(
      # direction=="pos" is group A (above axis), "neg" is group B (below axis).
      dominant_panel = if_else(freq_a >= freq_b, "pos", "neg"),
      # Position-specific coloring: color only if significant at this specific position AND dominant
      color = case_when(
        !is_sig                          ~ "#cccccc", # Not significant at this position
        direction == dominant_panel      ~ aa_colors[aa], # Significant and in dominant panel
        TRUE                             ~ "#cccccc"  # Significant but not in dominant panel
      ),
      color = replace_na(color, "#cccccc")
    )

  if (nrow(d_long) == 0) {
    message("No amino acids above min_freq threshold.")
    return(ggplot() + theme_void())
  }

  # Significant AAs (full color) and non-significant (grey).
  # Within each position and direction, sort by frequency (descending) 
  # so that the most frequent AAs are stacked closest to the axis.
  d_long <- d_long |>
    mutate(is_colored = color != "#cccccc") |>
    group_by(imgt, direction) |>
    arrange(desc(freq), .by_group = TRUE) |>
    mutate(
      cum_f  = cumsum(freq),
      prev_f = lag(cum_f, default = 0),
      ymin   = if_else(direction == "pos",  prev_f,  -cum_f),
      ymax   = if_else(direction == "pos",  cum_f,   -prev_f),
      ymid   = (ymin + ymax) / 2
    ) |>
    ungroup()

  y_max <- max(d_long$ymax, na.rm = TRUE)
  y_min <- min(d_long$ymin, na.rm = TRUE)

  ggplot(d_long, aes(x = imgt)) +
    # faint background bars (grey = non-significant or not-dominant panel)
    geom_rect(
      data = d_long |> filter(!is_colored),
      aes(xmin = as.integer(imgt) - 0.42,
          xmax = as.integer(imgt) + 0.42,
          ymin = ymin, ymax = ymax, fill = color),
      alpha = 0.07, colour = NA
    ) +
    # coloured background bars (significant + dominant panel)
    geom_rect(
      data = d_long |> filter(is_colored),
      aes(xmin = as.integer(imgt) - 0.42,
          xmax = as.integer(imgt) + 0.42,
          ymin = ymin, ymax = ymax, fill = color),
      alpha = 0.18, colour = NA
    ) +
    geom_text(
      aes(x = imgt, y = ymid, label = aa, colour = color, size = freq),
      fontface = "bold",
      family   = "Courier"
    ) +
    scale_size_continuous(range = c(1.5, 15), guide = "none") +
    scale_fill_identity() +
    scale_colour_identity() +
    # IMGT labels on x-axis (from factor levels)
    scale_x_discrete(labels = levels(d_long$imgt)) +
    geom_hline(yintercept = 0, linewidth = 0.7, colour = "grey20") +
    annotate(
      "text", x = ordered_pos[1], y = y_max,
      label = paste0("\u2191 ", group_a),
      hjust = 0, vjust = 1.2, size = 3.5, colour = "grey30", fontface = "italic"
    ) +
    annotate(
      "text", x = ordered_pos[1], y = y_min,
      label = paste0("\u2193 ", group_b),
      hjust = 0, vjust = -0.2, size = 3.5, colour = "grey30", fontface = "italic"
    ) +
    labs(
      x        = "IMGT position",
      y        = "Frequency",
      title    = paste0("AA composition: ", group_a, " (above) vs ", group_b, " (below)"),
      subtitle = paste0(
        "Position-specific coloring: significant (p \u2264 ", p_cutoff, ", |D| \u2265 ", d_cutoff,
        ") AAs colored only at positions where significant  |  non-significant: grey"
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(colour = "grey40", size = 10)
    )
}

# ---------------------------------------------------------------------------
# plot_plogo_like()
#
# p-logo-style visualisation: only amino acids passing significance cutoffs
# are shown at each IMGT position. Letter height encodes either -log10(p-value)
# or |Cohen's D|, chosen via the `height_by` argument.
#
#   Letters ABOVE x-axis  →  positive Cohen's D (enriched in group A)
#   Letters BELOW x-axis  →  negative Cohen's D (enriched in group B)
#   Letter height         ∝  -log10(p_value) OR |Cohen's D| (stacked)
#   Letter colour         =  aa_colors mapping
#   Non-significant AAs   =  omitted entirely
#
# Default column names in df (override any via the cols argument):
#   imgt             (cols$imgt)     - IMGT position
#   aa               (cols$aa)       - amino acid
#   cohens_d_oo      (cols$cohens_d) - Cohen's D (sign sets direction)
#   p_value_oo       (cols$p_value)  - p-value (used for height and/or filtering)
#   label_group_a_oo (cols$label_a)  - group A label
#   label_group_b_oo (cols$label_b)  - group B label
#   (cols$mean_a, cols$mean_b accepted but not required)
#
# Args:
#   df           : data frame as described above
#   p_cutoff     : max p-value to include                  (default 0.05)
#   d_cutoff     : min |Cohen's D| to include              (default 0.5)
#   height_by    : "neg_log10_p" or "cohens_d"             (default "neg_log10_p")
#   positions    : character vector of imgt values; NULL = all
#   aa_colors    : named character vector; NULL uses built-in scheme
#   cols         : named list mapping canonical names to actual column names;
#                  partial lists merged with defaults; only imgt, aa, cohens_d,
#                  p_value, label_a, label_b are required
#   min_freq     : suppress AAs whose max frequency across groups is below this
#                  value; requires cols$mean_a and cols$mean_b             (default 0.005)
#   add_to_title : extra string appended to the plot title
#
# Returns: a ggplot object
# ---------------------------------------------------------------------------

plot_plogo_like <- function(
    df,
    p_cutoff     = 0.05,
    d_cutoff     = 0.5,
    height_by    = c("neg_log10_p", "cohens_d"),
    positions    = NULL,
    aa_colors    = NULL,
    min_freq     = 0.005,
    cols         = list(),
    add_to_title = ""
) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  height_by <- match.arg(height_by)
  if (is.null(aa_colors)) aa_colors <- .default_aa_colors
  if (add_to_title != "") add_to_title <- paste0("; ", add_to_title)

  # ---- resolve column names ----
  required_cols <- c("imgt", "aa", "cohens_d", "p_value", "label_a", "label_b")
  missing <- setdiff(required_cols, names(cols))
  if (length(missing) > 0) {
    stop(paste("The following required mapping(s) are missing from 'cols':",
               paste(missing, collapse = ", ")))
  }

  extra_pval_cols <- setdiff(names(cols)[grepl("^p_value", names(cols))], "p_value")
  extra_d_cols    <- setdiff(names(cols)[grepl("^cohens_d", names(cols))], "cohens_d")

  rename_map <- list(
    imgt             = cols$imgt,
    aa               = cols$aa,
    cohens_d_oo      = cols$cohens_d,
    p_value_oo       = cols$p_value,
    label_group_a_oo = cols$label_a,
    label_group_b_oo = cols$label_b
  )
  # optional frequency columns (accepted but not required)
  if (!is.null(cols$mean_a)) rename_map$mean_group_a_oo <- cols$mean_a
  if (!is.null(cols$mean_b)) rename_map$mean_group_b_oo <- cols$mean_b

  for (col in extra_pval_cols) rename_map[[paste0(col, "_oo")]] <- cols[[col]]
  for (col in extra_d_cols)    rename_map[[paste0(col, "_oo")]] <- cols[[col]]

  df <- .rename_cols(df, rename_map)

  # ---- filter positions ----
  if (!is.null(positions)) df <- df |> filter(imgt %in% positions)
  if (nrow(df) == 0) {
    message("No data for the specified positions.")
    return(ggplot() + theme_void())
  }

  # ---- canonical position order ----
  present_pos <- unique(df$imgt)
  ordered_pos <- .imgt_order[.imgt_order %in% present_pos]
  extra_pos   <- setdiff(present_pos, ordered_pos)
  ordered_pos <- c(ordered_pos, extra_pos)

  # ---- build list of all p_value / cohens_d column pairs ----
  pval_cols <- c("p_value_oo", paste0(extra_pval_cols, "_oo"))
  d_cols    <- c("cohens_d_oo", paste0(extra_d_cols, "_oo"))

  # ---- keep only rows passing ALL filter criteria ----
  keep <- rep(TRUE, nrow(df))
  for (i in seq_along(pval_cols)) {
    pc <- pval_cols[i]; dc <- d_cols[i]
    if (pc %in% colnames(df) && dc %in% colnames(df)) {
      keep <- keep & (df[[pc]] <= p_cutoff & abs(df[[dc]]) >= d_cutoff)
    }
  }
  df_sig <- df[keep, ]

  if (nrow(df_sig) == 0) {
    message("No amino acids pass the significance thresholds.")
    return(ggplot() + theme_void())
  }

  # ---- apply min_freq filter if frequency columns are available ----
  has_freq <- all(c("mean_group_a_oo", "mean_group_b_oo") %in% colnames(df_sig))
  if (!is.null(min_freq) && has_freq) {
    df_sig <- df_sig |>
      filter(pmax(mean_group_a_oo, mean_group_b_oo, na.rm = TRUE) >= min_freq)
    if (nrow(df_sig) == 0) {
      message("No amino acids pass the min_freq threshold.")
      return(ggplot() + theme_void())
    }
  }

  # ---- one row per (imgt, aa) ----
  df_agg <- df_sig |>
    rename(cohens_d = cohens_d_oo, p_value = p_value_oo) |>
    mutate(
      imgt      = factor(imgt, levels = ordered_pos),
      direction = if_else(cohens_d >= 0, "pos", "neg"),
      height    = if (height_by == "neg_log10_p") -log10(p_value) else abs(cohens_d)
    )

  group_a <- unique(df$label_group_a_oo)[1]
  group_b <- unique(df$label_group_b_oo)[1]

  # ---- stack letters within each (imgt, direction): tallest closest to axis ----
  df_agg <- df_agg |>
    group_by(imgt, direction) |>
    arrange(desc(height), .by_group = TRUE) |>
    mutate(
      cum_h  = cumsum(height),
      prev_h = lag(cum_h, default = 0),
      ymin   = if_else(direction == "pos",  prev_h,  -cum_h),
      ymax   = if_else(direction == "pos",  cum_h,   -prev_h),
      ymid   = (ymin + ymax) / 2
    ) |>
    ungroup() |>
    mutate(color = aa_colors[aa],
           color = if_else(is.na(color), "#999999", color))

  # symmetric y-limits so both panels are equal height
  y_range <- max(abs(df_agg$ymax), abs(df_agg$ymin), na.rm = TRUE)

  ggplot(df_agg, aes(x = imgt)) +
    geom_text(
      aes(x = imgt, y = ymid, label = aa, colour = color, size = height),
      fontface = "bold",
      family   = "sans"
    ) +
    scale_size_continuous(range = c(3, 22), guide = "none") +
    scale_colour_identity() +
    scale_x_discrete() +
    geom_hline(yintercept = 0, linewidth = 0.7, colour = "grey20") +
    { if (height_by == "cohens_d") {
        h_min    <- min(df_agg$height, na.rm = TRUE)
        h_max    <- max(df_agg$height, na.rm = TRUE)
        size_ref <- if (h_max > h_min)
          3 + (1 - h_min) / (h_max - h_min) * (22 - 3)
        else 12
        size_ref <- max(3, min(22, size_ref))
        list(
          annotate("text", x = Inf, y =  0.5, label = "A",
                   size = size_ref, colour = "grey75", hjust = 1.2,
                   fontface = "bold", family = "sans"),
          annotate("text", x = Inf, y = -0.5, label = "A",
                   size = size_ref, colour = "grey75", hjust = 1.2,
                   fontface = "bold", family = "sans"),
          annotate("text", x = Inf, y =  0,   label = "|d|=1",
                   size = 2.8, colour = "grey55", hjust = 1.1,
                   vjust = 0.5)
        )
      }
    } +
    annotate(
      "text", x = ordered_pos[1], y =  y_range,
      label = group_a,
      hjust = 0, vjust = 1.2, size = 3.5, colour = "grey30", fontface = "italic"
    ) +
    annotate(
      "text", x = ordered_pos[1], y = -y_range,
      label = group_b,
      hjust = 0, vjust = -0.2, size = 3.5, colour = "grey30", fontface = "italic"
    ) +
    coord_cartesian(ylim = c(-y_range, y_range)) +
    labs(
      x        = "IMGT position",
      y        = NULL,
      title    = paste0("AA enrichment logo: ", group_a, " (above) vs ", group_b, " (below)", add_to_title),
      subtitle = paste0(
        "Letter size = ", if (height_by == "neg_log10_p") "-log10(p)" else "|Cohen's D|",
        "  |  shown: p \u2264 ", format_pval(p_cutoff), ", |D| \u2265 ", d_cutoff,
        "  |  direction = sign(Cohen's D)"
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(colour = "grey40", size = 10)
    )
}

# ---------------------------------------------------------------------------
# plot_property_logo()
#
# Line plot showing per-position mean values of a biochemical property
# (e.g. a Kidera factor) for two groups across IMGT positions.
#
# For each property (prop) one facet panel is produced; optionally a separate
# row per cell type (cells). Group A and group B means are drawn as coloured
# lines with points. Positions passing both q_cutoff and |Cohen's D| cutoffs
# are marked with vertical asterisks (*, **, ***) above the higher mean,
# coloured by the sign of Cohen's D (positive D = color_a, negative D = color_b).
#
# Default column names in df (override any via the cols argument):
#   imgt             (cols$imgt)     - IMGT position label
#   cells            (cols$cells)    - cell type (facet rows; can be constant)
#   prop             (cols$prop)     - property name (e.g. "KF1", "KF2")
#   mean_a_oo        (cols$mean_a)   - mean property value for group A
#   mean_b_oo        (cols$mean_b)   - mean property value for group B
#   cohens_d_oo      (cols$cohens_d) - Cohen's D
#   q_value_oo       (cols$q_value)  - FDR-adjusted q-value
#   label_group_a_oo (cols$label_a)  - group A label (optional; used in legend)
#   label_group_b_oo (cols$label_b)  - group B label (optional; used in legend)
#
# Args:
#   df             : data frame as described above
#   q_cutoff       : max q-value to mark as significant    (default 0.05)
#   d_cutoff       : min |Cohen's D| to mark as significant (default 0.5)
#   positions      : character vector of imgt values; NULL = all
#   properties     : character vector of prop values to include; NULL = all
#   cell_types     : character vector of cells values to include; NULL = all
#   color_a        : line colour for group A               (default "#e41a1c")
#   color_b        : line colour for group B               (default "#377eb8")
#   facet_by_cells : if TRUE and >1 cell type, facet rows by cells (default TRUE)
#   cols           : named list mapping canonical internal names to actual
#                    column names in df; partial lists merged with defaults
#
# Returns: a ggplot object faceted by prop (and optionally cells)
# ---------------------------------------------------------------------------

plot_property_logo <- function(
  df,
  q_cutoff       = 0.05,
  d_cutoff       = 0.5,
  positions      = NULL,
  properties     = NULL,
  cell_types     = NULL,
  color_a        = "#e41a1c",
  color_b        = "#377eb8",
  color_c        = "#4daf4a", # default color for third group
  facet_by_cells = TRUE,
  cols = list(),
  add_to_title = ""
) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  # ---- resolve column names (require all from cols) ----
  required_cols <- c("imgt", "cells", "prop", "mean_a", "mean_b", "cohens_d", "q_value", "label_a", "label_b")
  # Optional third group: mean_c, label_c
  has_group_c <- "mean_c" %in% names(cols) && !is.null(cols$mean_c)
  if (has_group_c) {
    required_cols <- c(required_cols, "mean_c", "label_c")
  }
  missing <- setdiff(required_cols, names(cols))
  if (length(missing) > 0) {
    stop(paste("The following required mapping(s) are missing from 'cols':", 
               paste(missing, collapse = ", ")))
  }
  wrong_names <-setdiff(cols, colnames(df))
  if (length(wrong_names) > 0) {
    stop(paste("The following columns are missing (wrongly named?) from 'cols':", 
               paste(wrong_names, collapse = ", ")))
  }

  # Identify additional q_value and cohens_d columns
  extra_qval_cols <- setdiff(names(cols)[grepl("^q_value", names(cols))], "q_value")
  extra_d_cols    <- setdiff(names(cols)[grepl("^cohens_d", names(cols))], "cohens_d")

  # Build rename mapping for all columns
  rename_map <- list(
    imgt             = cols$imgt,
    cells            = cols$cells,
    prop             = cols$prop,
    mean_a_oo        = cols$mean_a,
    mean_b_oo        = cols$mean_b,
    cohens_d_oo      = cols$cohens_d,
    q_value_oo       = cols$q_value,
    label_group_a_oo = cols$label_a,
    label_group_b_oo = cols$label_b
  )
  if (has_group_c) {
    rename_map$mean_c_oo <- cols$mean_c
    rename_map$label_group_c_oo <- cols$label_c
  }
  for (col in extra_qval_cols) {
    rename_map[[paste0(col, "_oo")]] <- cols[[col]]
  }
  for (col in extra_d_cols) {
    rename_map[[paste0(col, "_oo")]] <- cols[[col]]
  }

  df <- .rename_cols(df, rename_map)

    if(add_to_title!=""){add_to_title <- paste0("; ",add_to_title)}
  # ---- optional subsetting ----
  if (!is.null(positions))  df <- df |> filter(imgt %in% positions)
  #.imgt_order <-intersect(factor(.imgt_order, levels=.imgt_order), as.factor(positions))
  if (!is.null(properties)) df <- df |> filter(prop %in% properties)
  if (!is.null(cell_types)) df <- df |> filter(cells %in% cell_types)

  if (nrow(df) == 0) {
    message("No data after filtering.")
    return(ggplot() + theme_void())
  }

  # ---- canonical IMGT position order ----
  present_pos <- unique(df$imgt)
  ordered_pos <- .imgt_order[.imgt_order %in% present_pos]
  extra_pos   <- setdiff(present_pos, ordered_pos)
  ordered_pos <- c(ordered_pos, extra_pos)

  df <- df |>
    mutate(imgt = factor(imgt, levels = ordered_pos))
#print(cols)
  # ---- group labels ----
  label_a <- ifelse ("label_group_a_oo" %in% names(df),
    unique(na.omit(df$label_group_a_oo))[1], ifelse("label_a" %in% names(cols), cols$label_a,"Group A"))
  label_b <- ifelse ("label_group_b_oo" %in% names(df),
    unique(na.omit(df$label_group_b_oo))[1], ifelse("label_b" %in% names(cols), cols$label_b,"Group B"))

  # ---- significance flags per (imgt, prop, cells) ----
  # Require ALL filtering criteria to be fulfilled simultaneously
  qval_cols <- c("q_value_oo", paste0(extra_qval_cols, "_oo"))
  d_cols    <- c("cohens_d_oo", paste0(extra_d_cols, "_oo"))
  filter_expr <- rep(TRUE, nrow(df))
  for (i in seq_along(qval_cols)) {
    qcol <- qval_cols[i]
    dcol <- d_cols[i]
    if (qcol %in% colnames(df) && dcol %in% colnames(df)) {
      filter_expr <- filter_expr & (df[[qcol]] <= q_cutoff & abs(df[[dcol]]) >= d_cutoff)
    }
  }
  df <- df |>
    mutate(is_sig = filter_expr)

  # ---- reshape to long for lines ----
  d_long <- bind_rows(
    df |> mutate(group = label_a, mean_val = mean_a_oo),
    df |> mutate(group = label_b, mean_val = mean_b_oo)
  )
  if (has_group_c) {
    d_long <- bind_rows(d_long, df |> mutate(group = label_group_c_oo, mean_val = mean_c_oo))
    d_long <- d_long |> mutate(group = factor(group, levels = c(label_a, label_b, df$label_group_c_oo[1])))
  } else {
    d_long <- d_long |> mutate(group = factor(group, levels = c(label_a, label_b)))
  }

  # ---- asterisk positions and labels ----
  # *** q < 0.005,  ** q < 0.01,  * q < 0.05
  # Star colour: color_a for positive Cohen's D, color_b for negative
  sig_marks <- df |>
    filter(is_sig) |>
    rowwise() |>
    mutate(
      y_ast = if (has_group_c) {
        max(mean_a_oo, mean_b_oo, mean_c_oo, na.rm = TRUE)
      } else {
        max(mean_a_oo, mean_b_oo, na.rm = TRUE)
      },
      # Use most conservative (maximum) q-value across all tested q-value columns
      max_q = max(dplyr::c_across(dplyr::any_of(qval_cols)), na.rm = TRUE),
      star_label = dplyr::case_when(
        max_q < 0.005 ~ "***",
        max_q < 0.01 ~ "**",
        max_q < 0.05  ~ "*",
        TRUE          ~ ""
      )
    ) |>
    ungroup()
  print(sig_marks)
  # ---- build plot ----
  color_vals <- c(color_a, color_b)
  color_names <- c(label_a, label_b)
  #manual fix for title when I am using pval, not qval for plotting
  what_is_actually_plotted <- if(grepl(pat="q_value",cols$q_value)){
    'q'}else{
      'p'
  }
  
  if (has_group_c) {
    color_vals <- c(color_vals, color_c)
    color_names <- c(color_names, df$label_group_c_oo[1])
  }
  p <- ggplot(d_long,
              aes(x = imgt, y = mean_val,
                  colour = group, group = group)) +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    scale_colour_manual(
      values = setNames(color_vals, color_names),
      name   = NULL
    ) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey60",
               linetype = "dashed") +
    labs(
      x        = "IMGT position",
      y        = "Mean property value",
      title    = paste0("Property profiles: ", paste(color_names, collapse = " vs "), add_to_title),
      subtitle = paste0(
        "* ",what_is_actually_plotted, "<0.05  ** ",what_is_actually_plotted, "<0.01  *** ",what_is_actually_plotted, "<0.005  |  |D| \u2265 ", d_cutoff
      )
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
      legend.position    = "top",
      strip.text         = element_text(face = "bold", size = 9),
      plot.title         = element_text(face = "bold"),
      plot.subtitle      = element_text(colour = "grey40", size = 9)
    )

  # ---- add asterisks (vertical, star-scaled) ----
  if (nrow(sig_marks) > 0) {
      p <- p + geom_text(
        data        = sig_marks,
        aes(x = imgt, y = y_ast, label = star_label),
        inherit.aes = FALSE,
        colour      = "black",
        angle       = 90,
        hjust       = -0.1,
        size        = 5,
        fontface    = "bold"
      )

  }

  # ---- faceting ----
  if (facet_by_cells && length(unique(df$cells)) > 1) {
    p <- p + facet_grid(cells ~ prop, scales = "free_y")
  } else {
    p <- p + facet_wrap(~ prop, scales = "free_y")
  }

  p
}


###################
## Aminoacid colors for color logo plot


####original version






library(ggplot2)
library(dplyr)
library(tidyr)

# Define the color scheme and properties
aa_data <- data.frame(
  aa = c("A", "V", "I", "L", "M", "F", "W", "Y", "G", "P", "S", "T", "C", "N", "Q", "K", "R", "H", "D", "E"),
  hex = c("#80b1d3", "#80b1d3", "#80b1d3", "#80b1d3", "#80b1d3", 
          "#bc80bd", "#bc80bd", "#fdb462", "#b3de69", "#b3de69", 
          "#fb8072", "#fb8072", "#fdb462", "#bebada", "#bebada", 
          "#1f78b4", "#1f78b4", "#8dd3c7", "#ff7f00", "#ff7f00"),
  property = c(rep("aliphatic / nonpolar", 5), 
               rep("aromatic", 2), "aromatic/polar", 
               rep("small / special", 2), 
               rep("polar uncharged", 2), "polar/sulfur", 
               rep("polar uncharged", 2), 
               rep("charged positive", 2), "charged positive/weak", 
               rep("charged negative", 2))
)


aa_data <- aa_data %>%
  mutate(group = case_when(
    aa %in% c("A", "V", "I", "L", "M") ~ "aliphatic / nonpolar",
    aa %in% c("F", "W", "Y")           ~ "aromatic",
    aa %in% c("G", "P")                ~ "small / special",
    aa %in% c("S", "T", "C", "N", "Q") ~ "polar uncharged",
    aa %in% c("K", "R", "H")           ~ "charged positive",
    aa %in% c("D", "E")                ~ "charged negative"
  )) %>%
  mutate(group = factor(group, levels = c(
    "aliphatic / nonpolar", "aromatic", "small / special", 
    "polar uncharged", "charged positive", "charged negative"
  ))) %>%
  group_by(group) %>%
  mutate(x = row_number()) %>%
  ungroup()

# Plot
p <- ggplot(aa_data, aes(x = x, y = group, fill = hex)) +
  geom_tile(color = "white", linewidth = 2, width = 0.8, height = 0.8) +
  geom_text(aes(label = aa), color = "black", fontface = "bold", size = 8) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev(levels(aa_data$group))) +
  labs(
    title = "Amino Acid Color Scheme & Properties",
    subtitle = "Standard visual coding for TRB CDR3 analysis",
    x = "", 
    y = "Physicochemical Properties"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.y = element_text(face = "bold", size = 12)
  )

# Save the plot

# Print confirmation
#print("Plot saved to assets/aa_color_legend.png")
