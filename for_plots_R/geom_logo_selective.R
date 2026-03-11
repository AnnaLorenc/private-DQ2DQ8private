library(ggplot2)
library(ggseqlogo)

#' Selective Logo Plot with Position-Specific Coloring
#' 
#' Creates a sequence logo plot where only specific amino acids at specific 
#' positions are colored according to the color scheme, while all others 
#' appear in grey. Automatically labels x-axis positions using column names
#' from matrix data or custom position labels.
#' 
#' @param data Sequence data (same format as ggseqlogo::geom_logo)
#' @param color_df Data frame with two columns: 'position' and 'aa' (amino acid letter).
#'                 For matrix data with column names, 'position' should match column names.
#'                 For sequence strings, 'position' should be sequential (1, 2, 3, ...).
#' @param position_labels Optional character vector of position labels. If provided,
#'                        positions 1, 2, 3... will be mapped to these labels for 
#'                        matching with color_df$position and x-axis labeling.
#'                        If NULL and data is a matrix with column names, column names
#'                        are used automatically for x-axis labels.
#' @param method Method for logo calculation ("bits", "probability", or "custom")
#' @param seq_type Sequence type ("auto", "aa", "dna", "rna", or "other")
#' @param namespace Custom namespace for letters
#' @param font Font to use for letters
#' @param stack_width Width of letter stacks (0-1)
#' @param rev_stack_order Reverse the order of letter stacking
#' @param col_scheme Color scheme to use (only for colored letters)
#' @param grey_col Color to use for non-highlighted letters
#' @param low_col Low color for gradient schemes
#' @param high_col High color for gradient schemes
#' @param plot Whether to return a plot or just the data
#' @param ... Additional parameters passed to the polygon geom
#' 
#' @details 
#' Position Labeling:
#' - For matrix data with column names (e.g., "111", "111.1", "112"), these names 
#'   automatically appear as x-axis labels
#' - For sequence strings, provide position_labels to specify custom x-axis labels
#' - X-axis labels are rotated 45 degrees to prevent overlap
#' 
#' @return A ggplot layer or list of layers
#' @export
geom_logo_selective <- function(data = NULL, color_df = NULL, position_labels = NULL, 
                               method = "bits", seq_type = "auto", namespace = NULL, 
                               font = "roboto_medium", stack_width = 0.95, 
                               rev_stack_order = FALSE, col_scheme = "auto", 
                               grey_col = "grey70", low_col = "black", 
                               high_col = "yellow", plot = TRUE, ...) {
  
  if (stack_width > 1 | stack_width <= 0) 
    stop("\"stack_width\" must be between 0 and 1")
  if (is.null(data)) 
    stop("Missing \"data\" parameter!")
  if (is.null(color_df))
    stop("Missing \"color_df\" parameter! Provide a data frame with 'position' and 'aa' columns.")
  
  # Validate color_df
  if (!is.data.frame(color_df) || !all(c("position", "aa") %in% names(color_df))) {
    stop("color_df must be a data frame with columns 'position' and 'aa'")
  }
  
  # Handle position labels/mapping
  if (is.matrix(data) && !is.null(colnames(data))) {
    # For matrix data with column names, use those as position labels
    position_labels <- colnames(data)
  } else if (is.list(data) && !is.null(names(data))) {
    # For list data, check if the first element is a matrix with column names
    first_data <- data[[1]]
    if (is.matrix(first_data) && !is.null(colnames(first_data))) {
      position_labels <- colnames(first_data)
    }
  }
  
  if (!is.null(namespace)) 
    seq_type = "other"
  
  all_methods = c("bits", "probability", "custom")
  pind = pmatch(method, all_methods)
  method = all_methods[pind]
  if (is.na(method)) 
    stop("method must be one of 'bits', 'probability', or 'custom'")
  
  # Process data (same as original geom_logo)
  if (is.character(data) | is.matrix(data)) 
    data = list(`1` = data)
  
  if (is.list(data)) {
    if (is.null(names(data))) 
      names(data) = seq_along(data)
    lvls = names(data)
    data_sp = lapply(names(data), function(n) {
      curr_seqs = data[[n]]
      ggseqlogo:::logo_data(seqs = curr_seqs, method = method, 
                           stack_width = stack_width, 
                           rev_stack_order = rev_stack_order, 
                           seq_group = n, seq_type = seq_type, 
                           font = font, namespace = namespace)
    })
    data = do.call(rbind, data_sp)
    data$seq_group = factor(data$seq_group, levels = lvls)
  }
  
  if (!plot) 
    return(data)
  
  seq_type = attr(data, "seq_type")
  
  # Get the original color scheme
  if (is.character(col_scheme) && length(col_scheme) == 1) {
    # Standard ggseqlogo color scheme name
    cs = ggseqlogo:::get_col_scheme(col_scheme, seq_type)
  } else if (is.list(col_scheme) || (is.character(col_scheme) && length(col_scheme) > 1)) {
    # Custom color scheme - create from named vector or list
    if (is.character(col_scheme) && !is.null(names(col_scheme))) {
      # Named character vector (like aa_colors)
      aa_order <- c("A","C","D","E","F","G","H","I","K","L",
                    "M","N","P","Q","R","S","T","V","W","Y")
      cs_colors <- sapply(aa_order, function(aa) {
        if (aa %in% names(col_scheme)) {
          col_scheme[[aa]]
        } else {
          grey_col
        }
      })
      cs <- data.frame(
        letter = aa_order,
        group = aa_order,
        col = unname(cs_colors),
        stringsAsFactors = FALSE
      )
    } else {
      # Fall back to auto
      cs = ggseqlogo:::get_col_scheme("auto", seq_type)
    }
  } else {
    # Fall back to auto
    cs = ggseqlogo:::get_col_scheme("auto", seq_type)
  }
  
  legend_title = attr(cs, "cs_label")
  if (is.null(legend_title)) legend_title = "Amino Acid"
  
  # Create a selective color scheme
  # First, merge with original color scheme to get all color information
  data = merge(data, cs, by = "letter", all.x = TRUE)
  
  # Handle position mapping for color_df
  if (!is.null(position_labels)) {
    # Create mapping from sequential positions to actual labels
    pos_mapping <- data.frame(
      seq_pos = 1:length(position_labels),
      actual_pos = position_labels,
      stringsAsFactors = FALSE
    )
    
    # Convert color_df positions to sequential positions
    if (nrow(color_df) > 0) {
      color_df_mapped <- merge(color_df, pos_mapping, 
                              by.x = "position", by.y = "actual_pos", 
                              all.x = TRUE)
      
      if (any(is.na(color_df_mapped$seq_pos))) {
        missing_pos <- color_df$position[is.na(color_df_mapped$seq_pos)]
        warning(paste("The following positions in color_df were not found in the data:", 
                     paste(missing_pos, collapse = ", ")))
        color_df_mapped <- color_df_mapped[!is.na(color_df_mapped$seq_pos), ]
      }
      
      # Use sequential positions for lookup
      color_lookup_df <- data.frame(
        position = color_df_mapped$seq_pos,
        aa = color_df_mapped$aa,
        should_color = TRUE,
        stringsAsFactors = FALSE
      )
    } else {
      # Empty color_df case
      color_lookup_df <- data.frame(
        position = integer(0),
        aa = character(0),
        should_color = logical(0),
        stringsAsFactors = FALSE
      )
    }
  } else {
    # Use positions as-is (assume they're already sequential)
    if (nrow(color_df) > 0) {
      color_lookup_df <- data.frame(
        position = color_df$position,
        aa = color_df$aa,
        should_color = TRUE,
        stringsAsFactors = FALSE
      )
    } else {
      # Empty color_df case
      color_lookup_df <- data.frame(
        position = integer(0),
        aa = character(0),
        should_color = logical(0),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Create a lookup for positions and amino acids that should be colored
  color_lookup <- merge(data.frame(position = data$position, letter = data$letter), 
                       color_lookup_df, 
                       by.x = c("position", "letter"), 
                       by.y = c("position", "aa"), 
                       all.x = TRUE)
  color_lookup$should_color[is.na(color_lookup$should_color)] <- FALSE
  
  # Add the should_color information back to data
  data = merge(data, color_lookup[, c("position", "letter", "should_color")], 
              by = c("position", "letter"), all.x = TRUE)
  
  # Modify colors: use grey for letters that shouldn't be colored
  data$original_col = data$col
  data$original_group = data$group
  data$col = ifelse(data$should_color, data$col, grey_col)
  data$group = ifelse(data$should_color, data$group, "grey")
  
  data = data[order(data$order), ]
  
  # Handle color scale
  colscale_gradient = is.numeric(cs$group)
  colscale_opts = NULL
  
  if (colscale_gradient) {
    # For gradient schemes, we need to handle the grey separately
    unique_groups = unique(data$group[data$should_color])
    if (length(unique_groups) > 0) {
      colscale_opts = scale_fill_gradient(name = legend_title, 
                                         low = low_col, high = high_col, 
                                         na.value = grey_col)
    } else {
      colscale_opts = scale_fill_manual(values = c("grey" = grey_col), 
                                       name = legend_title)
    }
  } else {
    # For discrete schemes
    tmp = cs[!duplicated(cs$group) & !is.na(cs$group), ]
    col_map = unlist(split(tmp$col, tmp$group))
    
    # Add grey color
    col_map = c(col_map, "grey" = grey_col)
    
    colscale_opts = scale_fill_manual(values = col_map, name = legend_title, 
                                     na.value = grey_col)
  }
  
  guides_opts = NULL
  if (identical(cs$letter, cs$group)) 
    guides_opts = guides(fill = FALSE)
  
  y_lim = NULL
  extra_opts = NULL
  
  # Set y-axis label based on method
  if (method == "custom") {
    y_lab = ""
  } else {
    y_lab = method
    substr(y_lab, 1, 1) = toupper(substr(y_lab, 1, 1))
  }
  
  data$group_by = with(data, interaction(seq_group, letter, position))
  data$x = data$x
  
  logo_layer = layer(stat = "identity", data = data, 
                    mapping = aes_string(x = "x", y = "y", fill = "group", 
                                        group = "group_by"), 
                    geom = "polygon", position = "identity", 
                    show.legend = NA, inherit.aes = FALSE, 
                    params = list(na.rm = TRUE, ...))
  
  breaks_fun = function(lim) {
    1:floor(lim[2]/1.05)
  }
  
  # Create custom x-axis labels if position_labels are available
  x_scale_opts = NULL
  if (!is.null(position_labels) && length(position_labels) > 0) {
    # Create breaks and labels for positions present in the data
    max_pos = max(data$position, na.rm = TRUE)
    breaks_vec = 1:max_pos
    labels_vec = position_labels[1:max_pos]
    
    x_scale_opts = scale_x_continuous(
      breaks = breaks_vec,
      labels = labels_vec,
      name = "Position"
    )
  } else {
    x_scale_opts = scale_x_continuous(breaks = breaks_fun, labels = identity)
  }
  
  list(logo_layer, 
       x_scale_opts,
       ylab(y_lab), 
       xlab(""), 
       colscale_opts, 
       guides_opts, 
       coord_cartesian(ylim = y_lim), 
       extra_opts,
       # Add theme to rotate x-axis labels if they are long
       theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  )
}

#' Selective Logo Plot Function (wrapper similar to ggseqlogo)
#' 
#' Creates a sequence logo plot with selective amino acid coloring and automatic
#' position labeling using column names from input data.
#' 
#' @param data Sequence data
#' @param color_df Data frame with 'position' and 'aa' columns for selective coloring
#' @param position_labels Optional character vector for custom position labels
#' @param facet Faceting option ("wrap" or "grid")
#' @param scales Scale option for facets
#' @param ncol Number of columns for facet_wrap
#' @param nrow Number of rows for facet_wrap
#' @param ... Additional parameters passed to geom_logo_selective
#' 
#' @details
#' This function automatically labels the x-axis with actual position names:
#' - Matrix data: uses column names (e.g., "111", "111.1", "112")  
#' - Sequence strings: uses position_labels parameter
#' - Labels are rotated 45° to prevent overlap
#' 
#' @return A ggplot object
#' @export
#' @param data Sequence data
#' @param color_df Data frame with 'position' and 'aa' columns for selective coloring
#' @param facet Faceting option ("wrap" or "grid")
#' @param scales Scale option for facets
#' @param ncol Number of columns for facet_wrap
#' @param nrow Number of rows for facet_wrap
#' @param ... Additional parameters passed to geom_logo_selective
#' 
#' @return A ggplot object
#' @export
ggseqlogo_selective <- function(data, color_df = NULL, position_labels = NULL, 
                               facet = "wrap", scales = "free_x", ncol = NULL, 
                               nrow = NULL, ...) {
  
  p = ggplot() + geom_logo_selective(data = data, color_df = color_df, 
                                    position_labels = position_labels, ...) + 
      ggseqlogo::theme_logo()
  
  if (!"list" %in% class(data)) 
    return(p)
  
  facet_opts = c("grid", "wrap")
  pind = pmatch(facet, facet_opts)
  facet = facet_opts[pind]
  if (is.na(facet)) 
    stop("facet option must be set to 'wrap' or 'grid'")
  
  if (facet == "grid") {
    p = p + facet_grid(~seq_group, scales = scales)
  } else if (facet == "wrap") {
    p = p + facet_wrap(~seq_group, scales = scales, nrow = nrow, ncol = ncol)
  }
  
  return(p)
}

# Example usage:
# 
# # Method 1: Matrix data with column names (IMGT positions)
# # Your matrix should have columns named like "111", "111.1", "112", etc.
# sequences_matrix <- your_sequence_matrix  # columns are IMGT positions
# color_data <- data.frame(
#   position = c("111", "111", "111.1"),  # Use actual IMGT position names
#   aa = c("A", "C", "R")
# )
# p1 <- ggseqlogo_selective(sequences_matrix, color_df = color_data)
#
# # Method 2: Sequence strings with explicit position labels
# sequences <- c("ACDEFGHIKLMNPQRSTVWY", "ACDAFGHIKLMNPQRSTVWY")
# position_labels <- c("111", "111.1", "111.2", "111.3", "111.4", "111.5",
#                      "112.5", "112.4", "112.3", "112.2", "112.1", "112",
#                      "113", "114", "115", "116", "117", "118", "119", "120")
# color_data <- data.frame(
#   position = c("111", "111.1", "112"),  # Use actual IMGT position names  
#   aa = c("A", "C", "E")
# )
# p2 <- ggseqlogo_selective(sequences, color_df = color_data, 
#                          position_labels = position_labels)
#
# print(p1)
# print(p2)