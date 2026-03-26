library(Peptides)
data(AAdata)
library(scico)
library(scales)


source("~/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/geom_logo_selective.R")
source("~/Documents/DQ2DQ8/pipeline/DQ2DQ8/for_plots_R/plot_aa_logo.R")


color_scale_properties_property <- function(AAdata_factors, fact_to_scale, palette=NA){
  #AAdata_factors for example AAdata$kideraFactors
  if(is.na(palette)){
    palette <-  scico(256, palette = "vik")
  } 
  
  df <- AAdata_factors%>%bind_rows(.id="param")
  
  df_colors <- df %>%
    pivot_longer( -param, names_to = "aa", values_to = "value") %>%
    group_by(param) %>%
    mutate(
      z = (value - mean(value)) / sd(value),   # row-wise z-score
      z = pmax(pmin(z, 2), -2),                # clip to fixed dynamic range
      color = col_numeric(
        palette = palette,
        domain = c(-2, 2)
      )(z)
    ) %>%
    ungroup() %>%
    select(param, aa, color) 
  
  df_colors|>filter(param==fact_to_scale)|>
    select(-param)|> with(setNames(color, aa))
  
  
}

change_imgt_format <- function(f) {
  #to change from 11A into 11.1 etc
  is_letter <- grepl("[A-Z]$", f)
  prefix <- ifelse(is_letter, substr(f, 1, nchar(f)-1), f)
  suffix <- ifelse(is_letter, substr(f, nchar(f), nchar(f)), "")
  suffix_num <- ifelse(is_letter, match(suffix, LETTERS), "")
  paste0(prefix, ifelse(is_letter, paste0(".", suffix_num), ""))
}


plot_canvas_shade_all_lengths <-function(shade_col="red"){
  shade_col <- alpha(shade_col,0.2)
  #plots an empty canvas for AAs in CDR3. lengths 8:22 only, pos 109:114; shades the unused part
  imgt_all <- factor(c(109:111,
                       paste0("111.",1:4),
                       paste0("112.",5:1),
                       112:114),
                     levels=c(109:111,paste0("111.",1:4),paste0("112.",5:1), 112:114))
  
  
  avail_aa <- sapply( 1:7, function(z){
    c(rep(1,z),rep(0, 15-2*z),rep(1,z),
      rep(1,z),rep(0,15-2*z-1),rep(1,z+1) )
  })%>%as.vector() 
  
  toplot_canvas <-tibble(imgt=rep(imgt_all, 14),
                         avail_aa =avail_aa,
                         length=rep(9:22,each=15))%>%
    mutate(imgt_numeric=as.numeric(factor(imgt, levels=imgt_all) ))
  
  for_x_labels <-toplot_canvas%>%
    select(imgt, length, imgt_numeric)%>%filter(length==length(imgt_all))
  
  p <- toplot_canvas%>%
    ggplot(aes(x=imgt_numeric, y=length))+
    geom_tile(aes(fill=as.factor(avail_aa)), color="darkgray",lwd=0.5)+
    scale_x_continuous(  breaks =  for_x_labels$imgt_numeric,
                         labels =  for_x_labels$imgt)+
    scale_y_continuous(breaks=9:22,labels=9:22)+
    scale_fill_manual(values = c(shade_col,"white"))+
    guides(fill="none")+theme_classic()+
    xlab("IMGT position")
}   




plot_enrichment_depletion3 <- function(ready_extracted, length_col="length", estim_col="estim_hom_vs_DQ2DQ8",
                                       estim_expr_for_plotting="estim_hom_vs_DQ2DQ8<0", 
                                       aminoacid_mapping=kidera_hydro_aacids, inset_colour="red"){
  #newest version - using position
  imgt_all <- factor(c(109:111,
                       paste0("111.",1:4),
                       paste0("112.",5:1),
                       112:114),
                     levels=c(109:111,paste0("111.",1:4),paste0("112.",5:1), 112:114))
  ##assumption: up to 5 letters per imgt position
  
  #colour aminoacid by their properties
  #Requires splitting aminoacids from ready_extracted, to add each one in a different colour
  ready_extracted_splitted <- ready_extracted%>%
    filter( eval(parse(text=estim_expr_for_plotting)), imgt %in%imgt_all)%>%
    group_by(length, cells, imgt)%>%
    summarise(aa=paste0(aa, collapse=""), avail_aa=1)%>%
    ungroup%>%
    mutate(aa_split = strsplit(aa, "")) %>%       # split string into letters
    unnest(aa_split) %>%                          # one row per letter
    group_by(length, cells, imgt) %>%            # keep original grouping
    mutate(position = row_number(),               # letter position in string
           color = aminoacid_mapping[aa_split]) %>%      # map color
    ungroup()%>%
    mutate(imgt_numeric=as.numeric(factor(imgt, levels=imgt_all)),
           imgt_numeric_pos=imgt_numeric + (position-1)/nchar(aa)- 0.25)
  
  
  
  p <- plot_canvas_shade_all_lengths(shade_col=inset_colour )
  p <-p+ geom_text(data=ready_extracted_splitted, aes(y=length, x=imgt_numeric_pos, label=aa_split, col=color),family = "Courier", position="identity", size=5)+ 
    scale_color_identity()
  p
}

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

