This is a set of R functions for plotting. Please rewrite them as a python function. Inputs should be inputs of the function plot_enrichment_depletion2.

ready_extracted is a dataframe from a csv, like:
in ~/sangermac/sr_projects/DR3DR4/analysis_rusted/processed_data/aa_pos_freq/all_together.csv.gz


kidera_hydro_aacids is a named factor vector, like:
kidera_hydro_aacids
     A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V 
 green    red orange orange   blue    red    red  green orange  green   blue    red  green   blue  green orange orange   blue  green  green 
Levels: blue green orange red
 

library(tidyverse)
plot_canvas_shade <-function(type="length", shade_col="red"){
  shade_col <- alpha(shade_col,0.2)
  #plots an empty canvas for AAs in CDR3. lengths 8:22 only, pos 109:114; shades the unused part
  imgt <- factor(c(109:111,
                   paste0("111.",1:4),
                   paste0("112.",5:1),
                   112:114),
                 levels=c(109:111,paste0("111.",1:4),paste0("112.",5:1), 112:114))
  
  if(type=="length"){
    
    avail_aa <- sapply( 1:7, function(z){
      c(rep(1,z),rep(0, 15-2*z),rep(1,z),
        rep(1,z),rep(0,15-2*z-1),rep(1,z+1) )
    })%>%as.vector() 
    
    toplot_canvas <-tibble(imgt=rep(imgt,14),
                           avail_aa =avail_aa,
                           length=factor(rep(paste0("l",9:22),each=15),levels=paste0("l",9:22) ))
    
    p <- toplot_canvas%>%
      ggplot(aes(x=imgt,y=length))+
      geom_tile(aes(fill=as.factor(avail_aa)), color="darkgray",lwd=1)+scale_fill_manual(values = c(shade_col,"white"))+guides(fill="none")+theme_classic()}else{
        
      }
  if(type=="len_bin"){
    avail_aa =c(rep(1,2),rep(0,11),rep(1,2),
                rep(1,3),rep(0,9),rep(1,3),
                rep(1,4),rep(0,7),rep(1,4),
                rep(1,5),rep(0,5),rep(1,5),
                rep(1,15))
    toplot_canvas_bin <-tibble(imgt=rep(imgt,5),
                               avail_aa,
                               length=factor(rep(c("(0,11]",  "(11,13]" ,"(13,15]", "(15,17]", "(17,28]"), each=15), levels=c("(0,11]",  "(11,13]" ,"(13,15]", "(15,17]" ,"(17,28]")))
    
    
    p <- toplot_canvas_bin%>%
      ggplot(aes(x=imgt,y=length,fill=as.factor(avail_aa)))+
      geom_tile(color="darkgray",lwd=1)+scale_fill_manual(values = c(shade_col,"white"))+guides(fill="none")+theme_classic()
    
  }
  return(p)
}


plot_enrichment_depletion2 <- function(ready_extracted, length_col="length", estim_col="estim_hom_vs_DQ2DQ8",
                                       estim_expr_for_plotting="estim_hom_vs_DQ2DQ8<0", 
                                       aminoacid_mapping=kidera_hydro_aacids, inset_colour="red"){
  #newer version
 
  #colour aminoacid by their properties
  #Requires splitting aminoacids from ready_extracted
  ready_extracted_splitted <- ready_extracted%>%
    filter( eval(parse(text=estim_expr_for_plotting)))%>%
    group_by(length, cells, imgt)%>%
    summarise(aa=paste0(aa, collapse=""),avail_aa=1)%>%
    ungroup%>%
    mutate(split_aa=strsplit(aa, split=""))%>%
    unnest(cols=c(split_aa))%>%
    mutate(aa_colour=aminoacid_mapping[split_aa])
  
  layers <- ready_extracted_splitted$aa_colour %>%unique()%>%as.character
  
  all_layers <- sapply(layers, function(color){
    aas_to_replace_by_spaces <- names(aminoacid_mapping[aminoacid_mapping!=color])
    ready_extracted%>%
      filter( eval(parse(text=estim_expr_for_plotting)))%>%
      mutate(aa=gsub(pat=paste(aas_to_replace_by_spaces,collapse ="|"), rep="  ",x=aa),
             aa_col=!!color)%>%
      group_by(length, cells, imgt)%>%
      summarise(aa=paste0(aa, collapse="",sep=""),avail_aa=1, aa_col=unique(aa_col))
  }, USE.NAMES = TRUE, simplify = FALSE)
  
  colors=names(all_layers)
  names(colors)<- names(all_layers)
    
  p <- plot_canvas_shade(type=length_col,shade_col=inset_colour )
  #how many layers do I need
 
  for(i in 1:length(all_layers)){
    
  
  p=p + geom_text(data=all_layers[[i]], aes(label=aa, y=length, x=imgt, col=aa_col),  size=5, fontface='bold')+ scale_colour_manual(values=colors)
}

p
}