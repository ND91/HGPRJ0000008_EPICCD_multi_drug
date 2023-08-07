dmr_effect_plot <- function(region_of_interest, tophits_gr, enhancers = F, ylim = c(-0.25, 0.25), significance = T, smooth = T){
  if(enhancers){
    cpgs <- tophits_gr[grep(paste0("(^|;)", as.character(region_of_interest$gene), "($|;)"), paste0(tophits_gr$UCSC_RefGene_Name, ";", tophits_gr$enhancer_gene)),] 
  } else{
    cpgs <- tophits_gr[grep(paste0("(^|;)", as.character(region_of_interest$gene), "($|;)"), tophits_gr$UCSC_RefGene_Name),] 
  }
  
  mdiff_track <- ggplot(cpgs, aes(x = start, y = Beta)) +
    geom_alignment() + 
    geom_hline(yintercept = 0) +
    #geom_line(alpha = 0.5) +
    ylim(ylim) +
    theme_bw() +
    theme(legend.pos = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 12))
  
  if(smooth){
    mdiff_track <- mdiff_track + geom_smooth(method = "loess", alpha = 0.5)
  } else{
    mdiff_track <- mdiff_track + geom_line(alpha = 0.5)
  }
  
  
  if(significance){
    cpgs$significance <- cpgs$adj.P.Val<0.05
    mdiff_track <- mdiff_track + geom_point(aes(size = -log10(P.Value), fill = significance), shape = 21)
  } else{
    mdiff_track <- mdiff_track + geom_point(aes(size = -log10(P.Value))) 
  }
  
  gene_track <- ggplot() + 
    geom_alignment(Homo.sapiens, which = range(cpgs, ignore.strand = T)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 12))

  plotobj <- tracks("Gene" = gene_track, 
                    "Methylation difference" = mdiff_track, 
                    heights = c(3, 4))
}
