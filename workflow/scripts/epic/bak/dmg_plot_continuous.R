dmg_plot_continuous <- function(gene_of_interest, tophits_gr, xlim = NULL, ylim = NULL, smooth = T, stat = "reduce", pval_colname){
  require(ggbio)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  data(genesymbol, package = "biovizBase")
  
  cpgs <- tophits_gr[grep(paste0("(^|;)", as.character(gene_of_interest), "($|;)"), tophits_gr$UCSC_RefGene_Name),] 
  
  #methylation difference
  mdiff_track <- ggplot(cpgs, aes(x = start, y = Beta)) +
    #geom_alignment() + 
    geom_hline(yintercept = 0) +
    #geom_line(alpha = 0.5) +
    theme_bw() +
    theme(legend.pos = "top",
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 12),
          panel.border = element_blank())
  
  if(!is.null(xlim)) mdiff_track <- mdiff_track + xlim(xlim)
  if(!is.null(ylim)) mdiff_track <- mdiff_track + ylim(ylim)
  
  if(smooth){
    mdiff_track <- mdiff_track + geom_smooth(method = "loess")
  } else{
    mdiff_track <- mdiff_track + geom_line()
  }
  
  mdiff_track <- mdiff_track + 
    geom_point(aes(size = -log10(eval(parse(text = pval_colname)))), alpha = 0.5) +
    scale_size_continuous(name = expression(paste(-log[10], "(p-value)")))
    
  gene_track <- ggplot() + 
    #geom_alignment(Homo.sapiens, which = range(cpgs, ignore.strand = T), columns = "GENEID", stat = "reduce") + 
    geom_alignment(TxDb.Hsapiens.UCSC.hg19.knownGene, which = genesymbol[gene_of_interest], stat = stat) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          text = element_text(size = 12),
          panel.border = element_blank())

  if(!is.null(xlim)) gene_track <- gene_track + xlim(xlim)
  
  plotobj <- tracks("Gene" = gene_track, 
                    "% Methylation difference" = mdiff_track, 
                    heights = c(1, 4), 
                    main = gene_of_interest)
  return(plotobj)
}
