cpg_timeplot <- function(cgid, anno, betas, response, timepoint, donor, response_cols, enlarged = F){
  #cgid: CpG identifier
  #anno: annotation object as obtained through minfi::getAnnotation()
  #betas: beta object as obtained through minfi::getBeta()
  #response: response variable per sample
  #timepoint: timeooint variable per sample
  #donor: donor variable per sample
  #response_cols: named color vector
  
  beta_df <- data.frame(Beta = betas[cgid,], Response = response, Timepoint = timepoint, Donor = donor)
  
  #Annotation
  hgnc <- paste(unique(strsplit(anno[cgid, "UCSC_RefGene_Name"], ";")[[1]]), collapse = ";")
  coords <- paste0(anno[cgid, "chr"], ":", anno[cgid, "pos"])
  
  if(hgnc != ""){
    stitle <- paste0(coords, " (", hgnc, ")")
  } else{
    stitle <- paste0(coords)
  }
  
  #Visualization
  gplot <- ggplot(beta_df, aes(x = Timepoint, y = Beta, col = Response)) +
    stat_summary(fun.data = mean_se, geom = "crossbar",  aes(color = Response)) +
    geom_point(aes(group = Donor, shape = Response), alpha = 0.25, show.legend = F, position = position_dodge(0.1)) +
    #geom_jitter(aes(group = Donor, shape = Response), width = 0.1, alpha = 0.25, show.legend = F) +
    geom_line(aes(group = Donor, alpha = 0.25), linetype = "dotted", show.legend = F, position = position_dodge(0.1)) +
    #geom_smooth(method = "lm", aes(group = Response), se = F) +
    scale_color_manual(values = response_cols) +
    xlab("Timepoint") +
    ylab("% Methylation") +
    labs(title = cgid,
         subtitle = stitle) +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,3,1)) +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          text = element_text(size = 12))
  
  if(enlarged == F) gplot <- gplot + coord_cartesian(ylim = c(0, 100))
  
  return(gplot)
}

dmg_plot <- function(gene_of_interest, beta_column, pval_column, significance_column, tophits_gr, xlim = NULL, ylim = NULL, significance = T, smooth = T, stat = "reduce"){
  require(ggbio)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  data(genesymbol, package = "biovizBase")
  
  cpgs <- tophits_gr[grep(paste0("(^|;)", as.character(gene_of_interest), "($|;)"), tophits_gr$UCSC_RefGene_Name),] 
  
  cpgs$significant <- "NS"
  #cpgs$significant[data.frame(cpgs)[,pval_column]<0.05] <- "Significant"
  cpgs$significant[data.frame(cpgs)[,significance_column] == "Significant"] <- "Significant"
  
  #methylation difference
  mdiff_track <- ggplot(cpgs, aes(x = start, y = !!sym(beta_column))) +
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
  
  if(significance){
    mdiff_track <- mdiff_track + geom_point(aes(fill = significant, size = -log10(!!sym(pval_column))), shape = 21) +
      labs(size = "-log10(pvalue)")
  } else{
    mdiff_track <- mdiff_track + geom_point(aes(size = -log10(!!sym(pval_column)))) +
      labs(size = "-log10(pvalue)")
  }
  
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
                    "Methylation difference" = mdiff_track, 
                    heights = c(1, 4))
}

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

dmr_plot <- function(dmr_gr, meth_data, anno_gr, meth_groups, groupline = c("group", "sample"), important_cpgs = NULL, flanks = NULL, title = NULL, legend = T, highlight = NULL, rug = T, gridlines = F, smooth = F){
  if(class(dmr_gr) != "GRanges") stop("dmr_gr must be a GRanges object")
  if(is.null(flanks)) flanks <- width(dmr_gr)/2
  if(!is.null(highlight)){
    if(length(highlight) == 1){
      print("Highlight drawn around a single point")
      highlight <- c(highlight, highlight)
    } else if(length(highlight) > 2){
      stop("Highlight can be one or two values")
    }
    highlight <- sort(highlight)
    highlight <- data.frame(x0 = highlight[1], x1 = highlight[2])
  }
  groupline <- match.arg(groupline)
  
  plotrange <- range(dmr_gr + flanks)
  
  cg_ids <- queryHits(findOverlaps(query = anno_gr, subject = dmr_gr))
  cg_ids <- c(min(cg_ids)-1, cg_ids, max(cg_ids)+1)
  anno_sub <- anno_gr[cg_ids,]
  
  meth_df <- data.frame(pos = start(anno_sub), 
                        meth_data[names(anno_sub),],
                        label = names(anno_sub))
  
  meth_df_melt <- melt(meth_df, id = c("pos", "label"), variable.name = "sample", value.name = "methylation")
  meth_df_melt$Group <- rep(meth_groups, each = length(cg_ids))
  
  if(!is.null(important_cpgs)){
    #Dilation factor of the highlighted region
    dil_factor <- ceiling(width(plotrange)/1000)
    if(all(!important_cpgs %in% names(anno_gr))) stop("important_cpgs cannot be found in anno_gr")
    if(all(!important_cpgs %in% meth_df_melt$label)) stop("important_cpgs cannot be found in dmr_gr")
    important_cpgs_coords <- data.frame(imp_cpgs_start = start(anno_gr[important_cpgs,])-dil_factor,
                                        imp_cpgs_end = start(anno_gr[important_cpgs,])+dil_factor)
  } 
  
  if(is.null(title)) title <- "Methylation"
  
  #Plot object
  gplot_obj <- ggplot(meth_df_melt, aes(x = pos, y = methylation))
  
  if(!is.null(highlight)){
    #Dilation factor of the highlighted region
    dil_factor <- ceiling(width(plotrange)/150)
    
    gplot_obj <- gplot_obj + 
      geom_rect(data = highlight,
                inherit.aes = FALSE,
                aes(xmin = x0 - dil_factor, xmax = x1 + dil_factor, ymin = -Inf, ymax = Inf), alpha = 0.15)
  }
  
  gplot_obj <- gplot_obj + 
    geom_point(aes(col = Group, group = Group), alpha = 0.1) +
    coord_cartesian(xlim = c(start(plotrange), end(plotrange))) +
    ylim(0,1) +
    ylab("Beta") +
    labs(title = title,
         subtitle = as.character(dmr_gr)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12))
  
  
  if(groupline == "group"){
      gplot_obj <- gplot_obj + geom_smooth(aes(col = Group, group = Group), method = "loess")
  } else if(groupline == "sample"){
      gplot_obj <- gplot_obj + geom_smooth(aes(col = Group, group = sample), method = "loess", se = F, alpha = 0.2)
  }
  
  # if(groupline == "group"){
  #   if(smooth){
  #     gplot_obj <- gplot_obj + geom_smooth(aes(col = Group, group = Group), method = "loess")
  #   } else{
  #     gplot_obj <- gplot_obj + stat_summary(aes_string(col = "Group", group = Group), fun = mean, geom = "smooth")
  #   }
  # } else if(groupline == "sample"){
  #   if(smooth){
  #     gplot_obj <- gplot_obj + geom_smooth(aes(col = Group, group = sample), method = "loess", se = F, alpha = 0.2)
  #   } else{
  #     gplot_obj <- gplot_obj + stat_summary(aes_string(col = "Group", group = sample), fun = mean, geom = "smooth", alpha = 0.2)
  #   }
  # }
    
  #Label CpGs of interest
  if(!is.null(important_cpgs)){
    gplot_obj <- gplot_obj + 
      geom_rect(data = important_cpgs_coords,
                inherit.aes = F,
                mapping = aes(xmin = imp_cpgs_start, 
                              xmax = imp_cpgs_end, 
                              ymin = 0, 
                              ymax = 1),
                alpha = 0.5)
  } 
  #Add rugplot at the bottom
  if(rug) gplot_obj <- gplot_obj + geom_rug(sides = "b")
  #Add gridline
  if(!gridlines) gplot_obj <- gplot_obj + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #Remove legend
  if(!legend) gplot_obj <- gplot_obj + theme(legend.position = "none")
  
  return(gplot_obj)
}

