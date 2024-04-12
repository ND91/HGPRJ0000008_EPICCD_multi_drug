cpg_summaryplot <- function(cgid, anno, betas, response, timepoint, donor, dmp_results, icc_results, response_cols, enlarged = F){
  #cgid: CpG identifier
  #anno: annotation object as obtained through minfi::getAnnotation()
  #betas: beta object as obtained through minfi::getBeta()
  #response: response variable per sample
  #timepoint: timeooint variable per sample
  #donor: donor variable per sample
  #dmp_results: results from the DMP analysis. This is the output from the script workflow/epic/dmp_analyses.R
  #icc_results: results from the ICC analysis. This is the output from the script workflow/epic/icc_analyses.R
  #response_cols: named color vector
    
  dmp_subset <- rbind(data.frame(dmp_results[match(cgid, dmp_results$CGID), c("CGID", "Betadiff_T1RvNR", "P.Value_T1RvNR", "adj.P.Val_T1RvNR")], Timepoint = "T1") %>%
                        dplyr::rename(Betadiff = Betadiff_T1RvNR,
                                      P.Value = P.Value_T1RvNR,
                                      adj.P.Val = adj.P.Val_T1RvNR),
                      data.frame(dmp_results[match(cgid, dmp_results$CGID), c("CGID", "Betadiff_T2RvNR", "P.Value_T2RvNR", "adj.P.Val_T2RvNR")], Timepoint = "T2") %>%
                        dplyr::rename(Betadiff = Betadiff_T2RvNR,
                                      P.Value = P.Value_T2RvNR,
                                      adj.P.Val = adj.P.Val_T2RvNR))
  
  beta_df <- data.frame(Beta = betas[cgid,], 
                        Response = response, 
                        Timepoint = timepoint, 
                        Donor = donor) %>%
    dplyr::left_join(dmp_subset, by = "Timepoint") %>%
    dplyr::left_join(icc_results, by = "CGID") %>%
    #dplyr::mutate(label = paste0(Timepoint, "\np-value = ", round(P.Value, 3)))
    dplyr::mutate(label = ifelse(P.Value<0.01, paste0(Timepoint, "\np-value = ", formatC(P.Value, format = "e", digits = 3)), paste0(Timepoint, "\np-value = ", round(P.Value, 3))))
  
  #Annotation
  hgnc <- paste(unique(strsplit(anno[cgid, "UCSC_RefGene_Name"], ";")[[1]]), collapse = ";")
  coords <- paste0(anno[cgid, "chr"], ":", anno[cgid, "pos"])
  
  if(hgnc != ""){
    stitle <- paste0(coords, " (", hgnc, ")\nICC = ", round(icc_results[match(cgid, dmp_results$CGID),]$icc, 2), " [", round(icc_results[match(cgid, dmp_results$CGID),]$icc_lbound, 2), "; ", round(icc_results[match(cgid, dmp_results$CGID),]$icc_rbound, 2), "]")
  } else{
    stitle <- paste0(coords, "\nICC = ", round(icc_results[match(cgid, dmp_results$CGID),]$icc, 2), " [", round(icc_results[match(cgid, dmp_results$CGID),]$icc_lbound, 2), "; ", round(icc_results[match(cgid, dmp_results$CGID),]$icc_rbound, 2), "]")
  }
  
  #Visualization
  tplot <- ggplot(beta_df, aes(x = Timepoint, y = Beta*100)) +
    geom_point(aes(group = Donor, shape = Response, col = Response), show.legend = F, position = position_dodge(0.1)) +
    geom_line(aes(group = Donor, alpha = 0.25), linetype = "dotted", show.legend = F, position = position_dodge(0.1)) +
    xlab("Timepoint") +
    ylab("% Methylation") +
    labs(title = cgid,
         subtitle = stitle) +
    scale_color_manual(values = response_cols) + 
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          text = element_text(size = 12))
  
  if(enlarged == F) tplot <- tplot + coord_cartesian(ylim = c(0, 100))
  
  bplot <- ggplot(beta_df, aes(x = Response, y = Beta*100, fill = Response)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(size = 2, shape = 21) +
    facet_wrap(~label, nrow = 1) +
    labs(x = "Response",
         y = "% Methylation") +
    scale_fill_manual(values = response_cols) + 
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          text = element_text(size = 12))
  
  if(enlarged == F) bplot <- bplot + coord_cartesian(ylim = c(0, 100))
  
  gplot <- ggarrange(tplot, bplot, nrow = 2, ncol = 1, align = "hv")
  
  return(gplot)
}