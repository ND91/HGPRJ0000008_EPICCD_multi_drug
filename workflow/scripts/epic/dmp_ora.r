#!/usr/bin/env R
# Perform gene set overrepresentation analysis

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

library(minfi)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(missMethyl)

dmps_csv <- args[1]
# ora_kegg_csv <- args[2]
ora_go_csv <- args[2]
# ora_promhyper_kegg_csv <- args[4]
ora_promhyper_go_csv <- args[3]
# ora_promhypo_kegg_csv <- args[6]
ora_promhypo_go_csv <- args[4]
comparison <- args[5]

dmps <- read.csv(dmps_csv)[,-1]

epic_anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

dmps_comparison <- dmps[,c(1, 
                           which(colnames(dmps) %in% paste0("P.Value_", comparison)), 
                           which(colnames(dmps) %in% paste0("adj.P.Val_", comparison)),
                           which(colnames(dmps) %in% paste0("Betadiff_", comparison)))]
dmps_sig <- dmps_comparison %>%
  dplyr::rename(pvalue = 2,
                padj = 3) %>%
  dplyr::filter(padj<0.05)

if(nrow(dmps_sig) > 0){
  # All
  ora_go <- tryCatch({
    go_results <- missMethyl::gometh(sig.cpg = dmps_sig$CGID,
                       all.cpg = dmps_comparison$CGID,
                       collection = c("GO"),
                       array.type = c("EPIC"),
                       plot.bias = FALSE,
                       prior.prob = TRUE,
                       anno = epic_anno,
                       equiv.cpg = TRUE,
                       fract.counts = TRUE,
                       genomic.features = c("ALL"),
                       sig.genes = FALSE)
    go_results <- go_results[order(go_results$P.DE),]
    return(go_results)
  }, error = function(e) {
    NULL
  })
  
  if(!is.null(ora_go)){
    write.csv(ora_go, ora_go_csv)
  } else{
    file.create(ora_go_csv)
  }
  
  # ora_kegg <- missMethyl::gometh(sig.cpg = dmps_sig$CGID,
  #                                all.cpg = dmps_comparison$CGID,
  #                                collection = c("KEGG"),
  #                                array.type = c("EPIC"),
  #                                plot.bias = FALSE,
  #                                prior.prob = TRUE,
  #                                anno = epic_anno,
  #                                equiv.cpg = TRUE,
  #                                fract.counts = TRUE,
  #                                genomic.features = c("ALL"),
  #                                sig.genes = FALSE)
  # ora_kegg <- ora_kegg[order(ora_kegg$P.DE),]
  # write.csv(ora_kegg, ora_kegg_csv)
  
  # Promoter hypermethylation
  dmps_sig_hyper <- dmps_sig %>%
    dplyr::filter(!!sym(paste0("Betadiff_", comparison)) > 0)
  
  if(nrow(dmps_sig_hyper) > 0){
    ora_promhyper_go <- tryCatch({
      go_results <- missMethyl::gometh(sig.cpg = dmps_sig_hyper$CGID,
                                       all.cpg = dmps_comparison$CGID,
                                       collection = c("GO"),
                                       array.type = c("EPIC"),
                                       plot.bias = FALSE,
                                       prior.prob = TRUE,
                                       anno = epic_anno,
                                       equiv.cpg = TRUE,
                                       fract.counts = TRUE,
                                       genomic.features = c("TSS200", "TSS1500", "1stExon"),
                                       sig.genes = FALSE)
      go_results <- go_results[order(go_results$P.DE),]
      return(go_results)
    }, error = function(e) {
      NULL
    })
    
    if(!is.null(ora_promhyper_go)){
      write.csv(ora_promhyper_go, ora_promhyper_go_csv)
    } else{
      file.create(ora_promhyper_go_csv)
    }
    
    # ora_promhyper_kegg <- missMethyl::gometh(sig.cpg = dmps_sig_hyper$CGID,
    #                                          all.cpg = dmps_comparison$CGID,
    #                                          collection = c("KEGG"),
    #                                          array.type = c("EPIC"),
    #                                          plot.bias = FALSE,
    #                                          prior.prob = TRUE,
    #                                          anno = epic_anno,
    #                                          equiv.cpg = TRUE,
    #                                          fract.counts = TRUE,
    #                                          genomic.features = c("TSS200", "TSS1500", "1stExon"),
    #                                          sig.genes = FALSE)
    # ora_promhyper_kegg <- ora_promhyper_kegg[order(ora_promhyper_kegg$P.DE),]
    # write.csv(ora_promhyper_kegg, ora_promhyper_kegg_csv)
  } else{
    file.create(ora_promhyper_go_csv)
    # file.create(ora_promhyper_kegg_csv)
  }
  
  # Promoter hypomethylation
  dmps_sig_hypo <- dmps_sig %>%
    dplyr::filter(!!sym(paste0("Betadiff_", comparison)) < 0)
  
  if(nrow(dmps_sig_hypo) > 0){
    ora_promhypo_go <- tryCatch({
      go_results <- missMethyl::gometh(sig.cpg = dmps_sig_hypo$CGID,
                                       all.cpg = dmps_comparison$CGID,
                                       collection = c("GO"),
                                       array.type = c("EPIC"),
                                       plot.bias = FALSE,
                                       prior.prob = TRUE,
                                       anno = epic_anno,
                                       equiv.cpg = TRUE,
                                       fract.counts = TRUE,
                                       genomic.features = c("TSS200", "TSS1500", "1stExon"),
                                       sig.genes = FALSE)
      go_results <- go_results[order(go_results$P.DE),]
      return(go_results)
    }, error = function(e) {
      NULL
    })
    
    if(!is.null(ora_promhypo_go)){
      write.csv(ora_promhypo_go, ora_promhypo_go_csv)
    } else{
      file.create(ora_promhypo_go_csv)
    }
        
    # ora_promhypo_kegg <- missMethyl::gometh(sig.cpg = dmps_sig_hypo$CGID,
    #                                         all.cpg = dmps_comparison$CGID,
    #                                         collection = c("KEGG"),
    #                                         array.type = c("EPIC"),
    #                                         plot.bias = FALSE,
    #                                         prior.prob = TRUE,
    #                                         anno = epic_anno,
    #                                         equiv.cpg = TRUE,
    #                                         fract.counts = TRUE,
    #                                         genomic.features = c("TSS200", "TSS1500", "1stExon"),
    #                                         sig.genes = FALSE)
    # ora_promhypo_kegg <- ora_promhypo_kegg[order(ora_promhypo_kegg$P.DE),]
    # write.csv(ora_promhypo_kegg, ora_promhypo_kegg_csv)
    } else{
      file.create(ora_promhypo_go_csv)
      # file.create(ora_promhypo_kegg_csv)
    }
} else{
  file.create(ora_go_csv)
  file.create(ora_promhyper_go_csv)
  file.create(ora_promhypo_go_csv)
  # file.create(ora_kegg_csv)
}

sessionInfo()