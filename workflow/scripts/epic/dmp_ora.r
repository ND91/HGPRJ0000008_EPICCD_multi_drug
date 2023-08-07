#!/usr/bin/env R
# Perform gene set overrepresentation analysis

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste0("Script needs 8 arguments. Current input is:", args))
}

# if(!"missMethyl" %in% rownames(installed.packages())){
#   BiocManager::install("missMethyl")
#   require(missMethyl)
# } else{
#   require(missMethyl)
# }

suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

dmps_csv_path <- args[1]
ora_kegg_csv_path <- args[2]
ora_go_csv_path <- args[3]
ora_promhyper_kegg_csv_path <- args[4]
ora_promhyper_go_csv_path <- args[5]
ora_promhypo_kegg_csv_path <- args[6]
ora_promhypo_go_csv_path <- args[7]
comparison <- args[8]

dmps <- read.csv(dmps_csv_path)[,-1]

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
  ora_go <- missMethyl::gometh(sig.cpg = dmps_sig$CGID,
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
  ora_go <- ora_go[order(ora_go$P.DE),]
  write.csv(ora_go, ora_go_csv_path)
  
  ora_kegg <- missMethyl::gometh(sig.cpg = dmps_sig$CGID,
                                 all.cpg = dmps_comparison$CGID,
                                 collection = c("KEGG"),
                                 array.type = c("EPIC"),
                                 plot.bias = FALSE,
                                 prior.prob = TRUE,
                                 anno = epic_anno,
                                 equiv.cpg = TRUE,
                                 fract.counts = TRUE,
                                 genomic.features = c("ALL"),
                                 sig.genes = FALSE)
  ora_kegg <- ora_kegg[order(ora_kegg$P.DE),]
  write.csv(ora_kegg, ora_kegg_csv_path)
  
  # Promoter hypermethylation
  dmps_sig_hyper <- dmps_sig %>%
    dplyr::filter(Betadiff_RT2vT1 > 0)
  
  if(nrow(dmps_sig_hyper) > 0){
    ora_promhyper_go <- missMethyl::gometh(sig.cpg = dmps_sig_hyper$CGID,
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
    ora_promhyper_go <- ora_promhyper_go[order(ora_promhyper_go$P.DE),]
    write.csv(ora_promhyper_go, ora_promhyper_go_csv_path)
    
    ora_promhyper_kegg <- missMethyl::gometh(sig.cpg = dmps_sig_hyper$CGID,
                                             all.cpg = dmps_comparison$CGID,
                                             collection = c("KEGG"),
                                             array.type = c("EPIC"),
                                             plot.bias = FALSE,
                                             prior.prob = TRUE,
                                             anno = epic_anno,
                                             equiv.cpg = TRUE,
                                             fract.counts = TRUE,
                                             genomic.features = c("TSS200", "TSS1500", "1stExon"),
                                             sig.genes = FALSE)
    ora_promhyper_kegg <- ora_promhyper_kegg[order(ora_promhyper_kegg$P.DE),]
    write.csv(ora_promhyper_kegg, ora_promhyper_kegg_csv_path)
  } else{
    file.create(ora_promhyper_go_csv_path)
    file.create(ora_promhyper_kegg_csv_path)
  }
  
  # Promoter hypomethylation
  dmps_sig_hypo <- dmps_sig %>%
    dplyr::filter(Betadiff_RT2vT1 < 0)
  
  if(nrow(dmps_sig_hypo) < 0){
    ora_promhypo_go <- missMethyl::gometh(sig.cpg = dmps_sig_hypo$CGID,
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
    ora_promhypo_go <- ora_promhypo_go[order(ora_promhypo_go$P.DE),]
    write.csv(ora_promhypo_go, ora_promhypo_go_csv_path)
    
    ora_promhypo_kegg <- missMethyl::gometh(sig.cpg = dmps_sig_hypo$CGID,
                                            all.cpg = dmps_comparison$CGID,
                                            collection = c("KEGG"),
                                            array.type = c("EPIC"),
                                            plot.bias = FALSE,
                                            prior.prob = TRUE,
                                            anno = epic_anno,
                                            equiv.cpg = TRUE,
                                            fract.counts = TRUE,
                                            genomic.features = c("TSS200", "TSS1500", "1stExon"),
                                            sig.genes = FALSE)
    ora_promhypo_kegg <- ora_promhypo_kegg[order(ora_promhypo_kegg$P.DE),]
    write.csv(ora_promhypo_kegg, ora_promhypo_kegg_csv_path)
    } else{
      file.create(ora_promhypo_go_csv_path)
      file.create(ora_promhypo_kegg_csv_path)
    }
} else{
  file.create(ora_go_csv_path)
  file.create(ora_kegg_csv_path)
}

sessionInfo()