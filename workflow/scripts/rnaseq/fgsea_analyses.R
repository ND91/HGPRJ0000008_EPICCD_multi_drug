#!/usr/bin/env Rscript
# The goal of this script is to perform gene set enrichment analyses.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

require(dplyr)
require(DESeq2)
require(fgsea)
require(msigdbr)

degs_csv <- args[1] 
fgsea_go_list_rds <- args[2]
fgsea_go_csv <- args[3]
fgsea_kegg_list_rds <- args[4]
fgsea_kegg_csv <- args[5]

degs <- read.csv(degs_csv)

comparisons <- c("T1RvNR", "T2RvNR", "RvNR", "RT2vT1", "NRT2vT1", "T2vT1", "RvNRvT2vT1")

# GO
go_gs <- msigdbr(species = "human", category = "C5") %>%
  dplyr::filter(gs_subcat %in% c("GO:BP", "GO:CC", "GO:MF"))
go_gs_list <- lapply(split(go_gs, go_gs$gs_name), function(gs){
  gs$ensembl_gene
})

fgsea_go_list <- lapply(comparisons, function(comparison){
  
  stat_comparison_column <- paste0("stat_", comparison) 
  
  degs_comparison <- degs %>%
    dplyr::filter(!is.na(!!sym(stat_comparison_column)))
  stat_comparison <- degs_comparison[,stat_comparison_column]
  names(stat_comparison) <- degs_comparison$ENSG
  
  fgsea_go_comparison <- fgsea(pathways = go_gs_list, 
                               stats = stat_comparison, 
                               minSize = 15)
  fgsea_go_comparison <- fgsea_go_comparison[order(fgsea_go_comparison$pval),]
  return(fgsea_go_comparison)
})
names(fgsea_go_list) <- comparisons

fgsea_go_df <- fgsea_go_list$T1RvNR[,c(1:7)] %>%
  dplyr::rename_with( ~ paste0(.x, "_T1RvNR")) %>%
  dplyr::rename(pathway = 1) %>%
  dplyr::full_join(fgsea_go_list$T2RvNR[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_T2RvNR")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_go_list$RvNR[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_RvNR")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_go_list$RT2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_RT2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_go_list$NRT2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_NRT2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_go_list$T2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_T2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_go_list$RvNRvT2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_RvNRvT2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway")%>%
  data.frame()

# KEGG

kegg_gs <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
kegg_gs_list <- lapply(split(kegg_gs, kegg_gs$gs_name), function(gs){
  gs$ensembl_gene
})

fgsea_kegg_list <- lapply(comparisons, function(comparison){
  
  stat_comparison_column <- paste0("stat_", comparison) 
  
  degs_comparison <- degs %>%
    dplyr::filter(!is.na(!!sym(stat_comparison_column)))
  stat_comparison <- degs_comparison[,stat_comparison_column]
  names(stat_comparison) <- degs_comparison$ENSG
  
  fgsea_kegg_comparison <- fgsea(pathways = kegg_gs_list, 
                               stats = stat_comparison, 
                               minSize = 15)
  fgsea_kegg_comparison <- fgsea_kegg_comparison[order(fgsea_kegg_comparison$pval),]
  return(fgsea_kegg_comparison)
})
names(fgsea_kegg_list) <- comparisons

fgsea_kegg_df <- fgsea_kegg_list$T1RvNR[,c(1:7)] %>%
  dplyr::rename_with( ~ paste0(.x, "_T1RvNR")) %>%
  dplyr::rename(pathway = 1) %>%
  dplyr::full_join(fgsea_kegg_list$T2RvNR[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_T2RvNR")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_kegg_list$RvNR[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_RvNR")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_kegg_list$RT2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_RT2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_kegg_list$NRT2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_NRT2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_kegg_list$T2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_T2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway") %>%
  dplyr::full_join(fgsea_kegg_list$RvNRvT2vT1[,c(1:7)] %>%
                     dplyr::rename_with( ~ paste0(.x, "_RvNRvT2vT1")) %>%
                     dplyr::rename(pathway = 1),
                   by = "pathway")%>%
  data.frame()

#Save
write.csv(fgsea_go_df, fgsea_go_csv)
write.csv(fgsea_kegg_df, fgsea_kegg_csv)
saveRDS(fgsea_go_list, fgsea_go_list_rds, compress = "gzip")
saveRDS(fgsea_kegg_list, fgsea_kegg_list_rds, compress = "gzip")

sessionInfo()