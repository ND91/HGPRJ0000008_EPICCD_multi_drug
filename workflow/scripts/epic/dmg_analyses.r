#!/usr/bin/env R
# Perform DMG analysis

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(EmpiricalBrownsMethod))

dmps_csv_path <- args[1]
gmset_path <- args[2]
dmgs_csv_path <- args[3]
comparison <- args[4]

dmps <- read.csv(dmps_csv_path)[,-1]
pvals <- dmps[,c(1, grep("P.Value", colnames(dmps)))]

gmset <- readRDS(gmset_path)
betas <- minfi::getBeta(gmset)

dmps$genes <- paste0(dmps$UCSC_RefGene_Name, ";", dmps$GencodeBasicV12_NAME, ";", dmps$Enhancer_gene)
dmp_genes <- lapply(strsplit(dmps$genes, ";"), unique)
names(dmp_genes) <- dmps$CGID

# Remove unannotated CpGs
dmp_genes <- lapply(dmp_genes, function(cgid){
  cgid[-which(cgid %in% c("", "NA"))]
})
dmp_genes_nz <- dmp_genes[which(lapply(dmp_genes, length) != 0)]

dmp_genes_nz_df <- data.frame(CpG = gsub("(^cg[0-9]{8}).+$", "\\1", names(unlist(dmp_genes_nz))), 
                              Gene = unlist(dmp_genes_nz)) %>%
  dplyr::left_join(pvals, by = c("CpG" = "CGID"))

dmps_comparison <- dmp_genes_nz_df[,c(1, 2, which(colnames(dmps) %in% paste0("P.Value_", comparison)))] %>%
  dplyr::rename(pvalue = 3)

dmps_comparison_list <- split(dmps_comparison, dmps_comparison$Gene)
dmps_comparison_list <- dmps_comparison_list[unlist(lapply(dmps_comparison_list, nrow)) != 1]

dmp_genes_nz_fisher <- lapply(X = dmps_comparison_list, FUN = function(gene){
  unlist(EmpiricalBrownsMethod::empiricalBrownsMethod(data_matrix = betas[which(rownames(betas) %in% gene$CpG),],
                                                      p_values = gene$pvalue,
                                                      extra_info = T))
})

rm(dmps_comparison_list)

dmp_ma <- data.frame(Gene = names(dmp_genes_nz_fisher), do.call(rbind, dmp_genes_nz_fisher))
colnames(dmp_ma) <- c("Gene", "Brown_pval", "Fisher_pval", "Scale_Factor_C", "DF")
dmp_ma$Brown_padj <- p.adjust(dmp_ma$Brown_pval, method = "BH")
dmp_ma <- dmp_ma[order(dmp_ma$Brown_pval),]

write.csv(dmp_ma, dmgs_csv_path)

sessionInfo()
