#!/usr/bin/env Rscript
# The goal of this script is to import and preprocess the RNA-sequencing data counts.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

require(dplyr)
require(DEseq2)
require(ggplot2)

dds_path <- args[1] #"output/rnaseq/deseq2/dds.Rds"
rlog_path <- args[2] #"output/rnaseq/deseq2/rld.Rds"
degs_path <- args[3] #"output/rnaseq/deseq2/degs.csv"

dds <- readRDS(dds_path)
rld <- readRDS(rlog_path)

colData(dds)$Smoking_status_baseline <- gsub("( |\\-)", "", colData(dds)$Smoking_status_baseline)

colData(dds)$Response <- factor(colData(dds)$Response, levels = c("NR", "R"))
colData(dds)$Sex <- factor(colData(dds)$Sex)
colData(dds)$Smoking_status_baseline <- factor(colData(dds)$Smoking_status_baseline, levels = c("Neversmoked", "Exsmoker", "Currentsmoker"))

design(dds) <- formula(~Response+Age+Sex+Smoking_status_baseline)

dds <- DESeq(dds)
degs <- results(dds, name="Response_R_vs_NR") %>%
  data.frame(., ENSGv = rownames(.)) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::left_join(data.frame(rowData(dds)[,1:9]), by = "ENSGv")

write.csv(degs, degs_path)

degs_sig <- degs %>%
  dplyr::filter(padj<0.05)

for(i in 1:length(degs_sig)){
  ggplotobj <- data.frame(Expr = assay(rld)[degs_sig[i,"ENSGv"],],
                          Response = colData(rld)$Response,
                          Timepoint = colData(rld)$Timepoint,
                          Donor = colData(rld)$DonorID) %>%
    ggplot(aes(x = Timepoint, y = Expr, col = Response)) +
    stat_summary(fun.data = mean_se, geom = "crossbar",  aes(color = Response)) +
    geom_point(aes(group = Donor, shape = Response), alpha = 0.45, show.legend = F, position = position_dodge(0.1)) +
    geom_line(aes(group = Donor, alpha = 0.45), linetype = "dotted", show.legend = F, position = position_dodge(0.1)) +
    labs(title = ifelse(is.na(degs_sig[i,"HGNC"]), "", degs_sig[i,"HGNC"]),
         subtitle = paste0(degs_sig[i,"ENSGv"], "\np-value = ", formatC(x = degs_sig[i,"pvalue"], digits = 3, format = "e")), y = bquote(log[2]~"(counts+1)")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          text = element_text(size = 12))
}

sessionInfo()
