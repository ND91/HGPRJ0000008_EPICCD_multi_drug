#!/usr/bin/env R
# Perform gene set overrepresentation analysis on horaizon probes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

library(minfi)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(missMethyl)

horaizon_predictor_probes_xlsx <- args[1]
ora_go_csv <- args[2]
treatment <- args[3]

horaizon_predictor_probes <- readxl::read_excel(horaizon_predictor_probes_xlsx)

epic_anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

horaizon_predictor_probes_comparison <- horaizon_predictor_probes %>%
  dplyr::filter(Treatment == treatment)

ora_go <- missMethyl::gometh(sig.cpg = horaizon_predictor_probes_comparison$CpG,
                             all.cpg = epic_anno$CGID,
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

write.csv(ora_go, ora_go_csv)

sessionInfo()