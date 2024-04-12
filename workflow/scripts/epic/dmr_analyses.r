#!/usr/bin/env R
# Perform DMR analysis

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

library(minfi)
library(limma)
library(DMRcate)
library(dplyr)

gmset_path <- args[1]
dmrcate_rds_path <- args[2]
treatment <- args[3]
comparison <- args[4]
dmr_csv_path <- args[5]

gmset <- readRDS(gmset_path)

prefix <- case_when(
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")

sample_metadata <- pData(gmset) %>%
  data.frame() %>%
  dplyr::mutate(resptime = paste0(!!sym(response_column), "_", !!sym(timepoint_column)),
                DNA_plate = gsub("(-| )", "_", DNA_plate),
                DNA_plate = gsub("\\[\\]", "", DNA_plate))

design_mat <- model.matrix(~0 + resptime + Sex + Age + DNA_plate, data = sample_metadata)

colnames(design_mat)[1:4] <- c("T1_NR", "T2_NR", "T1_R", "T2_R")

contrast_mat <- makeContrasts(
  T1RvNR = T1_R-T1_NR,
  T2RvNR = T2_R-T2_NR,
  RvNR = (T1_R+T2_R)/2-(T1_NR+T2_NR)/2,
  RT2vT1 = T2_R-T1_R,
  NRT2vT1 = T2_NR-T1_NR,
  T2vT1 = (T2_R+T2_NR)/2-(T1_R+T1_NR)/2,
  RvNRvT2vT1 = (T2_R-T1_R)-(T2_NR-T1_NR),
  levels = design_mat
)

set.seed(8657564)
dmrs <- cpg.annotate(object = gmset,
                     datatype = "array",
                     arraytype = "EPIC",
                     what = "M",
                     analysis.type = "differential",
                     design = design_mat,
                     contrasts = T,
                     cont.matrix = contrast_mat,
                     coef = comparison)
saveRDS(dmrs, dmrcate_rds_path)

tryCatch(expr = {
  dmrs <- dmrcate(dmrs, lambda = 1000, C = 2, min.cpgs = 3, pcutoff = 0.05)
  dmrs_gr <- extractRanges(dmrs, genome = "hg19")
  write.csv(as.data.frame(dmrs_gr), dmr_csv_path)
},
error = function(cond) {
  message("Failed with message:")
  message(cond)
})

sessionInfo()