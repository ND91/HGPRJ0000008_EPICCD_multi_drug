#!/usr/bin/env Rscript
# The goal of this script is to import and preprocess the RNA-sequencing data counts.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

require(dplyr)
require(DESeq2)

dds_rds <- args[1] #"output/rnaseq/deseq2/dds.Rds"
treatment <- args[2]
dds_subset_rds <- args[3] #"output/rnaseq/deseq2/dds_subset.Rds"
rlog_subset_rds <- args[4] #"output/rnaseq/deseq2/rld_subset.Rds"

# Import full dds object
dds <- readRDS(dds_rds)

prefix <- case_when(
  treatment == "Adalimumab" ~ "ADA_",
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
  treatment == "Infliximab" ~ "IFX_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")

selected_samples <- colData(dds) %>%
  data.frame() %>%
  dplyr::filter(!!sym(timepoint_column) %in% c("T1", "T2")) %>%
  dplyr::pull(SampleID)

# t1t2_donors <- colData(dds) %>%
#   data.frame() %>%
#   dplyr::filter(!!sym(timepoint_column) %in% c("T1", "T2")) %>%
#   dplyr::group_by(DonorID) %>%
#   summarize(n = n()) %>%
#   dplyr::filter(n == 2,
#                 !is.na(DonorID)) %>%
#   dplyr::pull(DonorID)
# 
# t1t2_samples <- colData(dds) %>%
#   data.frame() %>%
#   dplyr::filter(DonorID %in% t1t2_donors,
#                 !!sym(timepoint_column) %in% c("T1", "T2")) %>%
#   tibble::rownames_to_column(var = "arrayname") %>%
#   dplyr::pull(arrayname)

dds_subset <- dds[,colnames(dds) %in% selected_samples]

saveRDS(object = dds_subset, file = dds_subset_rds, compress = "gzip")

rld_subset <- rlog(dds_subset)

saveRDS(object = rld_subset, file = rlog_subset_rds, compress = "gzip")

sessionInfo()