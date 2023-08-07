#!/usr/bin/env R
# Subset and normalize data

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

library(minfi)
library(dplyr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

getwd()
source("workflow/scripts/epic/functions.r")

rgset_qc_path <- args[1]
gmset_path <- args[2]
treatment <- args[3]

rgset <- readRDS(rgset_qc_path)

prefix <- case_when(
  treatment == "Adalimumab" ~ "ADA_",
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
  treatment == "Infliximab" ~ "IFX_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")

t1t2_donors <- pData(rgset) %>%
  data.frame() %>%
  dplyr::filter(!!sym(timepoint_column) %in% c("T1", "T2")) %>%
  dplyr::group_by(DonorID) %>%
  summarize(n = n()) %>%
  dplyr::filter(n == 2,
                !is.na(DonorID)) %>%
  dplyr::pull(DonorID)

t1t2_samples <- pData(rgset) %>%
  data.frame() %>%
  dplyr::filter(DonorID %in% t1t2_donors,
                !!sym(timepoint_column) %in% c("T1", "T2")) %>%
  tibble::rownames_to_column(var = "arrayname") %>%
  dplyr::pull(arrayname)

rgset_t1t2 <- rgset[,colnames(rgset) %in% t1t2_samples]

gmset <- minfi::preprocessFunnorm(rgset_t1t2)

saveRDS(gmset, gmset_path, compress = "gzip")

sessionInfo()