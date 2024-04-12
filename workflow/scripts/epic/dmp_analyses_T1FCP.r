#!/usr/bin/env R
# Perform DMP analyses comparing FCP at T1.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

library(statmod)
library(minfi)
library(limma)
library(dplyr)

gmset_path <- args[1]
dmps_path <- args[2]
treatment <- args[3]

gmset <- readRDS(gmset_path)

prefix <- case_when(
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")

betas <- minfi::getBeta(gmset)
mvals <- minfi::getM(gmset)

sample_metadata_t1 <- pData(gmset) %>%
  data.frame() %>%
  dplyr::mutate(DNA_plate = gsub("(-| )", "_", DNA_plate),
                DNA_plate = gsub("\\[\\]", "", DNA_plate),
                FCP = as.numeric(Calprotectin)) %>%
  dplyr::filter(!!sym(timepoint_column) == "T1")

gmset_t1 <- gmset[,rownames(sample_metadata_t1)]

betas_t1 <- minfi::getBeta(gmset_t1)
mvals_t1 <- minfi::getM(gmset_t1)

design_mat_t1 <- model.matrix(~FCP, data = sample_metadata_t1)

mvals_t1_lmfit <- lmFit(mvals_t1[,rownames(design_mat_t1)], 
                        design = design_mat_t1)
betas_t1_lmfit <- lmFit(betas_t1[,rownames(design_mat_t1)], 
                        design = design_mat_t1)

mvals_ebayes <- eBayes(mvals_t1_lmfit)

betas_ebayes <- eBayes(betas_t1_lmfit)

dmps_T1FCP <- data.frame(topTable(mvals_ebayes, coef = "FCP", number = "Inf", adjust.method = "BH", sort.by = "p")) %>%
  dplyr::rename(Mdiff = logFC) %>%
  tibble::rownames_to_column(var = "CGID") %>%
  dplyr::left_join(data.frame(coefficients(betas_ebayes)) %>%
                     dplyr::rename_with(function(cname){paste0("Betadiff_", cname)}) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID")

write.csv(dmps_T1FCP, dmps_path)

sessionInfo()