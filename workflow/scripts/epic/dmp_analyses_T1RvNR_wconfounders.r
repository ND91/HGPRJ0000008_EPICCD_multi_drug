#!/usr/bin/env R
# Perform DMP analyses comparing R with NR at T1 correcting for age, sex, blood cell composition, and plate.

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

sample_metadata_t1 <- pData(gmset) %>%
  data.frame() %>%
  dplyr::mutate(DNA_plate = gsub("(-| )", "_", DNA_plate),
                DNA_plate = gsub("\\[\\]", "", DNA_plate),
                Response = !!sym(response_column),
                Smoking_status = case_when(
                  Smoking_status_baseline == "Active smoker" ~ "Active",
                  Smoking_status_baseline == "Ex-smoker" ~ "Ex",
                  Smoking_status_baseline == "Never smoked" ~ "Never")) %>%
  dplyr::filter(!!sym(timepoint_column) == "T1")

gmset_t1 <- gmset[,rownames(sample_metadata_t1)]

betas_t1 <- minfi::getBeta(gmset_t1)
mvals_t1 <- minfi::getM(gmset_t1)

design_mat_t1 <- model.matrix(~0 + Response + DNA_plate + Age + Sex + Smoking_status + CD8T + CD4T + NK + Bcell + Mono + Neu, data = sample_metadata_t1)

colnames(design_mat_t1)[1:2] <- c("NR", "R")

mvals_t1_lmfit <- lmFit(mvals_t1[,rownames(design_mat_t1)], 
                        design = design_mat_t1)
betas_t1_lmfit <- lmFit(betas_t1[,rownames(design_mat_t1)], 
                        design = design_mat_t1)

contrast_mat_t1_wconfounders <- makeContrasts(
  T1RvNR = R-NR,
  levels = design_mat_t1
)

mvals_contrastfit <- contrasts.fit(mvals_t1_lmfit, contrast_mat_t1_wconfounders)
mvals_ebayes <- eBayes(mvals_contrastfit)

betas_contrastfit <- contrasts.fit(betas_t1_lmfit, contrast_mat_t1_wconfounders)
betas_ebayes <- eBayes(betas_contrastfit)

dmps_T1RvNR <- data.frame(topTable(mvals_ebayes, coef = "T1RvNR", number = "Inf", adjust.method = "BH", sort.by = "p")) %>%
  dplyr::rename(Mdiff = logFC) %>%
  tibble::rownames_to_column(var = "CGID") %>%
  dplyr::left_join(data.frame(coefficients(betas_ebayes)) %>%
                     dplyr::rename_with(function(cname){paste0("Betadiff_", cname)}) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID")

write.csv(dmps_T1RvNR, dmps_path)

sessionInfo()