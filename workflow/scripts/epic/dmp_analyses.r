#!/usr/bin/env R
# Perform DMP analyses

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
  treatment == "Adalimumab" ~ "ADA_",
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
  treatment == "Infliximab" ~ "IFX_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")

betas <- minfi::getBeta(gmset)
mvals <- minfi::getM(gmset)

sample_metadata <- pData(gmset) %>%
  data.frame() %>%
  dplyr::mutate(resptime = paste0(!!sym(response_column), "_", !!sym(timepoint_column)),
                DNA_plate = gsub("(-| )", "_", DNA_plate),
                DNA_plate = gsub("\\[\\]", "", DNA_plate))

design_mat <- model.matrix(~0 + resptime + Sex + Age + DNA_plate, data = sample_metadata)

colnames(design_mat)[1:4] <- c("T1_NR", "T2_NR", "T1_R", "T2_R")

dupcor <- duplicateCorrelation(mvals, 
                               design = design_mat, 
                               block = sample_metadata$DonorID)

mvals_lmfit <- lmFit(mvals, 
                     design = design_mat, 
                     block = sample_metadata$DonorID, 
                     correlation = dupcor$consensus)
betas_lmfit <- lmFit(betas, 
                     design = design_mat)

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

mvals_contrastfit <- contrasts.fit(mvals_lmfit, contrast_mat)
mvals_ebayes <- eBayes(mvals_contrastfit)

betas_contrastfit <- contrasts.fit(betas_lmfit, contrast_mat)
betas_ebayes <- eBayes(betas_contrastfit)

dmps <- data.frame(topTable(mvals_ebayes, coef = "T1RvNR", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
  dplyr::rename(Mdiff = logFC) %>%
  dplyr::rename_with(function(cname){paste0(cname, "_T1RvNR")}) %>%
  tibble::rownames_to_column(var = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_ebayes, coef = "T2RvNR", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_T2RvNR")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_ebayes, coef = "RvNR", number = "Inf", adjust.method = "BH")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_RvNR")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_ebayes, coef = "RvNRvT2vT1", number = "Inf", adjust.method = "BH")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_RvNRvT2vT1")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_ebayes, coef = "RT2vT1", number = "Inf", adjust.method = "BH")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_RT2vT1")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(data.frame(coefficients(betas_ebayes)) %>%
                     dplyr::rename_with(function(cname){paste0("Betadiff_", cname)}) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID")

write.csv(dmps, dmps_path)

sessionInfo()