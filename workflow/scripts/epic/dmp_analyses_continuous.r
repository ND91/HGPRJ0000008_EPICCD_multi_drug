#!/usr/bin/env R
# Perform DMP analyses on the decrease in inflammatory metrics continuous

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

# Paired samples

sample_metadata_inflammationscale <- pData(gmset) %>%
  data.frame() %>%
  dplyr::select(DonorID, !!sym(timepoint_column), CRP, Calprotectin, SES, HBI) %>%
  dplyr::mutate(CRP = as.numeric(CRP),
                Calprotectin = as.numeric(Calprotectin),
                SES = as.numeric(SES),
                HBI = as.numeric(HBI)) %>%
  tidyr::pivot_wider(names_from = !!sym(timepoint_column), values_from = c(CRP, Calprotectin, SES, HBI)) %>%
  dplyr::mutate(CRP_delta = CRP_T2-CRP_T1,
                Calprotectin_delta = Calprotectin_T2-Calprotectin_T1,
                SES_delta = SES_T2-SES_T1,
                HBI_delta = HBI_T2-HBI_T1)

samples_t1 <- pData(gmset) %>%
  data.frame() %>%
  dplyr::filter(!!sym(timepoint_column) == "T1")

gmset_t1 <- gmset[,rownames(samples_t1)]

betas_t1 <- getBeta(gmset_t1)
mvals_t1 <- getM(gmset_t1)

sample_metadata_t1 <- pData(gmset_t1) %>%
  data.frame() %>%
  dplyr::left_join(sample_metadata_inflammationscale, by = "DonorID") %>%
  dplyr::mutate(logCRP_delta = log10(CRP_delta+abs(min(CRP_delta, na.rm = T)+0.1)), #We do this as the CRP levels range from a scale of -150 to 50 with a strong cluster at 0.
                logCalprotectin_delta = log10(Calprotectin_delta+abs(min(Calprotectin_delta, na.rm = T)+0.1)), #We do this as the CRP levels range from a scale of -6000 to 4500 with a strong cluster at 0.
                resptime = paste0(!!sym(response_column), "_", !!sym(timepoint_column)),
                DNA_plate = gsub("(-| )", "_", DNA_plate),
                DNA_plate = gsub("\\[\\]", "", DNA_plate))

# SES

sample_metadata_t1_ses <- sample_metadata_t1 %>%
  dplyr::filter(!is.na(SES_delta))

design_mat_sesdelta <- model.matrix(~SES_delta + DNA_plate, data = sample_metadata_t1_ses)

mvals_t1_lmfit_sesdelta <- lmFit(mvals_t1[,sample_metadata_t1_ses$SXSPOS], 
                                 design = design_mat_sesdelta)
mvals_t1_ebayes_sesdelta <- eBayes(mvals_t1_lmfit_sesdelta)

betas_t1_lmfit_sesdelta <- lmFit(betas_t1[,sample_metadata_t1_ses$SXSPOS], 
                                 design = design_mat_sesdelta)
betas_t1_ebayes_sesdelta <- eBayes(betas_t1_lmfit_sesdelta)

# HBI

sample_metadata_t1_hbi <- sample_metadata_t1 %>%
  dplyr::filter(!is.na(HBI_delta))

design_mat_hbidelta <- model.matrix(~HBI_delta + DNA_plate, data = sample_metadata_t1_hbi)

mvals_t1_lmfit_hbidelta <- lmFit(mvals_t1[,sample_metadata_t1_hbi$SXSPOS], 
                                 design = design_mat_hbidelta)
mvals_t1_ebayes_hbidelta <- eBayes(mvals_t1_lmfit_hbidelta)

betas_t1_lmfit_hbidelta <- lmFit(betas_t1[,sample_metadata_t1_hbi$SXSPOS], 
                                 design = design_mat_hbidelta)
betas_t1_ebayes_hbidelta <- eBayes(betas_t1_lmfit_hbidelta)

# CRP

sample_metadata_t1_crp <- sample_metadata_t1 %>%
  dplyr::filter(!is.na(CRP_delta))

design_mat_crpdelta <- model.matrix(~CRP_delta + DNA_plate, data = sample_metadata_t1_crp)

mvals_t1_lmfit_crpdelta <- lmFit(mvals_t1[,sample_metadata_t1_crp$SXSPOS], 
                                 design = design_mat_crpdelta)
mvals_t1_ebayes_crpdelta <- eBayes(mvals_t1_lmfit_crpdelta)

betas_t1_lmfit_crpdelta <- lmFit(betas_t1[,sample_metadata_t1_crp$SXSPOS], 
                                 design = design_mat_crpdelta)
betas_t1_ebayes_crpdelta <- eBayes(betas_t1_lmfit_crpdelta)

# logCRP

sample_metadata_t1_logcrp <- sample_metadata_t1 %>%
  dplyr::filter(!is.na(logCRP_delta))

design_mat_logcrpdelta <- model.matrix(~logCRP_delta + DNA_plate, data = sample_metadata_t1_logcrp)

mvals_t1_lmfit_logcrpdelta <- lmFit(mvals_t1[,sample_metadata_t1_logcrp$SXSPOS], 
                                    design = design_mat_logcrpdelta)
mvals_t1_ebayes_logcrpdelta <- eBayes(mvals_t1_lmfit_logcrpdelta)

betas_t1_lmfit_logcrpdelta <- lmFit(betas_t1[,sample_metadata_t1_logcrp$SXSPOS], 
                                    design = design_mat_logcrpdelta)
betas_t1_ebayes_logcrpdelta <- eBayes(betas_t1_lmfit_logcrpdelta)

# Calprotectin

sample_metadata_t1_calprotectin <- sample_metadata_t1 %>%
  dplyr::filter(!is.na(Calprotectin_delta))

design_mat_calprotectindelta <- model.matrix(~Calprotectin_delta + DNA_plate, data = sample_metadata_t1_calprotectin)

mvals_t1_lmfit_calprotectindelta <- lmFit(mvals_t1[,sample_metadata_t1_calprotectin$SXSPOS], 
                                          design = design_mat_calprotectindelta)
mvals_t1_ebayes_calprotectindelta <- eBayes(mvals_t1_lmfit_calprotectindelta)

betas_t1_lmfit_calprotectindelta <- lmFit(betas_t1[,sample_metadata_t1_calprotectin$SXSPOS], 
                                          design = design_mat_calprotectindelta)
betas_t1_ebayes_calprotectindelta <- eBayes(betas_t1_lmfit_calprotectindelta)

# log Calprotectin

sample_metadata_t1_logcalprotectin <- sample_metadata_t1 %>%
  dplyr::filter(!is.na(logCalprotectin_delta))

design_mat_logcalprotectindelta <- model.matrix(~logCalprotectin_delta + DNA_plate, data = sample_metadata_t1_logcalprotectin)

mvals_t1_lmfit_logcalprotectindelta <- lmFit(mvals_t1[,sample_metadata_t1_logcalprotectin$SXSPOS], 
                                             design = design_mat_logcalprotectindelta)
mvals_t1_ebayes_logcalprotectindelta <- eBayes(mvals_t1_lmfit_logcalprotectindelta)

betas_t1_lmfit_logcalprotectindelta <- lmFit(betas_t1[,sample_metadata_t1_logcalprotectin$SXSPOS], 
                                             design = design_mat_logcalprotectindelta)
betas_t1_ebayes_logcalprotectindelta <- eBayes(betas_t1_lmfit_logcalprotectindelta)

# Collect results

dmps <- data.frame(topTable(mvals_t1_ebayes_sesdelta, coef = "SES_delta", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
  dplyr::rename(Mdiff = logFC) %>%
  dplyr::rename_with(function(cname){paste0(cname, "_dSES")}) %>%
  tibble::rownames_to_column(var = "CGID") %>%
  dplyr::left_join(coefficients(betas_t1_ebayes_sesdelta) %>%
                     data.frame() %>%
                     dplyr::select(SES_delta) %>%
                     dplyr::rename(Betadiff_dSES = SES_delta) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_t1_ebayes_hbidelta, coef = "HBI_delta", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_dHBI")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(coefficients(betas_t1_ebayes_hbidelta) %>%
                     data.frame() %>%
                     dplyr::select(HBI_delta) %>%
                     dplyr::rename(Betadiff_dHBI = HBI_delta) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_t1_ebayes_crpdelta, coef = "CRP_delta", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_dCRP")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(coefficients(betas_t1_ebayes_crpdelta) %>%
                     data.frame() %>%
                     dplyr::select(CRP_delta) %>%
                     dplyr::rename(Betadiff_dCRP = CRP_delta) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_t1_ebayes_logcrpdelta, coef = "logCRP_delta", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_dlogCRP")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(coefficients(betas_t1_ebayes_logcrpdelta) %>%
                     data.frame() %>%
                     dplyr::select(logCRP_delta) %>%
                     dplyr::rename(Betadiff_dlogCRP = logCRP_delta) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_t1_ebayes_calprotectindelta, coef = "Calprotectin_delta", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_dCalprotectin")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(coefficients(betas_t1_ebayes_calprotectindelta) %>%
                     data.frame() %>%
                     dplyr::select(Calprotectin_delta) %>%
                     dplyr::rename(Betadiff_dCalprotectin = Calprotectin_delta) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID") %>%
  dplyr::left_join(data.frame(topTable(mvals_t1_ebayes_logcalprotectindelta, coef = "logCalprotectin_delta", number = "Inf", adjust.method = "BH", sort.by = "none")) %>%
                     dplyr::rename(Mdiff = logFC) %>%
                     dplyr::rename_with(function(cname){paste0(cname, "_dlogCalprotectin")}) %>%
                     tibble::rownames_to_column(var = "CGID"),
                   by = "CGID") %>%
  dplyr::left_join(coefficients(betas_t1_ebayes_logcalprotectindelta) %>%
                     data.frame() %>%
                     dplyr::select(logCalprotectin_delta) %>%
                     dplyr::rename(Betadiff_dlogCalprotectin = logCalprotectin_delta) %>%
                     tibble::rownames_to_column(var = "CGID"), 
                   by = "CGID")

write.csv(dmps, dmps_path)

sessionInfo()