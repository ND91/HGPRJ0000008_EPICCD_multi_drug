#!/usr/bin/env R
# Prepare the discovery + validation data regressed for the confounders for validation at HorAIzon. 
# The goal is to train on the discovery and generate predictions on the validation data. 
# The confounders have been regressed out using sva::combat for categorical variables and limma:removeBatchEffects for continuous variables.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

library(minfi)
library(sva)
library(dplyr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

rgset_path <- args[1]
treatment <- args[2]
X_path <- args[3]
y_path <- args[4]

rgset <- readRDS(rgset_path)

prefix <- case_when(
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")

t1_discval_samples <- pData(rgset) %>%
  data.frame() %>%
  dplyr::filter(!!sym(timepoint_column) %in% c("T1"),
                !!sym(cohort_column) %in% c("EPIC-CD Discovery", "EPIC-CD Validation"),
                !is.na(!!sym(response_column)),
                Disease == "CD") %>%
  tibble::rownames_to_column(var = "arrayname") %>%
  dplyr::pull(arrayname)

rgset_t1_discval <- rgset[,colnames(rgset) %in% t1_discval_samples]

# Normalization

## Functional normalization

gmset <- preprocessFunnorm(rgSet = rgset_t1_discval)

probe_annotations <- minfi::getAnnotation(gmset)

betas <- minfi::getBeta(gmset)

## ComBat batch effect removal

betas_combat_c <- ComBat(dat = betas, 
                         batch = paste0(pData(gmset)$Center_source, "_", pData(gmset)[,cohort_column]), 
                         mod = NULL)
betas_combat_cs <- ComBat(dat = betas_combat_c, 
                          batch = gsub("([0-9]+)_R[0-9]{2}C[0-9]{2}", "\\1", pData(gmset)$SXSPOS), 
                          mod = NULL)
betas_combat_csp <- ComBat(dat = betas_combat_cs, 
                           batch = gsub("[0-9]+_(R[0-9]{2}C[0-9]{2})", "\\1", pData(gmset)$SXSPOS), 
                           mod = NULL)

### ComBat removal of categorical variables sex and smoker status

betas_combat_csp_nosex <- ComBat(dat = betas_combat_csp, 
                                 batch = pData(gmset)$predictedSex, 
                                 mod = NULL)

pData(gmset)$Smoking_status_baseline[is.na(pData(gmset)$Smoking_status_baseline)] <- "Unknown"

betas_combat_csp_nosex_nosmoke <- ComBat(dat = betas_combat_csp_nosex, 
                                         batch = pData(gmset)$Smoking_status_baseline, 
                                         mod = NULL)

## Regress sex, age, smoker status, and blood cell distribution out

betas_combat_csp_nosex_nosmoke_regressedconfounders <- limma::removeBatchEffect(betas_combat_csp_nosex_nosmoke, 
                                                                                covariates = data.frame(Age = pData(gmset)$Age,
                                                                                                        CD8T = pData(gmset)$CD8T,
                                                                                                        CD4T = pData(gmset)$CD4T,
                                                                                                        NK = pData(gmset)$NK,
                                                                                                        Bcell = pData(gmset)$Bcell,
                                                                                                        Mono = pData(gmset)$Mono,
                                                                                                        Neu = pData(gmset)$Neu))

# Prepare annotations

## Feature annotations

genenames <- unlist(lapply(lapply(strsplit(minfi::getAnnotation(gmset)$UCSC_RefGene_Name, ";"), unique),
                           function(genes){
                             paste(genes, collapse = ".")
                           })
)

featurenames <- unlist(lapply(lapply(strsplit(minfi::getAnnotation(gmset)$UCSC_RefGene_Group, ";"), unique),
                              function(gfeats){
                                paste(gfeats, collapse = ".")
                              })
)

coordinates <- with(getAnnotation(gmset), paste0(chr, ".", pos))

cpg_annotation <- paste0(rownames(gmset), "_", coordinates, "_", genenames, "_", featurenames)
cpg_annotation <- gsub("(^.+?)_+$", "\\1", cpg_annotation)

betas_combat_csp_nosex_nosmoke_regressedconfounders_t <- t(betas_combat_csp_nosex_nosmoke_regressedconfounders)
colnames(betas_combat_csp_nosex_nosmoke_regressedconfounders_t) <- cpg_annotation
rownames(betas_combat_csp_nosex_nosmoke_regressedconfounders_t) <- colnames(betas_combat_csp_nosex_nosmoke_regressedconfounders)

## Sample annotation

sample_metadata <- pData(gmset) %>%
  data.frame() %>%
  dplyr::select(SampleID,
                !!sym(response_column),
                !!sym(cohort_column),
                !!sym(timepoint_column),
                Center_source,
                Disease) %>%
  dplyr::mutate(Treatment = treatment) %>%
  dplyr::rename(Response = !!sym(response_column),
                Timepoint = !!sym(timepoint_column),
                Cohort = !!sym(cohort_column))
# Save data
data.table::fwrite(data.frame(betas_combat_csp_nosex_nosmoke_regressedconfounders_t), X_path, row.names = T)
write.csv(sample_metadata, y_path)