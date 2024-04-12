#!/usr/bin/env R
# Prepare the validation data machine learning at HorAIzon. 
# The goal is to train on t1 from AmsterdamUMC (same as the training set) and generate predictions on Oxford (and Exeter) validation data.
# Normalize the data using functional normalization.
# Remove the batch effects (run, slide, and batch).

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

## ComBat normalization

betas_combat_c <- ComBat(dat = betas, 
                         batch = paste0(pData(gmset)$Center_source, "_", pData(gmset)[,cohort_column]), 
                         mod = NULL)
betas_combat_cs <- ComBat(dat = betas_combat_c, 
                          batch = gsub("([0-9]+)_R[0-9]{2}C[0-9]{2}", "\\1", pData(gmset)$SXSPOS), 
                          mod = NULL)
betas_combat_csp <- ComBat(dat = betas_combat_cs, 
                           batch = gsub("[0-9]+_(R[0-9]{2}C[0-9]{2})", "\\1", pData(gmset)$SXSPOS), 
                           mod = NULL)

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

betas_combat_csp_t <- t(betas_combat_csp)
colnames(betas_combat_csp_t) <- cpg_annotation
rownames(betas_combat_csp_t) <- colnames(betas_combat_csp)

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
data.table::fwrite(data.frame(betas_combat_csp_t), X_path, row.names = T)
write.csv(sample_metadata, y_path)

sessionInfo()