#!/usr/bin/env R
# The goal is to train a model using DNA methylation data from the T1 discovery data (collected at the AmsterdamUMC) and to predict response in the T1 ancillary validation cohort (collected at the AmsterdamUMC and John Radcliffe Hospital). 
# Prepare the validation data for HorAIzon. 
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

rgset_rds <- args[1]
treatment <- args[2]
X_path <- args[3]
y_path <- args[4]

rgset <- readRDS(rgset_rds)

d1_prefix <- case_when(
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
)

d2_prefix <- case_when(
  treatment == "Vedolizumab" ~ "UST_",
  treatment == "Ustekinumab" ~ "VDZ_",
)

d1_cohort_column <- paste0(d1_prefix, "cohort")
d1_timepoint_column <- paste0(d1_prefix, "timepoint")
d1_response_column <- paste0(d1_prefix, "response")

d2_cohort_column <- paste0(d2_prefix, "cohort")
d2_timepoint_column <- paste0(d2_prefix, "timepoint")
d2_response_column <- paste0(d2_prefix, "response")

t1_d1disc_samples <- pData(rgset) %>%
  data.frame() %>%
  dplyr::filter(!!sym(d1_timepoint_column) %in% c("T1"),
                !!sym(d1_cohort_column) %in% c("EPIC-CD Discovery", "EPIC-CD Validation"),
                !is.na(!!sym(d1_response_column)),
                Disease == "CD",
                !(!!sym(d1_cohort_column) == "EPIC-CD Validation" & Center_source == "John Radcliffe Hospital")) %>%
  tibble::rownames_to_column(var = "arrayname")

t1_d2disc_samples <- pData(rgset) %>%
  data.frame() %>%
  dplyr::filter(!!sym(d2_timepoint_column) %in% c("T1"),
                !!sym(d2_cohort_column) %in% c("EPIC-CD Discovery", "EPIC-CD Validation"),
                !is.na(!!sym(d2_response_column)),
                Disease == "CD") %>%
  tibble::rownames_to_column(var = "arrayname")

t1_d1discd2discval_samples <- c(t1_d1disc_samples$arrayname, t1_d2disc_samples$arrayname)

rgset_t1_d1discd2discval <- rgset[,colnames(rgset) %in% t1_d1discd2discval_samples]

pData(rgset_t1_d1discd2discval)$Cohort <- c(paste0(d1_prefix, t1_d1disc_samples[,d1_cohort_column]),
                                            paste0(d2_prefix, t1_d2disc_samples[,d2_cohort_column]))
pData(rgset_t1_d1discd2discval)$Timepoint <- c(paste0(d1_prefix, t1_d1disc_samples[,d1_timepoint_column]),
                                               paste0(d2_prefix, t1_d2disc_samples[,d2_timepoint_column]))
pData(rgset_t1_d1discd2discval)$Response <- c(paste0(d1_prefix, t1_d1disc_samples[,d1_response_column]),
                                              paste0(d2_prefix, t1_d2disc_samples[,d2_response_column]))

# Normalization

## Functional normalization

gmset <- preprocessFunnorm(rgSet = rgset_t1_d1discd2discval)

probe_annotations <- minfi::getAnnotation(gmset)

betas <- minfi::getBeta(gmset)

## ComBat normalization

betas_combat_c <- ComBat(dat = betas, 
                         batch = paste0(pData(gmset)$Center_source, "_", pData(gmset)$Cohort), 
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
                Response,
                Cohort,
                Timepoint,
                Center_source,
                Disease) %>%
  dplyr::mutate(Treatment = treatment)

# Save data
data.table::fwrite(data.frame(betas_combat_csp_t), X_path, row.names = T)
write.csv(sample_metadata, y_path)