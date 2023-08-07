#!/usr/bin/env R
# Prepare the data for machine learning at HorAIzon. 
# Normalize the data using functional normalization.
# Remove the batch effects (run, slide, and batch).
# Subset the data to T1, remove the SNPs (both catalogued as well as predicted), allosomes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

library(minfi)
library(sva)
library(dplyr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

rgset_qc_path <- args[1]
treatment <- args[2]
X_path <- args[3]
y_path <- args[4]

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

t1_discval_samples <- pData(rgset) %>%
  data.frame() %>%
  dplyr::filter(!!sym(timepoint_column) %in% c("T1"),
                !!sym(cohort_column) %in% c("EPIC-CD Discovery", "EPIC-CD Validation"),
                Center_source == "AmsterdamUMC") %>%
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
                         batch = gsub(" ", "_", paste0(pData(gmset)$Center_source, "_", pData(gmset)[,cohort_column])), 
                         mod = NULL)
betas_combat_cs <- ComBat(dat = betas_combat_c, 
                          batch = gsub("([0-9]+)_R[0-9]{2}C[0-9]{2}", "\\1", pData(gmset)$SXSPOS), 
                          mod = NULL)
betas_combat_csp <- ComBat(dat = betas_combat_cs, 
                           batch = gsub("[0-9]+_(R[0-9]{2}C[0-9]{2})", "\\1", pData(gmset)$SXSPOS), 
                           mod = NULL)

# Remove probes

## Identify catalogued GVs (only in CpG)

snp_cpgs <- minfi::getSnpInfo(gmset) %>%
  data.frame() %>%
  dplyr::filter(CpG_maf > 0 | SBE_maf > 0)

## Identify predicted GVs

gmset_gh <- minfi::gaphunter(gmset, threshold = 0.25)
gap_cpgs <- gmset_gh$proberesults

## Identify allosomal probes

allosome_cpgs <- probe_annotations %>%
  data.frame() %>%
  dplyr::filter(chr %in% c("chrX", "chrY"))

## Filter out aforementioned probes

cpgs_remove <- unique(c(rownames(snp_cpgs), rownames(gap_cpgs), rownames(allosome_cpgs)))

betas_combat_csp_clean <- betas_combat_csp[-which(rownames(betas_combat_csp) %in% cpgs_remove),]
gmset_clean <- gmset[rownames(gmset) %in% rownames(betas_combat_csp_clean),]

# Prepare annotations

## Feature annotations

genenames <- unlist(lapply(lapply(strsplit(minfi::getAnnotation(gmset_clean)$UCSC_RefGene_Name, ";"), unique),
                           function(genes){
                             paste(genes, collapse = ".")
                           })
)

featurenames <- unlist(lapply(lapply(strsplit(minfi::getAnnotation(gmset_clean)$UCSC_RefGene_Group, ";"), unique),
                              function(gfeats){
                                paste(gfeats, collapse = ".")
                              })
)

coordinates <- with(getAnnotation(gmset_clean), paste0(chr, ".", pos))

cpg_annotation <- paste0(rownames(gmset_clean), "_", coordinates, "_", genenames, "_", featurenames)
cpg_annotation <- gsub("(^.+?)_+$", "\\1", cpg_annotation)

betas_combat_csp_clean_t <- t(betas_combat_csp_clean)
colnames(betas_combat_csp_clean_t) <- cpg_annotation
rownames(betas_combat_csp_clean_t) <- colnames(betas_combat_csp_clean)

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
data.table::fwrite(data.frame(betas_combat_csp_clean_t), X_path, row.names = T)
write.csv(sample_metadata, y_path)