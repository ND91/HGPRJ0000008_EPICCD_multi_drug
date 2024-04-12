#!/usr/bin/env Rscript
# The goal of this script is to perform differential expression analyses.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

require(dplyr)
require(DESeq2)

dds_rds <- args[1] 
treatment <- args[2]
degs_csv <- args[3]

dds <- readRDS(dds_rds)

prefix <- case_when(
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")

colData(dds)$Smoking_status_baseline <- gsub("( |\\-)", "", colData(dds)$Smoking_status_baseline)
colData(dds)$Timepoint <- factor(colData(dds)[,timepoint_column], levels = c("T1", "T2"))
colData(dds)$Response <- factor(colData(dds)[,response_column], levels = c("NR", "R"))
colData(dds)$DonorID <- factor(colData(dds)$DonorID)
colData(dds)$Sex <- factor(colData(dds)$Sex)
colData(dds)$Smoking_status_baseline <- factor(colData(dds)$Smoking_status_baseline, levels = c("Neversmoked", "Exsmoker", "Currentsmoker"))

# Trick from https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#group-specific-condition-effects-individuals-nested-within-groups
rvnr_donorindex <- do.call(rbind, lapply(colData(dds) %>%
                                           data.frame() %>%
                                           dplyr::select(DonorID, Response) %>%
                                           split(., .$Response), function(response_group){
                                             response_group %>%
                                               dplyr::mutate(DonorIndex = paste0("D", as.numeric(as.factor(as.character(DonorID))))) %>%
                                               dplyr::select(DonorID, DonorIndex)
                                           })) %>%
  unique()


colData(dds)$DonorIndex <- factor(rvnr_donorindex[match(colData(dds)$DonorID, rvnr_donorindex$DonorID),"DonorIndex"])

# T1RvNR

dds_t1rvnr <- dds[,which(colData(dds)[,timepoint_column] == "T1")]

design(dds_t1rvnr) <- formula(~Response)

dds_t1rvnr <- DESeq(dds_t1rvnr)
degs_t1rvnr <- results(dds_t1rvnr, name="Response_R_vs_NR", independentFiltering = F) %>%
  data.frame() %>%
  dplyr::rename_with( ~ paste0(.x, "_T1RvNR")) %>%
  dplyr::mutate(ENSGv = rownames(.))

# T2RvNR

dds_t2rvnr <- dds[,which(colData(dds)[,timepoint_column] == "T2")]

design(dds_t2rvnr) <- formula(~Response)

dds_t2rvnr <- DESeq(dds_t2rvnr)
degs_t2rvnr <- results(dds_t2rvnr, name="Response_R_vs_NR", independentFiltering = F) %>%
  data.frame() %>%
  dplyr::rename_with( ~ paste0(.x, "_T2RvNR")) %>%
  dplyr::mutate(ENSGv = rownames(.))

# RvNR

dds_rvnr <- dds

design_mat <- model.matrix(~Response + Response:DonorIndex + Response:Timepoint, colData(dds_rvnr))
design_mat <- design_mat[,colSums(design_mat) > 0]

dds_rvnr <- DESeq(dds_rvnr, full = design_mat)
degs_rvnr <- results(dds_rvnr, name="ResponseR", independentFiltering = F) %>%
  data.frame() %>%
  dplyr::rename_with( ~ paste0(.x, "_RvNR")) %>%
  dplyr::mutate(ENSGv = rownames(.))

# RT2vT1

degs_rt2vt1 <- results(dds_rvnr, name="ResponseR.TimepointT2", independentFiltering = F) %>%
  data.frame() %>%
  dplyr::rename_with( ~ paste0(.x, "_RT2vT1")) %>%
  dplyr::mutate(ENSGv = rownames(.))

# NRT2vT1

degs_nrt2vt1 <- results(dds_rvnr, name="ResponseNR.TimepointT2", independentFiltering = F) %>%
  data.frame() %>%
  dplyr::rename_with( ~ paste0(.x, "_NRT2vT1")) %>%
  dplyr::mutate(ENSGv = rownames(.))

# T2vT1

dds_t2vt1 <- dds

design(dds_t2vt1) <- formula(~Response+Timepoint+Response:Timepoint)

dds_t2vt1 <- DESeq(dds_t2vt1)
degs_t2vt1 <- results(dds_t2vt1, name="Timepoint_T2_vs_T1", independentFiltering = F) %>%
  data.frame() %>%
  dplyr::rename_with( ~ paste0(.x, "_T2vT1")) %>%
  dplyr::mutate(ENSGv = rownames(.))

# RvNRvT2vT1

degs_rvnrt2vt1 <- results(dds_t2vt1, name="ResponseR.TimepointT2", independentFiltering = F) %>%
  data.frame() %>%
  dplyr::rename_with( ~ paste0(.x, "_RvNRvT2vT1")) %>%
  dplyr::mutate(ENSGv = rownames(.))

# Combine

degs <- degs_t1rvnr %>%
  dplyr::full_join(degs_t2rvnr, by = "ENSGv") %>%
  dplyr::full_join(degs_rvnr, by = "ENSGv") %>%
  dplyr::full_join(degs_rt2vt1, by = "ENSGv") %>%
  dplyr::full_join(degs_nrt2vt1, by = "ENSGv") %>%
  dplyr::full_join(degs_t2vt1, by = "ENSGv") %>%
  dplyr::full_join(degs_rvnrt2vt1, by = "ENSGv") %>%
  dplyr::left_join(data.frame(rowData(dds)[,c(1,7:9)]), by = "ENSGv")

write.csv(degs, degs_csv)

sessionInfo()
