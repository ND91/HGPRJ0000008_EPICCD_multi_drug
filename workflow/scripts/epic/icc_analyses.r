#!/usr/bin/env R
# Perform ICC analyses

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

library(minfi)
library(limma)
library(dplyr)
library(irr)

gmset_path <- args[1]
icc_path <- args[2]
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

t1t2_data <- data.frame(Timepoint = pData(gmset)[,timepoint_column],
                        DonorID = pData(gmset)$DonorID,
                        t(getBeta(gmset)))

icc_t1t2 <- foreach(i = 3:ncol(t1t2_data), .combine='rbind') %do%                     
  {
    iccscore <- t1t2_data[,c(1,2,i)] %>%
      tidyr::pivot_wider(names_from = Timepoint, values_from = 3) %>% 
      dplyr::select(T1, T2) %>%
      irr::icc(., 
               model = "twoway", 
               type = "consistency")
    
    data.frame(icc = iccscore$value,
               icc_lbound = iccscore$lbound,
               icc_rbound = iccscore$ubound,
               fvalue = iccscore$Fvalue,
               pvalue = iccscore$p.value,
               row.names = colnames(t1t2_data)[i])
  }

icc_t1t2$padj <- p.adjust(icc_t1t2$pvalue, method = "BH")

write.csv(icc_t1t2, icc_path)

sessionInfo()