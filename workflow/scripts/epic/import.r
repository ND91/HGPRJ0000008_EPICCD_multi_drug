#!/usr/bin/env R
# Import data

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

# Reinstall preprocessCore as otherwise preprocessFunnorm will not work.
BiocManager::install("preprocessCore", configure.args="--disable-threading", force = T)

library(minfi)
library(dplyr)

source("workflow/scripts/epic/functions.r")

sample_metadata_xlsx <- args[1]
epic_files_xlsx <- args[2]
rgset_rds <- args[3]

#import samples
sample_metadata <- readxl::read_excel(sample_metadata_xlsx) %>%
  dplyr::inner_join(readxl::read_excel(epic_files_xlsx), by = "SampleID") %>%
  dplyr::rename(Basename = Filepath)

#prepare rgset
rgset <- read.metharray.exp(targets = sample_metadata, force = T)

saveRDS(rgset, rgset_rds, compress = "gzip")

sessionInfo()