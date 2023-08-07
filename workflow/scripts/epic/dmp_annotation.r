#!/usr/bin/env R
# Append annotations to DMPs

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0("Script needs 3 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(dplyr))

dmps_csv_path <- args[1]
epic_annotation_csv <- args[2]
dmps_anno_csv_path <- args[3]

dmps <- read.csv(dmps_csv_path)[,-1]
epic_annotation <- read.csv(epic_annotation_csv)[,-1]

dmps_anno <- dmps %>%
  dplyr::left_join(epic_annotation,
                   by = c("CGID" = "Name"))

write.csv(dmps_anno, dmps_anno_csv_path)

sessionInfo()