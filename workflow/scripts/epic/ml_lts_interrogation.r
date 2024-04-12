#!/usr/bin/env R
# Extract the ICC values for the predictor CpGs from Joustra et al. 2022 (DOI: 10.1016/j.jcmgh.2022.12.011).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

library(ggplot2)
library(dplyr)

horaizon_markers_path <- args[1]
cdpbmethlts_path <- args[2]
horaizon_markers_cdpbmethlts_ICC_csv <- args[3]
horaizon_markers_cdpbmethlts_ICC_pdf <- args[4]

horaizon_markers <- readxl::read_excel(horaizon_markers_path)

cdpbmethlts <- read.csv(cdpbmethlts_path)

horaizon_markers_icc <- horaizon_markers %>% 
  dplyr::filter(Disease == "CD",
                Cohort %in% c("Discovery", "Discovery + Validation")) %>%
  dplyr::left_join(cdpbmethlts, by = c("CGID" = "Name")) 

pdf(width = 8, height = 8, file = horaizon_markers_cdpbmethlts_ICC_pdf)
ggplot(horaizon_markers_icc, aes(x = ICC, y = Biological)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.1) +
  geom_point() +
  geom_vline(xintercept = c(0.5, 0.75, 0.9)) +
  facet_wrap(~Cohort, nrow = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())
dev.off()

# Save data
write.csv(data.frame(horaizon_markers_icc), horaizon_markers_cdpbmethlts_ICC_csv, row.names = T)