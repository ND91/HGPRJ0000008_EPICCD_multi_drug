#!/usr/bin/env R
# Scatterplot of the t-statistic when comparing responders with non-responders at T1 on the x-axis and responders with non-responders at T2 on the y-axis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

library(tidyverse)
library(ggplot2)
library(ggrastr)

dmp_csv <- args[1]
horaizon_markers_xlsx <- args[2]
treatment <- args[3]
scatterplot_T1RvNR_T2RvNR_pdf <- args[4]
scatterplot_T1RvNR_T2RvNR_png <- args[5]

dmps <- read.csv(dmp_csv)

predictor_cpgs <- readxl::read_excel(horaizon_markers_xlsx) %>%
  dplyr::filter(Treatment == treatment)

dmps <- dmps %>%
  dplyr::mutate(predictor = ifelse(CGID %in% predictor_cpgs$CpG, CGID, NA))

predictor_cpg_stats <- dmps %>%
  dplyr::filter(!is.na(predictor))

predictor_cpg_stats_cor <- cor.test(predictor_cpg_stats$t_T1RvNR, predictor_cpg_stats$t_T2RvNR, method = "spearman")

plot_boundary <- ceiling(max(abs(dmps$t_T1RvNR), abs(dmps$t_T2RvNR)))

scatterplot_T1RvNR_T2RvNR_ggplotobj <- dmps %>%
  ggplot(aes(x = t_T1RvNR, y = t_T2RvNR)) +
  # annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10, alpha = .2, fill = "#AEC6CF") +
  # annotate("rect", xmin = 0, xmax = -10, ymin = 0, ymax = 10, alpha = .2, fill = "#77DD77") +
  # annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = -10, alpha = .2, fill = "#FFFAA0") +
  # annotate("rect", xmin = 0, xmax = -10, ymin = 0, ymax = -10, alpha = .2, fill = "#FAA0A0") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point_rast(col = "#808080") +
  geom_point_rast(data = dmps %>%
                    dplyr::filter(!is.na(predictor)),
                  col = "#FF0000") +
  geom_smooth(data = . %>%
                dplyr::filter(!is.na(predictor)),
              method=lm) +
  theme_bw() +
  coord_cartesian(xlim = c(-plot_boundary, plot_boundary),
                  ylim = c(-plot_boundary, plot_boundary)) +
  labs(title = paste0(treatment, " response-associated CpGs at T1 and T2"),
       subtitle = paste0("Spearman rho = ", round(predictor_cpg_stats_cor$estimate, 3)),
       x = "T1",
       y = "T2")

pdf(width = 7, height = 7, file = scatterplot_T1RvNR_T2RvNR_pdf)
print(scatterplot_T1RvNR_T2RvNR_ggplotobj)
dev.off()

png(width = 450, height = 450, file = scatterplot_T1RvNR_T2RvNR_png)
print(scatterplot_T1RvNR_T2RvNR_ggplotobj)
dev.off()

sessionInfo()