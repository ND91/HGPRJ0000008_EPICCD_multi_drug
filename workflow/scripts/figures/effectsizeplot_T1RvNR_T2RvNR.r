#!/usr/bin/env R
# Scatterplot of the t-statistic when comparing responders with non-responders at T1 on the x-axis and responders with non-responders at T2 on the y-axis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

dmp_csv <- args[1]
horaizon_markers <- args[2]
treatment <- args[3]
scatterplot_T1RvNR_T2RvNR_pdf <- args[4]

dmps <- read.csv(dmp_csv)

dmps <- dmps %>%
  dplyr::mutate(predictor = ifelse(CGID %in% predictor_cpgs$CGID, CGID, NA))

predictor_cpg_stats <- tofa_dmps_RvNR_int %>%
  dplyr::filter(!is.na(predictor))

predictor_cpg_stats_cor <- cor.test(predictor_cpg_stats$t_RvNRT1, predictor_cpg_stats$t_RvNRT2, method = "spearman")

pdf(width = 7, height = 7, file = scatterplot_T1RvNR_T2RvNR_pdf)
ggplot(dmps, aes(x = t_RvNRT1, y = t_RvNRT2)) +
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
  geom_smooth(data = tofa_dmps_RvNR_int %>%
                dplyr::filter(!is.na(predictor)),
              method=lm) +
  theme_bw() +
  coord_cartesian(xlim = c(-6,6),
                  ylim = c(-6,6)) +
  labs(title = "Temporal stability response-associated CpGs between w0 and w8",
       subtitle = paste0("Spearman rho = ", round(predictor_cpg_stats_cor$estimate, 3)),
       x = "T1",
       y = "T2")
dev.off()
