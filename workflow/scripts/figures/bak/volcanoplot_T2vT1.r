#!/usr/bin/env R
# Volcanoplot of each treatment comparing responders over time

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

library(tidyverse)
library(ggplot2)
library(ggrastr)

dmp_csv <- args[1]
treatment <- args[2]
volcanoplot_RT2vT1_pdf <- args[3]
volcanoplot_NRT2vT1_pdf <- args[4]
scatterplot_RT2vT1_pdf <- args[5]
effectsizeplot_RT2vT1_NRT2vT1_pdf <- args[6]

dmps <- read.csv(dmps_csv)

# RT2vT1 volcanoplot

dmps_RT2vT1 <- dmps %>%
  dplyr::select(CGID, AveExpr_RT2vT1, Mdiff_RT2vT1, Betadiff_RT2vT1, t_RT2vT1, B_RT2vT1, P.Value_RT2vT1, adj.P.Val_RT2vT1) %>%
  dplyr::mutate(Significance = ifelse(adj.P.Val_RT2vT1<0.05, "Significant", "NS")) %>%
  dplyr::arrange(P.Value_RT2vT1)

volcanoplot_RT2vT1_ggplotobj <- dmps_RT2vT1 %>%
  ggplot(aes(x = Betadiff_RT2vT1, y = -log10(P.Value_RT2vT1), col = Significance)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Responders: T2 vs T1",
       x = "%Methylation difference",
       y = bquote('-'~log[10]~'(p-value)')) +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none")

dmps_RT2vT1 <- dmps %>%
  dplyr::select(CGID, AveExpr_RT2vT1, Mdiff_RT2vT1, Betadiff_RT2vT1, t_RT2vT1, B_RT2vT1, P.Value_RT2vT1, adj.P.Val_RT2vT1) %>%
  dplyr::mutate(Significance = ifelse(adj.P.Val_RT2vT1<0.05, "Significant", "NS")) %>%
  dplyr::arrange(P.Value_RT2vT1)

volcanoplot_RT2vT1_ggplotobj <- dmps_RT2vT1 %>%
  ggplot(aes(x = Betadiff_RT2vT1, y = -log10(P.Value_RT2vT1), col = Significance)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Responders: T2 vs T1",
       x = "%Methylation difference",
       y = bquote('-'~log[10]~'(p-value)')) +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none")

# NRT2vT1 volcanoplot

dmps_NRT2vT1 <- dmps %>%
  dplyr::select(CGID, Mdiff_NRT2vT1, Betadiff_NRT2vT1, t_NRT2vT1, B_NRT2vT1, P.Value_NRT2vT1, adj.P.Val_NRT2vT1) %>%
  dplyr::mutate(Significance = ifelse(adj.P.Val_NRT2vT1<0.05, "Significant", "NS")) %>%
  dplyr::arrange(P.Value_NRT2vT1)

volcanoplot_NRT2vT1_ggplotobj <- dmps_NRT2vT1 %>%
  ggplot(aes(x = Betadiff_NRT2vT1, y = -log10(P.Value_NRT2vT1), col = Significance)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Responders: T2 vs T1",
       x = "%Methylation difference",
       y = bquote('-'~log[10]~'(p-value)')) +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none")

dmps_NRT2vT1 <- dmps %>%
  dplyr::select(CGID, AveExpr_NRT2vT1, Mdiff_NRT2vT1, Betadiff_NRT2vT1, t_NRT2vT1, B_NRT2vT1, P.Value_NRT2vT1, adj.P.Val_NRT2vT1) %>%
  dplyr::mutate(Significance = ifelse(adj.P.Val_NRT2vT1<0.05, "Significant", "NS")) %>%
  dplyr::arrange(P.Value_NRT2vT1)

volcanoplot_NRT2vT1_ggplotobj <- dmps_NRT2vT1 %>%
  ggplot(aes(x = Betadiff_NRT2vT1, y = -log10(P.Value_NRT2vT1), col = Significance)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Non-responders: T2 vs T1",
       x = "%Methylation difference",
       y = bquote('-'~log[10]~'(p-value)')) +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none")

png(fig_volcanoplot_png_path, width = 7.5, height = 7.5, units = "in", res = 300)
print(vplot_all)
dev.off()

# RT2vT1 vs NRT2vT1 effectsize plot

dmps_RT2vT1_NRT2vT1 <- dmps %>%
  dplyr::mutate(Significance = dplyr::case_when(
    adj.P.Val_RvNRvT2vT1<0.05 ~ "Significant",
    .default = "NS")) %>%
  dplyr::arrange(P.Value_RvNRvT2vT1)

plot_boundary <- ceiling(max(abs(dmps_RT2vT1_NRT2vT1$t_T1RvNR), abs(dmps_RT2vT1_NRT2vT1$t_T2RvNR)))

scatterplot_T1RvNR_T2RvNR_ggplotobj <- dmps_RT2vT1_NRT2vT1 %>%
  ggplot(aes(x = t_RT2vT1, y = t_NRT2vT1)) +
  # annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 10, alpha = .2, fill = "#AEC6CF") +
  # annotate("rect", xmin = 0, xmax = -10, ymin = 0, ymax = 10, alpha = .2, fill = "#77DD77") +
  # annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = -10, alpha = .2, fill = "#FFFAA0") +
  # annotate("rect", xmin = 0, xmax = -10, ymin = 0, ymax = -10, alpha = .2, fill = "#FAA0A0") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point_rast(aes(col = Significance)) +
  theme_bw() +
  coord_cartesian(xlim = c(-plot_boundary, plot_boundary),
                  ylim = c(-plot_boundary, plot_boundary)) +
  labs(title = paste0(treatment, " time-associated CpGs"),
       x = "R",
       y = "NR")

png(fig_pairplot_png_path, width = 15, height = 15, units = "in", res = 300)
print(scatterplot_T1RvNR_T2RvNR_ggplotobj)
dev.off()
