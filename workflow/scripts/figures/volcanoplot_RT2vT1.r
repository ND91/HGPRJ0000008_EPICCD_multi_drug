#!/usr/bin/env R
# Volcanoplot of each treatment comparing responders over time

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggpubr))

dmps_ada_path <- args[1]
dmps_ifx_path <- args[2]
dmps_vdz_path <- args[3]
dmps_ust_path <- args[4]
fig_volcanoplot_png_path <- args[5]
fig_pairplot_png_path <- args[6]

dmps_RT2vT1_ada <- read.csv(dmps_ada_path) %>%
  dplyr::select(CGID, t_RT2vT1, P.Value_RT2vT1, adj.P.Val_RT2vT1, Betadiff_RT2vT1) %>%
  dplyr::mutate(Treatment = "Adalimumab")
dmps_RT2vT1_ifx <- read.csv(dmps_ifx_path) %>%
  dplyr::select(CGID, t_RT2vT1, P.Value_RT2vT1, adj.P.Val_RT2vT1, Betadiff_RT2vT1) %>%
  dplyr::mutate(Treatment = "Infliximab")
dmps_RT2vT1_vdz <- read.csv(dmps_vdz_path) %>%
  dplyr::select(CGID, t_RT2vT1, P.Value_RT2vT1, adj.P.Val_RT2vT1, Betadiff_RT2vT1) %>%
  dplyr::mutate(Treatment = "Vedolizumab")
dmps_RT2vT1_ust <- read.csv(dmps_ust_path) %>%
  dplyr::select(CGID, t_RT2vT1, P.Value_RT2vT1, adj.P.Val_RT2vT1, Betadiff_RT2vT1) %>%
  dplyr::mutate(Treatment = "Ustekinumab")

dmps_RT2vT1 <- rbind(dmps_RT2vT1_ada, dmps_RT2vT1_ifx, dmps_RT2vT1_vdz, dmps_RT2vT1_ust)

dmps_RT2vT1$Significant <- ifelse(dmps_RT2vT1$adj.P.Val_RT2vT1<0.05, "Significant", "NS")

# Volcanoplot

vplot_all <- dmps_RT2vT1 %>%
  ggplot(aes(x = Betadiff_RT2vT1, y = -log10(P.Value_RT2vT1), col = Significant)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Responders: T2 vs T1",
       x = "%Methylation difference",
       y = bquote('-'~log[10]~'(p-value)')) +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  facet_wrap(~Treatment, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "none")

png(fig_volcanoplot_png_path, width = 7.5, height = 7.5, units = "in", res = 300)
print(vplot_all)
dev.off()

# Pairplot

vplot_ada <- dmps_RT2vT1 %>%
  dplyr::filter(Treatment == "Adalimumab") %>%
  ggplot(aes(x = Betadiff_RT2vT1, y = -log10(P.Value_RT2vT1), col = Significant)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Adalimumab",
       y = "Adalimumab") +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "none")

vplot_ifx <- dmps_RT2vT1 %>%
  dplyr::filter(Treatment == "Infliximab") %>%
  ggplot(aes(x = Betadiff_RT2vT1, y = -log10(P.Value_RT2vT1), col = Significant)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Infliximab",
       x = "%Methylation difference (T2vT1)",
       y = bquote('-'~log[10]~'(p-value)')) +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "none")

vplot_vdz <- dmps_RT2vT1 %>%
  dplyr::filter(Treatment == "Vedolizumab") %>%
  ggplot(aes(x = Betadiff_RT2vT1, y = -log10(P.Value_RT2vT1), col = Significant)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Vedolizumab",
       x = "%Methylation difference (T2vT1)",
       y = bquote('-'~log[10]~'(p-value)')) +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "none")

vplot_ust <- dmps_RT2vT1 %>%
  dplyr::filter(Treatment == "Ustekinumab") %>%
  ggplot(aes(x = Betadiff_RT2vT1, y = -log10(P.Value_RT2vT1), col = Significant)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs(title = "Ustekinumab",
       x = "Ustekinumab") +
  scale_color_manual(values = list(Significant = "#000000", NS = "#d3d3d3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "none")

dmps_RT2vT1_wide <- dmps_RT2vT1_ada %>%
  dplyr::select(-Treatment) %>%
  dplyr::rename_with(~ paste0(., "_ADA"), -1) %>%
  dplyr::left_join(dmps_RT2vT1_ifx %>%
                     dplyr::select(-Treatment) %>%
                     dplyr::rename_with(~ paste0(., "_IFX"), -1), 
                   by = "CGID") %>%
  dplyr::left_join(dmps_RT2vT1_vdz %>%
                     dplyr::select(-Treatment) %>%
                     dplyr::rename_with(~ paste0(., "_VDZ"), -1), 
                   by = "CGID") %>%
  dplyr::left_join(dmps_RT2vT1_ust %>%
                     dplyr::select(-Treatment) %>%
                     dplyr::rename_with(~ paste0(., "_UST"), -1), 
                   by = "CGID")



eplot_ada_ifx <- ggplot(dmps_RT2vT1_wide, aes(x = Betadiff_RT2vT1_ADA, y = Betadiff_RT2vT1_IFX)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  labs() +
  theme_bw() +
  labs(y = "Infliximab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

eplot_ada_vdz <- ggplot(dmps_RT2vT1_wide, aes(x = Betadiff_RT2vT1_ADA, y = Betadiff_RT2vT1_VDZ)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  theme_bw() +
  labs(y = "Vedolizumab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

eplot_ada_ust <- ggplot(dmps_RT2vT1_wide, aes(x = Betadiff_RT2vT1_ADA, y = Betadiff_RT2vT1_UST)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  theme_bw() +
  labs(x = "Adalimumab",
       y = "Ustekinumab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

eplot_ifx_vdz <- ggplot(dmps_RT2vT1_wide, aes(x = Betadiff_RT2vT1_IFX, y = Betadiff_RT2vT1_VDZ)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

eplot_ifx_ust <- ggplot(dmps_RT2vT1_wide, aes(x = Betadiff_RT2vT1_IFX, y = Betadiff_RT2vT1_UST)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  theme_bw() +
  labs(x = "Infliximab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())

eplot_vdz_ust <- ggplot(dmps_RT2vT1_wide, aes(x = Betadiff_RT2vT1_VDZ, y = Betadiff_RT2vT1_UST)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast() +
  theme_bw() +
  labs(x = "Vedolizumab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())

pplot <- ggarrange(vplot_ada, ggplot() + theme_void(), ggplot() + theme_void(), ggplot() + theme_void(),
                   eplot_ada_ifx, vplot_ifx, ggplot() + theme_void(), ggplot() + theme_void(),
                   eplot_ada_vdz, eplot_ifx_vdz, vplot_vdz, ggplot() + theme_void(),
                   eplot_ada_ust, eplot_ifx_ust, eplot_vdz_ust, vplot_ust, 
                   align = "hv")

png(fig_pairplot_png_path, width = 15, height = 15, units = "in", res = 300)
print(pplot)
dev.off()
