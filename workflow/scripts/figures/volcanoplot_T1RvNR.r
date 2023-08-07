#!/usr/bin/env R
# Volcanoplot of each treatment comparing responders with non-responders at T1.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste0("Script needs 6 arguments. Current input is:", args))
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggarrange))
suppressPackageStartupMessages(library(readxl))

dmps_ada_path <- args[1]
dmps_ifx_path <- args[2]
dmps_vdz_path <- args[3]
dmps_ust_path <- args[4]
predictor_cpgs_path <- args[5]
fig_volcanoplot_png_path <- args[6]
fig_pairplot_png_path <- args[7]

predictor_cpgs <- readxl::read_excel(predictor_cpgs_path) %>%
  dplyr::mutate(Type = "Predictor")

dmps_T1RvNR_ada <- read.csv(dmps_ada_path) %>%
  dplyr::select(CGID, t_T1RvNR, P.Value_T1RvNR, adj.P.Val_T1RvNR, Betadiff_T1RvNR) %>%
  dplyr::mutate(Treatment = "Adalimumab")
dmps_T1RvNR_ifx <- read.csv(dmps_ifx_path) %>%
  dplyr::select(CGID, t_T1RvNR, P.Value_T1RvNR, adj.P.Val_T1RvNR, Betadiff_T1RvNR) %>%
  dplyr::mutate(Treatment = "Infliximab")
dmps_T1RvNR_vdz <- read.csv(dmps_vdz_path) %>%
  dplyr::select(CGID, t_T1RvNR, P.Value_T1RvNR, adj.P.Val_T1RvNR, Betadiff_T1RvNR) %>%
  dplyr::mutate(Treatment = "Vedolizumab")
dmps_T1RvNR_ust <- read.csv(dmps_ust_path) %>%
  dplyr::select(CGID, t_T1RvNR, P.Value_T1RvNR, adj.P.Val_T1RvNR, Betadiff_T1RvNR) %>%
  dplyr::mutate(Treatment = "Ustekinumab")

dmps_T1RvNR <- rbind(dmps_T1RvNR_ada, dmps_T1RvNR_ifx, dmps_T1RvNR_vdz, dmps_T1RvNR_ust)

#dmps_T1RvNR$Significant <- ifelse(dmps_T1RvNR$adj.P.Val_T1RvNR<0.05, "Significant", "NS")

dmps_T1RvNR <- dmps_T1RvNR %>%
  dplyr::left_join(predictor_cpgs, by = c("CGID" = "CpG", "Treatment"))

# Pairplot

vplot_ada <- dmps_T1RvNR %>%
  dplyr::filter(Treatment == "Adalimumab") %>%
  ggplot(aes(x = Betadiff_T1RvNR, y = -log10(P.Value_T1RvNR))) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR %>% 
                    dplyr::filter(Treatment == "Adalimumab" & Type == "Predictor"),
                  col = "#000000") +
  labs(title = "Adalimumab",
       x = "%Methylation difference (T1:RvNR)",
       y = "Adalimumab") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        legend.pos = "none")

vplot_ifx <- dmps_T1RvNR %>%
  dplyr::filter(Treatment == "Infliximab") %>%
  ggplot(aes(x = Betadiff_T1RvNR, y = -log10(P.Value_T1RvNR))) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR %>% 
                    dplyr::filter(Treatment == "Infliximab" & Type == "Predictor"),
                  col = "#000000") +
  labs(title = "Infliximab",
       x = "%Methylation difference (T1:RvNR)",
       y = bquote('-'~log[10]~'(p-value)')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "none")

vplot_vdz <- dmps_T1RvNR %>%
  dplyr::filter(Treatment == "Vedolizumab") %>%
  ggplot(aes(x = Betadiff_T1RvNR, y = -log10(P.Value_T1RvNR))) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR %>% 
                    dplyr::filter(Treatment == "Vedolizumab" & Type == "Predictor"),
                  col = "#000000") +
  labs(title = "Vedolizumab",
       x = "%Methylation difference (T1:RvNR)",
       y = bquote('-'~log[10]~'(p-value)')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "none")

vplot_ust <- dmps_T1RvNR %>%
  dplyr::filter(Treatment == "Ustekinumab") %>%
  ggplot(aes(x = Betadiff_T1RvNR, y = -log10(P.Value_T1RvNR))) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR %>% 
                    dplyr::filter(Treatment == "Ustekinumab" & Type == "Predictor"),
                  col = "#000000") +
  labs(title = "Ustekinumab",
       x = "Ustekinumab") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "none")

dmps_T1RvNR_wide <- dmps_T1RvNR_ada %>%
  dplyr::select(-Treatment) %>%
  dplyr::rename_with(~ paste0(., "_ADA"), -1) %>%
  dplyr::left_join(dmps_T1RvNR_ifx %>%
                     dplyr::select(-Treatment) %>%
                     dplyr::rename_with(~ paste0(., "_IFX"), -1), 
                   by = "CGID") %>%
  dplyr::left_join(dmps_T1RvNR_vdz %>%
                     dplyr::select(-Treatment) %>%
                     dplyr::rename_with(~ paste0(., "_VDZ"), -1), 
                   by = "CGID") %>%
  dplyr::left_join(dmps_T1RvNR_ust %>%
                     dplyr::select(-Treatment) %>%
                     dplyr::rename_with(~ paste0(., "_UST"), -1), 
                   by = "CGID") %>%
  dplyr::left_join(predictor_cpgs, by = c("CGID" = "CpG"))

eplot_ada_ifx <- ggplot(dmps_T1RvNR_wide, aes(x = Betadiff_T1RvNR_ADA, y = Betadiff_T1RvNR_IFX)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR_wide %>% 
                    dplyr::filter(Type == "Predictor"),
                  aes(col = Treatment)) +
  labs() +
  theme_bw() +
  labs(y = "Infliximab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

eplot_ada_vdz <- ggplot(dmps_T1RvNR_wide, aes(x = Betadiff_T1RvNR_ADA, y = Betadiff_T1RvNR_VDZ)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR_wide %>% 
                    dplyr::filter(Type == "Predictor"),
                  aes(col = Treatment)) +
  theme_bw() +
  labs(y = "Vedolizumab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

eplot_ada_ust <- ggplot(dmps_T1RvNR_wide, aes(x = Betadiff_T1RvNR_ADA, y = Betadiff_T1RvNR_UST)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR_wide %>% 
                    dplyr::filter(Type == "Predictor"),
                  aes(col = Treatment)) +
  theme_bw() +
  labs(x = "Adalimumab",
       y = "Ustekinumab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

eplot_ifx_vdz <- ggplot(dmps_T1RvNR_wide, aes(x = Betadiff_T1RvNR_IFX, y = Betadiff_T1RvNR_VDZ)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR_wide %>% 
                    dplyr::filter(Type == "Predictor"),
                  aes(col = Treatment)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

eplot_ifx_ust <- ggplot(dmps_T1RvNR_wide, aes(x = Betadiff_T1RvNR_IFX, y = Betadiff_T1RvNR_UST)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR_wide %>% 
                    dplyr::filter(Type == "Predictor"),
                  aes(col = Treatment)) +
  theme_bw() +
  labs(x = "Infliximab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())

eplot_vdz_ust <- ggplot(dmps_T1RvNR_wide, aes(x = Betadiff_T1RvNR_VDZ, y = Betadiff_T1RvNR_UST)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR_wide %>% 
                    dplyr::filter(Type == "Predictor"),
                  aes(col = Treatment)) +
  theme_bw() +
  labs(x = "Vedolizumab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank())

pplot <- ggarrange(vplot_ada, ggplot() + theme_void(), ggplot() + theme_void(), ggplot() + theme_void(),
                   eplot_ada_ifx, vplot_ifx, ggplot() + theme_void(), ggplot() + theme_void(),
                   eplot_ada_vdz, eplot_ifx_vdz, vplot_vdz, ggplot() + theme_void(),
                   eplot_ada_ust, eplot_ifx_ust, eplot_vdz_ust, vplot_ust, 
                   align = "hv", common.legend = T, legend = "bottom")

png(fig_pairplot_png_path, width = 15, height = 16, units = "in", res = 300)
print(pplot)
dev.off()

overlap_predictor_ada_ifx <- ggplot(dmps_T1RvNR_wide, aes(x = Betadiff_T1RvNR_ADA, y = Betadiff_T1RvNR_IFX)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_T1RvNR_wide %>% 
                    dplyr::filter(Treatment %in% c("Adalimumab", "Infliximab") & Type == "Predictor"),
                  aes(col = Treatment)) +
  labs() +
  theme_bw() +
  labs(title = "T1: RvNR",
       subtitle = "Pearson correlation: -0.07",
       x = "Adalimumab",
       y = "Infliximab") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")

png(fig_pairplot_png_path, width = 5, height = 6, units = "in", res = 300)
print(overlap_predictor_ada_ifx)
dev.off()
