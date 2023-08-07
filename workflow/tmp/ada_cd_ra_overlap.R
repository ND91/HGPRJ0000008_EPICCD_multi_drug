require(aggregation)

cd_ada_predictor_cpgs <- readxl::read_excel("config/horaizon/predictor_cpgs.xlsx") %>%
  dplyr::filter(Treatment == "Adalimumab")
ra_ada_predictor_cpgs <- readxl::read_excel("~/lkgres/projects/PRJ0000043_RAADARESP/output/01_DM_analyses/211115/horaizon/output/220704_Selected_markers_annotated.xlsx") %>%
  dplyr::filter(Cohort == "Secondary_split")

ada_predictor_cpgs <- rbind(data.frame(CGID = cd_ada_predictor_cpgs$CpG, Treatment = "Adalimumab", Disease = "CD"),
                            data.frame(CGID = ra_ada_predictor_cpgs$CGID, Treatment = "Adalimumab", Disease = "RA"))

anno_epic <- read.csv("/media/ssd1/andrewliyim/multi_drug/resources/epic_annotations.csv")
dmps_cd_ada_path <- "/media/ssd1/andrewliyim/multi_drug/output/epic/dmp/dmp_Adalimumab_annotated.csv"
dmps_ra_ada_path <- "~/lkgres/projects/PRJ0000043_RAADARESP/output/01_DM_analyses/220816/dmps/rvnr.csv"

dmps_cd_ada <- read.csv(dmps_cd_ada_path)[,-1]
dmps_ra_ada <- read.csv(dmps_ra_ada_path)[,-1]

dmps_cd_ra_ada <- dmps_cd_ada %>%
  dplyr::select(Betadiff_T1RvNR, P.Value_T1RvNR, adj.P.Val_T1RvNR, CGID) %>%
  dplyr::rename(Betadiff_T1RvNR_CD = Betadiff_T1RvNR,
                P.Value_T1RvNR_CD = P.Value_T1RvNR,
                adj.P.Val_T1RvNR_CD = adj.P.Val_T1RvNR) %>%
  dplyr::inner_join(dmps_ra_ada %>%
                     dplyr::select(Beta, P.Value, adj.P.Val, Name) %>%
                     dplyr::rename(Betadiff_T1RvNR_RA = Beta,
                                   P.Value_T1RvNR_RA = P.Value,
                                   adj.P.Val_T1RvNR_RA = adj.P.Val),
                   by = c("CGID" = "Name")) %>%
  dplyr::left_join(ada_predictor_cpgs, by = "CGID")

dmps_cd_ra_ada_predictor_pvals <- dmps_cd_ra_ada %>%
  dplyr::filter(Treatment %in% c("Adalimumab")) %>%
  dplyr::select(CGID, P.Value_T1RvNR_CD, P.Value_T1RvNR_RA) %>%
  tidyr::pivot_longer(-CGID, names_to = "Comparison", values_to = "pvalue") %>%
  split(., .$CGID)

dmps_cd_ra_ada_predictor_pvals_aggregated <- unlist(lapply(dmps_cd_ra_ada_predictor_pvals, function(cgid){
  aggregation::lancaster(cgid$pvalue, c(27, 93))
}))

dmps_cd_ra_ada <- dmps_cd_ra_ada %>%
  dplyr::left_join(data.frame(CGID = names(dmps_cd_ra_ada_predictor_pvals_aggregated),
                              pvalues = dmps_cd_ra_ada_predictor_pvals_aggregated) %>%
                     dplyr::mutate(padj = p.adjust(pvalues, method = "bonferroni")), 
                   by = "CGID")

dmps_cd_ra_ada_sig <- dmps_cd_ra_ada %>%
  dplyr::filter(Treatment %in% c("Adalimumab"),
                padj<0.05)

png("ada_ra_cd.png", width = 5, height = 5.5, units = "in", res = 300)
ggplot(dmps_cd_ra_ada, aes(x = Betadiff_T1RvNR_CD, y = Betadiff_T1RvNR_RA)) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_vline(xintercept = 0, col = "grey") +
  geom_point_rast(col = "#d3d3d3") +
  geom_point_rast(data = dmps_cd_ra_ada %>%
                    dplyr::filter(Treatment %in% c("Adalimumab")),
                  aes(col = Disease)) +
  theme_bw() +
  labs(title = "Adalimumab",
       x = "CD",
       y = "RA") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.pos = "bottom")
dev.off()

