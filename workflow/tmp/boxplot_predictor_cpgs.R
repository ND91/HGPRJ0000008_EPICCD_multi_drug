predictor_cpgs <- readxl::read_excel("config/horaizon/predictor_cpgs.xlsx")
anno_epic <- read.csv("/media/ssd1/andrewliyim/multi_drug/resources/epic_annotations.csv")

rgset <- readRDS("/media/ssd1/andrewliyim/multi_drug/output/epic/rgset_qc.Rds")
gmset <- minfi::preprocessFunnorm(rgset)
betas <- getBeta(gmset)

beta_df <- data.frame(pData(gmset)[,c("Timepoint", "Response")], t(betas[predictor_cpgs$CpG,])) %>%
  tidyr::pivot_longer(-c(Timepoint, Response), names_to = "CGID", values_to = "Methylation") %>%
  dplyr::left_join(predictor_cpgs, by = c("CGID" = "CpG"))%>%
  dplyr::filter(Timepoint == "T1")


ada_ifx_predictors <- beta_df %>%
  dplyr::filter(Treatment %in% c("Adalimumab", "Infliximab")) %>%
  ggplot(aes(x = forcats::fct_reorder(ID, -Methylation), y = Methylation, col = Response)) +
  geom_jitter(alpha = 0.25, position=position_dodge(0.8)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, position=position_dodge(0.8)) +
  facet_grid(~Treatment, space = "free_x", scales = "free_x") +
  labs(title = "Predictor CpGs",
       y = "%Methylation") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.pos = "bottom")

vdz_ust_predictors <- beta_df %>%
  dplyr::filter(Treatment %in% c("Vedolizumab", "Ustekinumab")) %>%
  ggplot(aes(x = forcats::fct_reorder(ID, -Methylation), y = Methylation, col = Response)) +
  geom_jitter(alpha = 0.25, position=position_dodge(0.8)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, position=position_dodge(0.8)) +
  facet_grid(~Treatment, space = "free_x", scales = "free_x") +
  labs(y = "%Methylation") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.pos = "bottom")

predictor_plot <- ggarrange(ada_ifx_predictors, vdz_ust_predictors, nrow = 2, common.legend = T, legend = "bottom")

png("boxjitterplot_ada_ifx_predictors.png", width = 20, height = 7.5, units = "in", res = 300)
predictor_plot
dev.off()

predictor_cpgs_anno <- predictor_cpgs %>%
  dplyr::left_join(anno_epic, by = c("CpG" = "Name")) %>%
  dplyr::mutate(gene = paste0(UCSC_RefGene_Name, ";", GencodeCompV12_NAME, ";", Enhancer_gene))

predictor_cpgs_anno_list <- split(predictor_cpgs_anno, predictor_cpgs_anno$Treatment)

lapply(predictor_cpgs_anno_list, function(treatment){
  genes <- unique(unlist(strsplit(treatment$gene, ";")))
  genes[!genes %in% c("", "NA")]
})
