ora_promhypo_kegg[ora_promhypo_kegg$FDR<0.05,]
ora_promhyper_kegg[ora_promhyper_kegg$FDR<0.05,]

png("barplot_ada_ora_kegg.png", width = 7, height = 7, units = "in", res = 300)
rbind(data.frame(ora_promhypo_kegg[1:10,], Direction = "Hypomethylated"),
      data.frame(ora_promhyper_kegg[1:10,], Direction = "Hypermethylated")) %>%
  dplyr::mutate(Significance = ifelse(FDR<0.05, "Significant", "NS")) %>%
  ggplot(aes(x = -log10(P.DE), y = forcats::fct_reorder(Description, -log10(P.DE)), fill = Significance)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Significant" = "#000000", "NS" = "#d3d3d3")) +
  facet_wrap(~Direction, scales = "free_y", nrow = 2) +
  labs(title = "Adalimumab",
       subtitle = "Responders: treated vs pretreated",
       x = bquote('-'~log[10]~'(p-value)')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        legend.pos = "bottom")
dev.off()
