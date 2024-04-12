#!/usr/bin/env R
# ROC plot of the samples collected at the AmsterdamUMC and the John Radcliffe Hospital at timepoint 1.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste0("Script needs 4 arguments. Current input is:", args))
}

require(ggplot2)
require(readxl)
require(dplyr)
require(plotROC)
require(Metrics)

horaizon_predictions_xlsx <- args[1]
rocplot_vdz_t1_pdf <- args[2]
rocplot_ust_t1_pdf <- args[3]
rocplot_vdz_ust_t1_pdf <- args[4]
rocplot_vdz_ust_t1_t2_pdf <- args[5]

# T1: AUMC & JRH

horaizon_t1_predictions <- readxl::read_excel(horaizon_predictions_xlsx) %>%
  dplyr::filter(Timepoint == "T1",
                Cohort %in% c("AmsterdamUMC", "John Radcliffe Hospital")) %>%
  dplyr::mutate(Response_coded = ifelse(Response == "R", 1, 0))

horaizon_t1_metrics <- horaizon_t1_predictions %>%
  dplyr::group_by(Cohort, Drug) %>%
  dplyr::summarize(auroc = Metrics::auc(Response_coded, Prediction),
                   recall = Metrics::recall(Response_coded, Prediction),
                   tpr = Metrics::precision(Response_coded, Prediction),
                   f1 = Metrics::f1(Response_coded, Prediction))

rocplot_vdz_t1_ggplotobj <- horaizon_t1_predictions %>% 
  dplyr::filter(Drug == "Vedolizumab") %>%
  ggplot(aes(d = Response_coded, m = Prediction, col = Cohort)) +
  geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
  geom_roc(n.cuts = 0) +
  labs(title = "Vedolizumab",
       y = "TPR",
       x = "FPR") +
  theme_bw() +
  theme(legend.pos = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 7, height = 7.5, file = rocplot_vdz_t1_pdf)
print(rocplot_vdz_t1_ggplotobj)
dev.off()

rocplot_ust_t1_ggplotobj <- horaizon_t1_predictions %>% 
  dplyr::filter(Drug == "Ustekinumab") %>%
  ggplot(aes(d = Response_coded, m = Prediction, col = Cohort)) +
  geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
  geom_roc(n.cuts = 0) +
  labs(title = "Ustekinumab",
       y = "TPR",
       x = "FPR") +
  theme_bw() +
  theme(legend.pos = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 7, height = 7.5, file = rocplot_ust_t1_pdf)
print(rocplot_ust_t1_ggplotobj)
dev.off()

rocplot_ust_vdz_t1_ggplotobj <- horaizon_t1_predictions %>% 
  ggplot(aes(d = Response_coded, m = Prediction, col = Cohort)) +
  geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
  geom_roc(n.cuts = 0) +
  labs(y = "TPR",
       x = "FPR") +
  facet_wrap(~Drug, ncol = 2) +
  theme_bw() +
  theme(legend.pos = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 14, height = 7.5, file = rocplot_vdz_ust_t1_pdf)
print(rocplot_ust_vdz_t1_ggplotobj)
dev.off()

# T1 & T2: AUMC

horaizon_aumc_t1_t2_predictions <- readxl::read_excel(horaizon_predictions_xlsx) %>%
  dplyr::filter(Cohort %in% c("AmsterdamUMC")) %>%
  dplyr::mutate(Response_coded = ifelse(Response == "R", 1, 0))

horaizon_t1_t2_metrics <- horaizon_aumc_t1_t2_predictions %>%
  dplyr::group_by(Timepoint, Drug) %>%
  dplyr::summarize(auroc = Metrics::auc(Response_coded, Prediction),
                   recall = Metrics::recall(Response_coded, Prediction),
                   tpr = Metrics::precision(Response_coded, Prediction),
                   f1 = Metrics::f1(Response_coded, Prediction)) %>%
  dplyr::arrange(Drug)

rocplot_vdz_ust_t1_t2_ggplotobj <- horaizon_aumc_t1_t2_predictions %>% 
  ggplot(aes(d = Response_coded, m = Prediction, col = Timepoint)) +
  geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
  geom_roc(n.cuts = 0) +
  labs(y = "TPR",
       x = "FPR") +
  facet_wrap(~Drug, ncol = 2) +
  scale_color_manual(values = c("T1" = "#d3d3d3", "T2" = "#000000")) +
  theme_bw() +
  theme(legend.pos = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 14, height = 7.5, file = rocplot_vdz_ust_t1_t2_pdf)
print(rocplot_vdz_ust_t1_t2_ggplotobj)
dev.off()

sessionInfo()