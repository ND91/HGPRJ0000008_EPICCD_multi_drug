#!/usr/bin/env R
# ROC plot of the samples collected at the AmsterdamUMC at timepoint 1.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste0("Script needs 2 arguments. Current input is:", args))
}

require(ggplot2)
require(readxl)
require(dplyr)
require(plotROC)
require(Metrics)

horaizon_predictions_xlsx <- args[1]
rocplot_pdf <- args[2]

horaizon_predictions <- readxl::read_excel(horaizon_predictions_xlsx) %>%
  dplyr::filter(Timepoint == "T1",
                Center_source == "John Radcliffe Hospital") %>%
  dplyr::mutate(Response_coded = ifelse(Response == "R", 1, 0))

horaizon_metrics <- horaizon_predictions %>%
  dplyr::group_by(Cohort) %>%
  dplyr::summarize(auroc = Metrics::auc(Response_coded, Prediction),
                   recall = Metrics::recall(Response_coded, Prediction),
                   tpr = Metrics::precision(Response_coded, Prediction),
                   f1 = Metrics::f1(Response_coded, Prediction))

rocplot_ggplotobj <- ggplot(horaizon_predictions, aes(d = Response_coded, m = Prediction, col = Cohort)) +
  geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
  geom_roc(n.cuts = 0) +
  labs(title = "John Radcliffe Hospital cohort pretreatment",
       #subtitle = "AUROC-VDZ = 0.87\nAUROC-UST = 0.89",
       y = "TPR",
       x = "FPR") +
  theme_bw() +
  theme(legend.pos = "bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(width = 7, height = 7.5, file = rocplot_pdf)
print(rocplot_ggplotobj)
dev.off()
