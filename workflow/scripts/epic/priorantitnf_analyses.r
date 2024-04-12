#!/usr/bin/env R
# The goal is to train a model using prior aTNF exposure to predict response (Reviewer 2). 
# We will train on the discovery samples and predict on the validation samples.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

library(minfi)
library(sva)
library(dplyr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

rgset_rds <- args[1]
treatment <- args[2]
glmmodel_rds <- args[3]
summary_glmmodel_txt <- args[4]
validation_predictions_csv <- args[5]

rgset <- readRDS(rgset_rds) # This is not memory-efficient as we are only using the metadata and not the methylation data.

prefix <- case_when(
  treatment == "Vedolizumab" ~ "VDZ_",
  treatment == "Ustekinumab" ~ "UST_",
)

cohort_column <- paste0(prefix, "cohort")
timepoint_column <- paste0(prefix, "timepoint")
response_column <- paste0(prefix, "response")
response_type_column <- paste0(prefix, "response_type")

t1_discval_samples <- pData(rgset) %>%
  data.frame() %>%
  dplyr::filter(!!sym(timepoint_column) %in% c("T1"),
                !!sym(cohort_column) %in% c("EPIC-CD Discovery", "EPIC-CD Validation"),
                !is.na(!!sym(response_column)),
                Disease == "CD") %>%
  dplyr::filter(!is.na(Prior_aTNF_exposure))

t1_disc_samples <- t1_discval_samples %>%
  dplyr::filter(Center_source == "AmsterdamUMC") %>%
  dplyr::mutate(Prior_aTNF_exposure_recoded = as.numeric(ifelse(Prior_aTNF_exposure == "Exposed", 1, 0)),
                Response_recoded = as.numeric(ifelse(!!sym(response_column) == "R", 1, 0)))

t1_val_samples <- t1_discval_samples %>%
  dplyr::filter(Center_source == "John Radcliffe Hospital") %>%
  dplyr::mutate(Prior_aTNF_exposure_recoded = as.numeric(ifelse(Prior_aTNF_exposure == "Exposed", 1, 0)),
                Response_recoded = as.numeric(ifelse(!!sym(response_column) == "R", 1, 0)))

#Train model

antitnf_model <- glm(Response_recoded ~ Prior_aTNF_exposure, family = 'binomial', data = t1_disc_samples)
saveRDS(antitnf_model, glmmodel_rds, compress = "gzip")

sink(summary_glmmodel_txt)
summary(antitnf_model)
sink()

#Validate model

t1_val_samples <- t1_val_samples %>%
  dplyr::mutate(Predictions = predict.glm(antitnf_model, newdata = t1_val_samples, type = "response"))

write.csv(t1_val_samples, validation_predictions_csv)

sessionInfo()

# t1_val_samples %>%
#   ggplot(aes(d = Response_recoded, m = Predictions)) +
#   geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
#   geom_roc(n.cuts = 0) +
#   labs(y = "TPR",
#        x = "FPR") +
#   # facet_wrap(~Drug, nrow = 2) +
#   # scale_color_manual(values = center_colors) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())


# vdz_antitnf_predictions <- data.frame(Response = prior_antitnf_vdz_val$Response, 
#                                       Predicted = predict.glm(vdz_antitnf_model, newdata = prior_antitnf_vdz_val, type = "response")) %>%
#   dplyr::mutate(predictions = (Predicted > 0.52)*1) %>%
#   dplyr::mutate(Success = Response == predictions)
# 
# 
# prior_antitnf <- readxl::read_excel("config/horaizon/probabilities_VDZ_UST.xlsx")
# prior_antitnf_vdz <- prior_antitnf |>
#   dplyr::filter(Drug == "Vedolizumab") |>
#   dplyr::mutate(aTNF_exposed = as.numeric(ifelse(aTNF_exposed == "Exposed", 1, 0)),
#                 Response = as.numeric(ifelse(Response == "R", 1, 0)))
# prior_antitnf_vdz_disc <- prior_antitnf_vdz |>
#   dplyr::filter(Cohort == "AmsterdamUMC")
# prior_antitnf_vdz_val <- prior_antitnf_vdz |>
#   dplyr::filter(Cohort == "John Radcliffe Hospital")
# 
# vdz_antitnf_model <- glm(Response ~ aTNF_exposed, family = 'binomial', data = prior_antitnf_vdz_val)
# vdz_antitnf_predictions <- data.frame(Response = prior_antitnf_vdz_val$Response, 
#                           Predicted = predict.glm(vdz_antitnf_model, newdata = prior_antitnf_vdz_val, type = "response")) %>%
#   dplyr::mutate(predictions = (Predicted > 0.52)*1) %>%
#   dplyr::mutate(Success = Response == predictions)
# 
# vdz_antitnf_predictions %>% 
#   ggplot(aes(d = Response, m = Predicted)) +
#   geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
#   geom_roc(n.cuts = 0) +
#   labs(y = "TPR",
#        x = "FPR") +
#   # facet_wrap(~Drug, nrow = 2) +
#   # scale_color_manual(values = center_colors) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())
# 
# vdz_antitnf_predictions %>% 
#   dplyr::summarize(auroc = Metrics::auc(Response, Predicted),
#                    recall = Metrics::recall(Response, Predicted),
#                    tpr = Metrics::precision(Response, Predicted),
#                    f1 = Metrics::f1(Response, Predicted))
# 
# # Ustekinumab
# 
# prior_antitnf_ust <- prior_antitnf |>
#   dplyr::filter(Drug == "Ustekinumab") |>
#   dplyr::mutate(aTNF_exposed = as.numeric(ifelse(aTNF_exposed == "Exposed", 1, 0)),
#                 Response = as.numeric(ifelse(Response == "R", 1, 0)))
# prior_antitnf_ust_disc <- prior_antitnf_ust |>
#   dplyr::filter(Cohort == "AmsterdamUMC")
# prior_antitnf_ust_val <- prior_antitnf_ust |>
#   dplyr::filter(Cohort == "John Radcliffe Hospital")
# 
# ust_antitnf_model <- glm(Response ~ aTNF_exposed, family = 'binomial', data = prior_antitnf_ust_val)
# ust_antitnf_predictions <- data.frame(Response = prior_antitnf_ust_val$Response, 
#                                       Predicted = predict.glm(ust_antitnf_model, newdata = prior_antitnf_ust_val, type = "response")) %>%
#   dplyr::mutate(predictions = (Predicted > 0.52)*1) %>%
#   dplyr::mutate(Success = Response == predictions)
# 
# ust_antitnf_predictions %>% 
#   ggplot(aes(d = Response, m = Predicted)) +
#   geom_abline(slope = 1, intercept = 0, col = "#808080", linetype = 2) +
#   geom_roc(n.cuts = 0) +
#   labs(y = "TPR",
#        x = "FPR") +
#   # facet_wrap(~Drug, nrow = 2) +
#   # scale_color_manual(values = center_colors) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())
# 
# ust_antitnf_predictions %>% 
#   dplyr::summarize(auroc = Metrics::auc(Response, Predicted),
#                    recall = Metrics::recall(Response, Predicted),
#                    tpr = Metrics::precision(Response, Predicted),
#                    f1 = Metrics::f1(Response, Predicted))
           