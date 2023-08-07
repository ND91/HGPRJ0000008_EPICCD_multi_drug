sample_metadata <- readxl::read_excel("config/samples/sample_metadata.xlsx")
files_metadata <- readxl::read_excel("config/samples/epic_files.xlsx")

#VDZ
vdz_sample_metadata <- sample_metadata %>%
  dplyr::filter(!is.na(VDZ_cohort))

vdz_train_t1_annotated <- read.csv("config/horaizon/VEDO_TRAIN_avg_pred.csv") %>%
  dplyr::mutate(DonorID = unlist(lapply(strsplit(vdz_train$SampleID, "_"), function(i)i[1])),
                Timepoint = "T1",
                Cohort = "Vedolizumab",
                Criteria = "Full") %>%
  dplyr::left_join(vdz_sample_metadata %>%
                     dplyr::filter(VDZ_timepoint == "T1"), 
                   by = c("DonorID")) %>% 
  dplyr::mutate(Database = ifelse(VDZ_cohort == "EPIC-CD Discovery", "Discovery cohort", "Validation cohort")) %>%
  dplyr::left_join(readxl::read_excel("config/samples/TNF_experienced_UST_VDZ.xlsx", sheet = 1) %>%
                     dplyr::mutate(CASTORID = as.character(CASTORID)),
                   by = c("ParticipantID" = "CASTORID",
                          "Database" = "Cohort AMC")) %>% 
  dplyr::rename(aTNF_exposed = `Previous TNF`) %>%
  dplyr::select(SampleID.y, DonorID, ParticipantID, Center_source, VDZ_timepoint, Cohort, VDZ_response, aTNF_exposed, Criteria, y_pred) %>%
  dplyr::rename(SampleID = SampleID.y, 
                Timepoint = VDZ_timepoint,
                Response = VDZ_response,
                Prediction = y_pred) %>%
  dplyr::mutate(aTNF_exposed = ifelse(aTNF_exposed == "yes", "Exposed", "Non-exposed"))

vdz_validation_t1_annotated <- read.csv("config/horaizon/VEDO_VAL_avg_pred.csv") %>%
  dplyr::mutate(Timepoint = "T1",
                Cohort = "Vedolizumab",
                Criteria = ifelse(SampleID %in% c("EPIC_OX_07", "EPIC_OX_09", "EPIC_OX_17", "EPIC_OX_33", "EPIC_OX_40", "EPIC_OX_48", "EPIC_OX_65", "EPIC_OX_102"), "Full", "Partial"),
                aTNF_exposed = ifelse(SampleID %in% c("EPIC_OX_65", "EPIC_OX_90", "EPIC_OX_09", "EPIC_OX_32", "EPIC_OX_17", "EPIC_OX_24", "EPIC_OX_113", "EPIC_OX_69", "EPIC_OX_105", "EPIC_OX_16", "EPIC_OX_64", "EPIC_OX_76"), "Exposed", "Non-exposed")) %>%
  dplyr::left_join(vdz_sample_metadata %>%
                     dplyr::filter(VDZ_timepoint == "T1"), 
                   by = c("SampleID")) %>%
  dplyr::select(SampleID, DonorID, Center_source, VDZ_timepoint, Cohort, VDZ_response, aTNF_exposed, Criteria, y_pred) %>%
  dplyr::rename(Timepoint = VDZ_timepoint,
                Response = VDZ_response,
                Prediction = y_pred)

vdz_validation_t2_annotated <- read.csv("config/horaizon/VEDO_T2_VAL_avg_pred.csv") %>%
  dplyr::left_join(files_metadata,
                   by = c("SampleID" = "SXSPOS")) %>%
  dplyr::left_join(sample_metadata,
                   by = c("SampleID.y" = "SampleID")) %>%
  dplyr::mutate(Cohort = "Vedolizumab",
                Criteria = "Full") %>%
  dplyr::mutate(Database = ifelse(VDZ_cohort == "EPIC-CD Discovery", "Discovery cohort", "Validation cohort")) %>%
  dplyr::left_join(readxl::read_excel("config/samples/TNF_experienced_UST_VDZ.xlsx", sheet = 1) %>%
                     dplyr::mutate(CASTORID = as.character(CASTORID)),
                   by = c("ParticipantID" = "CASTORID",
                          "Database" = "Cohort AMC")) %>% 
  dplyr::rename(aTNF_exposed = `Previous TNF`) %>%
  dplyr::select(SampleID.y, DonorID, ParticipantID, Center_source, VDZ_timepoint, Cohort, VDZ_response, aTNF_exposed, Criteria, y_pred) %>%
  dplyr::rename(SampleID = SampleID.y, 
                Timepoint = VDZ_timepoint,
                Response = VDZ_response,
                Prediction = y_pred) %>%
  dplyr::mutate(aTNF_exposed = ifelse(aTNF_exposed == "yes", "Exposed", "Non-exposed"))

#UST
ust_sample_metadata <- sample_metadata %>%
  dplyr::filter(!is.na(UST_cohort))

ust_train_t1_annotated <- read.csv("config/horaizon/USTE_TRAIN_avg_pred.csv") %>%
  dplyr::mutate(DonorID = unlist(lapply(strsplit(SampleID, "_"), function(i)i[1])),
                Timepoint = "T1",
                Cohort = "Ustekinumab",
                Criteria = "Full") %>%
  dplyr::left_join(ust_sample_metadata %>%
                     dplyr::filter(UST_timepoint == "T1"), 
                   by = c("DonorID")) %>%
  dplyr::mutate(Database = ifelse(UST_cohort == "EPIC-CD Discovery", "Discovery cohort", "Validation cohort")) %>%
  dplyr::left_join(readxl::read_excel("config/samples/TNF_experienced_UST_VDZ.xlsx", sheet = 2) %>%
                     dplyr::mutate(CASTORID = as.character(CASTORID)),
                   by = c("ParticipantID" = "CASTORID",
                          "Database" = "Cohort AMC")) %>% 
  dplyr::rename(aTNF_exposed = `Previous TNF`) %>%
  dplyr::select(SampleID.y, DonorID, Center_source, UST_timepoint, Cohort, UST_response, aTNF_exposed, Criteria, y_pred) %>%
  dplyr::rename(SampleID = SampleID.y, 
                Timepoint = UST_timepoint,
                Response = UST_response,
                Prediction = y_pred) %>%
  dplyr::mutate(aTNF_exposed = ifelse(aTNF_exposed == "yes", "Exposed", "Non-exposed"))

ust_validation_t1_annotated <- read.csv("config/horaizon/USTE_VAL_avg_pred.csv") %>%
  dplyr::mutate(Timepoint = "T1",
                Cohort = "Ustekinumab",
                Criteria = ifelse(SampleID %in% c("EPIC_OX_14", "EPIC_OX_23", "EPIC_OX_26", "EPIC_OX_46", "EPIC_OX_78", "EPIC_OX_81", "EPIC_OX_82", "EPIC_OX_86"), "Full", "Partial"),
                aTNF_exposed = ifelse(SampleID %in% c("EPIC_OX_81", "EPIC_OX_70", "EPIC_OX_82", "EPIC_OX_21", "EPIC_OX_56", "EPIC_OX_62", "EPIC_OX_61", "EPIC_OX_78", "EPIC_OX_14", "EPIC_OX_23", "EPIC_OX_66", "EPIC_OX_18", "EPIC_OX_58", "EPIC_OX_46", "EPIC_OX_27", "EPIC_OX_22", "EPIC_OX_28", "EPIC_OX_50", "EPIC_OX_55", "EPIC_OX_63", "EPIC_OX_71"), "Exposed", "Non-exposed")) %>%
  dplyr::left_join(ust_sample_metadata %>%
                     dplyr::filter(UST_timepoint == "T1"), 
                   by = c("SampleID")) %>%
  dplyr::select(SampleID, DonorID, Center_source, UST_timepoint, Cohort, UST_response, aTNF_exposed, Criteria, y_pred) %>%
  dplyr::rename(Timepoint = UST_timepoint,
                Response = UST_response,
                Prediction = y_pred)

ust_validation_t2_annotated <- read.csv("config/horaizon/USTE_T2_VAL_avg_pred.csv") %>%
  dplyr::left_join(files_metadata,
                   by = c("SampleID" = "SXSPOS")) %>%
  dplyr::left_join(sample_metadata,
                   by = c("SampleID.y" = "SampleID")) %>%
  dplyr::mutate(Cohort = "Ustekinumab",
                Criteria = "Full") %>%
  dplyr::mutate(Database = ifelse(VDZ_cohort == "EPIC-CD Discovery", "Discovery cohort", "Validation cohort")) %>%
  dplyr::left_join(readxl::read_excel("config/samples/TNF_experienced_UST_VDZ.xlsx", sheet = 1) %>%
                     dplyr::mutate(CASTORID = as.character(CASTORID)),
                   by = c("ParticipantID" = "CASTORID",
                          "Database" = "Cohort AMC")) %>% 
  dplyr::rename(aTNF_exposed = `Previous TNF`) %>%
  dplyr::select(SampleID.y, DonorID, Center_source, UST_timepoint, Cohort, UST_response, aTNF_exposed, Criteria, y_pred) %>%
  dplyr::rename(SampleID = SampleID.y, 
                Timepoint = UST_timepoint,
                Response = UST_response,
                Prediction = y_pred) %>%
  dplyr::mutate(aTNF_exposed = ifelse(aTNF_exposed == "yes", "Exposed", "Non-exposed"))

probabilities_VDZ_UST <- vdz_train_t1_annotated %>%
  dplyr::rows_append(vdz_validation_t1_annotated) %>%
  dplyr::rows_append(vdz_validation_t2_annotated) %>%
  dplyr::rows_append(ust_train_t1_annotated) %>%
  dplyr::rows_append(ust_validation_t1_annotated) %>%
  dplyr::rows_append(ust_validation_t2_annotated)
  
write.csv(probabilities_VDZ_UST, "config/horaizon/probabilities_VDZ_UST.csv")
