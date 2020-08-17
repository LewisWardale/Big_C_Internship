breast_data <- read.csv("/Users/lewiswardale/Desktop/TCGA-BRCA/clinical_data_TCGA-BRCA_eridanus.csv")

library(tidyverse)
library(survminer)
library(survival)
library(DataExplorer)
library(magrittr)
library(plyr)
# Filter code -------------------------------------------------------------

Breast_Summary_Stages_Filtered_Blanks_Age_Race <- breast_data %>%
  drop_columns(
    c("tumor_tissue_site_other", 
      "days_to_last_known_alive",
      "pos_finding_metastatic_breast_carcinoma_estrogen_receptor_other_measuremenet_scale_text", 
      "metastatic_breast_carcinoma_estrogen_receptor_detection_method_text", 
      "metastatic_breast_carcinoma_pos_finding_progesterone_receptor_other_measure_scale_text", 
      "metastatic_breast_carcinoma_her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text",
      "metastatic_breast_carcinoma_pos_finding_other_scale_measurement_text",
      "metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_ratio_value",
      "her2_neu_and_centromere_17_copy_number_metastatic_breast_carcinoma_analysis_input_total_number_count",
      "metastatic_breast_carcinoma_fluorescence_in_situ_hybridization_diagnostic_proc_centromere_17_signal_result_range",
      "her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count",
      "metastatic_breast_carcinoma_progesterone_receptor_detection_method_text",
      "metastatic_breast_carcinoma_pos_finding_her2_erbb2_other_measure_scale_text",
      "metastatic_breast_carcinoma_her2_erbb_method_calculation_method_text",
      "her2_neu_breast_carcinoma_copy_analysis_input_total_number",
      "stage_event_clinical_stage",
      "stage_event_psa",
      "stage_event_gleason_grading",
      "stage_event_ann_arbor",
      "stage_event_serum_markers",
      "stage_event_igcccg_stage",
      "stage_event_masaoka_stage",
      "her2_neu_metastatic_breast_carcinoma_copy_analysis_input_total_number")) %>%
  mutate(stage_event_pathologic_stage = as.character(stage_event_pathologic_stage)) %>%
  mutate(stage_event_pathologic_stage = case_when(
    stage_event_pathologic_stage == "Stage IA" ~ "Stage I",
    stage_event_pathologic_stage == "Stage IB" ~ "Stage I",
    stage_event_pathologic_stage == "Stage IIA" ~ "Stage II",
    stage_event_pathologic_stage == "Stage IIB" ~ "Stage II",
    stage_event_pathologic_stage == "Stage IIIA" ~ "Stage III",
    stage_event_pathologic_stage == "Stage IIIB" ~ "Stage III",
    stage_event_pathologic_stage == "Stage IIIC" ~ "Stage III",
    TRUE ~ stage_event_pathologic_stage
  )) %>%
  mutate(Years_to_death = days_to_death / 365) %>%
  mutate(Age = abs(days_to_birth / 365)) %>%
  mutate(Years_to_last_followup = abs(days_to_last_followup /365)) %>%
  dplyr::filter(!stage_event_pathologic_stage %in% c("", "Stage X")) %>%
  mutate(vital_status_TRUE_FALSE = ifelse(vital_status == "Alive", FALSE, TRUE)) %>%
  mutate(Years_to_death_and_last_followup = ifelse(is.na(Years_to_death), Years_to_last_followup,Years_to_death)) %>%
  mutate(age_group = case_when(Age >= 49.07 & Age < 58.60 ~ "Lower_Quartile", 
                               Age >= 58.60 & Age <= 67.66 ~ "Upper_Quartile", 
                               Age > 67.66 ~ "Above_Upper_Q",
                               Age < 49.07 ~ "Below_Lower_Q")) 

# Mutate Stage ------------------------------------------------------------


g <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  dplyr::filter(stage_event_pathologic_stage %in% c("Stage I", "Stage II", "Stage III", "Stage IV"))

g_factor <- g %>%
  mutate(Stage_group = as.factor(case_when(stage_event_pathologic_stage == "Stage I" ~ 0, 
                                 stage_event_pathologic_stage == "Stage II" ~ 0, 
                                 stage_event_pathologic_stage == "Stage III" ~ 1, 
                                 stage_event_pathologic_stage == "Stage IV" ~ 1))) 

 m <- lm(Age ~ Stage_group , g_factor)
summary(m)$r.squared

ggplot(g_factor, aes(x = Stage_group , y = Age)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE) +  
  geom_text(label = summary(m)$r.squared, x = 0.3, y = 1.0, show.legend = FALSE)

