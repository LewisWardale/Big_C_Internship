# Importing Data ----------------------------------------------------------

Breast_data_w_LPD <- read.csv("/Users/lewiswardale/Desktop/clinical_BRCA_assigned.csv")

library(survminer)
library(survival)
library(lubridate)
library(dplyr)
library(MASS)
library(tidyr)
library(DataExplorer)
library(magrittr)
library(ggplot2)

Breast_Summary_Stages_Filtered_Blanks_Age_Race <- Breast_data_w_LPD %>%
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
  dplyr::filter(!race_list == "" ) %>%
  mutate(vital_status_TRUE_FALSE = ifelse(vital_status == "Alive", FALSE, TRUE)) %>%
  mutate(Years_to_death_and_last_followup = ifelse(is.na(Years_to_death), Years_to_last_followup,Years_to_death)) %>%
  mutate(age_group = case_when(Age >= 49.07 & Age < 58.60 ~ "Lower_Quartile", 
                               Age >= 58.60 & Age <= 67.66 ~ "Upper_Quartile", 
                               Age > 67.66 ~ "Above_Upper_Q",
                               Age < 49.07 ~ "Below_Lower_Q"))

#Removed columns containing no observations, 
#Summarised each stage e.g. IA = I,
#Changed, days to death, days to birth and days to last follow up to years,
#Removed all "" observations from Stage and Race columns and also Stage X,as they are uncharacterised tumour samples
#Created new column that characterised vital status to TRUE or FALSE, Alive = False, Dead = TRUE,
#Combined the alive patients last follow up (in Years) with the years to death. 
#Created new column "age group" which puts age into quartiles

# Box Plot Age vs Stage ---------------------------------------------------

Age_v_LPD <- ggplot(Breast_Summary_Stages_Filtered_Blanks_Age_Race, aes(x = Max_LPD, y = Age, fill = Max_LPD)) + 
  geom_boxplot(outlier.colour = "red")

Age_v_LPD + ggtitle("Box Plot of Age and LPD") + labs(y = "Age (Years)", x = "Max LPD") + labs(fill = "Max LPD")

Age_v_LPD_Race <- Age_v_LPD <- ggplot(Breast_Summary_Stages_Filtered_Blanks_Age_Race, aes(x = Max_LPD, y = Age, fill = race_list )) + 
  geom_boxplot(outlier.colour = "red")

Age_v_LPD_Race + ggtitle("Box Plot of Age and LPD and Race") + labs(y = "Age (Years)", x = "Max LPD") + labs(fill = "Race")

Breast_Summary_Stages_Filtered_Blanks_Age_Race <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  dplyr::filter(!is.na(stage_event_pathologic_stage))

Age_v_LPD_Stage <- Age_v_LPD <- ggplot(Breast_Summary_Stages_Filtered_Blanks_Age_Race, aes(x = Max_LPD, y = Age, fill = stage_event_pathologic_stage )) + 
  geom_boxplot(outlier.colour = "red")

Age_v_LPD_Stage + ggtitle("Box Plot of Age and LPD and Stage") + labs(y = "Age (Years)", x = "Max LPD") + labs(fill = "Stage")

# Pairwise LPD ----------------------------------------------------------

Pairwise_LPD <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ Max_LPD,
                          data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
Pairwise_LPD

# Filtering for LPD 2 and 5 --------------------------------------------------

Breast_unknown_LPD_2_5 <-  Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  dplyr::filter(Max_LPD %in% c("LPD_2", "LPD_5"))

Breast_unknown_LPD_2_5$Max_LPD

# KM LPD 2 vs 5 --------------------------------------------------------------

surv_object <- Surv(time = Breast_unknown_LPD_2_5$Years_to_death_and_last_followup, event =Breast_unknown_LPD_2_5$vital_status_TRUE_FALSE)
surv_object
fit1 <- survfit(surv_object ~ Max_LPD, data = Breast_unknown_LPD_2_5) 
fit1

LPD_2_5 <- ggsurvplot(
  fit1,
  data = Breast_unknown_LPD_2_5, 
  pval = TRUE, legend.title = "Max LPD", 
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs =
    c("LPD 2","LPD 5"),
  ggtheme = theme_bw(),
  break.time.by = 2,
  tables.height = 0.4
) + xlab("Years")  

LPD_2_5$plot <- LPD_2_5$plot + labs(
  title = "KM Max LPD plot for TCGA-BRCA, 2 and 5")

LPD_2_5 <- ggpar(
  LPD_2_5,
  font.title = c(16, "bold", "darkblue"),
  font.x = c(14, "plain", "black"),
  font.y = c(14, "plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab =  c(12, "plain", "black"),
  legend = "bottom"
)

LPD_2_5

# Filtering for LPD 5 and 6 --------------------------------------------------

Breast_unknown_LPD_5_6 <-  Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  dplyr::filter(Max_LPD %in% c("LPD_5", "LPD_6"))

Breast_unknown_LPD_5_6$Max_LPD

# KM 5 vs 6 --------------------------------------------------------------

surv_object <- Surv(time = Breast_unknown_LPD_5_6$Years_to_death_and_last_followup, event =Breast_unknown_LPD_5_6$vital_status_TRUE_FALSE)
surv_object
fit1 <- survfit(surv_object ~ Max_LPD, data = Breast_unknown_LPD_5_6) 
fit1

LPD_5_6 <- ggsurvplot(
  fit1,
  data = Breast_unknown_LPD_5_6, 
  pval = TRUE, legend.title = "Max LPD", 
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs =
    c("LPD 5","LPD 6"),
  ggtheme = theme_bw(),
  break.time.by = 2,
  tables.height = 0.4
) + xlab("Years")  

LPD_5_6$plot <- LPD_5_6$plot + labs(
  title = "KM Max LPD plot for TCGA-BRCA, 5 and 6")

LPD_5_6 <- ggpar(
  LPD_5_6,
  font.title = c(16, "bold", "darkblue"),
  font.x = c(14, "plain", "black"),
  font.y = c(14, "plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab =  c(12, "plain", "black"),
  legend = "bottom"
)

LPD_5_6

# Race and LPD Tally -----------------------------------------------

Race_LPD <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  group_by(Max_LPD, race_list) %>%
  tally() %>%
  ungroup() %>%
  spread(key = Max_LPD, value = n, fill = 0) %>%
  as.data.frame(Race_LPD) %>%
  set_colnames(c("Race_list","LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7" ))

Race_LPD

# Chi Squared Race LPD ----------------------------------------------------

#The above data frame is the reference for each table made below

RACE_LPD_1_AIoAN <-matrix(c(0,1,125,882),ncol=2,byrow=TRUE)
rownames(RACE_LPD_1_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_LPD_1_AIoAN)<-c("LPD_1","Not_LPD_1")
RACE_LPD_1_AIoAN <- as.table(RACE_LPD_1_AIoAN)

chisq.test(RACE_LPD_1_AIoAN)

#p-value = 1

RACE_LPD_1_ASIAN <-matrix(c(9,52,116,831),ncol=2,byrow=TRUE)
rownames(RACE_LPD_1_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_LPD_1_ASIAN)<-c("LPD_1","Not_LPD_1")
RACE_LPD_1_ASIAN <- as.table(RACE_LPD_1_ASIAN)

chisq.test(RACE_LPD_1_ASIAN)

#p-value = 0.7077

RACE_LPD_1_BoAA <-matrix(c(17,161,108,722),ncol=2,byrow=TRUE)
rownames(RACE_LPD_1_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_LPD_1_BoAA)<-c("LPD_1","Not_LPD_1")
RACE_LPD_1_BoAA <- as.table(RACE_LPD_1_BoAA)

chisq.test(RACE_LPD_1_BoAA)

#p-value = 0.2517

RACE_LPD_1_W <-matrix(c(99,669,26,214),ncol=2,byrow=TRUE)
rownames(RACE_LPD_1_W)<-c("White", "All other races")
colnames(RACE_LPD_1_W)<-c("LPD_1","LPD_1_Total")
RACE_LPD_1_W <- as.table(RACE_LPD_1_W)

chisq.test(RACE_LPD_1_W)

#p-value = 0.4642

pvalues <- c(1, 0.7077, 0.2517, 0.4642)

p.adjust(pvalues, method= "BH")

#1.0000 0.9436 0.9284 0.9284

#The above data frame is the reference for each table made below

RACE_LPD_2_AIoAN <-matrix(c(0,1,144,863),ncol=2,byrow=TRUE)
rownames(RACE_LPD_2_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_LPD_2_AIoAN)<-c("LPD_2","Not_LPD_2")
RACE_LPD_2_AIoAN <- as.table(RACE_LPD_2_AIoAN)

chisq.test(RACE_LPD_2_AIoAN)

#p-value = 1

RACE_LPD_2_ASIAN <-matrix(c(4,57,140,807),ncol=2,byrow=TRUE)
rownames(RACE_LPD_2_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_LPD_2_ASIAN)<-c("LPD_2","Not_LPD_2")
RACE_LPD_2_ASIAN <- as.table(RACE_LPD_2_ASIAN)

chisq.test(RACE_LPD_2_ASIAN)

#p-value = 0.1116

RACE_LPD_2_BoAA <-matrix(c(11,167,133,697),ncol=2,byrow=TRUE)
rownames(RACE_LPD_2_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_LPD_2_BoAA)<-c("LPD_2","Not_LPD_2")
RACE_LPD_2_BoAA <- as.table(RACE_LPD_2_BoAA)

chisq.test(RACE_LPD_2_BoAA)

#p-value = 0.00101

RACE_LPD_2_W <-matrix(c(129,639,15,225),ncol=2,byrow=TRUE)
rownames(RACE_LPD_2_W)<-c("White", "All other races")
colnames(RACE_LPD_2_W)<-c("LPD_2","LPD_2_Total")
RACE_LPD_2_W <- as.table(RACE_LPD_2_W)

chisq.test(RACE_LPD_2_W)

#p-value = 7.186e-05

pvalues <- c(1, 0.1116, 0.00101, 7.186e-05)

p.adjust(pvalues, method= "BH")

#1.00000000 0.14880000 0.00202000 0.00028744

#The above data frame is the reference for each table made below

RACE_LPD_3_AIoAN <-matrix(c(0,1,141,866),ncol=2,byrow=TRUE)
rownames(RACE_LPD_3_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_LPD_3_AIoAN)<-c("LPD_3","Not_LPD_3")
RACE_LPD_3_AIoAN <- as.table(RACE_LPD_3_AIoAN)

chisq.test(RACE_LPD_3_AIoAN)

#p-value = 1

RACE_LPD_3_ASIAN <-matrix(c(9,52,132,815),ncol=2,byrow=TRUE)
rownames(RACE_LPD_3_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_LPD_3_ASIAN)<-c("LPD_3","Not_LPD_3")
RACE_LPD_3_ASIAN <- as.table(RACE_LPD_3_ASIAN)

chisq.test(RACE_LPD_3_ASIAN)

#p-value = 1

RACE_LPD_3_BoAA <-matrix(c(20,158,121,709),ncol=2,byrow=TRUE)
rownames(RACE_LPD_3_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_LPD_3_BoAA)<-c("LPD_3","Not_LPD_3")
RACE_LPD_3_BoAA <- as.table(RACE_LPD_3_BoAA)

chisq.test(RACE_LPD_3_BoAA)

#p-value = 0.2949

RACE_LPD_3_W <-matrix(c(112,656,29,211),ncol=2,byrow=TRUE)
rownames(RACE_LPD_3_W)<-c("White", "All other races")
colnames(RACE_LPD_3_W)<-c("LPD_3","LPD_3_Total")
RACE_LPD_3_W <- as.table(RACE_LPD_3_W)

chisq.test(RACE_LPD_3_W)

#p-value = 0.3854

pvalues <- c(1, 1, 0.2949, 0.3854)

p.adjust(pvalues, method= "BH")

#1.0000 1.0000 0.7708 0.7708

#The above data frame is the reference for each table made below

RACE_LPD_4_AIoAN <-matrix(c(1,0,223,784),ncol=2,byrow=TRUE)
rownames(RACE_LPD_4_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_LPD_4_AIoAN)<-c("LPD_4","Not_LPD_4")
RACE_LPD_4_AIoAN <- as.table(RACE_LPD_4_AIoAN)

chisq.test(RACE_LPD_4_AIoAN)

#p-value = 0.5038

RACE_LPD_4_ASIAN <-matrix(c(17,44,207,740),ncol=2,byrow=TRUE)
rownames(RACE_LPD_4_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_LPD_4_ASIAN)<-c("LPD_4","Not_LPD_4")
RACE_LPD_4_ASIAN <- as.table(RACE_LPD_4_ASIAN)

chisq.test(RACE_LPD_4_ASIAN)

#p-value = 0.3495

RACE_LPD_4_BoAA <-matrix(c(73,105,119,711),ncol=2,byrow=TRUE)
rownames(RACE_LPD_4_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_LPD_4_BoAA)<-c("LPD_4","Not_LPD_4")
RACE_LPD_4_BoAA <- as.table(RACE_LPD_4_BoAA)

chisq.test(RACE_LPD_4_BoAA)

#p-value = 4.717e-16

RACE_LPD_4_W <-matrix(c(133,635,91,149),ncol=2,byrow=TRUE)
rownames(RACE_LPD_4_W)<-c("White", "All other races")
colnames(RACE_LPD_4_W)<-c("LPD_4","LPD_4_Total")
RACE_LPD_4_W <- as.table(RACE_LPD_4_W)

chisq.test(RACE_LPD_4_W)

#p-value = 3.814e-11

pvalues <- c(0.5038, 0.3495, 4.717e-16, 3.814e-11)

p.adjust(pvalues, method= "BH")

#5.0380e-01 4.6600e-01 1.8868e-15 7.6280e-11


#The above data frame is the reference for each table made below

RACE_LPD_5_AIoAN <-matrix(c(0,1,30,977),ncol=2,byrow=TRUE)
rownames(RACE_LPD_5_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_LPD_5_AIoAN)<-c("LPD_5","Not_LPD_5")
RACE_LPD_5_AIoAN <- as.table(RACE_LPD_5_AIoAN)

chisq.test(RACE_LPD_5_AIoAN)

#p-value = 1

RACE_LPD_5_ASIAN <-matrix(c(0,61,30,917),ncol=2,byrow=TRUE)
rownames(RACE_LPD_5_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_LPD_5_ASIAN)<-c("LPD_5","Not_LPD_5")
RACE_LPD_5_ASIAN <- as.table(RACE_LPD_5_ASIAN)

chisq.test(RACE_LPD_5_ASIAN)

#p-value = 0.3065

RACE_LPD_5_BoAA <-matrix(c(3,175,27,803),ncol=2,byrow=TRUE)
rownames(RACE_LPD_5_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_LPD_5_BoAA)<-c("LPD_5","Not_LPD_5")
RACE_LPD_5_BoAA <- as.table(RACE_LPD_5_BoAA)

chisq.test(RACE_LPD_5_BoAA)

#p-value = 0.3822

RACE_LPD_5_W <-matrix(c(27,741,3,237),ncol=2,byrow=TRUE)
rownames(RACE_LPD_5_W)<-c("White", "All other races")
colnames(RACE_LPD_5_W)<-c("LPD_5","LPD_5_Total")
RACE_LPD_5_W <- as.table(RACE_LPD_5_W)

chisq.test(RACE_LPD_5_W)

#p-value =0.1129

pvalues <- c(1, 0.3065, 0.3822, 0.1129)

p.adjust(pvalues, method= "BH")

# 1.0000 0.5096 0.5096 0.4516

#The above data frame is the reference for each table made below

RACE_LPD_6_AIoAN <-matrix(c(0,1,253,754),ncol=2,byrow=TRUE)
rownames(RACE_LPD_6_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_LPD_6_AIoAN)<-c("LPD_6","Not_LPD_6")
RACE_LPD_6_AIoAN <- as.table(RACE_LPD_6_AIoAN)

chisq.test(RACE_LPD_6_AIoAN)

#p-value = 1

RACE_LPD_6_ASIAN <-matrix(c(16,45,237,710),ncol=2,byrow=TRUE)
rownames(RACE_LPD_6_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_LPD_6_ASIAN)<-c("LPD_6","Not_LPD_6")
RACE_LPD_6_ASIAN <- as.table(RACE_LPD_6_ASIAN)

chisq.test(RACE_LPD_6_ASIAN)

#p-value = 0.954

RACE_LPD_6_BoAA <-matrix(c(35,143,218,612),ncol=2,byrow=TRUE)
rownames(RACE_LPD_6_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_LPD_6_BoAA)<-c("LPD_6","Not_LPD_6")
RACE_LPD_6_BoAA <- as.table(RACE_LPD_6_BoAA)

chisq.test(RACE_LPD_6_BoAA)

#p-value = 0.08043

RACE_LPD_6_W <-matrix(c(202,566,51,189),ncol=2,byrow=TRUE)
rownames(RACE_LPD_6_W)<-c("White", "All other races")
colnames(RACE_LPD_6_W)<-c("LPD_6","LPD_6_Total")
RACE_LPD_6_W <- as.table(RACE_LPD_6_W)

chisq.test(RACE_LPD_6_W)

#p-value =0.1361

pvalues <- c(1, 0.954, 0.08043, 0.1361)

p.adjust(pvalues, method= "BH")

#  1.0000 1.0000 0.2722 0.2722

#The above data frame is the reference for each table made below

RACE_LPD_7_AIoAN <-matrix(c(0,1,91,916),ncol=2,byrow=TRUE)
rownames(RACE_LPD_7_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_LPD_7_AIoAN)<-c("LPD_7","Not_LPD_7")
RACE_LPD_7_AIoAN <- as.table(RACE_LPD_7_AIoAN)

chisq.test(RACE_LPD_7_AIoAN)

#p-value = 1

RACE_LPD_7_ASIAN <-matrix(c(6,55,85,862),ncol=2,byrow=TRUE)
rownames(RACE_LPD_7_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_LPD_7_ASIAN)<-c("LPD_7","Not_LPD_7")
RACE_LPD_7_ASIAN <- as.table(RACE_LPD_7_ASIAN)

chisq.test(RACE_LPD_7_ASIAN)

#p-value = 1

RACE_LPD_7_BoAA <-matrix(c(19,159,72,758),ncol=2,byrow=TRUE)
rownames(RACE_LPD_7_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_LPD_7_BoAA)<-c("LPD_7","Not_LPD_7")
RACE_LPD_7_BoAA <- as.table(RACE_LPD_7_BoAA)

chisq.test(RACE_LPD_7_BoAA)

#p-value = 1

RACE_LPD_7_W <-matrix(c(66,702,25,215),ncol=2,byrow=TRUE)
rownames(RACE_LPD_7_W)<-c("White", "All other races")
colnames(RACE_LPD_7_W)<-c("LPD_7","LPD_7_Total")
RACE_LPD_7_W <- as.table(RACE_LPD_7_W)

chisq.test(RACE_LPD_7_W)

#p-value =0.4647

pvalues <- c(1,1,1,0.4647)

p.adjust(pvalues, method= "BH")

# 1 1 1 1

# Stage and LPD Tally -----------------------------------------------

Stage_LPD <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  group_by(Max_LPD, stage_event_pathologic_stage) %>%
  tally() %>%
  ungroup() %>%
  spread(key = Max_LPD, value = n, fill = 0) %>%
  as.data.frame(Stage_LPD) %>%
  set_colnames(c("Stage","LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7" ))

Stage_LPD


# Chi Squared For Stage I v LPD ---------------------------------------------

Stage_I_LPD_1 <-matrix(c(17,162,107,712),ncol=2,byrow=TRUE)
rownames(Stage_I_LPD_1)<-c("Stage_I", "All_other_Stages")
colnames(Stage_I_LPD_1)<-c("LPD_1","Not_LPD_1")
Stage_I_LPD_1 <- as.table(Stage_I_LPD_1)

chisq.test(Stage_I_LPD_1)

#0.2357

Stage_II_LPD_1 <-matrix(c(84,488,40, 386),ncol=2,byrow=TRUE)
rownames(Stage_II_LPD_1)<-c("Stage_II", "All_other_Stages")
colnames(Stage_II_LPD_1)<-c("LPD_1","Not_LPD_1")
Stage_II_LPD_1 <- as.table(Stage_II_LPD_1)

chisq.test(Stage_II_LPD_1)

#0.01589

Stage_III_LPD_1 <-matrix(c(23,207,101,667 ),ncol=2,byrow=TRUE)
rownames(Stage_III_LPD_1)<-c("Stage_III", "All_other_Stages")
colnames(Stage_III_LPD_1)<-c("LPD_1","Not_LPD_1")
Stage_III_LPD_1 <- as.table(Stage_III_LPD_1)

chisq.test(Stage_III_LPD_1)

#0.2473

Stage_IV_LPD_1 <-matrix(c(0,17,124, 857 ),ncol=2,byrow=TRUE)
rownames(Stage_IV_LPD_1)<-c("Stage_IV", "All_other_Stages")
colnames(Stage_IV_LPD_1)<-c("LPD_1","Not_LPD_1")
Stage_IV_LPD_1 <- as.table(Stage_IV_LPD_1)

chisq.test(Stage_IV_LPD_1)

#0.2318

pvalues <- c(0.2357, 0.2318, 0.2473, 0.01589)

p.adjust(pvalues, method= "BH")

# 0.24730 0.24730 0.24730 0.06356

# Chi Squared For Stage II v LPD ---------------------------------------------

Stage_I_LPD_2 <-matrix(c(39,140,104,715),ncol=2,byrow=TRUE)
rownames(Stage_I_LPD_2)<-c("Stage_I", "All_other_Stages")
colnames(Stage_I_LPD_2)<-c("LPD_2","Not_LPD_2")
Stage_I_LPD_2 <- as.table(Stage_I_LPD_2)

chisq.test(Stage_I_LPD_2)

#0.002474

Stage_II_LPD_2 <-matrix(c(69,503,74,352),ncol=2,byrow=TRUE)
rownames(Stage_II_LPD_2)<-c("Stage_II", "All_other_Stages")
colnames(Stage_II_LPD_2)<-c("LPD_2","Not_LPD_2")
Stage_II_LPD_2 <- as.table(Stage_II_LPD_2)

chisq.test(Stage_II_LPD_2)

#0.02285

Stage_III_LPD_2 <-matrix(c(34,196,109,659),ncol=2,byrow=TRUE)
rownames(Stage_III_LPD_2)<-c("Stage_III", "All_other_Stages")
colnames(Stage_III_LPD_2)<-c("LPD_2","Not_LPD_2")
Stage_III_LPD_2 <- as.table(Stage_III_LPD_2)

chisq.test(Stage_III_LPD_2)

#0.9071

Stage_IV_LPD_2 <-matrix(c(1,16,142,839 ),ncol=2,byrow=TRUE)
rownames(Stage_IV_LPD_2)<-c("Stage_IV", "All_other_Stages")
colnames(Stage_IV_LPD_2)<-c("LPD_2","Not_LPD_2")
Stage_IV_LPD_2 <- as.table(Stage_IV_LPD_2)

chisq.test(Stage_IV_LPD_2)

#0.5135

pvalues <- c(0.002474, 0.02285, 0.9071, 0.5135)

p.adjust(pvalues, method= "BH")

# 0.0098960 0.0457000 0.9071000 0.6846667

# Chi Squared For Stage v LPD 3 ---------------------------------------------

Stage_I_LPD_3 <-matrix(c(26,153,114,705),ncol=2,byrow=TRUE)
rownames(Stage_I_LPD_3)<-c("Stage_I", "All_other_Stages")
colnames(Stage_I_LPD_3)<-c("LPD_3","Not_LPD_3")
Stage_I_LPD_3 <- as.table(Stage_I_LPD_3)

chisq.test(Stage_I_LPD_3)

#0.9262

Stage_II_LPD_3 <-matrix(c(72,500,68,358),ncol=2,byrow=TRUE)
rownames(Stage_II_LPD_3)<-c("Stage_II", "All_other_Stages")
colnames(Stage_II_LPD_3)<-c("LPD_3","Not_LPD_3")
Stage_II_LPD_3 <- as.table(Stage_II_LPD_3)

chisq.test(Stage_II_LPD_3)

#0.1537

Stage_III_LPD_3 <-matrix(c(40,190,100,668),ncol=2,byrow=TRUE)
rownames(Stage_III_LPD_3)<-c("Stage_III", "All_other_Stages")
colnames(Stage_III_LPD_3)<-c("LPD_3","Not_LPD_3")
Stage_III_LPD_3 <- as.table(Stage_III_LPD_3)

chisq.test(Stage_III_LPD_3)

#0.1173

Stage_IV_LPD_3 <-matrix(c(2,15,138,843 ),ncol=2,byrow=TRUE)
rownames(Stage_IV_LPD_3)<-c("Stage_IV", "All_other_Stages")
colnames(Stage_IV_LPD_3)<-c("LPD_3","Not_LPD_3")
Stage_IV_LPD_3 <- as.table(Stage_IV_LPD_3)

chisq.test(Stage_IV_LPD_3)

#1

pvalues <- c(0.9262, 0.1537, 0.1173, 1)

p.adjust(pvalues, method= "BH")

#  1.0000 0.3074 0.3074 1.0000

# Chi Squared For Stage v LPD 4 ---------------------------------------------

Stage_I_LPD_4 <-matrix(c(28,151,192,627),ncol=2,byrow=TRUE)
rownames(Stage_I_LPD_4)<-c("Stage_I", "All_other_Stages")
colnames(Stage_I_LPD_4)<-c("LPD_4","Not_LPD_4")
Stage_I_LPD_4 <- as.table(Stage_I_LPD_4)

chisq.test(Stage_I_LPD_4)

#0.02917

Stage_II_LPD_4 <-matrix(c(153,419,67,359),ncol=2,byrow=TRUE)
rownames(Stage_II_LPD_4)<-c("Stage_II", "All_other_Stages")
colnames(Stage_II_LPD_4)<-c("LPD_4","Not_LPD_4")
Stage_II_LPD_4 <- as.table(Stage_II_LPD_4)

chisq.test(Stage_II_LPD_4)

# 4.565e-05

Stage_III_LPD_4 <-matrix(c(34,196,186,582),ncol=2,byrow=TRUE)
rownames(Stage_III_LPD_4)<-c("Stage_III", "All_other_Stages")
colnames(Stage_III_LPD_4)<-c("LPD_4","Not_LPD_4")
Stage_III_LPD_4 <- as.table(Stage_III_LPD_4)

chisq.test(Stage_III_LPD_4)

#0.003307

Stage_IV_LPD_4 <-matrix(c(5,12,215,766),ncol=2,byrow=TRUE)
rownames(Stage_IV_LPD_4)<-c("Stage_IV", "All_other_Stages")
colnames(Stage_IV_LPD_4)<-c("LPD_4","Not_LPD_4")
Stage_IV_LPD_4 <- as.table(Stage_IV_LPD_4)

chisq.test(Stage_IV_LPD_4)

#0.657

pvalues <- c(0.02917, 4.565e-05, 0.003307, 0.657)

p.adjust(pvalues, method= "BH")

# 0.03889333 0.00018260 0.00661400 0.65700000

# Chi Squared For Stage v LPD 5 ---------------------------------------------

Stage_I_LPD_5 <-matrix(c(5,172,24,797),ncol=2,byrow=TRUE)
rownames(Stage_I_LPD_5)<-c("Stage_I", "All_other_Stages")
colnames(Stage_I_LPD_5)<-c("LPD_5","Not_LPD_5")
Stage_I_LPD_5 <- as.table(Stage_I_LPD_5)

chisq.test(Stage_I_LPD_5)

#1

Stage_II_LPD_5 <-matrix(c(10,562,19,407),ncol=2,byrow=TRUE)
rownames(Stage_II_LPD_5)<-c("Stage_II", "All_other_Stages")
colnames(Stage_II_LPD_5)<-c("LPD_5","Not_LPD_5")
Stage_II_LPD_5 <- as.table(Stage_II_LPD_5)

chisq.test(Stage_II_LPD_5)

# 0.01969

Stage_III_LPD_5 <-matrix(c(12,218,17,751),ncol=2,byrow=TRUE)
rownames(Stage_III_LPD_5)<-c("Stage_III", "All_other_Stages")
colnames(Stage_III_LPD_5)<-c("LPD_5","Not_LPD_5")
Stage_III_LPD_5 <- as.table(Stage_III_LPD_5)

chisq.test(Stage_III_LPD_5)

#0.03113

Stage_IV_LPD_5 <-matrix(c(2,15,27,954),ncol=2,byrow=TRUE)
rownames(Stage_IV_LPD_5)<-c("Stage_IV", "All_other_Stages")
colnames(Stage_IV_LPD_5)<-c("LPD_5","Not_LPD_5")
Stage_IV_LPD_5 <- as.table(Stage_IV_LPD_5)

chisq.test(Stage_IV_LPD_5)

#0.1429

pvalues <- c(1, 0.01969, 0.03113, 0.1429)

p.adjust(pvalues, method= "BH")

# 1.0000000 0.0622600 0.0622600 0.1905333

# Chi Squared For Stage v LPD 6 ---------------------------------------------

Stage_I_LPD_6 <-matrix(c(55,124,197,622),ncol=2,byrow=TRUE)
rownames(Stage_I_LPD_6)<-c("Stage_I", "All_other_Stages")
colnames(Stage_I_LPD_6)<-c("LPD_6","Not_LPD_6")
Stage_I_LPD_6 <- as.table(Stage_I_LPD_6)

chisq.test(Stage_I_LPD_6)

#0.07731

Stage_II_LPD_6 <-matrix(c(130,442,122,304),ncol=2,byrow=TRUE)
rownames(Stage_II_LPD_6)<-c("Stage_II", "All_other_Stages")
colnames(Stage_II_LPD_6)<-c("LPD_6","Not_LPD_6")
Stage_II_LPD_6 <- as.table(Stage_II_LPD_6)

chisq.test(Stage_II_LPD_6)

# 0.04013

Stage_III_LPD_6 <-matrix(c(66,164,186,582),ncol=2,byrow=TRUE)
rownames(Stage_III_LPD_6)<-c("Stage_III", "All_other_Stages")
colnames(Stage_III_LPD_6)<-c("LPD_6","Not_LPD_6")
Stage_III_LPD_5 <- as.table(Stage_III_LPD_6)

chisq.test(Stage_III_LPD_6)

#0.199

Stage_IV_LPD_6 <-matrix(c(1,16,251,730),ncol=2,byrow=TRUE)
rownames(Stage_IV_LPD_6)<-c("Stage_IV", "All_other_Stages")
colnames(Stage_IV_LPD_6)<-c("LPD_6","Not_LPD_5")
Stage_IV_LPD_5 <- as.table(Stage_IV_LPD_6)

chisq.test(Stage_IV_LPD_6)

#0.1158

pvalues <- c(0.07731, 0.04013, 0.119, 0.1158)

p.adjust(pvalues, method= "BH")

# 0.119 0.119 0.119 0.119

# Chi Squared For Stage v LPD 7 ---------------------------------------------

Stage_I_LPD_7 <-matrix(c(9,170,81,738),ncol=2,byrow=TRUE)
rownames(Stage_I_LPD_7)<-c("Stage_I", "All_other_Stages")
colnames(Stage_I_LPD_7)<-c("LPD_7","Not_LPD_7")
Stage_I_LPD_7 <- as.table(Stage_I_LPD_7)

chisq.test(Stage_I_LPD_7)

#0.05571

Stage_II_LPD_7 <-matrix(c(54,518,36,390),ncol=2,byrow=TRUE)
rownames(Stage_II_LPD_7)<-c("Stage_II", "All_other_Stages")
colnames(Stage_II_LPD_7)<-c("LPD_7","Not_LPD_7")
Stage_II_LPD_7 <- as.table(Stage_II_LPD_7)

chisq.test(Stage_II_LPD_7)

# 0.6685

Stage_III_LPD_7 <-matrix(c(21,209,69,699),ncol=2,byrow=TRUE)
rownames(Stage_III_LPD_7)<-c("Stage_III", "All_other_Stages")
colnames(Stage_III_LPD_7)<-c("LPD_7","Not_LPD_7")
Stage_III_LPD_7 <- as.table(Stage_III_LPD_7)

chisq.test(Stage_III_LPD_7)

#1

Stage_IV_LPD_7 <-matrix(c(6,11,84,897),ncol=2,byrow=TRUE)
rownames(Stage_IV_LPD_7)<-c("Stage_IV", "All_other_Stages")
colnames(Stage_IV_LPD_7)<-c("LPD_7","Not_LPD_7")
Stage_IV_LPD_7 <- as.table(Stage_IV_LPD_7)

chisq.test(Stage_IV_LPD_7)

#0.0007044

pvalues <- c(0.05571, 0.6685, 1, 0.0007044)

p.adjust(pvalues, method= "BH")

# 0.1114200 0.8913333 1.0000000 0.0028176

# Filtering for Subtypes ----------------------------------------

Breast_unknown_ER_HER2_PR <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_estrogen_receptor_status %in% c("Positive", "Negative")) %>%
  filter(lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Positive", "Negative")) %>%
  filter(breast_carcinoma_progesterone_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(Known_Subtypes = case_when(
    breast_carcinoma_estrogen_receptor_status == "Positive" | breast_carcinoma_progesterone_receptor_status == "Positive" & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative" ~ "Luminal A",
    breast_carcinoma_estrogen_receptor_status == "Negative" & breast_carcinoma_progesterone_receptor_status == "Negative" & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive" ~ "HER2_Over_Expression",
    breast_carcinoma_estrogen_receptor_status == "Negative" & breast_carcinoma_progesterone_receptor_status == "Negative" & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative" ~ "Basal",
    breast_carcinoma_estrogen_receptor_status == "Positive" | breast_carcinoma_progesterone_receptor_status == "Positive" & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive" ~ "Luminal B", 
    TRUE ~ "Unknown"))

Breast_unknown_ER_HER2_PR$Known_Subtypes


# Subtypes and LPD Table --------------------------------------------------

Subtype_LPD <- Breast_unknown_ER_HER2_PR %>%
  group_by(Max_LPD, Known_Subtypes) %>%
  tally() %>%
  ungroup() %>%
  spread(key = Max_LPD, value = n, fill = 0) %>%
  as.data.frame(Subtype_LPD) %>% 
set_colnames(c("Known_Subtypes","LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7" ))

Subtype_LPD


# Chi Squared -------------------------------------------------------------

Luminal_LPD_1 <-matrix(c(68,456,0,112),ncol=2,byrow=TRUE)
rownames(Luminal_LPD_1)<-c("Luminal", "Basal")
colnames(Luminal_LPD_1)<-c("LPD_1","Not_LPD_1")
Luminal_LPD_1 <- as.table(Luminal_LPD_1)

chisq.test(Luminal_LPD_1)

pvalues <- c(0.0001108, 0.0001108)

p.adjust(pvalues, method= "BH")

#0.0001108 0.0001108

Luminal_LPD_2 <-matrix(c(96,428,2,110),ncol=2,byrow=TRUE)
rownames(Luminal_LPD_2)<-c("Luminal", "Basal")
colnames(Luminal_LPD_2)<-c("LPD_2","Not_LPD_2")
Luminal_LPD_2 <- as.table(Luminal_LPD_2)

chisq.test(Basal_LPD_2)

pvalues <- c(2.088e-05, 2.088e-05)

p.adjust(pvalues, method= "BH")

#2.088e-05 2.088e-05

Luminal_LPD_3 <-matrix(c(96,428,2,110),ncol=2,byrow=TRUE)
rownames(Luminal_LPD_3)<-c("Luminal", "Basal")
colnames(Luminal_LPD_3)<-c("LPD_3","Not_LPD_3")
Luminal_LPD_3 <- as.table(Luminal_LPD_3)

pvalues <- c(2.088e-05, 2.088e-05)

p.adjust(pvalues, method= "BH")

#2.088e-05 2.088e-05

Luminal_LPD_4 <-matrix(c(49,475,100,12),ncol=2,byrow=TRUE)
rownames(Luminal_LPD_4)<-c("Luminal", "Basal")
colnames(Luminal_LPD_4)<-c("LPD_4","Not_LPD_4")
Luminal_LPD_4 <- as.table(Luminal_LPD_4)

chisq.test(Luminal_LPD_4)

pvalues <- c(2.2e-16,2.2e-16)

p.adjust(pvalues, method= "BH")

#2.2e-16 2.2e-16

Luminal_LPD_5 <-matrix(c(12,512,4,108),ncol=2,byrow=TRUE)
rownames(Luminal_LPD_5)<-c("Luminal", "Basal")
colnames(Luminal_LPD_5)<-c("LPD_5","Not_LPD_5")
Luminal_LPD_5 <- as.table(Luminal_LPD_5)

chisq.test(Luminal_LPD_5)

pvalues <- c(0.6501, 0.6501)

p.adjust(pvalues, method= "BH")

#0.6501, 0.6501

Luminal_LPD_6 <-matrix(c(153,371,2,110),ncol=2,byrow=TRUE)
rownames(Luminal_LPD_6)<-c("Luminal", "Basal")
colnames(Luminal_LPD_6)<-c("LPD_6","Not_LPD_6")
Luminal_LPD_6 <- as.table(Luminal_LPD_6)

chisq.test(Luminal_LPD_6)

pvalues <- c(1.828e-09, 1.828e-09)

p.adjust(pvalues, method= "BH")

#1.828e-09 1.828e-09

Luminal_LPD_7 <-matrix(c(50,474,2,110),ncol=2,byrow=TRUE)
rownames(Luminal_LPD_7)<-c("Luminal", "Basal")
colnames(Luminal_LPD_7)<-c("LPD_7","Not_LPD_7")
Luminal_LPD_7 <- as.table(Luminal_LPD_7)

chisq.test(Luminal_LPD_7)

pvalues <- c( 0.01143,  0.01143)

p.adjust(pvalues, method= "BH")

#0.01143 0.01143




