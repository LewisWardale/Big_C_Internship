
# Importing Data ----------------------------------------------------------

breast_data <- read.csv("/Users/lewiswardale/Desktop/TCGA-BRCA/clinical_data_TCGA-BRCA_eridanus.csv")
names(breast_data)
# install.packages("survival")
# install.packages("survminer")
# install.packages("lubridate")
# install.packages("MASS")
# install.packages("DataExplorer")
library(survminer)
library(survival)
library(lubridate)
library(dplyr)
library(MASS)
library(tidyr)
library(DataExplorer)
library(magrittr)
library(ggplot2)


# Cleaning all data ----------------------------------------------------------------

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
  dplyr::filter(!race_list %in% c("", "AMERICAN INDIAN OR ALASKA NATIVE")) %>%
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

Age_v_Stage <- ggplot(Breast_Summary_Stages_Filtered_Blanks_Age_Race, aes(x = stage_event_pathologic_stage, y = Age, fill = stage_event_pathologic_stage)) + 
  geom_boxplot(outlier.colour = "red")

Age_v_Stage + ggtitle("Box Plot of Age and Stage") + labs(y = "Age (Years)", x = "Pathological Stage") + labs(fill = "Pathological Stage") 

Age_v_Stage_Race <- Age_v_Stage <- ggplot(Breast_Summary_Stages_Filtered_Blanks_Age_Race, aes(x = stage_event_pathologic_stage, y = Age, fill = race_list )) + 
  geom_boxplot(outlier.colour = "red")

Age_v_Stage_Race + ggtitle("Box Plot of Age and Stage and Race") + labs(y = "Age (Years)", x = "Pathological Stage") + labs(fill = "Race")


# Box plot Age vs Receptor Status -----------------------------------------

Breast_unknown_PR <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_progesterone_receptor_status %in% c("Positive", "Negative")) 

Age_v_PR <- ggplot(Breast_unknown_PR, aes(x = breast_carcinoma_progesterone_receptor_status, y = Age, fill = breast_carcinoma_progesterone_receptor_status )) + 
  geom_boxplot(outlier.colour = "red")

Age_v_PR + ggtitle("Box Plot of Age and PR Status") + labs(y = "Age (Years)", x = "PR Status") + labs(fill = "PR Status") 

Breast_unknown_ER <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_estrogen_receptor_status %in% c("Positive", "Negative")) 

Age_v_ER <- ggplot(Breast_unknown_ER, aes(x = breast_carcinoma_estrogen_receptor_status, y = Age, fill = breast_carcinoma_estrogen_receptor_status )) + 
  geom_boxplot(outlier.colour = "red")

Age_v_ER + ggtitle("Box Plot of Age and ER Status") + labs(y = "Age (Years)", x = "ER Status") + labs(fill = "ER Status") 

Breast_unknown_HER2 <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Positive", "Negative"))

Age_v_HER2 <- ggplot(Breast_unknown_HER2, aes(x = lab_proc_her2_neu_immunohistochemistry_receptor_status, y = Age, fill = lab_proc_her2_neu_immunohistochemistry_receptor_status )) + 
  geom_boxplot(outlier.colour = "red")

Age_v_HER2 + ggtitle("Box Plot of Age and HER2 Status") + labs(y = "Age (Years)", x = "HER2 Status") + labs(fill = "HER2 Status") 


# Insigificant Pairwise ---------------------------------------------------


# Pairwise Race -----------------------------------------------------------

res <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ race_list,
                         data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
res

                       #AMERICAN INDIAN OR ALASKA NATIVE ASIAN BLACK OR AFRICAN AMERICAN
#ASIAN                     1                                -     -                        
#BLACK OR AFRICAN AMERICAN 1                                1     -                        
#WHITE                     1                                1     1              


# Pairwise Age ------------------------------------------------------------

res2 <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ age_group,
                          data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
res2

summary(Breast_Summary_Stages_Filtered_Blanks_Age_Race$Age)

#   Min.   1st Qu.  Median  Mean   3rd Qu.  Max.    NA's 
#  26.59   49.07   58.60   58.61   67.66   90.06    15 

              #Above_Upper_Q Below_Lower_Q Lower_Quartile
#Below_Lower_Q  0.0111        -             -             
#Lower_Quartile 0.0017        0.5231        -             
#Upper_Quartile 0.0124        0.9238        0.5231  

# Pairwise for PR ---------------------------------------------------------

res3 <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ breast_carcinoma_progesterone_receptor_status,
                          data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
res3

          # Negative
#Positive   0.057  


# Pairwise loop for Stage v Race ------------------------------------------

RaceStageVector <- as.character(unique(Breast_Summary_Stages_Filtered_Blanks_Age_Race$stage_event_pathologic_stage))

for (Stage in RaceStageVector) {
  
#Filters for each Stage
  StagevRace <- Breast_Summary_Stages_Filtered_Blanks_Age_Race  %>%
    dplyr::filter(stage_event_pathologic_stage == Stage)

  #Pairwise for each stage vs Race
Pairwise <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ race_list,
                               data = StagevRace)

  #Symbolic number coding for P-Value
if(Pairwise$p.value > 0.05) {
  next()
}

print(Pairwise)

surv_object <- Surv(time = StagevRace$Years_to_death_and_last_followup, event = StagevRace$vital_status_TRUE_FALSE)
surv_object
fit1 <- survfit(surv_object ~ race_list, data = StagevRace) 
fit1

Stage_vs_Race <- ggsurvplot(
  fit1,
  data = StagevRace, 
  pval = TRUE, legend.title = paste0(" ", Stage, " vs race"), 
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs =
    c("Black or African American", "White"),
  ggtheme = theme_bw(),
  break.time.by = 2
) + xlab("Years")  

Stage_vs_Race$plot <- Stage_vs_Race$plot + labs(
  title = paste0 ("KM overall", Stage, " vs Race status plot for TCGA-BRCA"))

Stage_vs_Race <- ggpar(
  Stage_vs_Race,
  font.title = c(16, "bold", "darkblue"),
  font.x = c(14, "plain", "black"),
  font.y = c(14, "plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab =  c(12, "plain", "black"),
  legend = "bottom"
)

Stage_vs_Race

pdf(paste0("StagevRace/KM plot of", Stage, " and_Race.pdf"),
            onefile = FALSE)
print(Stage_vs_Race)
dev.off()

}

# Significant Pairwise ----------------------------------------------------

# Pairwise ER Status ------------------------------------------------------

res5 <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ breast_carcinoma_estrogen_receptor_status,
                          data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
res5

                #Negative
#Positive       0.0068  

# KM plot for ER Status ---------------------------------------------------

Breast_unknown_ER <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_estrogen_receptor_status %in% c("Positive", "Negative")) 


surv_object <- Surv(time = Breast_unknown_ER$Years_to_death_and_last_followup, event = Breast_unknown_ER$vital_status_TRUE_FALSE)
surv_object
fit1 <- survfit(surv_object ~ breast_carcinoma_estrogen_receptor_status, data = Breast_unknown_ER) 
fit1

ER_status <- ggsurvplot(
  fit1,
  data = Breast_unknown_ER, 
  pval = TRUE, legend.title = "Breast Carcinoma Oestrogen Receptor Status", 
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs =
    c("Negative", "Positive"),
  ggtheme = theme_bw(),
  break.time.by = 2
) + xlab("Years")  

ER_status$plot <- ER_status$plot + labs(
  title = "KM overall breast carcinoma oestrogen receptor status plot for TCGA-BRCA")

ER_status <- ggpar(
  ER_status,
  font.title = c(16, "bold", "darkblue"),
  font.x = c(14, "plain", "black"),
  font.y = c(14, "plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab =  c(12, "plain", "black"),
  legend = "bottom"
)

ER_status

# Pairwise HER2 -----------------------------------------------------------

res5 <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ lab_proc_her2_neu_immunohistochemistry_receptor_status,
                          data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
res5

            # Negative
# Positive    0.0098  

# KM HER2 -----------------------------------------------------------------

Breast_unknown_HER2 <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Positive", "Negative")) 

surv_object <- Surv(time = Breast_unknown_HER2$Years_to_death_and_last_followup, event = Breast_unknown_HER2$vital_status_TRUE_FALSE)
surv_object
fit1 <- survfit(surv_object ~ lab_proc_her2_neu_immunohistochemistry_receptor_status, data = Breast_unknown_HER2) 
fit1

Positive_v_Negative <- ggsurvplot(
  fit1,
  data = Breast_unknown_HER2, 
  pval = TRUE, legend.title = "Human Epidermal Growth Factor Receptor Status", 
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs =
    c("Negative", "Positive"),
  ggtheme = theme_bw(),
  break.time.by = 2
) + xlab("Years")  

Positive_v_Negative$plot <- Positive_v_Negative$plot + labs(
  title = "KM human epidermal growth factor receptor status plot for TCGA-BRCA, positive and negative")

Positive_v_Negative <- ggpar(
  Positive_v_Negative,
  font.title = c(16, "bold", "darkblue"),
  font.x = c(14, "plain", "black"),
  font.y = c(14, "plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab =  c(12, "plain", "black"),
  legend = "bottom"
)

Positive_v_Negative

# Pairwise Stage ----------------------------------------------------------

StageVector <- as.character(unique(Breast_Summary_Stages_Filtered_Blanks_Age_Race$stage_event_pathologic_stage))

for (Stage in StageVector) {
  
  #Pairwise for each stage vs Race
  Pairwise <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ stage_event_pathologic_stage,
                                data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
  
  #Symbolic number coding for P-Value
  if(Pairwise$p.value > 0.05) {
    next()
  }
  
  colvector <- colnames(Pairwise)
  rowvector <- rownames(Pairwise)
  
  surv_object <- Surv(time = StagevRace$Years_to_death_and_last_followup, event = StagevRace$vital_status_TRUE_FALSE)
  surv_object
  fit1 <- survfit(surv_object ~ race_list, data = StagevRace) 
  fit1
  
  Stage_vs_Race <- ggsurvplot(
    fit1,
    data = StagevRace, 
    pval = TRUE, legend.title = paste0(" ", Stage, " vs race"), 
    risk.table = TRUE,
    risk.table.col = "strata",
    legend.labs =
      c("Black or African American", "White"),
    ggtheme = theme_bw(),
    break.time.by = 2
  ) + xlab("Years")  
  
  Stage_vs_Race$plot <- Stage_vs_Race$plot + labs(
    title = paste0 ("KM overall", Stage, " vs Race status plot for TCGA-BRCA"))
  
  Stage_vs_Race <- ggpar(
    Stage_vs_Race,
    font.title = c(16, "bold", "darkblue"),
    font.x = c(14, "plain", "black"),
    font.y = c(14, "plain", "black"),
    font.xtickslab = c(12, "plain", "black"),
    font.ytickslab =  c(12, "plain", "black"),
    legend = "bottom"
  )
  
  Stage_vs_Race
  
  pdf(paste0("StagevRace/KM plot of", Stage, " and_Race.pdf"),
      onefile = FALSE)
  print(Stage_vs_Race)
  dev.off()
  
}

res4 <- pairwise_survdiff(Surv(Years_to_death_and_last_followup, event = vital_status_TRUE_FALSE) ~ stage_event_pathologic_stage,
                          data = Breast_Summary_Stages_Filtered_Blanks_Age_Race)
res4

            #Stage I  Stage II Stage III
#Stage II  0.5010  -        -        
#Stage III 0.0098  0.0098   -        
#Stage IV  5.4e-07 1.1e-06  0.0098   

  
# KM for Stages -----------------------------------------------------------

 surv_object <- Surv(time = Breast_Summary_Stages_Filtered_Blanks_Age_Race$Years_to_death_and_last_followup, event = Breast_Summary_Stages_Filtered_Blanks_Age_Race$vital_status_TRUE_FALSE)
 surv_object
 fit1 <- survfit(surv_object ~ stage_event_pathologic_stage, data = Breast_Summary_Stages_Filtered_Blanks_Age_Race) 
 fit1
 
 Pathologic_stage <- ggsurvplot(
   fit1,
   data = Breast_Summary_Stages_Filtered_Blanks_Age_Race, 
   pval = TRUE, legend.title = "Overall Pathologic Stage", 
   risk.table = TRUE,
   risk.table.col = "strata",
   legend.labs =
     c("I","II","III","IV"),
   ggtheme = theme_bw(),
   break.time.by = 2,
   tables.height = 0.4
 ) + xlab("Years")  
 
 Pathologic_stage$plot <- Pathologic_stage$plot + labs(
   title = "KM overall pathologic stages plot for TCGA-BRCA")
 
 Pathologic_stage <- ggpar(
   Pathologic_stage,
   font.title = c(16, "bold", "darkblue"),
   font.x = c(14, "plain", "black"),
   font.y = c(14, "plain", "black"),
   font.xtickslab = c(12, "plain", "black"),
   font.ytickslab =  c(12, "plain", "black"),
   legend = "bottom"
 )
 
 Pathologic_stage
 

# Filtering for I and II --------------------------------------------------

 Breast_unknown_I_II <-  Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
   dplyr::filter(stage_event_pathologic_stage %in% c("Stage I", "Stage II"))
 
 Breast_unknown_I_II$stage_event_pathologic_stage

# KM I vs II --------------------------------------------------------------

 surv_object <- Surv(time = Breast_unknown_I_II$Years_to_death_and_last_followup, event =  Breast_unknown_I_II$vital_status_TRUE_FALSE)
 surv_object
 fit1 <- survfit(surv_object ~ stage_event_pathologic_stage, data = Breast_unknown_I_II) 
 fit1
 
 Pathologic_stage <- ggsurvplot(
   fit1,
   data = Breast_unknown_I_II, 
   pval = TRUE, legend.title = "Pathologic Stage", 
   risk.table = TRUE,
   risk.table.col = "strata",
   legend.labs =
     c("I","II"),
   ggtheme = theme_bw(),
   break.time.by = 2,
   tables.height = 0.4
 ) + xlab("Years")  
 
 Pathologic_stage$plot <- Pathologic_stage$plot + labs(
   title = "KM pathologic stages plot for TCGA-BRCA, I and II")
 
 Pathologic_stage <- ggpar(
   Pathologic_stage,
   font.title = c(16, "bold", "darkblue"),
   font.x = c(14, "plain", "black"),
   font.y = c(14, "plain", "black"),
   font.xtickslab = c(12, "plain", "black"),
   font.ytickslab =  c(12, "plain", "black"),
   legend = "bottom"
 )
 
 Pathologic_stage
 
 # Filtering for I and III --------------------------------------------------
 
 Breast_unknown_I_III <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
   dplyr::filter(stage_event_pathologic_stage %in% c("Stage I", "Stage III"))
 
 Breast_unknown_I_III$stage_event_pathologic_stage
 
 # KM I vs III --------------------------------------------------------------
 
 surv_object <- Surv(time = Breast_unknown_I_III$Years_to_death_and_last_followup, event = Breast_unknown_I_III$vital_status_TRUE_FALSE)
 surv_object
 fit1 <- survfit(surv_object ~ stage_event_pathologic_stage, data = Breast_unknown_I_III) 
 fit1
 
 Pathologic_stage <- ggsurvplot(
   fit1,
   data = Breast_unknown_I_III, 
   pval = TRUE, legend.title = "Pathologic Stage", 
   risk.table = TRUE,
   risk.table.col = "strata",
   legend.labs =
     c("I","III"),
   ggtheme = theme_bw(),
   break.time.by = 2,
   tables.height = 0.4
 ) + xlab("Years")  
 
 Pathologic_stage$plot <- Pathologic_stage$plot + labs(
   title = "KM pathologic stages plot for TCGA-BRCA, I and III")
 
 Pathologic_stage <- ggpar(
   Pathologic_stage,
   font.title = c(16, "bold", "darkblue"),
   font.x = c(14, "plain", "black"),
   font.y = c(14, "plain", "black"),
   font.xtickslab = c(12, "plain", "black"),
   font.ytickslab =  c(12, "plain", "black"),
   legend = "bottom"
 )
 
 Pathologic_stage
 
 # Filtering for I and IV --------------------------------------------------
 
 Breast_unknown_I_IV <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
   dplyr::filter(stage_event_pathologic_stage %in% c("Stage I", "Stage IV"))
 
 Breast_unknown_I_IV$stage_event_pathologic_stage
 
 # KM I vs IV --------------------------------------------------------------
 
 surv_object <- Surv(time = Breast_unknown_I_IV$Years_to_death_and_last_followup, event = Breast_unknown_I_IV$vital_status_TRUE_FALSE)
 surv_object
 fit1 <- survfit(surv_object ~ stage_event_pathologic_stage, data = Breast_unknown_I_IV) 
 fit1
 
 Pathologic_stage <- ggsurvplot(
   fit1,
   data = Breast_unknown_I_IV, 
   pval = TRUE, legend.title = "Pathologic Stage", 
   risk.table = TRUE,
   risk.table.col = "strata",
   legend.labs =
     c("I","IV"),
   ggtheme = theme_bw(),
   break.time.by = 2,
   tables.height = 0.4
 ) + xlab("Years")  
 
 Pathologic_stage$plot <- Pathologic_stage$plot + labs(
   title = "KM pathologic stages plot for TCGA-BRCA, I and IV")
 
 Pathologic_stage <- ggpar(
   Pathologic_stage,
   font.title = c(16, "bold", "darkblue"),
   font.x = c(14, "plain", "black"),
   font.y = c(14, "plain", "black"),
   font.xtickslab = c(12, "plain", "black"),
   font.ytickslab =  c(12, "plain", "black"),
   legend = "bottom"
 )
 
 Pathologic_stage
 
 # Filtering for II and III --------------------------------------------------
 
 Breast_unknown_II_III <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
   dplyr::filter(stage_event_pathologic_stage %in% c("Stage II", "Stage III"))
 
 Breast_unknown_II_III$stage_event_pathologic_stage
 
 # KM II vs III --------------------------------------------------------------
 
 surv_object <- Surv(time = Breast_unknown_II_III$Years_to_death_and_last_followup, event = Breast_unknown_II_III$vital_status_TRUE_FALSE)
 surv_object
 fit1 <- survfit(surv_object ~ stage_event_pathologic_stage, data = Breast_unknown_II_III) 
 fit1
 
 Pathologic_stage <- ggsurvplot(
   fit1,
   data = Breast_unknown_II_III, 
   pval = TRUE, legend.title = "Pathologic Stage", 
   risk.table = TRUE,
   risk.table.col = "strata",
   legend.labs =
     c("II","III"),
   ggtheme = theme_bw(),
   break.time.by = 2,
   tables.height = 0.4
 ) + xlab("Years")  
 
 Pathologic_stage$plot <- Pathologic_stage$plot + labs(
   title = "KM pathologic stages plot for TCGA-BRCA, II and III")
 
 Pathologic_stage <- ggpar(
   Pathologic_stage,
   font.title = c(16, "bold", "darkblue"),
   font.x = c(14, "plain", "black"),
   font.y = c(14, "plain", "black"),
   font.xtickslab = c(12, "plain", "black"),
   font.ytickslab =  c(12, "plain", "black"),
   legend = "bottom"
 )
 
 Pathologic_stage
 
 # Filtering for II and IV --------------------------------------------------
 
 Breast_unknown_II_IV <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
   dplyr::filter(stage_event_pathologic_stage %in% c("Stage II", "Stage IV"))
 
 Breast_unknown_II_IV$stage_event_pathologic_stage
 
 # KM II vs IV --------------------------------------------------------------
 
 surv_object <- Surv(time = Breast_unknown_II_IV$Years_to_death_and_last_followup, event = Breast_unknown_II_IV$vital_status_TRUE_FALSE)
 surv_object
 fit1 <- survfit(surv_object ~ stage_event_pathologic_stage, data = Breast_unknown_II_IV) 
 fit1
 
 Pathologic_stage <- ggsurvplot(
   fit1,
   data = Breast_unknown_II_IV, 
   pval = TRUE, legend.title = "Pathologic Stage", 
   risk.table = TRUE,
   risk.table.col = "strata",
   legend.labs =
     c("II","IV"),
   ggtheme = theme_bw(),
   break.time.by = 2,
   tables.height = 0.4
 ) + xlab("Years")  
 
 Pathologic_stage$plot <- Pathologic_stage$plot + labs(
   title = "KM pathologic stages plot for TCGA-BRCA, II and IV")
 
 Pathologic_stage <- ggpar(
   Pathologic_stage,
   font.title = c(16, "bold", "darkblue"),
   font.x = c(14, "plain", "black"),
   font.y = c(14, "plain", "black"),
   font.xtickslab = c(12, "plain", "black"),
   font.ytickslab =  c(12, "plain", "black"),
   legend = "bottom"
 )
 
 Pathologic_stage
 
 # Filtering for IV and III --------------------------------------------------
 
 Breast_unknown_III_IV <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
   dplyr::filter(stage_event_pathologic_stage %in% c("Stage III", "Stage IV"))
 
 Breast_unknown_III_IV$stage_event_pathologic_stage
 
 # KM III vs IV --------------------------------------------------------------
 
 surv_object <- Surv(time = Breast_unknown_III_IV$Years_to_death_and_last_followup, event = Breast_unknown_III_IV$vital_status_TRUE_FALSE)
 surv_object
 fit1 <- survfit(surv_object ~ stage_event_pathologic_stage, data = Breast_unknown_III_V) 
 fit1
 
 Pathologic_stage <- ggsurvplot(
   fit1,
   data = Breast_unknown_III_IV, 
   pval = TRUE, legend.title = "Pathologic Stage", 
   risk.table = TRUE,
   risk.table.col = "strata",
   legend.labs =
     c("III","IV"),
   ggtheme = theme_bw(),
   break.time.by = 2,
   tables.height = 0.4
 ) + xlab("Years")  
 
 Pathologic_stage$plot <- Pathologic_stage$plot + labs(
   title = "KM pathologic stages plot for TCGA-BRCA, III and IV")
 
 Pathologic_stage <- ggpar(
   Pathologic_stage,
   font.title = c(16, "bold", "darkblue"),
   font.x = c(14, "plain", "black"),
   font.y = c(14, "plain", "black"),
   font.xtickslab = c(12, "plain", "black"),
   font.ytickslab =  c(12, "plain", "black"),
   legend = "bottom"
 )
 
 Pathologic_stage
 
# Chi Squared Stage I Races -----------------------------------------------

Race_Stage <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  group_by(stage_event_pathologic_stage, race_list) %>%
  tally() %>%
  ungroup() %>%
  spread(key = stage_event_pathologic_stage, value = n, fill = 0) %>%
  as.data.frame(Race_Stage) %>%
  set_colnames(c("Race_list","Stage_I", "Stage_II", "Stage_III", "Stage_IV"))

Race_Stage

#The above data frame is the reference for each table made below

RACE_STAGE_I_AIoAN <-matrix(c(0,1,176,802),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_I_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_STAGE_I_AIoAN)<-c("Stage_I","Not_Stage_1")
RACE_STAGE_I_AIoAN <- as.table(RACE_STAGE_I_AIoAN)

chisq.test(RACE_STAGE_I_AIoAN)

#p-value = 1

RACE_STAGE_I_ASIAN <-matrix(c(4,57,172,746),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_I_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_STAGE_I_ASIAN)<-c("Stage_I","Not_Stage_I")
RACE_STAGE_I_ASIAN <- as.table(RACE_STAGE_I_ASIAN)

chisq.test(RACE_STAGE_I_ASIAN)

#p-value = 0.02598

RACE_STAGE_I_BoAA <-matrix(c(32,145,144,658),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_I_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_STAGE_I_BoAA)<-c("Stage_I","Not_Stage_I")
RACE_STAGE_I_BoAA <- as.table(RACE_STAGE_I_BoAA)

chisq.test(RACE_STAGE_I_BoAA)

#p-value = 1

RACE_STAGE_I_W <-matrix(c(140,600,36,203),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_I_W)<-c("White", "All other races")
colnames(RACE_STAGE_I_W)<-c("Stage_I","Stage_I_Total")
RACE_STAGE_I_W <- as.table(RACE_STAGE_I_W)

chisq.test(RACE_STAGE_I_W)

#p-value = 0.2103

pvalues <- c(1, 0.02598, 1, 0.2103)

p.adjust(pvalues, method= "BH")

#1.00000 0.10392 1.00000 0.42060

# Chi Squared Stage II Races ----------------------------------------------------------

RACE_STAGE_II_AIoAN <-matrix(c(0,1,562,416),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_II_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_STAGE_II_AIoAN)<-c("Stage_II","Not_Stage_II")
RACE_STAGE_II_AIoAN <- as.table(RACE_STAGE_II_AIoAN)

chisq.test(RACE_STAGE_II_AIoAN)

#p-value = 0.8809

RACE_STAGE_II_ASIAN <-matrix(c(43,18,519,399),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_II_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_STAGE_II_ASIAN)<-c("Stage_II","Not_Stage_II")
RACE_STAGE_II_ASIAN <- as.table(RACE_STAGE_II_ASIAN)

chisq.test(RACE_STAGE_II_ASIAN)

#p-value = 0.04541

RACE_STAGE_II_BoAA <-matrix(c(106,71,456,346),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_II_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_STAGE_II_BoAA)<-c("Stage_II","Not_Stage_II")
RACE_STAGE_II_BoAA <- as.table(RACE_STAGE_II_BoAA)

chisq.test(RACE_STAGE_II_BoAA)

#p-value = 0.5133

RACE_STAGE_II_W <-matrix(c(413,327,149,90),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_II_W)<-c("White", "All other races")
colnames(RACE_STAGE_II_W)<-c("Stage_II","Not_Stage_II_Total")
RACE_STAGE_II_W <- as.table(RACE_STAGE_II_W)

chisq.test(RACE_STAGE_II_W)

#p-value = 0.08907

pvalues <- c(0.8809, 0.04541, 0.5133, 0.08907)

p.adjust(pvalues, method= "BH")

#0.88090 0.17814 0.68440 0.17814


# Chi Squared Stage III Races ---------------------------------------------------

RACE_STAGE_III_AIoAN <-matrix(c(1,0,224,754),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_III_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_STAGE_III_AIoAN)<-c("Stage_III","Not_Stage_III")
RACE_STAGE_III_AIoAN <- as.table(RACE_STAGE_III_AIoAN)

chisq.test(RACE_STAGE_III_AIoAN)

#p-value = 0.5206

RACE_STAGE_III_ASIAN <-matrix(c(14,47,210,708),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_III_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_STAGE_II_ASIAN)<-c("Stage_III","Not_Stage_III")
RACE_STAGE_III_ASIAN <- as.table(RACE_STAGE_III_ASIAN)

chisq.test(RACE_STAGE_III_ASIAN)

#p-value = 1

RACE_STAGE_III_BoAA <-matrix(c(34,143,191,611),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_III_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_STAGE_III_BoAA)<-c("Stage_III","Not_Stage_III")
RACE_STAGE_III_BoAA <- as.table(RACE_STAGE_III_BoAA)

chisq.test(RACE_STAGE_III_BoAA)

#p-value = 0.2226

RACE_STAGE_III_W <-matrix(c(176,564,49,190),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_III_W)<-c("White", "All other races")
colnames(RACE_STAGE_III_W)<-c("Stage_III","Not_Stage_III_Total")
RACE_STAGE_III_W <- as.table(RACE_STAGE_III_W)

chisq.test(RACE_STAGE_III_W)

#p-value = 0.3371

pvalues <- c(0.5206, 1, 0.2226, 0.3371)

p.adjust(pvalues, method= "BH")

#0.6941333 1.0000000 0.6742000 0.6742000

# Chi Squared Stage IV Race -------------------------------------------------------

RACE_STAGE_IV_AIoAN <-matrix(c(0,1,16,962),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_IV_AIoAN)<-c("AMERICAN INDIAN OR ALASKA NATIVE", "All other Races")
colnames(RACE_STAGE_IV_AIoAN)<-c("Stage_IV","Not_Stage_IV")
RACE_STAGE_IV_AIoAN <- as.table(RACE_STAGE_IV_AIoAN)

chisq.test(RACE_STAGE_IV_AIoAN)

#p-value = 1

RACE_STAGE_IV_ASIAN <-matrix(c(0,47,16,916),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_IV_ASIAN)<-c("ASIAN", "All other Races")
colnames(RACE_STAGE_IV_ASIAN)<-c("Stage_IV","Not_Stage_IV")
RACE_STAGE_IV_ASIAN <- as.table(RACE_STAGE_IV_ASIAN)

chisq.test(RACE_STAGE_IV_ASIAN)

#p-value = 0.7519

RACE_STAGE_IV_BoAA <-matrix(c(5,172,11,791),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_IV_BoAA)<-c("BoAA", "All other Races")
colnames(RACE_STAGE_IV_BoAA)<-c("Stage_IV","Not_Stage_IV")
RACE_STAGE_IV_BoAA <- as.table(RACE_STAGE_IV_BoAA)

chisq.test(RACE_STAGE_IV_BoAA)

#p-value = 0.2925

RACE_STAGE_IV_W <-matrix(c(11,729,5,234),ncol=2,byrow=TRUE)
rownames(RACE_STAGE_IV_W)<-c("White", "All other races")
colnames(RACE_STAGE_IV_W)<-c("Stage_IV","Not_Stage_IV_Total")
RACE_STAGE_IV_W <- as.table(RACE_STAGE_IV_W)

chisq.test(RACE_STAGE_IV_W)

#p-value = 0.7274

pvalues <- c(1, 0.7519, 0.2925, 0.7274)

p.adjust(pvalues, method= "BH")

#0.6941333 1.0000000 0.6742000 0.6742000

# ER-/PR- and ER+/PR+ cleaning ----------------------------------------------------------

Breast_unknown_ER_PR_NegNeg_PosPos <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  dplyr::filter(!breast_carcinoma_estrogen_receptor_status %in% c("","Indeterminate")) %>%
  dplyr::filter(!breast_carcinoma_progesterone_receptor_status %in% c("","Indeterminate")) %>%
  filter((breast_carcinoma_progesterone_receptor_status == "Negative" & breast_carcinoma_estrogen_receptor_status == "Negative")|
         (breast_carcinoma_progesterone_receptor_status == "Positive" & breast_carcinoma_estrogen_receptor_status == "Positive"))

# KM ER+/PR+ vs ER-/PR- ---------------------------------------------------

surv_object <- Surv(time = Breast_unknown_ER_PR_NegNeg_PosPos$Years_to_death_and_last_followup, event = Breast_unknown_ER_PR_NegNeg_PosPos$vital_status_TRUE_FALSE)
surv_object
fit1 <- survfit(surv_object ~ breast_carcinoma_estrogen_receptor_status + breast_carcinoma_progesterone_receptor_status , data = Breast_unknown_ER_PR_NegNeg_PosPos) 
fit1

ER_PR_NegNeg_PosPos <- ggsurvplot(
  fit1,
  data = Breast_unknown_ER_PR_NegNeg_PosPos, 
  pval = TRUE, legend.title = "Receptor Status", 
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs =
    c("ER-/PR-","ER+/PR+"),
  ggtheme = theme_bw(),
  break.time.by = 2,
  tables.height = 0.4
) + xlab("Years")  

ER_PR_NegNeg_PosPos$plot <- ER_PR_NegNeg_PosPos$plot + labs(
  title = "KM of ER+/PR+ vs ER-/PR-")

ER_PR_NegNeg_PosPos <- ggpar(
  ER_PR_NegNeg_PosPos,
  font.title = c(16, "bold", "darkblue"),
  font.x = c(14, "plain", "black"),
  font.y = c(14, "plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab =  c(12, "plain", "black"),
  legend = "bottom"
)

ER_PR_NegNeg_PosPos

# ER+/PR- and ER+/PR+ cleaning ----------------------------------------------------------

Breast_unknown_ER_PR_PosNeg_PosPos <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  dplyr::filter(!breast_carcinoma_estrogen_receptor_status %in% c("","Indeterminate")) %>%
  dplyr::filter(!breast_carcinoma_progesterone_receptor_status %in% c("","Indeterminate")) %>%
  filter((breast_carcinoma_progesterone_receptor_status == "Negative" & breast_carcinoma_estrogen_receptor_status == "Positive")|
           (breast_carcinoma_progesterone_receptor_status == "Positive" & breast_carcinoma_estrogen_receptor_status == "Positive"))

# KM ER+/PR+ vs ER+/PR- ---------------------------------------------------

surv_object <- Surv(time = Breast_unknown_ER_PR_PosNeg_PosPos$Years_to_death_and_last_followup, event = Breast_unknown_ER_PR_PosNeg_PosPos$vital_status_TRUE_FALSE)
surv_object
fit1 <- survfit(surv_object ~ breast_carcinoma_estrogen_receptor_status + breast_carcinoma_progesterone_receptor_status , data = Breast_unknown_ER_PR_PosNeg_PosPos) 
fit1

ER_PR_PosNeg_PosPos <- ggsurvplot(
  fit1,
  data = Breast_unknown_ER_PR_PosNeg_PosPos, 
  pval = TRUE, legend.title = "Receptor Status", 
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs =
    c("ER+/PR-","ER+/PR+"),
  ggtheme = theme_bw(),
  break.time.by = 2,
  tables.height = 0.4
) + xlab("Years")  

ER_PR_PosNeg_PosPos$plot <- ER_PR_PosNeg_PosPos$plot + labs(
  title = "KM of ER+/PR+ vs ER+/PR-")

ER_PR_PosNeg_PosPos <- ggpar(
  ER_PR_PosNeg_PosPos,
  font.title = c(16, "bold", "darkblue"),
  font.x = c(14, "plain", "black"),
  font.y = c(14, "plain", "black"),
  font.xtickslab = c(12, "plain", "black"),
  font.ytickslab =  c(12, "plain", "black"),
  legend = "bottom"
)

ER_PR_PosNeg_PosPos

# Testing for Proportional Hazard  ----------------------------------------

Breast_unknown_ER <- Breast_unknown_ER %>%
  mutate(breast_carcinoma_estrogen_receptor_status_n = case_when(
    breast_carcinoma_estrogen_receptor_status == "Positive" ~ 2,
    breast_carcinoma_estrogen_receptor_status == "Negative" ~ 1))

cox_ER <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ breast_carcinoma_estrogen_receptor_status_n, data = Breast_unknown_ER)
cox_ER

summary(cox_ER)

#The hazard is higher , and thus the prognosis is worse for patients with a Negative ER status than patients of ER positive. 
#Hence at any given point ER- patients have an instantaneous risk of death that is 1.92 x that of ER + patients.

Breast_unknown_PR <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  mutate(breast_carcinoma_progesterone_receptor_status_n = case_when(
    breast_carcinoma_progesterone_receptor_status == "Positive" ~ 2,
    breast_carcinoma_progesterone_receptor_status == "Negative" ~ 1))

cox_PR <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ breast_carcinoma_progesterone_receptor_status_n, data = Breast_unknown_PR)
summary(cox_PR)

#The hazard is higher , and thus the prognosis is worse for patients with a Negative PR status than patients of PR positive. 
#Hence at any given point PR- patients have an instantaneous risk of death that is 1.92 x that of PR + patients.

cox_Age <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ Age, data = Breast_unknown_ER)
summary(cox_Age)

# an additional year of age increases the Yearly hazard of death by a factor of eb2 = 1.026 on average—that is, by 2.6 percent.

Breast_unknown_HER2 <- Breast_unknown_HER2 %>%
  mutate(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = case_when(
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive" ~ 2,
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative" ~ 1))

cox_HER2 <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ lab_proc_her2_neu_immunohistochemistry_receptor_status_n, data = Breast_unknown_HER2)
summary(cox_HER2)

#The hazard is lower , and thus the prognosis is worse for patients with a Positive HER2 status than patients HER2 negative. 
#Hence being HER2 negative reduces the hazard by a factor of 0.41, or 59%. Being HER2 negative is associated with good prognosis.

Breast_unknown_ER_HER2 <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_estrogen_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(breast_carcinoma_estrogen_receptor_status_n = case_when(
    breast_carcinoma_estrogen_receptor_status == "Positive" ~ 2,
    breast_carcinoma_estrogen_receptor_status == "Negative" ~ 1)) %>%
  filter(lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = case_when(
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive" ~ 2,
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative" ~ 1)) %>%
  filter(breast_carcinoma_progesterone_receptor_status %in% c("Positive", "Negative"))
 

cox_ER_HER2 <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ lab_proc_her2_neu_immunohistochemistry_receptor_status_n + breast_carcinoma_estrogen_receptor_status_n , data = Breast_unknown_ER_HER2)
cox_ER_HER2
summary(cox_ER_HER2)

#The hazard is lower , and thus the prognosis is worse for patients with a Positive HER2 status than patients HER2 negative. 
#Hence being HER2 negative reduces the hazard by a factor of 0.42, or 58%. Being HER2 negative is associated with good prognosis.
#The hazard is higher , and thus the prognosis is worse for patients with a Negative ER status than patients of ER positive. 
#Hence at any given point ER- patients have an instantaneous risk of death that is 1.83 x that of ER + patients.

Breast_unknown_PR_HER2 <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_progesterone_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(breast_carcinoma_progesterone_receptor_status_n = case_when(
    breast_carcinoma_progesterone_receptor_status == "Positive" ~ 2,
    breast_carcinoma_progesterone_receptor_status == "Negative" ~ 1)) %>%
  filter(lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = case_when(
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive" ~ 2,
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative" ~ 1)) %>%
  filter(breast_carcinoma_progesterone_receptor_status %in% c("Positive", "Negative"))


cox_PR_HER2 <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ lab_proc_her2_neu_immunohistochemistry_receptor_status_n + breast_carcinoma_progesterone_receptor_status_n , data = Breast_unknown_PR_HER2)
cox_PR_HER2
summary(cox_PR_HER2)

#The hazard is lower , and thus the prognosis is worse for patients with a Positive HER2 status than patients HER2 negative. Holding ER status constant, being HER2 negative reduces the hazard by a factor of 0.39, or 61%. Being HER2 negative is associated with good prognosis.
#The p-value for PR is  p = 0.491. The hazard ratio HR = exp(coef) = 0.306, with a 95% confidence interval of 0.6778 to 2.247.  Because the confidence interval for HR includes 1, these results indicate that PR Status makes a smaller contribution to the difference in the HR after adjusting for the HER2 status, and only trend toward significance.

Breast_unknown_ER_PR <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_estrogen_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(breast_carcinoma_estrogen_receptor_status_n = case_when(
    breast_carcinoma_estrogen_receptor_status == "Positive" ~ 2,
    breast_carcinoma_estrogen_receptor_status == "Negative" ~ 1)) %>%
  filter(breast_carcinoma_progesterone_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(breast_carcinoma_progesterone_receptor_status_n = case_when(
    breast_carcinoma_progesterone_receptor_status == "Positive" ~ 2,
    breast_carcinoma_progesterone_receptor_status == "Negative" ~ 1)) 


cox_ER_PR <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ breast_carcinoma_progesterone_receptor_status_n + breast_carcinoma_estrogen_receptor_status_n , data = Breast_unknown_ER_PR)
cox_ER_PR
summary(cox_ER_PR)

#The p-value for PR is  p = 0.66. The hazard ratio HR = exp(coef) = 0.372, with a 95% confidence interval of 0.5683 to 2.440.  Because the confidence interval for HR includes 1, these results indicate that PR Status makes a smaller contribution to the difference in the HR after adjusting for the ER status, and only trend toward significance.
#The p-value for ER is  p = 0.208. The hazard ratio HR = exp(coef) = 0.382, with a 95% confidence interval of 0.7651 to 3.421.  Because the confidence interval for HR includes 1, these results indicate that ER Status makes a smaller contribution to the difference in the HR after adjusting for the PR status, and only trend toward significance.

Breast_unknown_ER_HER2_PR <- Breast_Summary_Stages_Filtered_Blanks_Age_Race %>%
  filter(breast_carcinoma_estrogen_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(ER_Factor = case_when(
    breast_carcinoma_estrogen_receptor_status == "Positive" ~ 2,
    breast_carcinoma_estrogen_receptor_status == "Negative" ~ 1)) %>%
  filter(lab_proc_her2_neu_immunohistochemistry_receptor_status %in% c("Positive", "Negative")) %>%
  mutate(HER2_Factor = case_when(
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive" ~ 2,
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative" ~ 1)) %>%
  filter(breast_carcinoma_progesterone_receptor_status %in% c("Positive", "Negative")) %>%
mutate(PR_Factor = case_when(
  breast_carcinoma_progesterone_receptor_status == "Positive" ~ 2,
  breast_carcinoma_progesterone_receptor_status == "Negative" ~ 1))

str(Breast_unknown_ER_HER2_PR)

cox_ER_HER2_PR <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ HER2_Factor + ER_Factor + PR_Factor, data = Breast_unknown_ER_HER2_PR)
cox_ER_HER2_PR

summary(cox_ER_HER2_PR)

ggforest(cox_ER_HER2_PR)

#The hazard is lower , and thus the prognosis is worse for patients with a Positive HER2 status than patients HER2 negative. 
#Hence being HER2 negative reduces the hazard by a factor of 0.39, or 61%. Being HER2 negative is associated with good prognosis.
#The hazard is higher , and thus the prognosis is worse for patients with a Negative ER status than patients of ER positive. 
#Indicating a strong relationship between the ER status value and increased risk of death.
#the p-value for PR is now p=0.22. The hazard ratio HR = exp(coef) = 0.43, with a 95% confidence interval of 0.114 to 1.647. 
#Because the confidence interval for HR includes 1, these results indicate that PR Status makes a smaller contribution to the difference in the HR after adjusting for the HER2 status and ER status, and only trend toward significance

HER2_df <- with(Breast_unknown_ER_HER2_PR,
                data.frame(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = c(1,2),
                           breast_carcinoma_progesterone_receptor_status_n = c(1,1),
                           breast_carcinoma_estrogen_receptor_status_n = c(1,1)))
HER2_df

fit <- survfit(cox_ER_HER2_PR, newdata = HER2_df)
fit

ggsurvplot(fit, data = Breast_unknown_ER_HER2_PR, pval = TRUE, conf.int = TRUE, title = "Cox Regression of HER2+/ER+/PR+ vs HER2-/ER+/PR+",legend.title = "Major Receptor Status",legend.labs=c("HER2+/ER+/PR+", "HER2-/ER+/PR+"),
           ggtheme = theme_minimal())


PR_df <- with(Breast_unknown_ER_HER2_PR,
                data.frame(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = c(1,1),
                           breast_carcinoma_progesterone_receptor_status_n = c(1,2),
                           breast_carcinoma_estrogen_receptor_status_n = c(1,1)))
PR_df

fit <- survfit(cox_ER_HER2_PR, newdata = PR_df)
fit

ggsurvplot(fit, data = Breast_unknown_ER_HER2_PR, conf.int = TRUE, title = "Cox Regression of HER2+/ER+/PR+ vs HER2+/ER+/PR-",legend.title = "Major Receptor Status", legend.labs=c("HER2+/ER+/PR+", "HER2+/ER+/PR-"),
           ggtheme = theme_minimal())

ER_df <- with(Breast_unknown_ER_HER2_PR,
              data.frame(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = c(1,1),
                         breast_carcinoma_progesterone_receptor_status_n = c(1,1),
                         breast_carcinoma_estrogen_receptor_status_n = c(1,2)))
ER_df

fit <- survfit(cox_ER_HER2_PR, newdata = ER_df)
fit

ggsurvplot(fit, data = Breast_unknown_ER_HER2_PR, conf.int = TRUE, title = "Cox Regression of HER2+/ER+/PR+ vs HER2+/ER-/PR+",legend.title = "Major Receptor Status", legend.labs=c("HER2+/ER+/PR+", "HER2+/ER-/PR+"),
           ggtheme = theme_minimal())

# Cox Negatives and average Age -----------------------------------------------------------

cox_ER_HER2_PR_Age <- coxph(Surv(Years_to_death_and_last_followup, vital_status_TRUE_FALSE) ~ HER2_Factor + ER_Factor + PR_Factor + Age, data = Breast_unknown_ER_HER2_PR)
cox_ER_HER2_PR_Age

summary(cox_ER_HER2_PR_Age)

ggforest(cox_ER_HER2_PR_Age)

ggadjustedcurves(cox_ER_HER2_PR_Age, data = Breast_unknown_ER_HER2_PR, variable  ="breast_carcinoma_progesterone_receptor_status_n", xlab = "Time (Years)") 
ggadjustedcurves(cox_ER_HER2_PR_Age, data = Breast_unknown_ER_HER2_PR, variable  ="breast_carcinoma_estrogen_receptor_status_n",  xlab = "Time (Years)")
ggadjustedcurves(cox_ER_HER2_PR_Age, data = Breast_unknown_ER_HER2_PR, variable  ="lab_proc_her2_neu_immunohistochemistry_receptor_status_n",  xlab = "Time (Years)") 

#The hazard is lower , and thus the prognosis is worse for patients with a Positive HER2 status than patients HER2 negative. 
#Hence being HER2 negative reduces the hazard by a factor of 0.34, or 66%. Being HER2 negative is associated with good prognosis.
#holding the other covariates constant, an additional year of age increases the Yearly hazard of death by a factor of eb2 = 1.044 on average—that is, by 4.4 percent.
#the p-value for PR is now p=0.24. The hazard ratio HR = exp(coef) = 0.44, with a 95% confidence interval of 0.11 to 1.73. 
#Because the confidence interval for HR includes 1, these results indicate that PR Status makes a smaller contribution to the difference in the HR after adjusting for the HER2 status and Age, and only trend toward significance
#the p-value for ER is now p=0.11. The hazard ratio HR = exp(coef) = 3.17, with a 95% confidence interval of 0.78 to 12.9. 
#Because the confidence interval for HR includes 1, these results indicate that ER Status makes a smaller contribution to the difference in the HR after adjusting for the HER2 status and Age, and only trend toward significance

test.ph_ER_HER2_PR2 <- cox.zph(cox_ER_HER2_PR2)
test.ph_ER_HER2_PR2

ggcoxzph(test.ph_ER_HER2_PR2)

HER2_df2 <- with(Breast_unknown_ER_HER2_PR,
                data.frame(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = c(1,2),
                           breast_carcinoma_progesterone_receptor_status_n = c(1,1),
                           breast_carcinoma_estrogen_receptor_status_n = c(1,1),
                           Age = rep(mean(Age, na.rm = TRUE), 2)))
HER2_df2

fit <- survfit(cox_ER_HER2_PR, newdata = HER2_df2)
fit

ggsurvplot(fit, data = Breast_unknown_ER_HER2_PR, conf.int = TRUE, title = "Cox Regression of HER2+/ER+/PR+ vs HER2-/ER+/PR+, Age = 58.69 (mean)",legend.title = "Major Receptor Status",legend.labs=c("HER2+/ER+/PR+", "HER2-/ER+/PR+"),
           ggtheme = theme_minimal())


PR_df2 <- with(Breast_unknown_ER_HER2_PR,
              data.frame(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = c(1,1),
                         breast_carcinoma_progesterone_receptor_status_n = c(1,2),
                         breast_carcinoma_estrogen_receptor_status_n = c(1,1),
                         Age = rep(mean(Age, na.rm = TRUE), 2)))
PR_df2

fit <- survfit(cox_ER_HER2_PR, newdata = PR_df2)
fit

ggsurvplot(fit, data = Breast_unknown_ER_HER2_PR2, conf.int = TRUE, title = "Cox Regression of HER2+/ER+/PR+ vs HER2+/ER+/PR-, Age = 58.69 (mean)",legend.title = "Major Receptor Status", legend.labs=c("HER2+/ER+/PR+", "HER2+/ER+/PR-"),
           ggtheme = theme_minimal())

ER_df2 <- with(Breast_unknown_ER_HER2_PR,
              data.frame(lab_proc_her2_neu_immunohistochemistry_receptor_status_n = c(1,1),
                         breast_carcinoma_progesterone_receptor_status_n = c(1,1),
                         breast_carcinoma_estrogen_receptor_status_n = c(1,2),
                         Age = rep(mean(Age, na.rm = TRUE), 2)))
ER_df2

fit <- survfit(cox_ER_HER2_PR, newdata = ER_df2)
fit

ggsurvplot(fit, data = Breast_unknown_ER_HER2_PR, pval = TRUE, conf.int = TRUE, title = "Cox Regression of HER2+/ER+/PR+ vs HER2+/ER-/PR+, Age = 58.69 (mean)",legend.title = "Major Receptor Status", legend.labs=c("HER2+/ER+/PR+", "HER2+/ER-/PR+"),
           ggtheme = theme_minimal())



