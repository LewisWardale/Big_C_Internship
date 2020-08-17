LPD_w_Gamma <- read.csv("/Users/lewiswardale/Desktop/TCGA-BRCA Assigned LPD/GammaValues_TCGA-BRCA_eridanus.csv")

library(tidyverse)
library(Rfast)
library(survminer)
library(survival)
library(lubridate)
library(dplyr)
library(MASS)
library(tidyr)
library(DataExplorer)
library(magrittr)
library(ggplot2)
library(reshape2)

# Most Likely LPD column --------------------------------------------------

LPD_w_Gamma_Max <- LPD_w_Gamma %>%
  mutate(LPD_Max_Gamma =
    ifelse(Max_LPD == "LPD_1", LPD_1,
    ifelse(Max_LPD == "LPD_2", LPD_2,
    ifelse(Max_LPD == "LPD_3", LPD_3,
    ifelse(Max_LPD == "LPD_4", LPD_4,
    ifelse(Max_LPD == "LPD_5", LPD_5,
    ifelse(Max_LPD == "LPD_6", LPD_6, LPD_7
  )))))))


#Box Plot of most likely LPD

Max_Gamma_v_LPD <- ggplot(LPD_w_Gamma_Max, aes(x = Max_LPD, y = LPD_Max_Gamma, fill = Max_LPD)) + 
  geom_boxplot(outlier.colour = "red")

Max_Gamma_v_LPD + ggtitle("Box Plot of Gamma and Most Likely LPD") + labs(y = "Gamma", x = "Max LPD") + labs(fill = "Max LPD")


# Box Plot for each Max LPD and their correspoding gamma ------------------

signaturesVector <- as.character(unique(LPD_w_Gamma$Max_LPD)) %>%
  as.character(unique(LPD_w_Gamma$Index))

allGammaList <- list()

for (LPD_n in signaturesVector) {
  
  #Filters for each max LPD
  allGamma <- LPD_w_Gamma %>%
    dplyr::filter(Max_LPD == LPD_n)
  
  if(nrow(allGamma) == 0){
    next()
  }
  
  #Combines all LPD gamma values into one column and runs a boxplot for these values
  LPD_w_Gamma_melt <- melt(allGamma, id.vars = "Index", measure.vars = c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7"))
  Boxplot <- ggplot(LPD_w_Gamma_melt, aes(x = variable, y = value, fill = variable)) + geom_boxplot(outlier.colour = "red") + ggtitle(paste0("Box Plot of each LPD Gamma for ", LPD_n)) + labs(y = "Gamma", x = "LPD") + labs(fill = "LPD")
 
  pdf(paste0("mySuperScatterplots/Boxplot_", LPD_n, " all gamma values.pdf"))
  print(Boxplot)
  dev.off()
  
  
  
  # Saves the mean gamma values
  allGammaList[[LPD_n]] <- allGamma
  
}

# Function to run Max LPD vs Second Highest Gamma -------------------------

# You want to make everything below this text into one function that returns the gammavalues in a nicely formatted 
# datafra,e like the one i just did


# Then, you can use this functions in your original script you showed 5 mins ago to have the gamma values, and you can then
# create a loop that itinerates to each row of the dataframe so in each itineration uses the proper gammavalues
#So itineration 1 (LPD-1) will use the values from the first row, itineration 2 the values from the second row and so on
# You will do great

Mean_gamma_table <- function(pathway){
  
  LPD_w_Gamma <- read.csv(pathway)
  
  signaturesVector <- as.character(unique(LPD_w_Gamma$Max_LPD))
  
  

  
}
  
  # Users/lewiswardale/Desktop/TCGA-BRCA Assigned LPD/GammaValues_TCGA-BRCA_eridanus.csv


meanGammaList <- list()

for (LPD_sign in signaturesVector) {
  
 # Creates a vector with the mean values for each LPD signature
  meanGamma <- LPD_w_Gamma %>%
    dplyr::filter(Max_LPD == LPD_sign) %>%
    summarise_at(c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7"), mean, na.rm = TRUE) %>%
    t() %>%
    as.data.frame() %>%
    pull()
  
  # Extracts the name of the second hgihest LPD signature for the sctatterplot
  names(meanGamma) = c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7")
  secondLPD <- names(sort(meanGamma, decreasing = TRUE)[2])
  
  # Creates an scatterplot with the corresponding LPD signatuer aand the second highest one
  superplot <- ggplot(LPD_w_Gamma, aes_string(x = LPD_sign, y = secondLPD)) +
    geom_point(shape = 1) + ggtitle(paste0("Scatter Plot of Max LPD = ", LPD_sign, " gamma values and the gamma values of ", secondLPD)) 
     
  pdf(paste0("mySuperScatterplots/Scatterplot_", LPD_sign, "_vs_", secondLPD, ".pdf"))
  print(superplot)
  dev.off()
  
  
  # Saves the mean gamma values
  meanGammaList[[LPD_sign]] <- meanGamma

 
}

meanGamma


