library(tidyverse)
library(survminer)
library(survival)
library(DataExplorer)
library(magrittr)
library(plyr)
library(grep)
# install.packages("janitor")
library(janitor)

Gene_count <- read_csv("/Users/lewiswardale/Desktop/Gene_count/clean_matrix_RNA-seq_TCGA-BRCA_eridanus.csv")
Gamma_samples <- read.csv("/Users/lewiswardale/Desktop/TCGA-BRCA Assigned LPD/GammaValues_TCGA-BRCA_eridanus.csv")


# Changing row names and melting ------------------------------------------

# Replaces . by -
colnames(Gene_count) <- gsub("\\.", "-", colnames(Gene_count))

# Take only the first 16 characters of the sample code
colnames(Gene_count) <- substr(colnames(Gene_count), 1, 16)

# merhmehmeh
for(process in 6:7){
  i = 0
  
  correlationResult <- tibble()
  
  for(selectedGene in unique(Gene_count$X1)){
    i = i + 1
    
    selectedGeneCount <- Gene_count %>%
      filter(X1 == selectedGene) %>%
      gather(key = "Sample", "Counts", -X1) %>%
      filter(Counts > 0)
      
      
    if(nrow(selectedGeneCount) < 2){
      next()
    }
    
    mergedGeneCount <- left_join(selectedGeneCount, Gamma_samples)
    
  # Calculates correlation
    result <- cor.test(mergedGeneCount$Counts, pull(mergedGeneCount[,process + 3]), method = "spearman")
    
    tempVector <- c(Gene = selectedGene, rho = result$estimate, pval = result$p.value)
    correlationResult <- bind_rows(correlationResult, tempVector)
  
  }
  
  write_csv(correlationResult, paste0("correlationResult_LPD_", process, ".csv"))
}




