library(tidyverse)
library(Rfast)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpmisc)
# install.packages("ggpmisc")

Gamma_Gene_count_LPD6 <- read.csv("/Users/lewiswardale/Desktop/Gene_count/Gene_count/correlationResult_LPD_6.csv") %>%
  mutate(Gene = gsub("(.*)\\..*", "\\1", Gene)) 

GeneSymbol = bitr(Gamma_Gene_count_LP6$Gene, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") %>%
  left_join(Gamma_Gene_count_LPD6, by = c("ENSEMBL"  = "Gene"))




Gene_count <- read_csv("/Users/lewiswardale/Desktop/Gene_count/clean_matrix_RNA-seq_TCGA-BRCA_eridanus.csv")
Gamma_samples <- read.csv("/Users/lewiswardale/Desktop/TCGA-BRCA Assigned LPD/GammaValues_TCGA-BRCA_eridanus.csv")


# Filter for Top 30 -------------------------------------------------------

Gene_symbol_final <- GeneSymbol %>%
  dplyr::filter(pval < 0.05) %>%
  mutate(abs_rho = abs(rho.rho)) %>%
  slice_max(abs_rho, n = 30) 

Gamma_samples <- dplyr::select(Gamma_samples, LPD_6, Sample)
# Joining -----------------------------------------------------------------

# Replaces . by -
colnames(Gene_count) <- gsub("\\.", "-", colnames(Gene_count))

# Take only the first 16 characters of the sample code
colnames(Gene_count) <- substr(colnames(Gene_count), 1, 16)

Gene_count <- Gene_count %>%
  mutate(Gene = gsub("(.*)\\..*", "\\1", X1))

#NEED TO SORT THIS CODE OUT 
Joined  <- left_join(Gene_symbol_final, Gene_count, by = c("ENSEMBL" = "Gene")) %>%
  gather(key = "Sample", value = "Count", 5:ncol(.)) %>%
  left_join(Gamma_samples) %>%
  drop_na() %>%
  mutate(Gene_Count = as.numeric(Count)) %>%
  filter(Count > 0)



superplot <- ggplot(Joined, aes(x = LPD_6, y = Gene_Count)) +
  geom_point(shape = 1) + 
  geom_smooth(method = "lm", colour = "orange", se = FALSE) + 
  ylim(0,40000) +
  facet_wrap( ~ SYMBOL) 


superplot + stat_poly_eq(formula = y ~ x,aes(label = paste(..rr.label.., sep = "~~~")), 
                         parse=TRUE,label.x.npc = "right")
