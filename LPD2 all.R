Gene_expression_LPD1 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_1_eridanus.csv")
Gene_expression_LPD2 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_2_eridanus.csv")
Gene_expression_LPD3 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_3_eridanus.csv")
Gene_expression_LPD4 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_4_eridanus.csv")
Gene_expression_LPD5 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_5_eridanus.csv")
Gene_expression_LPD6 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_6_eridanus.csv")
Gene_expression_LPD7 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_7_eridanus.csv")

# install.packages("doBy")
library(doBy)
library(tidyverse)
library(survminer)
library(survival)
library(DataExplorer)
library(magrittr)
library(plyr)

Gene_expression_LPD1 <- Gene_expression_LPD1 %>%
  select(Gene,log2FoldChange, padj)
Gene_expression_LPD2 <- Gene_expression_LPD2 %>%
  select(Gene,log2FoldChange,padj)
Gene_expression_LPD3 <- Gene_expression_LPD3 %>%
  select(Gene,log2FoldChange,padj)
Gene_expression_LPD4 <- Gene_expression_LPD4 %>%
  select(Gene,log2FoldChange,padj)
Gene_expression_LPD5 <- Gene_expression_LPD5 %>%
  select(Gene,log2FoldChange,padj)
Gene_expression_LPD6 <- Gene_expression_LPD6 %>%
  select(Gene,log2FoldChange,padj)
Gene_expression_LPD7 <- Gene_expression_LPD7 %>%
  select(Gene,log2FoldChange,padj)

Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD2, Gene_expression_LPD1)
Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD3, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD4, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD5, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD6, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD7, Join_LPD1)

Join_LPD2 <- rename(Join_LPD1, c("log2FoldChange" = "log2FoldChange_LPD7",
                                 "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                 "log2FoldChange.y" = "log2FoldChange_LPD5",
                                 "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                 "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                 "log2FoldChange.x" = "log2FoldChange_LPD2",
                                 "log2FoldChange.y.y.y" = "log2FoldChange_LPD1"))

Join_LPD2 <- Join_LPD2 %>%
  select(Gene,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)

Join_LPD2_melt <- melt(Join_LPD2, id.vars = c("Gene"))

Join_LPD2_melt_abs <- Join_LPD2_melt %>%
  mutate(abs_Log2FC = abs(value))

Join_LP2_melt_abs_meangamma <- Join_LPD2_melt_abs %>%
  mutate(meangamma = case_when(
    variable == "log2FoldChange_LPD7" ~ 0.0151820,
    variable == "log2FoldChange_LPD6" ~ 0.15546842,
    variable == "log2FoldChange_LPD5" ~ 0.06612514,
    variable == "log2FoldChange_LPD4" ~ 0.01685848,
    variable == "log2FoldChange_LPD3" ~ 0.23478199,
    variable == "log2FoldChange_LPD2" ~ 0.46744255,
    variable == "log2FoldChange_LPD1" ~ 0.04414136,
  ) )

Ordered_Gene_LPD2 <- Join_LP2_melt_abs_meangamma %>%
  arrange(Gene)

NA_Gene_LPD2 <- na.omit(Ordered_Gene_LPD2)

GeneVector <- as.character(unique(NA_Gene_LPD2$Gene))

for (Gene1 in GeneVector) {
  
  #Filters for each Gene
  Single_gene <- NA_Gene_LPD2 %>%
    dplyr::filter(Gene == Gene1)
  
  #Run cor.test for each gene
  result <- cor.test(Single_gene$abs_Log2FC, Single_gene$meangamma, method = "spearman")
  
  if(result$p.value > 0.05) {
    next()
  } else { 
    print(result)
    print(paste0("Significant ", Gene1, " Pvalue"))
  }
  
}

Sig.Genes <- NA_Gene_LPD2 %>%
  dplyr::filter(Gene %in% c("ARHGAP40",
                            "CD200R1L-AS1",
                            "IGFN1",
                            "TAC1"))

write.csv(Sig.Genes, "Sig_genes_LPD_2.csv")


GeneVector2 <- as.character(unique(Sig.Genes$Gene))

for (Gene1 in GeneVector2) {
  
  #Filters for each Gene
  Single_gene1 <- Sig.Genes %>%
    dplyr::filter(Gene == Gene1)
  
  superplot <- ggplot(Single_gene1, aes(x = meangamma, y = abs_Log2FC)) +
    geom_point(shape = 1) + ggtitle(paste0("Scatter Plot of LPD 1 for " , Gene1)) + geom_smooth(method = "lm")
  
  pdf(paste0("LPD_2/Scatterplot_LPD1_",Gene1,".pdf"))
  print(superplot)
  dev.off()
}

ggplot(Sig.Genes, aes(x = meangamma , y = abs_Log2FC, color = Gene)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE)

# Sig diff genes 2 ------------------------------------------------------------------


Sig_genes_LPD2 <- Gene_expression_LPD2 %>%
  dplyr::filter(Gene %in% c("ARHGAP40",
                            "CD200R1L-AS1",
                            "IGFN1",
                            "TAC1"))

Sig_diff_exp_LPD2 <- Sig_genes_LPD2 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  select(Gene,log2FoldChange,padj)

write.csv(Sig_diff_exp_LPD2,"Sig_diff_genes_LPD_2.csv")  

