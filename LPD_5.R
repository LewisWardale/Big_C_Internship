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
library(reshape2)

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

Join_LPD5 <- rename(Join_LPD1, c("log2FoldChange" = "log2FoldChange_LPD7",
                                 "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                 "log2FoldChange.y" = "log2FoldChange_LPD5",
                                 "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                 "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                 "log2FoldChange.x" = "log2FoldChange_LPD2",
                                 "log2FoldChange.y.y.y" = "log2FoldChange_LPD1"))

Join_LPD5 <- Join_LPD5 %>%
  select(Gene,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)

Join_LPD5_melt <- melt(Join_LPD5, id.vars = c("Gene"))

Join_LPD5_melt_abs <- Join_LPD5_melt %>%
  mutate(abs_Log2FC = abs(value))

Join_LP5_melt_abs_meangamma <- Join_LPD5_melt_abs %>%
  mutate(meangamma = case_when(
    variable == "log2FoldChange_LPD7" ~ 0.01415445,
    variable == "log2FoldChange_LPD6" ~ 0.02080138,
    variable == "log2FoldChange_LPD5" ~ 0.70248198,
    variable == "log2FoldChange_LPD4" ~ 0.03219549,
    variable == "log2FoldChange_LPD3" ~ 0.07965159,
    variable == "log2FoldChange_LPD2" ~ 0.13639636,
    variable == "log2FoldChange_LPD1" ~ 0.01431875,
  ) )

Ordered_Gene_LPD5 <- Join_LP5_melt_abs_meangamma %>%
  arrange(Gene)

NA_Gene_LPD5 <- na.omit(Ordered_Gene_LPD5)

GeneVector <- as.character(unique(NA_Gene_LPD5$Gene))

for (Gene1 in GeneVector) {
  
  #Filters for each Gene
  Single_gene <- NA_Gene_LPD5 %>%
    dplyr::filter(Gene == Gene1)
  
  #Run cor.test for each gene
  result <- cor.test(Single_gene$abs_Log2FC, Single_gene$meangamma, method = "spearman")
  
  if(result$p.value > 0.001) {
    next()
  } else { 
    print(result)
    print(paste0("Significant ", Gene1, " Pvalue"))
  }
  
}

Sig.Genes <- NA_Gene_LPD5 %>%
  dplyr::filter(Gene %in% c("API5P2",
                            "BRCC3",
                            "CA6",
                            "COL2A1",
                            "DM1-AS",
                            "DPP3P2",
                            "DSG1",
                            "GABRP",
                            "IGHV1-2",
                            "IGHV3-13",
                            "IGHV3OR16-11",
                            "IGHV7-27",
                            "IGHV8-51-1",
                            "IGKV1-27",
                            "IGKV1OR10-1",
                            "IGKV1OR2-11",
                            "IGKV1OR22-1",
                            "IGKV2-24",
                            "IGKV2-28",
                            "IGKV2D-26",
                            "IGKV2D-29",
                            "IGKV3D-15",
                            "IGKV6D-21",
                            "IGLV1-36",
                            "IGLV3-27",
                            "IGLV3-29",
                            "IGLV5-48",
                            "IGLV8-61",
                            "KLK5",
                            "LINC02188",
                            "MIA",
                            "MIR1273C",
                            "OSTCP1",
                            "ROPN1",
                            "RPL5P35",
                            "SLC6A14"
                            ))

write.csv(Sig.Genes, "Sig_genes_LPD_5.csv")

GeneVector2 <- as.character(unique(Sig.Genes$Gene))

for (Gene1 in GeneVector2) {
  
  #Filters for each Gene
  Single_gene1 <- Sig.Genes %>%
    dplyr::filter(Gene == Gene1)
  
  superplot <- ggplot(Single_gene1, aes(x = meangamma, y = abs_Log2FC)) +
    geom_point(shape = 1) + ggtitle(paste0("Scatter Plot of LPD 4 for " , Gene1)) + geom_smooth(method = "lm")
  
  pdf(paste0("LPD_5/Scatterplot_LPD5_",Gene1,".pdf"))
  print(superplot)
  dev.off()
}  

Sig.Genes2 <- Sig.Genes %>%
  dplyr::filter(Gene %in% c("BRCC3",
                            "COL2A1",
                            "DM1-AS",
                            "DPP32P",
                            "OSTCP1"))

m <- lm(abs_Log2FC ~ meangamma, Sig.Genes2)
summary(m)$r.squared

ggplot(Sig.Genes2, aes(x = meangamma , y = abs_Log2FC, color = Gene)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE) + 
  geom_text(label = summary(m)$r.squared, x = 0.3, y = 0.2, show.legend = FALSE)



ggplot(Sig.Genes2, aes(x = meangamma , y = abs_Log2FC, color = Gene)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE)

# Sig diff genes 2 ------------------------------------------------------------------


Sig_genes_LPD5 <- Gene_expression_LPD5 %>%
  dplyr::filter(Gene %in% c("API5P2",
                            "BRCC3",
                            "CA6",
                            "COL2A1",
                            "DM1-AS",
                            "DPP3P2",
                            "DSG1",
                            "GABRP",
                            "IGHV1-2",
                            "IGHV3-13",
                            "IGHV3OR16-11",
                            "IGHV7-27",
                            "IGHV8-51-1",
                            "IGKV1-27",
                            "IGKV1OR10-1",
                            "IGKV1OR2-11",
                            "IGKV1OR22-1",
                            "IGKV2-24",
                            "IGKV2-28",
                            "IGKV2D-26",
                            "IGKV2D-29",
                            "IGKV3D-15",
                            "IGKV6D-21",
                            "IGLV1-36",
                            "IGLV3-27",
                            "IGLV3-29",
                            "IGLV5-48",
                            "IGLV8-61",
                            "KLK5",
                            "LINC02188",
                            "MIA",
                            "MIR1273C",
                            "OSTCP1",
                            "ROPN1",
                            "RPL5P35",
                            "SLC6A14"
  ))


Sig_diff_exp_LPD5 <- Sig_genes_LPD5 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  select(Gene,log2FoldChange,padj)

write.csv(Sig_diff_exp_LPD5,"Sig_diff_genes_LPD_5.csv")  



