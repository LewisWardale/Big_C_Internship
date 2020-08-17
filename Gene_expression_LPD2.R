Gene_expression_LPD1 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_1_eridanus.csv")
Gene_expression_LPD2 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_2_eridanus.csv")
Gene_expression_LPD3 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_3_eridanus.csv")
Gene_expression_LPD4 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_4_eridanus.csv")
Gene_expression_LPD5 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_5_eridanus.csv")
Gene_expression_LPD6 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_6_eridanus.csv")
Gene_expression_LPD7 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_7_eridanus.csv")


library(tidyverse)
library(survminer)
library(survival)
library(DataExplorer)
library(magrittr)
library(plyr)
library(reshape2)


# Filtering the data for  Log2fold abs 0-1, padj P < 0.05 and ENGS Genes  -------------------------

gene_overexpression_LPD2 <- Gene_expression_LPD2 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  top_n(n = 50, wt = log2FoldChange) %>%
  select(Gene,log2FoldChange,padj)

Gene_expression_LPD1 <- Gene_expression_LPD1 %>%
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

Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD1,gene_overexpression_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD3, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD4, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD5, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD6, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD7, Join_LPD2)

Join_LPD2 <- rename(Join_LPD2, c("log2FoldChange" = "log2FoldChange_LPD7",
                                 "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                 "log2FoldChange.y" = "log2FoldChange_LPD5",
                                 "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                 "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                 "log2FoldChange.x" = "log2FoldChange_LPD1",
                                 "log2FoldChange.y.y.y" = "log2FoldChange_LPD2"))

write.csv(Join_LPD2, 'LPD2_Top50_OverExpressed.csv')

summary(Join_LPD2)

MeanGamma <- c(0.04414136, 0.46744255, 0.23478199, 0.01685848, 0.06612514, 0.15546842, 0.01518204)
LPD <- c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7")
MeanLog2FC <- c(-0.97193, 1.355, -0.1856, -2.2773, 0.3433, -0.58198, -1.5970)

LPD2_Overexpression_gamma <- data.frame(LPD, MeanGamma, MeanLog2FC)

ggplot(LPD2_Overexpression_gamma, aes(x = MeanGamma, y = MeanLog2FC)) +
  geom_point() + 
  ggtitle("Scatter Plot of Max LPD 2 mean gamma values and the mean Log2 Fold Change of Top 50 overexpressed genes") +
  labs(x = "Mean Gamma", y = "Mean Log2 Fold Change") +
  geom_smooth(method = "lm") +
  geom_text(aes(label= LPD), hjust = -0.1)


# Underexpression ---------------------------------------------------------

gene_underexpression_LPD2 <- Gene_expression_LPD2 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  top_n(n = -50, wt = log2FoldChange) %>%
  select(Gene,log2FoldChange,padj)

Gene_expression_LPD1 <- Gene_expression_LPD1 %>%
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

Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD1,gene_underexpression_LPD2)
Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD3, Join_LPD2_under)
Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD4, Join_LPD2_under)
Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD5, Join_LPD2_under)
Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD6, Join_LPD2_under)
Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD7, Join_LPD2_under)

Join_LPD2_under <- rename(Join_LPD2_under, c("log2FoldChange" = "log2FoldChange_LPD7",
                                             "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                             "log2FoldChange.y" = "log2FoldChange_LPD5",
                                             "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                             "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                             "log2FoldChange.x" = "log2FoldChange_LPD1",
                                             "log2FoldChange.y.y.y" = "log2FoldChange_LPD2"))

write.csv(Join_LPD2_under, 'LPD2_Top50_UnderExpressed.csv')

summary(Join_LPD2_under)


MeanGamma <- c(0.04414136, 0.46744255, 0.23478199, 0.01685848, 0.06612514, 0.15546842, 0.01518204)
LPD <- c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7")
MeanLog2FC <- c(-4.038, -4.972, -3.7296, 3.390, -3.111, -3.409, -2.943)

LPD2_Underexpression_gamma <- data.frame(LPD, MeanGamma, MeanLog2FC)

ggplot(LPD2_Underexpression_gamma, aes(x = MeanGamma, y = MeanLog2FC)) +
  geom_point() + 
  ggtitle("Scatter Plot of Max LPD 2 mean gamma values and the mean Log2 Fold Change of Top 50 Underexpressed genes") +
  labs(x = "Mean Gamma", y = "Mean Log2 Fold Change") +
  geom_smooth(method = "lm") +
  geom_text(aes(label= LPD), hjust = -0.1)


# OverExpression vs Expression ------------------------------------------------

gene_overexpression_LPD2 <- Gene_expression_LPD2 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  top_n(n = 50, wt = log2FoldChange) %>%
  select(Gene,log2FoldChange,padj)

Gene_expression_LPD1 <- Gene_expression_LPD1 %>%
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

Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD1,gene_overexpression_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD3, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD4, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD5, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD6, Join_LPD2)
Join_LPD2 <- right_join(by = "Gene", Gene_expression_LPD7, Join_LPD2)

Join_LPD2 <- rename(Join_LPD2, c("log2FoldChange" = "log2FoldChange_LPD7",
                                 "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                 "log2FoldChange.y" = "log2FoldChange_LPD5",
                                 "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                 "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                 "log2FoldChange.x" = "log2FoldChange_LPD1",
                                 "log2FoldChange.y.y.y" = "log2FoldChange_LPD2"))
Join_LPD2 <- Join_LPD2 %>%
  select(Gene,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)

Join_LPD2

Join_LPD2_melt <- melt(Join_LPD2, id.vars = c("log2FoldChange_LPD2", "Gene"))

 ggplot(Join_LPD2_melt, aes(x = log2FoldChange_LPD2 , y = value, color = variable)) +
  geom_point() + 
  ggtitle("Scatter Plot of Max LPD 2 Top 50 overexpressed genes") +
  labs(x = "Log2 Fold Change LPD2", y = "Log2 Fold Change") +
  geom_smooth(method = "lm") + labs(color = "Log2 Fold Change of LPD")
  # geom_text(aes(label= Gene), hjust = -0.1)


# Underexpression vs expression -------------------------------------------

 gene_underexpression_LPD2 <- Gene_expression_LPD2 %>%
   dplyr::filter(abs(log2FoldChange) >= 1) %>%
   dplyr::filter(padj <= 0.05) %>%
   dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
   top_n(n = -50, wt = log2FoldChange) %>%
   select(Gene,log2FoldChange,padj)
 
 Gene_expression_LPD1 <- Gene_expression_LPD1 %>%
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
 
 Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD1,gene_underexpression_LPD2)
 Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD3, Join_LPD2_under)
 Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD4, Join_LPD2_under)
 Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD5, Join_LPD2_under)
 Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD6, Join_LPD2_under)
 Join_LPD2_under <- right_join(by = "Gene", Gene_expression_LPD7, Join_LPD2_under)
 
 Join_LPD2_under <- rename(Join_LPD2_under, c("log2FoldChange" = "log2FoldChange_LPD7",
                                              "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                              "log2FoldChange.y" = "log2FoldChange_LPD5",
                                              "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                              "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                              "log2FoldChange.x" = "log2FoldChange_LPD1",
                                              "log2FoldChange.y.y.y" = "log2FoldChange_LPD2"))
 
 Join_LPD2_under <- Join_LPD2_under %>%
   select(Gene,
          log2FoldChange_LPD7,
          log2FoldChange_LPD6,
          log2FoldChange_LPD5,
          log2FoldChange_LPD4,
          log2FoldChange_LPD3,
          log2FoldChange_LPD2,
          log2FoldChange_LPD1)
 
 Join_LPD2_under
 
 Join_LPD2_under_melt <- melt(Join_LPD2_under, id.vars = c("log2FoldChange_LPD2", "Gene"))
 
 ggplot(Join_LPD2_under_melt, aes(x = log2FoldChange_LPD2 , y = value, color = variable)) +
   geom_point() + 
   ggtitle("Scatter Plot of Max LPD 2 Top 50 underexpressed genes") +
   labs(x = "Log2 Fold Change LPD2", y = "Log2 Fold Change") +
   geom_smooth(method = "lm") + labs(color = "Log2 Fold Change of LPD")
