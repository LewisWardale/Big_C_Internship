Meth_expression_LPD1 <- read.csv("/Users/lewiswardale/Desktop/DMGenes/DM_allGenes_TCGA-BRCA_LPD_1_eridanus.csv")
Meth_expression_LPD2 <- read.csv("/Users/lewiswardale/Desktop/DMGenes/DM_allGenes_TCGA-BRCA_LPD_2_eridanus.csv")
Meth_expression_LPD3 <- read.csv("/Users/lewiswardale/Desktop/DMGenes/DM_allGenes_TCGA-BRCA_LPD_3_eridanus.csv")
Meth_expression_LPD4 <- read.csv("/Users/lewiswardale/Desktop/DMGenes/DM_allGenes_TCGA-BRCA_LPD_4_eridanus.csv")
Meth_expression_LPD5 <- read.csv("/Users/lewiswardale/Desktop/DMGenes/DM_allGenes_TCGA-BRCA_LPD_5_eridanus.csv")
Meth_expression_LPD6 <- read.csv("/Users/lewiswardale/Desktop/DMGenes/DM_allGenes_TCGA-BRCA_LPD_6_eridanus.csv")
Meth_expression_LPD7 <- read.csv("/Users/lewiswardale/Desktop/DMGenes/DM_allGenes_TCGA-BRCA_LPD_7_eridanus.csv")


# install.packages("doBy")
library(doBy)
library(tidyverse)
library(survminer)
library(survival)
library(DataExplorer)
library(magrittr)
library(plyr)
library(reshape2)

Meth_expression_LPD1 <- Meth_expression_LPD1 %>%
  select(Gene_Symbol,meanlogFC, p.val)
Meth_expression_LPD2 <- Meth_expression_LPD2 %>%
  select(Gene_Symbol,meanlogFC,p.val)
Meth_expression_LPD3 <- Meth_expression_LPD3 %>%
  select(Gene_Symbol,meanlogFC,p.val)
Meth_expression_LPD4 <- Meth_expression_LPD4 %>%
  select(Gene_Symbol,meanlogFC,p.val)
Meth_expression_LPD5 <- Meth_expression_LPD5 %>%
  select(Gene_Symbol,meanlogFC,p.val)
Meth_expression_LPD6 <- Meth_expression_LPD6 %>%
  select(Gene_Symbol,meanlogFC,p.val)
Meth_expression_LPD7 <- Meth_expression_LPD7 %>%
  select(Gene_Symbol,meanlogFC,p.val)

Meth_expression_LPD1 <- Meth_expression_LPD1[order(Meth_expression_LPD1$Gene_Symbol, -abs(Meth_expression_LPD1$meanlogFC) ), ]
Meth_expression_LPD1 <- Meth_expression_LPD1[ !duplicated(Meth_expression_LPD1$Gene_Symbol), ]

Meth_expression_LPD2 <- Meth_expression_LPD2[order(Meth_expression_LPD2$Gene_Symbol, -abs(Meth_expression_LPD2$meanlogFC) ), ]
Meth_expression_LPD2 <- Meth_expression_LPD2[ !duplicated(Meth_expression_LPD2$Gene_Symbol), ]

Meth_expression_LPD3 <- Meth_expression_LPD3[order(Meth_expression_LPD3$Gene_Symbol, -abs(Meth_expression_LPD3$meanlogFC) ), ]
Meth_expression_LPD3 <- Meth_expression_LPD3[ !duplicated(Meth_expression_LPD3$Gene_Symbol), ]

Meth_expression_LPD4 <- Meth_expression_LPD4[order(Meth_expression_LPD4$Gene_Symbol, -abs(Meth_expression_LPD4$meanlogFC) ), ]
Meth_expression_LPD4 <- Meth_expression_LPD4[ !duplicated(Meth_expression_LPD4$Gene_Symbol), ]

Meth_expression_LPD5 <- Meth_expression_LPD5[order(Meth_expression_LPD5$Gene_Symbol, -abs(Meth_expression_LPD5$meanlogFC) ), ]
Meth_expression_LPD5 <- Meth_expression_LPD5[ !duplicated(Meth_expression_LPD5$Gene_Symbol), ]

Meth_expression_LPD6 <- Meth_expression_LPD6[order(Meth_expression_LPD6$Gene_Symbol, -abs(Meth_expression_LPD6$meanlogFC) ), ]
Meth_expression_LPD6 <- Meth_expression_LPD6[ !duplicated(Meth_expression_LPD6$Gene_Symbol), ]

Meth_expression_LPD7 <- Meth_expression_LPD7[order(Meth_expression_LPD7$Gene_Symbol, -abs(Meth_expression_LPD7$meanlogFC) ), ]
Meth_expression_LPD7 <- Meth_expression_LPD7[ !duplicated(Meth_expression_LPD7$Gene_Symbol), ]


Join_LPD1 <- right_join(by = "Gene_Symbol", Meth_expression_LPD2, Meth_expression_LPD1)
Join_LPD1 <- right_join(by = "Gene_Symbol", Meth_expression_LPD3, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene_Symbol", Meth_expression_LPD4, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene_Symbol", Meth_expression_LPD5, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene_Symbol", Meth_expression_LPD6, Join_LPD1)
Join_LPD1 <- right_join(by = "Gene_Symbol", Meth_expression_LPD7, Join_LPD1)

Join_LPD7 <- rename(Join_LPD1, c("meanlogFC" = "log2FoldChange_LPD7",
                                 "meanlogFC.x.x.x" = "log2FoldChange_LPD6",
                                 "meanlogFC.y" = "log2FoldChange_LPD5",
                                 "meanlogFC.x.x" = "log2FoldChange_LPD4",
                                 "meanlogFC.y.y" = "log2FoldChange_LPD3",
                                 "meanlogFC.x" = "log2FoldChange_LPD2",
                                 "meanlogFC.y.y.y" = "log2FoldChange_LPD1"))

Join_LPD7 <- Join_LPD7 %>%
  select(Gene_Symbol,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)

Join_LPD7_melt <- melt(Join_LPD7, id.vars = c("Gene_Symbol"))

Join_LPD7_melt_abs <- Join_LPD7_melt %>%
  mutate(abs_Log2FC = abs(value))


Join_LPD7_melt_abs_meangamma <- Join_LPD7_melt_abs %>%
  mutate(meangamma = case_when(
    variable == "log2FoldChange_LPD7" ~ 0.39131970,
    variable == "log2FoldChange_LPD6" ~ 0.13265366,
    variable == "log2FoldChange_LPD5" ~ 0.04095057,
    variable == "log2FoldChange_LPD4" ~ 0.08301207,
    variable == "log2FoldChange_LPD3" ~ 0.04825003,
    variable == "log2FoldChange_LPD2" ~ 0.07256134,
    variable == "log2FoldChange_LPD1" ~ 0.23125263,
  ) )


Ordered_Gene_LPD7 <- Join_LPD7_melt_abs_meangamma %>%
  arrange(Gene_Symbol)

NA_Gene_LPD7 <- na.omit(Ordered_Gene_LPD7)

GeneVector <- as.character(unique(NA_Gene_LPD7$Gene_Symbol))

for (Gene1 in GeneVector) {
  
  #Filters for each Gene
  Single_gene <- NA_Gene_LPD7 %>%
    dplyr::filter(Gene_Symbol == Gene1)
  
  #Run cor.test for each gene
  result <- cor.test(Single_gene$abs_Log2FC, Single_gene$meangamma, method = "spearman")
  
  if(result$p.value > 0.005) {
    next()
  } else { 
    print(result)
    print(paste0("Significant ", Gene1, " Pvalue"))
  }
  
}

Sig.Genes <- NA_Gene_LPD7 %>%
  dplyr::filter(Gene_Symbol %in% c("AC068312.1", "ABCC10", "ACOX3",
                                   "BCL3", "CD81", "CD81-AS1",
                                   "COPB1", "FBXL16", "IRF2BP1",
                                   "MGAT4B", "PIDD1", "POP4", 
                                   "RBFOX2", "SLC4A4", "THOP1",
                                   "TRMT44", "ZNF613" ))

write.csv(Sig.Genes, "Sig_genes_LPD_7_005.csv")


GeneVector2 <- as.character(unique(Sig.Genes$Gene_Symbol))

for (Gene1 in GeneVector2) {
  
  #Filters for each Gene
  Single_gene1 <- Sig.Genes %>%
    dplyr::filter(Gene_Symbol == Gene1)
  
  superplot <- ggplot(Single_gene1, aes(x = meangamma, y = abs_Log2FC)) +
    geom_point(shape = 1) + ggtitle(paste0("Scatter Plot of LPD 1 for " , Gene1)) + geom_smooth(method = "lm")
  
  pdf(paste0("LPD_7/Scatterplot_LPD7_",Gene1,"005.pdf"))
  print(superplot)
  dev.off()
}

m <- lm(abs_Log2FC ~ meangamma, Sig.Genes2)
summary(m)$r.squared

ggplot(Sig.Genes2, aes(x = meangamma , y = abs_Log2FC, color = Gene_Symbol)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE) +
  geom_text(label = summary(m)$r.squared, x = 0.3, y = 0.08, show.legend = FALSE)

Sig.Genes2 <- Sig.Genes %>%
  dplyr::filter(Gene_Symbol %in% c("BCL3", "MGAT4B"))

ggplot(Sig.Genes2, aes(x = meangamma , y = abs_Log2FC, color = Gene_Symbol)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE)
ggtitle("Scatter Plot of Max LPD 1 Top 50 overexpressed genes") +
  labs(x = "Log2 Fold Change LPD1", y = "Log2 Fold Change") +
  geom_smooth(method = "lm") + labs(color = "Log1 Fold Change of LPD")

# Sig diff genes 2 ------------------------------------------------------------------

Sig_genes_LPD7 <- Meth_expression_LPD7 %>%
  dplyr::filter(Gene_Symbol %in% c("AC068312.1", "ABCC10", "ACOX3",
                                   "BCL3", "CD81", "CD81-AS1",
                                   "COPB1", "FBXL16", "IRF2BP1",
                                   "MGAT4B", "PIDD1", "POP4", 
                                   "RBFOX2", "SLC4A4", "THOP1",
                                   "TRMT44", "ZNF613" ))

Sig_diff_exp_LPD7 <- Sig_genes_LPD7 %>%
  dplyr::filter(abs(meanlogFC) >= 1) %>%
  dplyr::filter(p.val <= 0.05) %>%
  select(Gene_Symbol,meanlogFC,p.val)

write.csv(Sig_diff_exp_LPD7,"Sig_diff_genes_LPD_7.csv") 

