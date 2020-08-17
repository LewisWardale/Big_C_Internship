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

Join_LPD5 <- rename(Join_LPD1, c("meanlogFC" = "log2FoldChange_LPD7",
                                 "meanlogFC.x.x.x" = "log2FoldChange_LPD6",
                                 "meanlogFC.y" = "log2FoldChange_LPD5",
                                 "meanlogFC.x.x" = "log2FoldChange_LPD4",
                                 "meanlogFC.y.y" = "log2FoldChange_LPD3",
                                 "meanlogFC.x" = "log2FoldChange_LPD2",
                                 "meanlogFC.y.y.y" = "log2FoldChange_LPD1"))

Join_LPD5 <- Join_LPD5 %>%
  select(Gene_Symbol,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)

Join_LPD5_melt <- melt(Join_LPD5, id.vars = c("Gene_Symbol"))

Join_LPD5_melt_abs <- Join_LPD5_melt %>%
  mutate(abs_Log2FC = abs(value))


Join_LPD5_melt_abs_meangamma <- Join_LPD5_melt_abs %>%
  mutate(meangamma = case_when(
    variable == "log2FoldChange_LPD7" ~ 0.01415445,
    variable == "log2FoldChange_LPD6" ~ 0.02080138,
    variable == "log2FoldChange_LPD5" ~ 0.70248198,
    variable == "log2FoldChange_LPD4" ~ 0.03219549,
    variable == "log2FoldChange_LPD3" ~ 0.07965159,
    variable == "log2FoldChange_LPD2" ~ 0.13639636,
    variable == "log2FoldChange_LPD1" ~ 0.01431875,
  ) )

Ordered_Gene_LPD5 <- Join_LPD5_melt_abs_meangamma %>%
  arrange(Gene_Symbol)

NA_Gene_LPD5 <- na.omit(Ordered_Gene_LPD5)

GeneVector <- as.character(unique(NA_Gene_LPD5$Gene_Symbol))

for (Gene1 in GeneVector) {
  
  #Filters for each Gene
  Single_gene <- NA_Gene_LPD5 %>%
    dplyr::filter(Gene_Symbol == Gene1)
  
  #Run cor.test for each gene
  result <- cor.test(Single_gene$abs_Log2FC, Single_gene$meangamma, method = "spearman")
  
  if(result$p.value > 0.01) {
    next()
  } else { 
    print(result)
    print(paste0("Significant ", Gene1, " Pvalue"))
  }
  
}

Sig.Genes <- NA_Gene_LPD5 %>%
  dplyr::filter(Gene_Symbol %in% c("ABCC10", "AC068312.1", "ACOX3",
                                   "ADRM1", "AGBL4", "AL078471.5",
                                   "BCL3", "BEND5", "CDC25A",
                                   "COL13A1", "EHMT2", "EIF4G1",
                                   "GNA14", "HEPH", "HSF2",
                                   "IRF2BP1", "KIF23", "MAPKAPK3",
                                   "NDUFA4", "NIPBL", "NUP37",
                                   "POP4", "PPP3R1", "RBFOX2",
                                   "RELT", "RMND5A", "RP11-253M7.1",
                                   "RP11-474G23.3", "RP4-616B8.4", "RP5-855F16.1",
                                   "SLC25A3", "SLC4A4", "TPTE",
                                   "TRMT44", "XPNPEP3", "ZNF613"))

write.csv(Sig.Genes, "Sig_genes_LPD_5.csv")


GeneVector2 <- as.character(unique(Sig.Genes$Gene_Symbol))

for (Gene1 in GeneVector2) {
  
  #Filters for each Gene
  Single_gene1 <- Sig.Genes %>%
    dplyr::filter(Gene_Symbol == Gene1)
  
  superplot <- ggplot(Single_gene1, aes(x = meangamma, y = abs_Log2FC)) +
    geom_point(shape = 1) + ggtitle(paste0("Scatter Plot of LPD 1 for " , Gene1)) + geom_smooth(method = "lm")
  
  pdf(paste0("LPD_5/Scatterplot_LPD5_",Gene1,"P001.pdf"))
  print(superplot)
  dev.off()
}


Sig.Genes2 <- Sig.Genes %>%
  dplyr::filter(Gene_Symbol %in% c("BCL3", "RBFOX2"))

m <- lm(abs_Log2FC ~ meangamma, Sig.Genes2)
summary(m)$r.squared

ggplot(Sig.Genes2, aes(x = meangamma , y = abs_Log2FC, color = Gene_Symbol)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE) +
  geom_text(label = summary(m)$r.squared, x = 0.3, y = 0.3, show.legend = FALSE)

# Sig diff genes 2 ------------------------------------------------------------------

Sig_genes_LPD5 <- Meth_expression_LPD5 %>%
  dplyr::filter(Gene_Symbol %in% c("CTD-2231E14.2", "MIR4444-1","PSMB9",
                                   "TPM4", "TTK"))

Sig_diff_exp_LPD5 <- Sig_genes_LPD5 %>%
  dplyr::filter(abs(meanlogFC) >= 1) %>%
  dplyr::filter(p.val <= 0.05) %>%
  select(Gene_Symbol,meanlogFC,p.val)

write.csv(Sig_diff_exp_LPD5,"Sig_diff_genes_LPD_5.csv") 

