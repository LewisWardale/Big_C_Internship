Gene_expression_LPD1 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_1_eridanus.csv")
Gene_expression_LPD2 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_2_eridanus.csv")
Gene_expression_LPD3 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_3_eridanus.csv")
Gene_expression_LPD4 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_4_eridanus.csv")
Gene_expression_LPD5 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_5_eridanus.csv")
Gene_expression_LPD6 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_6_eridanus.csv")
Gene_expression_LPD7 <- read.csv("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_7_eridanus.csv")

# install.packages("doBy")
library(doBy)

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

Join_LPD1 <- rename(Join_LPD1, c("log2FoldChange" = "log2FoldChange_LPD7",
                                             "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                             "log2FoldChange.y" = "log2FoldChange_LPD5",
                                             "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                             "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                             "log2FoldChange.x" = "log2FoldChange_LPD2",
                                             "log2FoldChange.y.y.y" = "log2FoldChange_LPD1"))

Join_LPD1_under <- Join_LPD1 %>%
  select(Gene,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)

Join_LPD1_melt <- melt(Join_LPD1_under, id.vars = c("Gene"))

Join_LPD1_melt_abs <- Join_LPD1_melt %>%
  mutate(abs_Log2FC = abs(value))

Join_LPD1_melt_abs_meangamma <- Join_LPD1_melt_abs %>%
  mutate(meangamma = case_when(
    variable == "log2FoldChange_LPD7" ~ 0.178,
    variable == "log2FoldChange_LPD6" ~ 0.226,
    variable == "log2FoldChange_LPD5" ~ 0.0152,
    variable == "log2FoldChange_LPD4" ~ 0.0388,
    variable == "log2FoldChange_LPD3" ~ 0.0719,
    variable == "log2FoldChange_LPD2" ~ 0.0991,
    variable == "log2FoldChange_LPD1" ~ 0.371,
  ) )

Ordered_Gene_LPD1 <- Join_LPD1_melt_abs_meangamma %>%
  arrange(Gene)

NA_Gene_LPD1 <- na.omit(Ordered_Gene_LPD1)

GeneVector <- as.character(unique(NA_Gene_LPD1$Gene))

for (Gene1 in GeneVector) {
  
  #Filters for each Gene
  Single_gene <- NA_Gene_LPD1 %>%
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

Sig.Genes <- NA_Gene_LPD1 %>%
  dplyr::filter(Gene %in% c("CLPB",
                         "C9orf78",
                         "CTDSPL",
                         "CTXN2",
                         "CUL4A",
                         "FAM189B",
                         "FAM222B",
                         "GYPA",
                         "LINC01094",
                         "MEPCE",
                         "NOL6",
                         "OLA1",
                         "OLA1P1",
                         "POGK",
                         "POLR1A",
                         "PPIF",
                         "PTCHD4",
                         "RPL26P27",
                         "RRP12",
                         "SFPQ",
                         "SH3GL2",
                         "SULT1A3",
                         "TAF6",
                         "TECPR2",
                         "TRIM72",
                         "TUBA3C",
                         "TYRO3P",
                         "UBQLN4",
                         "UQCRBP2"))

write.csv(Sig.Genes, "Sig_genes_LPD_1.csv")


GeneVector2 <- as.character(unique(Sig.Genes$Gene))

for (Gene1 in GeneVector2) {
  
  #Filters for each Gene
  Single_gene1 <- Sig.Genes %>%
    dplyr::filter(Gene == Gene1)
  
  superplot <- ggplot(Single_gene1, aes(x = meangamma, y = abs_Log2FC)) +
    geom_point(shape = 1) + ggtitle(paste0("Scatter Plot of LPD 1 for " , Gene1)) + geom_smooth(method = "lm")
  
  pdf(paste0("LPD_1/Scatterplot_LPD1_",Gene1,".pdf"))
  print(superplot)
  dev.off()
}



m <- lm(abs_Log2FC ~ meangamma, Sig_Genes2)
summary(m)$r.squared

ggplot(Sig_Genes2, aes(x = meangamma , y = abs_Log2FC, color = Gene)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE) + 
  geom_text(label = summary(m)$r.squared, x = 0.3, y = 1.0, show.legend = FALSE)


Sig_Genes2 <- Sig.Genes %>%
  dplyr::filter(Gene == "GYPA")

ggplot(Negative_cor, aes(x = meangamma , y = abs_Log2FC, color = Gene)) +
  geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE)
  ggtitle("Scatter Plot of Max LPD 1 Top 50 overexpressed genes") +
  labs(x = "Log2 Fold Change LPD1", y = "Log2 Fold Change") +
  geom_smooth(method = "lm") + labs(color = "Log1 Fold Change of LPD")


# Up ----------------------------------------------------------------------

  Positive_cor <- Sig.Genes %>%
    dplyr::filter(Gene %in% c("GYPA", "TYRO3P"))
  
  ggplot(Positive_cor, aes(x = meangamma , y = abs_Log2FC, color = Gene)) +
    geom_point() + geom_smooth(method = "lm", colour = "orange", se = FALSE)

# Sig diff genes 2 ------------------------------------------------------------------

  Sig_genes_LPD1 <- Gene_expression_LPD1 %>%
    dplyr::filter(Gene %in% c("CLPB",
                            "C9orf78",
                            "CTDSPL",
                            "CTXN2",
                            "CUL4A",
                            "FAM189B",
                            "FAM222B",
                            "GYPA",
                            "LINC01094",
                            "MEPCE",
                            "NOL6",
                            "OLA1",
                            "OLA1P1",
                            "POGK",
                            "POLR1A",
                            "PPIF",
                            "PTCHD4",
                            "RPL26P27",
                            "RRP12",
                            "SFPQ",
                            "SH3GL2",
                            "SULT1A3",
                            "TAF6",
                            "TECPR2",
                            "TRIM72",
                            "TUBA3C",
                            "TYRO3P",
                            "UBQLN4",
                            "UQCRBP2"))
  
  Sig_diff_exp_LPD1 <- Sig_genes_LPD1 %>%
    dplyr::filter(abs(log2FoldChange) >= 1) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
    select(Gene,log2FoldChange,padj)

  write.csv(Sig_diff_exp_LPD1,"Sig_diff_genes_LPD_1.csv")  
  