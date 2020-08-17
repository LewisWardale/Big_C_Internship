

gene_expressionList <- list()

for(number in 1:7){
  
  gene_expressionList[[number]] <- read.csv(paste0("/Users/lewiswardale/Desktop/DEGenes_BReast/DE-Genes_TCGA-BRCA_LPD_", number, "_eridanus.csv"))
  
  # Filtering the data for  Log2fold abs 0-1, padj P < 0.05 and ENGS Genes
  
  gene_overexpression <- gene_expressionList[[number]] %>%
    dplyr::filter(abs(log2FoldChange) >= 1) %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
    top_n(n = 50, wt = log2FoldChange) %>%
    select(Gene,log2FoldChange,padj)
  
  # Select all other LPD for Log2FC, padj, Gene
  
  gene_overexpression_other <- lapply(gene_expressionList[-number], function(x) dplyr::select(x, Gene,log2FoldChange,padj))
  
  
  tellMeSomethingBoy <- gene_overexpression_other %>%
    imap(.x = . , ~ set_names(.x, c("Gene", .y))) %>%
    purrr::reduce(left_join, by = "Gene") 
    right_join(gene_overexpression, by = "Gene")
  
  
  #
  Join_LPD1_overexpression <- right_join(by = "Gene", gene_expressionList[[!number]],gene_expressionList[[number]])
  
  Join_LPD1_overexpression <- rename(Join_LPD1_overexpression, c("log2FoldChange" = "log2FoldChange_LPD7",
                                   "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                   "log2FoldChange.y" = "log2FoldChange_LPD5",
                                   "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                   "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                   "log2FoldChange.x" = "log2FoldChange_LPD2",
                                   "log2FoldChange.y.y.y" = "log2FoldChange_LPD1"))
  
  write.csv(Join_LPD1_overexpression, 'LPD1_Top50_OverExpressed.csv')
  
  MeanGamma <- c(mean(log2FoldChange_LPD1), 
                 mean(log2FoldChange_LPD2),
                 mean(log2FoldChange_LPD3), 
                 mean(log2FoldChange_LPD4), 
                 mean(log2FoldChange_LPD5), 
                 mean(log2FoldChange_LPD6), 
                 mean(log2FoldChange_LPD7))

  LPD <- c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7")
  
  Mean_gamma_table(pathway = "/Users/lewiswardale/Desktop/TCGA-BRCA Assigned LPD/GammaValues_TCGA-BRCA_eridanus.csv")
  
  MeanLog2FC <- meanGammaList[1,]
  
  Overexpression_gamma <- data.frame(LPD, MeanGamma, MeanLog2FC)
  
  ggplot(Overexpression_gamma, aes(x = MeanGamma, y = MeanLog2FC)) +
    geom_point() + 
    ggtitle("Scatter Plot of Max LPD 1 mean gamma values and the mean Log2 Fold Change of Top 50 overexpressed genes") +
    labs(x = "Mean Gamma", y = "Mean Log2 Fold Change") +
    geom_smooth(method = "lm") +
    geom_text(aes(label= LPD), hjust = -0.1)
}




library(tidyverse)
library(survminer)
library(survival)
library(DataExplorer)
library(magrittr)
library(plyr)



# Filtering the data for  Log2fold abs 0-1, padj P < 0.05 and ENGS Genes  -------------------------


  
  gene_overexpression <- Gene_expression_LPD1 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  top_n(n = 50, wt = log2FoldChange) %>%
  select(Gene,log2FoldChange,padj)


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

 Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD2,gene_overexpression_LPD1)
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

write.csv(Join_LPD1, 'LPD1_Top50_OverExpressed.csv')

summary(Join_LPD1)
        #Log2FC  MeanGamma
#LPD1 = 2.812     0.37068121
#LPD2 = -1.9537   0.09910208
#LPD3 = -2.06246     0.07189306
#LPD4 = -3.031     0.03878323
#LPD5 = -2.7675    0.01572141
#LPD6 = -0.9919     0.22591530
#LPD7 = 0.63708     0.17790370

MeanGamma <- c(0.37068121, 0.09910208, 0.07189306, 0.03878323, 0.01572141, 0.22591530, 0.17790370)
LPD <- c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7")
MeanLog2FC <- c(2.812, -1.9537, -2.06246, -3.031,-2.7675, -0.9919, 0.63708)

LPD1_Overexpression_gamma <- data.frame(LPD, MeanGamma, MeanLog2FC)

ggplot(LPD1_Overexpression_gamma, aes(x = MeanGamma, y = MeanLog2FC)) +
  geom_point() + 
  ggtitle("Scatter Plot of Max LPD 1 mean gamma values and the mean Log2 Fold Change of Top 50 overexpressed genes") +
  labs(x = "Mean Gamma", y = "Mean Log2 Fold Change") +
  geom_smooth(method = "lm") +
  geom_text(aes(label= LPD), hjust = -0.1)


# Underexpression ---------------------------------------------------------

gene_underexpression_LPD1 <- Gene_expression_LPD1 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  top_n(n = -50, wt = log2FoldChange) %>%
  select(Gene,log2FoldChange,padj)

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

Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD2,gene_underexpression_LPD1)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD3, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD4, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD5, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD6, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD7, Join_LPD1_under)

Join_LPD1_under <- rename(Join_LPD1_under, c("log2FoldChange" = "log2FoldChange_LPD7",
                                             "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                             "log2FoldChange.y" = "log2FoldChange_LPD5",
                                             "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                             "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                             "log2FoldChange.x" = "log2FoldChange_LPD2",
                                             "log2FoldChange.y.y.y" = "log2FoldChange_LPD1"))

write.csv(Join_LPD1_under, 'LPD1_Top50_UnderExpressed.csv')

summary(Join_LPD1_under)


MeanGamma <- c(0.37068121, 0.09910208, 0.07189306, 0.03878323, 0.01572141, 0.22591530, 0.17790370)
LPD <- c("LPD_1", "LPD_2", "LPD_3", "LPD_4", "LPD_5", "LPD_6", "LPD_7")
MeanLog2FC <- c(-5.937, -3.1910, -3.4441, 3.757, -0.1844, -4.3019, -5.171)

LPD1_Underexpression_gamma <- data.frame(LPD, MeanGamma, MeanLog2FC)

ggplot(LPD1_Underexpression_gamma, aes(x = MeanGamma, y = MeanLog2FC)) +
  geom_point() + 
  ggtitle("Scatter Plot of Max LPD 1 mean gamma values and the mean Log2 Fold Change of Top 50 Underexpressed genes") +
  labs(x = "Mean Gamma", y = "Mean Log2 Fold Change") +
  geom_smooth(method = "lm") +
  geom_text(aes(label= LPD), hjust = -0.1)

# OverExpression vs Expression ------------------------------------------------

gene_overexpression <- Gene_expression_LPD1 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  top_n(n = 50, wt = log2FoldChange) %>%
  select(Gene,log2FoldChange,padj)


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

Join_LPD1 <- right_join(by = "Gene", Gene_expression_LPD2,gene_overexpression_LPD1)
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
Join_LPD1 <- Join_LPD1 %>%
  select(Gene,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)



Join_LPD1_melt <- melt(Join_LPD1, id.vars = c("Gene"))

Join_LPD1_melt_gamma <- Join_LPD1_melt %>%
  mutate(meangamma = case_when(
    variable == "log2FoldChange_LPD7" ~ 0.178,
    variable == "log2FoldChange_LPD6" ~ 0.226,
    variable == "log2FoldChange_LPD5" ~ 0.0152,
    variable == "log2FoldChange_LPD4" ~ 0.0388,
    variable == "log2FoldChange_LPD3" ~ 0.0719,
    variable == "log2FoldChange_LPD2" ~ 0.0991,
    variable == "log2FoldChange_LPD1" ~ 0.371,
  ) )

ggplot(Join_LPD1_melt, aes(x = log2FoldChange_LPD1 , y = value, color = variable)) +
  geom_point() + 
  ggtitle("Scatter Plot of Max LPD 1 Top 50 overexpressed genes") +
  labs(x = "Log2 Fold Change LPD1", y = "Log2 Fold Change") +
  geom_smooth(method = "lm") + labs(color = "Log1 Fold Change of LPD")
# geom_text(aes(label= Gene), hjust = -0.1)

Boxplot <- ggplot(Join_LPD1_melt_gamma, aes(x = meangamma, y = value, fill = variable)) + geom_boxplot(outlier.colour = "red") + labs(y = "Log2 Fold", x = "mean gamma") + labs(fill = "LPD")
Boxplot

# UnderExpression vs Expression ------------------------------------------------

gene_underexpression_LPD1 <- Gene_expression_LPD1 %>%
  dplyr::filter(abs(log2FoldChange) >= 1) %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::filter(!str_detect(Gene, "^ENSG")) %>%
  top_n(n = -50, wt = log2FoldChange) %>%
  select(Gene,log2FoldChange,padj)

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

Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD2,gene_underexpression_LPD1)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD3, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD4, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD5, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD6, Join_LPD1_under)
Join_LPD1_under <- right_join(by = "Gene", Gene_expression_LPD7, Join_LPD1_under)

Join_LPD1_under <- rename(Join_LPD1_under, c("log2FoldChange" = "log2FoldChange_LPD7",
                                             "log2FoldChange.x.x.x" = "log2FoldChange_LPD6",
                                             "log2FoldChange.y" = "log2FoldChange_LPD5",
                                             "log2FoldChange.x.x" = "log2FoldChange_LPD4",
                                             "log2FoldChange.y.y" = "log2FoldChange_LPD3",
                                             "log2FoldChange.x" = "log2FoldChange_LPD2",
                                             "log2FoldChange.y.y.y" = "log2FoldChange_LPD1"))

Join_LPD1_under <- Join_LPD1_under %>%
  select(Gene,
         log2FoldChange_LPD7,
         log2FoldChange_LPD6,
         log2FoldChange_LPD5,
         log2FoldChange_LPD4,
         log2FoldChange_LPD3,
         log2FoldChange_LPD2,
         log2FoldChange_LPD1)



Join_LPD1_under_melt <- melt(Join_LPD1_under, id.vars = c("Gene"))

Join_LPD1_under_melt_gamma <- Join_LPD1_under_melt %>%
  mutate(meangamma = case_when(
    variable == "log2FoldChange_LPD7" ~ 0.178,
    variable == "log2FoldChange_LPD6" ~ 0.226,
    variable == "log2FoldChange_LPD5" ~ 0.0152,
    variable == "log2FoldChange_LPD4" ~ 0.0388,
    variable == "log2FoldChange_LPD3" ~ 0.0719,
    variable == "log2FoldChange_LPD2" ~ 0.0991,
    variable == "log2FoldChange_LPD1" ~ 0.371,
  ) )


ggplot(Join_LPD1_under_melt_gamma, aes(x = meangamma , y = value, color = variable)) +
  geom_point() + 
  ggtitle("Scatter Plot of Max LPD 1 Top 50 underexpressed genes") +
  labs(x = "Mean Gamma", y = "Log2 Fold Change") +
  geom_smooth(method = "lm") + labs(color = "Log2 Fold Change of LPD") +
  geom_text(aes(label= Gene), hjust = -0.1)


Boxplot <- ggplot(Join_LPD1_under_melt_gamma, aes(x = meangamma, y = value, fill = variable)) + geom_boxplot(outlier.colour = "red") + labs(y = "Log2 Fold", x = "mean gamma") + labs(fill = "LPD")
Boxplot

Mean_gamma_table(pathway = "/Users/lewiswardale/Desktop/TCGA-BRCA Assigned LPD/GammaValues_TCGA-BRCA_eridanus.csv")

#LPD1 0.371
#LPD2 0.0991
#LPD3 0.0719
#LPD4 0.0388
#LPD5 0.0157
#LPD6 0.226
#LPD7 0.178
