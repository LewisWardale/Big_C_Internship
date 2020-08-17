
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")

library(clusterProfiler)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)


# LPD1 --------------------------------------------------------------------

LPD1_genes <- c("CLPB",
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
                "UQCRBP2")

geneNamesTransformed <- bitr(LPD1_genes, fromType = "SYMBOL",
                             toType=c("ENTREZID", "ENSEMBL"),
                             OrgDb="org.Hs.eg.db") %>%
  dplyr::filter(is.na(ENTREZID) == FALSE)

keggAnalysis <- enrichKEGG(gene = geneNamesTransformed$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cnetplot(keggAnalysis, showCategory = 10)

keggAnalysis

# LPD2 --------------------------------------------------------------------

LPD2_genes <- c("ARHGAP40",
                "CD200R1L-AS1",
                "IGFN1",
                "TAC1")

geneNamesTransformed <- bitr(LPD2_genes, fromType = "SYMBOL",
                             toType=c("ENTREZID", "ENSEMBL"),
                             OrgDb="org.Hs.eg.db") %>%
  dplyr::filter(is.na(ENTREZID) == FALSE)

keggAnalysis <- enrichKEGG(gene = geneNamesTransformed$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cnetplot(keggAnalysis, showCategory = 10)

# LPD3 --------------------------------------------------------------------

LPD3_genes <- c("SLC7A3",
                "RNU7-96P",
                "LINC02023",
                "COL17A1",
                "AKT3")

geneNamesTransformed <- bitr(LPD3_genes, fromType = "SYMBOL",
                             toType=c("ENTREZID", "ENSEMBL"),
                             OrgDb="org.Hs.eg.db") %>%
  dplyr::filter(is.na(ENTREZID) == FALSE)

keggAnalysis <- enrichKEGG(gene = geneNamesTransformed$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cnetplot(keggAnalysis, showCategory = 10)

# LPD4--------------------------------------------------------------------

LPD4_genes <- c("CUL7",
                "HNRNPA3P7",
                "LINC02023",
                "MAP6",
                "LINC02706",
                "MON1A",
                "POLR2M",
                "RNU1-30P",
                "RNU6-157P",
                "SCARNA6")

geneNamesTransformed <- bitr(LPD4_genes, fromType = "SYMBOL",
                             toType=c("ENTREZID", "ENSEMBL"),
                             OrgDb="org.Hs.eg.db") %>%
  dplyr::filter(is.na(ENTREZID) == FALSE)

keggAnalysis <- enrichKEGG(gene = geneNamesTransformed$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cnetplot(keggAnalysis, showCategory = 10)

# LPD5 --------------------------------------------------------------------

LPD5_genes <- c("API5P2",
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
)
geneNamesTransformed <- bitr(LPD5_genes, fromType = "SYMBOL",
                             toType=c("ENTREZID", "ENSEMBL"),
                             OrgDb="org.Hs.eg.db") %>%
  dplyr::filter(is.na(ENTREZID) == FALSE)

keggAnalysis <- enrichKEGG(gene = geneNamesTransformed$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cnetplot(keggAnalysis, showCategory = 10)

# LPD6 --------------------------------------------------------------------

LPD6_genes <- c("C11orf54",
               "CCDC28A",
               "COX20P1",
               "DRICH1",
               "CHURC1-FNTB",
               "ETFDH",
               "FNDC10",
               "GORASP1",
               "IFITM5",
               "KLF16",
               "KNOP1P4",
               "LILRB4",
               "LINC02809",
               "PAQR6",
               "PCYOX1",
               "PRICKLE3",
               "RBM10",
               "STAC3",
               "TMED1",
               "TMEM52B",
               "UPF1",
               "ZNHIT1",
               "ZNRD2-AS1")

geneNamesTransformed <- bitr(LPD6_genes, fromType = "SYMBOL",
                             toType=c("ENTREZID", "ENSEMBL"),
                             OrgDb="org.Hs.eg.db") %>%
  dplyr::filter(is.na(ENTREZID) == FALSE)

keggAnalysis <- enrichKEGG(gene = geneNamesTransformed$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cnetplot(keggAnalysis, showCategory = 10)

# LPD7 --------------------------------------------------------------------

LPD7_genes <- c("ARHGEF34P",
                "C19orf67",
                "DIS3",
                "FTLP1",
                "GTF2H1",
                "HSPE1P8",
                "IGKV1-33",
                "IGKV2D-28",
                "IGLV1-41",
                "IGLV3-30",
                "IGLV4-3",
                "IGLV4-60",
                "IGLV9-49",
                "RRP15",
                "SNORA65")

geneNamesTransformed <- bitr(LPD7_genes, fromType = "SYMBOL",
                             toType=c("ENTREZID", "ENSEMBL"),
                             OrgDb="org.Hs.eg.db") %>%
  dplyr::filter(is.na(ENTREZID) == FALSE)

keggAnalysis <- enrichKEGG(gene = geneNamesTransformed$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.01) %>%
  setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

cnetplot(keggAnalysis, showCategory = 10)

