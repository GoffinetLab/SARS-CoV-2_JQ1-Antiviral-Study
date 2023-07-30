#### LOAD LIBRARIES ####

suppressMessages(library(DESeq2))
library(pheatmap)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(fgsea)
library(msigdbr)
library(knitr)
library(ggplot2)
library(sjmisc)
library(org.Hs.eg.db)
library(biomaRt)
library(dorothea)
library(ggpubr)
library(ashr)
library(clusterProfiler)
## -------------------------------------------------------------------------------


#### 1. Heatmaps of genes in selected genesets ####

## assign colours

col_ifnisg = "#FF0000"
col_autoph = "#0000FF"
col_ccycle = "forestgreen"
col_nrf2 = "#FFA500"


# load csv with pathways of interest and associated genes
all_genesets_genes <- read.csv("C:/Users/postmusd/OneDrive - Charité - Universitätsmedizin Berlin/Dokumente/Push/SARS2_JQ1/publication_figs/genesets_allgenes.csv")

# read in NRF2 regulated genes

nrf2_reg <- c("GSR", "TUBA4A", "OSGIN1", "HMOX1", "GCLM",
              "KEAP1", "PSMC5", "ITPKC", "GADD45C", "TNXB",
              "MTHFR", "MAFG", "TUBA1C", "FTL", "TNXB",
              "NQO1", "FOS", "NRF1", "PRDX1", "MAFF", "NOTCH4",
              "GCLC", "EPHX1", "GSTP1", "TXNRD1", "SQSTM1",
              "ALDH3A1", "CYP1A1", "EPHX1", "GSTP1", "TXNRD1", "TXNRD3",
              "GSS", "GSTA1", "PPARG", "GPX3", "PSMD2", "GADD45B",
              "EGR1", "ENO1", "HMOX2")


## add nrf2 reg genes
nrf2_reg_df <- data.frame( "gene" = nrf2_reg, "pway" = "NRF2_regulated")
all_genesets_genes <- rbind(all_genesets_genes, nrf2_reg_df)


## retrieve variance stabilised expression scores

mat <- assay(vst_dds)
mat <- data.frame(mat) %>%
  rownames_to_column(var = "ensID") %>% 
  rowwise() %>% 
  mutate(var = var(c_across(2:13))) %>%
  arrange(desc(var))

## retrieve genes
mat$gene <- mapIds(EnsDb.Hsapiens.v86, keys= mat$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")
mat$gene <- ifelse(is.na(mat$gene), mat$ensID, mat$gene)

## create expression matrix for genes in genesets of interest
allgs_mat <- mat %>%
  drop_na() %>%
  dplyr::filter(gene %in% all_genesets_genes$gene) %>%
  rowwise() %>% 
  group_by(gene) %>% dplyr::filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(match(gene, all_genesets_genes)) %>%
  column_to_rownames(var = "gene") %>%
  dplyr::select(-c( var, ensID)) %>%
  mutate(gene_set=all_genesets_genes$pway[ match(rownames(.), all_genesets_genes$gene) ])

## cluster genes within each gene set and reset order of expression matrix for plotting in heatmap
allgs_mat2 <- map_dfr(unique(allgs_mat$gene_set), ~ {
  
  gs <- .
  .mat <- allgs_mat %>% dplyr::filter(gene_set == gs)
  ord <- hclust(dist(as.matrix(.mat %>% dplyr::select( -gene_set))))$order
  .mat[ord, ]
})

allgs_mat2
colnames(allgs_mat2) <- gsub("^([0-9]*)_(.*)$", "\\2 \\1", colnames(allgs_mat2))

## plot heatmap

colors <- colorRampPalette(c("#0000FF", "#FFFF00"))(30)
x <- pheatmap(allgs_mat2 %>% dplyr::select( -gene_set), 
              annotation_row = allgs_mat %>% dplyr::select(gene_set),
              annotation_colors = list(gene_set=c(
                kegg_autophagy_animal=col_autoph,
                reactome_ifn_signalling=col_ifnisg,
                reactome_cellcycle=col_ccycle,
                NRF2_regulated=col_nrf2)),
              #labels_row = annot_allgs$gene,
              cluster_rows = F, treeheight_row=0, 
              color = colors, 
              show_rownames = F, scale = "row")
