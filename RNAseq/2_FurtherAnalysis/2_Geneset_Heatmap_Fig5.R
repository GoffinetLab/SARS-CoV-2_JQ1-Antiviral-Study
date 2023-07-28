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






mat <- assay(vst_dds)
mat <- data.frame(mat) %>%
  rownames_to_column(var = "ensID") %>% 
  rowwise() %>% 
  mutate(var = var(c_across(2:13))) %>%
  arrange(desc(var))
mat

## retrieve genes
mat$gene <- mapIds(EnsDb.Hsapiens.v86, keys= mat$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")
mat$gene <- ifelse(is.na(mat$gene), mat$ensID, mat$gene)

mat
allgs_mat
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
allgs_mat

unique(allgs_mat$gene_set)

gs <- "reactome_ifn_signalling"

.mat <- allgs_mat %>% dplyr::filter(gene_set == gs)
.mat
ord


allgs_mat2 <- map_dfr(unique(allgs_mat$gene_set), ~ {
  gs <- .
  .mat <- allgs_mat %>% dplyr::filter(gene_set == gs)
  ord <- hclust(dist(as.matrix(.mat %>% dplyr::select( -gene_set))))$order
  .mat[ord, ]
})

allgs_mat2
colnames(allgs_mat2) <- gsub("^([0-9]*)_(.*)$", "\\2 \\1", colnames(allgs_mat2))

##############################

expr_df <- mat %>%
  as.data.frame() %>%
  rownames_to_column("PrimaryID") %>%
  pivot_longer(-PrimaryID) %>%
  left_join(covar %>% select(label, treatment, group), by=c(name="label")) %>%
  left_join(annot %>% select(PrimaryID, annotation, SYMBOL, ENSEMBL), by="PrimaryID")


all_genesets_genes <- read.csv("genesets_allgenes.csv")
colors <- colorRampPalette(c("#0000FF", "#FFFF00"))(30)

allgs_mat <- mat %>% 
  as.data.frame() %>%
  rownames_to_column("PrimaryID") %>%
  mutate(gene = annot$SYMBOL[ match(PrimaryID, annot$PrimaryID) ]) %>%
  mutate(gene_set=all_genesets_genes$pway[ match(gene, all_genesets_genes$gene) ]) %>%
  filter(gene %in% all_genesets_genes$gene) %>%
  arrange(gene_set)
rownames(allgs_mat) <- allgs_mat$PrimaryID

allgs_mat2 <- map_dfr(unique(allgs_mat$gene_set), ~ {
  gs <- .
  .mat <- allgs_mat %>% dplyr::filter(gene_set == gs)
  ord <- hclust(dist(as.matrix(.mat %>% dplyr::select(-PrimaryID, -gene, -gene_set))))$order
  .mat[ord, ]
})

allgs_mat$gene_set
allgs_mat %>% dplyr::select(gene_set)
allgs_mat
colnames(allgs_mat2) <- gsub("^([0-9]*)_(.*)$", "\\2 \\1", colnames(allgs_mat2))
annotation_row = allgs_mat %>% select(gene_set)


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
#ggsave("ATACseq_heatmap_allGenesets_allcontrasts.pdf", plot = x, width = 8, height = 15)

clust_hmap <- ggplotify::as.ggplot(x)

ggsave(filename = "RNAseq_hmap_allgenesets_annotated.svg", width = 8, height = 12)
ggsave(filename = "RNAseq_hmap_allgenesets_annotated.pdf", width = 8, height = 12)