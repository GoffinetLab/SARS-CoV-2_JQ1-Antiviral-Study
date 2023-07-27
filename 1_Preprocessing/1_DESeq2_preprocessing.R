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


#read in counts and coldata
counts <- read.table("sars2_jq1_counts.Rmatrix.txt",
                     sep = "\t",
                     row.names = 1, header = T)

coldata <- read.csv("sars2_jq1_coldata.csv",
                    row.names = 1, header = T)

## add column for rep covariate and treatment group

coldata$treatment_group <- paste0(coldata$treat, "_", coldata$infection)
coldata$rep <- gsub(pattern = ".*_([1-3]{1})", replacement = "\\1", x = rownames(coldata))
coldata$rep <- paste0("R", coldata$rep)

## create DESEQ2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ 0 + treatment_group + rep)
dds<- DESeq(dds)



##prefilter genes with at least 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## perform variance stabilising transformation
vst_dds <- vst(dds)

## create expression matrix
#mat <- counts(dds, normalized=TRUE)
mat <- assay(vst_dds)
expr_df <- mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  gather(treatment_group, expression, 2:16) %>%
  mutate(treat = gsub("^(.*)_.*_[1-3]", "\\1", treatment_group)) %>%
  mutate(infect = gsub("^.*_(.*)_[1-3]", "\\1", treatment_group)) %>%
  mutate(rep = gsub("^.*_.*_([1-3])", "\\1", treatment_group)) 

# assign HGNC IDs
expr_df$hgnc <- mapIds(EnsDb.Hsapiens.v86, keys= expr_df$gene, column="SYMBOL", keytype="GENEID", multiVals="first")
# ensure SARS2 genes are assigned in hgnc column
expr_df$hgnc <- ifelse(is.na(expr_df$hgnc), expr_df$gene, expr_df$hgnc)
