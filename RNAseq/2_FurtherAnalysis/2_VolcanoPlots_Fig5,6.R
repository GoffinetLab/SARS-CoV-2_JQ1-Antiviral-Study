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


#### 1. Volcano Plots Fig5, Fig6 ####

## Figure 5A
## VOLCANO PLOTS of signficantly differentially regulated genes

## add signifier of statistical significance
all_res.df$Significant <- ifelse(all_res.df$padj > 0.05 | abs(all_res.df$log2FoldChange) < 1,
                                 "NS",
                                 ifelse(all_res.df$log2FoldChange > 0, "UP", "DOWN")) 

## approximate adj P-values of 0 to 5e-324
all_res.df <-all_res.df %>%
  mutate(padj = ifelse(padj == 0, 5e-324, padj))

all_res.df$yvar <- -log10(all_res.df$padj) 

## plot volcano plots showing up and downregulated genes
volc <- ggplot(all_res.df, aes(x=log2FoldChange, y=yvar, color=Significant)) + 
  geom_point() + scale_color_manual(values=c(NS="darkgrey", UP="#FFFF00",
                                             DOWN="#0000FF")) +
  xlim(-20,20) +
  ylim(0,400)+
  theme_bw() +
  ylab("-log(padj)")+
  facet_wrap(~contrast)


## Figure 6B,C,D
## VOLCANO PLOTS highlighting expression of genes from different genesets
## Interferon stimulated / NRF2 regulated

## assign HGNC identifiers to genes
all_res.df$hgnc <- mapIds(EnsDb.Hsapiens.v86, keys= all_res.df$ens_id, column="SYMBOL", keytype="GENEID", multiVals="first")

## load list of interferon stimulated genes from REACTOME R-HSA-913531
reactome_ifn_sig <- read.csv("interferon_signalling.csv")
## identify genes as ISGs
all_res.df$isg <- ifelse(all_res.df$hgnc %in% reactome_ifn_sig$gene, "isg", "non_isg")

## load list of NRF2 regulated genes
nrf2_reg <- c("GSR", "TUBA4A", "OSGIN1", "HMOX1", "GCLM",
              "KEAP1", "PSMC5", "ITPKC", "GADD45C", "TNXB",
              "MTHFR", "MAFG", "TUBA1C", "FTL", "TNXB",
              "NQO1", "FOS", "NRF1", "PRDX1", "MAFF", "NOTCH4",
              "GCLC", "EPHX1", "GSTP1", "TXNRD1", "SQSTM1",
              "ALDH3A1", "CYP1A1", "EPHX1", "GSTP1", "TXNRD1", "TXNRD3",
              "GSS", "GSTA1", "PPARG", "GPX3", "PSMD2", "GADD45B",
              "EGR1", "ENO1", "HMOX2")

# identify genes as NRF2(NFE2L2) and NRF2 regulation status
all_res.df$nrf2_stat <- ifelse(all_res.df$hgnc == "NFE2L2", "NRF2", 
                               ifelse(all_res.df$hgnc %in% nrf2_reg, "NRF2_reg", "NRF2_nonreg"))
all_res.df$nrf2_stat <- ifelse(is.na(all_res.df$hgnc), "NRF2_nonreg", all_res.df$nrf2_stat)


## identify if genes are ISGs, NRF2 regulated, both or NRF2

all_res.df$isg_nrf2 <- ifelse(all_res.df$isg == "isg" & all_res.df$nrf2_stat == "NRF2_nonreg" & all_res.df$significant == "Significant", "isg", 
                              ifelse(all_res.df$isg == "non_isg" & all_res.df$nrf2_stat == "NRF2_reg" & all_res.df$significant == "Significant", "nrf2_reg",
                                     ifelse(all_res.df$isg == "isg" & all_res.df$nrf2_stat == "NRF2_reg" & all_res.df$significant == "Significant", "both",
                                            "none")))
## mark nrf2
all_res.df$isg_nrf2[all_res.df$hgnc == "NFE2L2"] <- "NRF2"

## plot volcano plot indicating ISG status/NRF2 regulation status/whether a gene is NRF2

volc_nrf2_isg <- ggplot() +
  geom_point(data = all_res.df[all_res.df$isg_nrf2 == "none",], aes(x=log2FoldChange, y=yvar,    stat="identity",  label = "none", colour = isg_nrf2)) +
  geom_point(data = all_res.df[all_res.df$isg_nrf2 == "isg",], aes(x=log2FoldChange, y=yvar,    stat="identity",  label = "isg", colour = isg_nrf2)) +
  geom_point(data = all_res.df[all_res.df$isg_nrf2 == "nrf2_reg",], aes(x=log2FoldChange, y=yvar,    stat="identity",  label = "nrf2_reg", colour = isg_nrf2)) +
  geom_point(data = all_res.df[all_res.df$isg_nrf2 == "both",], aes(x=log2FoldChange, y=yvar,    stat="identity",  label = "both", colour = isg_nrf2)) +
  geom_point(data = all_res.df[all_res.df$isg_nrf2 == "NRF2",], aes(x=log2FoldChange, y=yvar,    stat="identity",  label = "NRF2", colour = isg_nrf2)) +
  
  theme_bw() +
  labs(y = "-log adj P", x = "log2 FC") +
  scale_colour_manual(values = c(isg = "red" , nrf2_reg = "forestgreen", both = "purple",  NRF2 = "orange", none = "grey")) +
  facet_wrap(~contrast)+
  #geom_text(data=subset(all_res.df[all_res.df$nrf2_stat == "NRF2_reg",], abs(log2FoldChange) > 1 | hgnc %in% c("NQO1", "HMOX1")),
  #aes(x=log2FoldChange, y=yvar, label= hgnc), position = position_jitter(width = 4, height=80), size = 2) +
  xlim(-20,20)  +
  ylim(0,400)
