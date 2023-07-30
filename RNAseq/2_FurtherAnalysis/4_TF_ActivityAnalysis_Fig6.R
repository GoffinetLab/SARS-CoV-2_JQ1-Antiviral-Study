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


#### 1. TF activity anaylsis ####


## retrieve normalised counts
norm_counts <- data.frame(counts(dds, normalized=T)) %>% 
  rownames_to_column(var = "ensID") %>% 
  rowwise() %>% 
  mutate(var = var(c_across(2:13))) %>%
  arrange(desc(var))
norm_counts
## retrieve genes
norm_counts$gene <- mapIds(EnsDb.Hsapiens.v86, keys= norm_counts$ensID, column="SYMBOL", keytype="GENEID", multiVals="first")


norm_counts1 <- norm_counts %>%
  drop_na() %>%
  group_by(gene) %>% filter(row_number() == 1) %>%
  ungroup() %>%
  column_to_rownames(var = "gene") %>%
  dplyr::select(-c( var, ensID))

## retrieve regulons
regulons = dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

tf_activities <- run_viper(norm_counts1, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))



#### 2. TF activity heatmap Fig 6A ####

## get mean TF score per treatment condition

tfs <- data.frame(tf_activities) %>% 
  rownames_to_column(var = "TF") %>% 
  gather( key = "sample", value = "score", -c(TF)) %>%
  mutate(treatment = gsub("_[123]$", "", sample)) %>%
  group_by(treatment, TF) %>%
  summarise(mean_score = mean(score), var_score = var(score))


#### group by tf and calc var in mean score, then rank
## pull out most var TFs by variance (60)
high_var_tfs  <- tfs %>%
  group_by(TF) %>%
  summarise(var =var(mean_score)) %>%
  ungroup() %>%
  top_n(60, var) %>%
  arrange(desc(var))


tfs2 <-tfs %>%
  dplyr::select(-var_score) %>%
  semi_join(high_var_tfs, by ="TF") %>%
  spread(TF, mean_score) %>%
  data.frame(row.names = 1, check.names = FALSE)


## draw heatmap


cols = colorRampPalette(c("#0000FF", "#FFFF00"))(30)

pheatmap(t(tfs2), color = cols, scale = "row",
              fontsize=14, 
              fontsize_row = 10,
              main = "TF Activity", angle_col = 45,
              treeheight_col = 10)