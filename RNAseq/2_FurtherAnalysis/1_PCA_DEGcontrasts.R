
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
library(data.table)
## -------------------------------------------------------------------------------


#### 1. PCA Analysis ####

## perform vst and PCA
vst_dds <- vst(dds)

data <- plotPCA(vst_dds, intgroup = c( "treatment_group"), returnData=TRUE)

data$treat <- gsub("(.*)_.*", "\\1", data$group)
data$infection <- gsub(".*_(.*)", "\\1", data$group)

### get percentage var of each component
percentVar <- round(100 * attr(data, "percentVar"))


### plot
pca <- ggplot(data, aes(PC1, PC2, color=infection, shape=treat)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic()



#### 2. Create contrasts ####

## create results list res_list containing all contrasts
res_list <- list()


#1 DMSO_Uninfected vs Naive_Uninfected
res <- results(dds, contrast=c(0,1,0,0,-1,0,0))
res <- lfcShrink(dds, res = res, type = "ashr")
du_nu <- res
res_list[["dmso.uninfect_naive.uninfect"]] <- du_nu

#2 JQ1_Uninfected vs Naive_Uninfected
res <- results(dds, contrast=c(0,0,0,1,-1,0,0))
res <- lfcShrink(dds, res = res, type = "ashr")
ju_nu <- res
res_list[["jq1.uninfect_naive.uninfect"]] <- ju_nu

#3 DMSO_Infected vs DMSO_Uninfected
res <- results(dds, contrast=c(1,-1,0,0,0,0,0))
res <- lfcShrink(dds, res = res, type = "ashr")
di_du <- res
res_list[["dmso.infect_dmso.uninfect"]] <- di_du

#4 JQ1_Infected vs JQ1_Uninfected
res <- results(dds, contrast=c(0,0,1,-1,0,0,0))
res <- lfcShrink(dds, res = res, type = "ashr")
ji_ju <- res
res_list[["jq1.infect_jq1.uninfect"]] <- ji_ju

#5 JQ1_Uninfected vs DMSO_Uninfected
res <- results(dds, contrast=c(0,-1,0,1,0,0,0))
res <- lfcShrink(dds, res = res, type = "ashr")
ju_du <- res
res_list[["jq1.uninfect_dmso.uninfect"]] <- ju_du

#6 JQ1_Infected vs DMSO_Infected
res <- results(dds, contrast=c(-1,0,1,0,0,0,0))
res <- lfcShrink(dds, res = res, type = "ashr")
ji_di <- res
res_list[["jq1.infect_dmso.infect"]] <- ji_di

#7 JQ1_Infection interaction
res <- results(dds, contrast=c(-1,1,1,-1,0,0,0))
res <- lfcShrink(dds, res = res, type = "ashr")
interact <- res
res_list[["interact"]] <- interact

## create dataframe with differential expression results of all contrasts

test.list <- lapply(res_list, as.data.frame)

for (i in names(test.list)){
  test.list[[i]]$contrast <- i
  test.list[[i]]$ens_id <- rownames(test.list[[i]])
}

## dataframe with all contrasts
all_res.df <- rbindlist(test.list)





