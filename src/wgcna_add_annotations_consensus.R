library(WGCNA)
library(DESeq2)
library(tximport)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

#####
# load data
#####
cons <- readRDS('data/wgcna/consensus_wgcna_objs.R')
nodes <- read_delim('data/wgcna/consensus.nodes',delim='\t')
colnames(nodes) <- c('nodeName','altName','module')


heart_data <- readRDS('data/wgcna/heart_wgcna_objs.rds')
brain_data <- readRDS('data/wgcna/brain_wgcna_objs.rds')
liver_data <- readRDS('data/wgcna/liver_wgcna_objs.rds')
kidney_data <- readRDS('data/wgcna/kidney_wgcna_objs.rds')

#####
# iter through datasets
# add relevant gene expression data
#####

genes <- nodes$nodeName
temps <- c('4C','12C','20C','25C','30C','SummerActive')
contrasts <- list(c("SummerActive", "4C"), c("SummerActive", "12C"), c("SummerActive", "20C"),
                  c("SummerActive", "25C"), c("SummerActive", "30C"),
                  c("30C", "4C"), c("30C", "12C"), c("30C", "20C"),c("30C", "25C"), 
                  c("25C", "4C"), c("25C", "12C"), c("25C", "20C"),
                  c("20C", "4C"), c("20C", "12C"),
                  c("12C", "4C"))
tissues <- c('heart','brain','liver','kidney')
for (tiss in tissues){
  dname <- paste0(tiss, '_data')
  data <- get(dname)
  for (t in temps){ #add rel expr level
    samples <- data$meta %>% filter(temp==t) %>% pull(sample)
    samples <- samples[which(samples %in% rownames(data$mat))]
    vname <- paste0(tiss,'_rel_expr_',t)
    if (length(samples)>1){nodes <- nodes %>% mutate(!!vname := data$mat[samples, genes] %>% colMeans())}
    else{nodes <- nodes %>% mutate(!!vname := data$mat[samples, genes])}
  }
  for (c in contrasts){ #add contrast DESEQ results
    res <- results(data$dds, contrast=c('temp',c))
    res <- res[genes,]
    vname <- paste0(tiss,'_contrast_',c[1],'_',c[2],'_lfc')
    nodes <- nodes %>% mutate(!!vname := res$log2FoldChange)
    vname <- paste0(tiss,'_contrast_',c[1],'_',c[2],'_stat')
    nodes <- nodes %>% mutate(!!vname := res$stat)
  }
}

write_delim(nodes,'data/wgcna/consensus.nodes.1',delim='\t')
