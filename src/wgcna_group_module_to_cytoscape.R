library(WGCNA)
library(DESeq2)
library(tximport)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)

######
# import data from single wgcna runs
#####

heart_data <- readRDS('data/wgcna/heart_wgcna_objs.rds')
brain_data <- readRDS('data/wgcna/brain_wgcna_objs.rds')
liver_data <- readRDS('data/wgcna/liver_wgcna_objs.rds')
kidney_data <- readRDS('data/wgcna/kidney_wgcna_objs.rds')



#####
#choose a dataset and a module to map in the heatmap
#####

data <- heart_data

modules <- c('lightcyan1','green')
group <- 'transition_torpid'
tissue <- 'heart'
module <- 'darkturquoise'
mat <- data$mat
meta <- data$meta %>% filter(sample %in% rownames(mat))
moduleColors <- data$moduleColors
MEs <- data$MEs
modNames <- substring(names(MEs),3)
datTraits <- meta[match(rownames(mat), meta$sample),] %>% 
  mutate(temp=case_when(temp=='SummerActive'~'37C', .default = meta$temp)) %>% 
  select(-c(sample))
nums <- as.numeric(str_sub(datTraits$temp, end=-2))
datTraits$temp <- nums 

geneModuleMembership <- as.data.frame(cor(mat, MEs, use = 'p'))
geneTraitSignificance <- as.data.frame(cor(mat, datTraits$temp, use = "p"))


moduleGenes = moduleColors %in% modules

#init extra node attributes starting with module color
atr <- tibble(module=moduleColors[moduleGenes])

#fill with relative expression grouped by temp
temps <- c('4C','12C','20C','25C','30C','SummerActive')
for (t in temps){
  samples <- data$meta %>% filter(temp==t) %>% pull(sample)
  samples <- samples[which(samples %in% rownames(data$mat))]
  vname <- paste0('rel_expr_',t)
  atr <- atr %>% mutate(!!vname := data$mat[samples, moduleGenes] %>% colMeans())
  if (!('gene' %in% colnames(atr))){atr <- atr %>% mutate(gene = colnames(data$mat[samples, moduleGenes]))}
  }

#fill with contrast data
contrasts <- list(c("SummerActive", "4C"), c("SummerActive", "12C"), c("SummerActive", "20C"),
                  c("SummerActive", "25C"), c("SummerActive", "30C"),
                  c("30C", "4C"), c("30C", "12C"), c("30C", "20C"),c("30C", "25C"), 
                  c("25C", "4C"), c("25C", "12C"), c("25C", "20C"),
                  c("20C", "4C"), c("20C", "12C"),
                  c("12C", "4C"))
genes <- atr$gene
for (c in contrasts){
  res <- results(data$dds, contrast=c('temp',c))
  res <- res[genes,]
  vname <- paste0('contrast_',c[1],'_',c[2],'_lfc')
  atr <- atr %>% mutate(!!vname := res$log2FoldChange)
  vname <- paste0('contrast_',c[1],'_',c[2],'_stat')
  atr <- atr %>% mutate(!!vname := res$stat)
}

cyt <- exportNetworkToCytoscape(data$TOM[moduleGenes, moduleGenes],
                                edgeFile = paste0('data/wgcna/networks/',tissue,'_',group,'.edges'),
                                nodeFile = paste0('data/wgcna/networks/',tissue,'_',group,'.nodes'),
                                weighted = T,
                                threshold = 0.02,
                                nodeNames = genes,
                                nodeAttr = atr)
