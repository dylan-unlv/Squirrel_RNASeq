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


column = match(module, modNames)
moduleGenes = moduleColors==module



cyt <- exportNetworkToCytoscape(data$TOM[moduleGenes, moduleGenes],
                                edgeFile = paste0('data/wgcna/networks/',tissue,'_',module,'.edges'),
                                nodeFile = paste0('data/wgcna/networks/',tissue,'_',module,'.nodes'),
                                weighted = T,
                                threshold = 0.02,
                                nodeNames = colnames(mat[,moduleGenes]),
                                nodeAttr = moduleColors[moduleGenes])