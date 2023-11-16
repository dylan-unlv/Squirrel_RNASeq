library(tidyverse)
library(biomaRt)
library(WGCNA)
library(ComplexHeatmap)
library(DESeq2)


cons_data <- read_delim('data/cytoscape/consensus.nodes.analyzed.csv', delim=',')

#####
#First, KEGG pathway analysis of modules
#


#need to convert to entrez ids for this
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
altNames <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id"),
                  values=cons_data$nodeName,mart= mart,uniqueRows = F)


library(clusterProfiler) #needs to be import after the previous line... seriously.

#iter modules
mods <- cons_data$module %>% unique()
for (imodule in mods){
  
  #subset data, add entrez ids
  data <- cons_data %>% filter(module==imodule) %>% 
    rename('ensembl_gene_id'=nodeName) %>% 
    dplyr::left_join(., altNames, by='ensembl_gene_id')


  #KEGG / GO Terms
  result <- enrichKEGG(data$entrezgene_id, keyType = 'kegg',
            organism = 'hsa', pAdjustMethod = 'hochberg')

  #write out
  result@result %>% write_csv(paste0('data/kegg/',imodule,'_KEGG.csv'))
}



#####
#Correlate consensus modules with metadata
#
cons_R <- readRDS('data/wgcna/tri_consensus_wgcna_objs.rds')
setLabels <- cons_R$setLabels
nSets <- length(setLabels)

meta <- read_csv('data/temp_metadata.csv')
meta <- meta %>% mutate(temp=case_when(temp=='SummerActive'~37,
                                       temp=='30C'~30, 
                                       temp=='25C'~25, 
                                       temp=='20C'~20,
                                       temp=='12C'~12,
                                       temp=='4C'~4)) 

#get module eigengenes
for (setn in 1:nSets){
  MEs0 <- moduleEigengenes(cons_R$multiExpr[[setn]]$data, cons_R$moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  tmeta <- meta %>% mutate(sample=paste0('CL',sample,str_to_title(setLabels[setn]),'RNA')) 
  tmeta <- tmeta %>% filter(sample %in% rownames(MEs))
  moduleTraitCor <- cor(MEs, tmeta$temp, use ='p', method="spearman")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(tmeta))

  if (setn==1) {
    gmat <- moduleTraitCor
    pmat <- moduleTraitPvalue
  }
  else{
    gmat <- cbind(gmat, moduleTraitCor[match(rownames(gmat), rownames(moduleTraitCor))] )
    pmat <- cbind(pmat, moduleTraitPvalue[match(rownames(pmat), rownames(moduleTraitPvalue))])
  }
  
}

setLabelsordered<- setLabels

gmat <- data.frame(gmat)
pmat <- data.frame(pmat)
colnames(gmat) <- setLabels
colnames(pmat) <- setLabels
gmat <- as.matrix(gmat[setLabelsordered])
pmat <- as.matrix(pmat[setLabelsordered])

#look at correlations within the module eigengenes to temp
pdf(file = paste0('figs/consensus_eigengene_temp_corr.pdf'), wi=9, he=9)
# Will display correlations and their p-values
textMatrix = paste(signif(gmat, 2), "\n(",
                   signif(pmat, 1), ")", sep = "")
dim(textMatrix) = dim(gmat)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot (if many meta values, change xlabels)
labeledHeatmap(Matrix = gmat,
               xLabels = setLabelsordered,
               yLabels = rownames(gmat),
               ySymbols = rownames(gmat),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.y = 0.9,
               zlim = c(-1,1),
               main = paste("Consensus Module-Temp Relationships "))
dev.off()



#####
# module expression graphs
####

#import individual networks
heart_data <- readRDS('data/wgcna/heart_wgcna_objs.rds')
brain_data <- readRDS('data/wgcna/brain_wgcna_objs.rds')
liver_data <- readRDS('data/wgcna/liver_wgcna_objs.rds')
kidney_data <- readRDS('data/wgcna/kidney_wgcna_objs.rds')



#select module, filter data for only selected genes above certain thresh
imodule <- 'blue'
mgenes <- colnames(cons_R$multiExpr[[1]]$data[,which(cons_R$moduleColors==imodule)])

#thresh up
pthresh <- 0.05
lfcthresh <- 2
hgenes <- data.frame(results(heart_data$dds, contrast=c('temp','SummerActive','20C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange>lfcthresh) %>% 
  rownames()
bgenes <- data.frame(results(brain_data$dds, contrast=c('temp','SummerActive','25C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange>lfcthresh) %>% 
  rownames()
kgenes <- data.frame(results(kidney_data$dds, contrast=c('temp','SummerActive','20C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange>lfcthresh) %>% 
  rownames()
lgenes <- data.frame(results(liver_data$dds, contrast=c('temp','SummerActive','20C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange>lfcthresh) %>% 
  rownames()

tgenes <- unique(c(hgenes, bgenes, kgenes, lgenes))

#collect vst data
hexpr <- data.frame(assay(heart_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='heart') %>% 
  left_join(., heart_data$meta[,c('sample','temp')], by='sample')

bexpr <- data.frame(assay(brain_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='brain') %>% 
  left_join(., brain_data$meta, by='sample')

kexpr <- data.frame(assay(kidney_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='kidney') %>% 
  left_join(., kidney_data$meta, by='sample')

lexpr <- data.frame(assay(liver_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='liver') %>% 
  left_join(., liver_data$meta, by='sample')

texpr <- rbind(hexpr, bexpr, kexpr, lexpr) %>% 
  mutate(temp=case_when(temp=='SummerActive'~37,
                               temp=='30C'~30, 
                               temp=='25C'~25, 
                               temp=='20C'~20,
                               temp=='12C'~12,
                               temp=='4C'~4)) %>% 
                  group_by(gene)

# Graph results

gg<-ggplot(texpr, aes(x=temp, y=expr, color=tiss, group=gene, shape=gene))+
  stat_summary(geom='line', fun.y=mean)+
  stat_summary(geom='point', fun.y=mean)+
  scale_x_continuous(breaks=c(4,12,20,25,30,37))+
  scale_shape_manual(values=0:19)+
  facet_wrap(~tiss)+
  theme_bw()

ggsave(paste0('figs/',imodule,'_up_DEG.png'), dpi=300)


#thresh down
lfcthresh <- 0 - lfcthresh
hgenes <- data.frame(results(heart_data$dds, contrast=c('temp','SummerActive','20C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange<lfcthresh) %>% 
  rownames()
bgenes <- data.frame(results(brain_data$dds, contrast=c('temp','SummerActive','25C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange<lfcthresh) %>% 
  rownames()
kgenes <- data.frame(results(kidney_data$dds, contrast=c('temp','SummerActive','20C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange<lfcthresh) %>% 
  rownames()
lgenes <- data.frame(results(liver_data$dds, contrast=c('temp','SummerActive','20C')))[mgenes, ] %>% 
  filter(padj<pthresh, log2FoldChange<lfcthresh) %>% 
  rownames()

tgenes <- unique(c(hgenes, bgenes, kgenes, lgenes))

#collect vst data
hexpr <- data.frame(assay(heart_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='heart') %>% 
  left_join(., heart_data$meta[,c('sample','temp')], by='sample')

bexpr <- data.frame(assay(brain_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='brain') %>% 
  left_join(., brain_data$meta, by='sample')

kexpr <- data.frame(assay(kidney_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='kidney') %>% 
  left_join(., kidney_data$meta, by='sample')

lexpr <- data.frame(assay(liver_data$vst)[tgenes, ]) %>% 
  rownames_to_column('gene') %>% 
  pivot_longer(!gene, names_to = 'sample', values_to = 'expr') %>% 
  mutate(tiss='liver') %>% 
  left_join(., liver_data$meta, by='sample')

texpr <- rbind(hexpr, bexpr, kexpr, lexpr) %>% 
  mutate(temp=case_when(temp=='SummerActive'~37,
                        temp=='30C'~30, 
                        temp=='25C'~25, 
                        temp=='20C'~20,
                        temp=='12C'~12,
                        temp=='4C'~4)) %>% 
  group_by(gene)

# Graph results

gg<-ggplot(texpr, aes(x=temp, y=expr, color=tiss, group=gene, shape=gene))+
  stat_summary(geom='line', fun.y=mean)+
  stat_summary(geom='point', fun.y=mean)+
  scale_x_continuous(breaks=c(4,12,20,25,30,37))+
  scale_shape_manual(values=0:19)+
  facet_wrap(~tiss)+
  theme_bw()

ggsave(paste0('figs/',imodule,'_down_DEG.png'), dpi=300)
