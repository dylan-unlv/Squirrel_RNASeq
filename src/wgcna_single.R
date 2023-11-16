library(WGCNA)
library(DESeq2)
library(tximport)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)


#####
#read in metadata, tximport quants
#####
tissue <- 'kidney'
meta <- read_csv('data/temp_metadata.csv') %>% mutate(sample=paste0('CL',sample,str_to_title(tissue),'RNA'))

#tximport quants
files<- paste0('data/quants/',tissue,'/',list.files(path=paste0('data/quants/',tissue,'/'), pattern='*\\.sf', recursive = T))

#sample names from salmon
snames <- tibble(fnames=files) %>% separate(fnames, into=c('a','b','c','sample','e')) %>% pull(sample)
#fix the ones that are wrong
#liver
#snames[which(!(snames %in% meta$sample))] <- c("CL4LiverRNA", "CL278LiverRNA", "CL764CLiverRNA",  "CL76SALiverRNA")
#heart
#snames[which(!(snames %in% meta$sample))] <- c("CL278HeartRNA", "CL764CHeartRNA",  "CL76SAHeartRNA" )
#brain
#snames[which(!(snames %in% meta$sample))] <- c("CL4BrainRNA", "CL278BrainRNA", "CL764CBrainRNA", "CL76SABrainRNA"  )
#kidney
snames[which(!(snames %in% meta$sample))] <- c("CL13KidneyRNA", "CL14KidneyRNA", "CL15KidneyRNA", "CL17KidneyRNA", 
   "CL19KidneyRNA", "CL20KidneyRNA", "CL21KidneyRNA", "CL22KidneyRNA",
  "CL23KidneyRNA", "CL278KidneyRNA", "CL27KidneyRNA", "CL29KidneyRNA", 
  "CL32KidneyRNA", "CL40KidneyRNA", "CL43KidneyRNA", "CL45KidneyRNA",
  "CL46KidneyRNA", "CL47KidneyRNA", "CL50KidneyRNA", "CL59KidneyRNA", 
  "CL72KidneyRNA", "CL73KidneyRNA", "CL75KidneyRNA", "CL76SAKidneyRNA")

names(files) <- snames
translate<- read_tsv('data/gid_translation.tsv') %>% drop_na()
setnames(translate, old=c('tnames','gnames'), new=c('TXNAME','GENEID'))
txi.salmon <- tximport(files, type = "salmon", tx2gene = translate)


#####
# Normalized counts from DESeq
#####

#make sure meta and counts are in the same order
meta <- meta %>% filter(sample %in% snames)
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample
meta<-meta[colnames(txi.salmon$counts),]

#run deseq
dds<-DESeqDataSetFromTximport(txi.salmon,
                              colData=meta,
                              design = ~ temp)
dds <- DESeq(dds)
rld <- vst(dds)


#a decision has to be made, either
#a)  pick *all* differentially expressed genes, to maximize the 
#variation between temperatures
#or 
#b) pick all genes (filter out low counts and variance) and identify DEG within those networks
#the manual says that using only DEG in the network gen process invalidates the
#scale-free topology assumption, but the results look pretty good so idk how much it matters
#So I checked both ways and the scale free topology assumption is not violated in either 
#direction, and since the DE genes give us a better signal for finding hub genes, that's the direction
#note that if DE genes were determined with only 2 categories / contrasts, the scale free assumption
#likely would be violated, but since we have 6C2 tests we ran, it has much less of an effect

### situation a)
q<-0.01

contrasts <- list(c("SummerActive", "4C"), c("SummerActive", "12C"), c("SummerActive", "20C"),
                  c("SummerActive", "25C"), c("SummerActive", "30C"),
                  c("30C", "4C"), c("30C", "12C"), c("30C", "20C"),c("30C", "25C"), 
                  c("25C", "4C"), c("25C", "12C"), c("25C", "20C"),
                  c("20C", "4C"), c("20C", "12C"),
                  c("12C", "4C"))

sig_genes <- c()
for (i in contrasts){
  print(i)
  res <- results(dds, contrast=c('temp',i))
  sig_genes <- c(sig_genes, 
                 rownames(subset(res, padj<=q)))
  
}
sig_genes <- unique(sig_genes)


### OR


#situation b)
#filter out genes with low variance, mostly just noise
#after trying this, looking at scale-free independence I think
#the original DESeq genes are more useful
'
mat <- t(assay(rld))
vars <- colVars(mat, useNames=T)
sums <- colSums(mat)
v_filtered_genes <- names(vars[which(vars>0.00001)])
c_filtered_genes <- names(sums[sums>150])
t_genes <- colnames(mat)
filtered_genes <- t_genes[(t_genes %in% c_filtered_genes)&(t_genes %in% v_filtered_genes)]
mat <- mat[,filtered_genes]
'

#subset on these genes to get count matrix mat
mat <- t(assay(rld[sig_genes,]))

#####
# Run WGCNA
#####

#good to check no NA, but unlikely at this step
#goodSamplesGenes(mat, verbose=3)

#check no outlier samples
sampleTree <- hclust(dist(mat), method = "average")

#save the sample clusters
png(filename = paste0('figs/',tissue,'_sample_clusters.png'), width=600, height=900, units = 'px')

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h=80, col='red')
dev.off()

#remove outlier if any
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
mat <- mat[keepSamples, ]
nGenes = ncol(mat)
nSamples = nrow(mat)

#create df in same order as mat, save
datTraits <- meta[match(rownames(mat), meta$sample),] %>% select(-c(sample))
nums <- as.numeric(str_sub(datTraits$temp, end=-2))
nums[is.na(nums)] <- 37
datTraits$temp <- nums 
#save(mat, datTraits, file = "data/wgcna_savepoint1.RData")


#set soft-thresholding power for network generation
powers <- c(c(1:10), seq(from = 12, to=70, by=2))
sft <- pickSoftThreshold(mat, powerVector = powers, verbose = 5, networkType='unsigned')

#sft for each tissue:
#heart = 9  (10 @ 0.9; 9 @ 0.8 r^2)  
#liver = 10  ( ; 10 0.8 r^2)
#brain = 12  (20 @ 0.9 r^2 ; 12 @ 0.8 r^2)
#kidney = 14 (had to use 0.8 r^2)

#sft for each tissue fullsend:
#kidney = 54 (54 @ 0.9)

#save thresholding results 
png(filename = paste0('figs/',tissue,'_sft_soft_threshold.png'), width=700, height=700, units = 'px')

cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue");
# this line corresponds to using an R^2 cut-off of .8
abline(h=0.90,col="red")
dev.off()


#####
# Generate networks
#####

softPower <-14
adj <- adjacency(mat, power=softPower)
TOM <- TOMsimilarity(adj, TOMType = 'unsigned')
dissTOM <- 1-TOM

#module detection using dynamic tree cutting
geneTree <- hclust(as.dist(dissTOM), method = "average")
minModuleSize <- 30 
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamStage=T, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

#this function outputs a ton of modules, we have to merge coexpressed modules
dynamicColors <- labels2colors(dynamicMods)
MEList <- moduleEigengenes(mat, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

#based on the threshold below, merge modules by eigengenes
MEDissThres <- 0.25
merged <- mergeCloseModules(mat, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merged$colors
mergedMEs <- merged$newMEs


#plot merge
#sizeGrWindow(12, 9)
pdf(file = paste0('figs/',tissue,'_clusters.pdf'), wi=9, he=12)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#final merge
moduleColors <- mergedColors
colorOrder<-c("grey", standardColors(50));
moduleLabels<-match(moduleColors, colorOrder)-1;
MEs <- mergedMEs

#recalculate eigengenes with merged modules
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(mat, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p", method = 'spearman')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

#look at correlations within the module eigengenes to temp

### it may make more sense to analyze these as factors, but correlation unclear
#sizeGrWindow(5, 7)
pdf(file = paste0('figs/',tissue,'_eigengene_temp_corr.pdf'), wi=5, he=9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.y = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships "))
dev.off()


#relationship between genes and modules
#within important modules
#from above, in heart cyan is negatively correlated with temperature, and brown4 positively
modNames <- substring(names(MEs),3)
geneModuleMembership <- as.data.frame(cor(mat, MEs, use = 'p'))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")
geneTraitSignificance <- as.data.frame(cor(mat, datTraits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(datTraits), sep="")
names(GSPvalue) <- paste("p.GS.", names(datTraits), sep="")


#####
#save progress, or load from here
#####
saveRDS(list(TOM=TOM, mat=mat, MEs=MEs, datTraits=datTraits, dds=dds, meta=meta, vst=rld, moduleColors=moduleColors),
        file = paste0('data/wgcna/',tissue,'_wgcna_objs.rds'))


#explore ...
module <- 'purple'
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for temperature",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col='#32a8a4')#'#32a8a4'

#expression of module genes wrt temp
modgenes <- colnames(mat[,moduleColors==module])
deg <- results(dds, contrast=c('temp','SummerActive','20C')) %>% as.data.frame() %>% filter(padj<0.1, abs(log2FoldChange) > 2) %>% rownames()
exprdat <- as.data.frame(mat[,modgenes]) %>% tibble::rownames_to_column('sample') %>% 
  pivot_longer(!sample, names_to = 'gene', values_to = 'expr', ) %>% 
  left_join(., meta, by='sample') %>% 
  mutate(temp=case_when(temp=='SummerActive'~37,
                        temp=='30C'~30, 
                        temp=='25C'~25, 
                        temp=='20C'~20,
                        temp=='12C'~12,
                        temp=='4C'~4)) %>% 
  filter(gene %in% deg)

ggplot(exprdat %>% group_by(gene))+
  stat_summary(mapping=aes(x=temp, y=expr, group=gene), color=module, geom='line', fun=mean)+
  stat_summary(mapping=aes(x=temp, y=expr, group=gene), color=module, geom='point', fun=mean)


heat <- t(scale(mat[,modgenes]))
#fill color gradient
myCol <- colorRampPalette(c('#ea6bdd', 'black','#50eea7'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

#set up temp gradient
htemp<- meta$temp
htemp.col <- brewer.pal(6, 'Reds')
meta<-meta %>% mutate(color=case_when(temp=='SummerActive'~"#A63603",
                                      temp=='4C'~"#FEEDDE",
                                      temp=='12C'~"#FDD0A2",
                                      temp=='20C'~"#FDAE6B",
                                      temp=='25C'~"#FD8D3C",
                                      temp=='30C'~"#E6550D"))
meta <- meta %>% mutate(temp=factor(meta$temp,
                                    levels=c('SummerActive',
                                             '4C',
                                             '12C',
                                             '20C',
                                             '25C',
                                             '30C')))
mcolors<- list(temp=c('SummerActive'="#A50F15",
                      '4C'="#FEE5D9",
                      '12C'="#FCBBA1",
                      '20C'="#FC9272",
                      '25C'="#FB6A4A",
                      '30C'="#DE2D26"))

#create annotation on heatmap
colAnn <- HeatmapAnnotation(df= meta %>%  
                              filter(sample %in% rownames(mat)) %>% 
                              select(temp) %>% as.data.frame(),
                            which='col',
                            na_col='white',
                            col=mcolors) 
#expression hist
boxplotCol <- HeatmapAnnotation(
  boxplot = anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'left')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'col')


boxplotRow <- HeatmapAnnotation(
  boxplot = row_anno_boxplot(
    heat,
    border = FALSE,
    gp = gpar(fill = '#CCCCCC'),
    pch = '.',
    size = unit(2, 'mm'),
    axis = TRUE,
    axis_param = list(
      gp = gpar(fontsize = 12),
      side = 'top')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'row')

#in this code snippet, we ‘step through’ the rownames and only retain every (nsamp)th successive label
nsamp <- 2
genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = seq(1, nrow(heat), nsamp),
    labels = rownames(heat)[seq(1, nrow(heat), nsamp)],
    labels_gp = gpar(fontsize = 6, fontface = 'bold'),
    padding = 0.75),
  width = unit(2.0, 'cm') +
    
    max_text_width(
      rownames(heat)[seq(1, nrow(heat), nsamp)],
      gp = gpar(fontsize = 6,  fontface = 'bold')))

#do a heirarchical cluster
pamClusters <- cluster::pam(heat, k = 3) # pre-select k = 3 centers
pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 2', 'Cluster 1', 'Cluster 3', 'Cluster 4'))

#get custom order of columns
mmeta <- meta
rownames(mmeta) <- c()
sorder <- mmeta %>% mutate(temp=factor(meta$temp,
                                       levels=c('SummerActive',
                                                '4C',
                                                '12C',
                                                '20C',
                                                '25C',
                                                '30C')))  %>% 
  arrange(temp) %>% pull(sample)

hmap <- Heatmap(heat,
                
                # split the genes / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = FALSE,
                
                name = 'Gene\nZ-\nscore',
                
                col = colorRamp2(myBreaks, myCol),
                
                # parameters for the colour-bar that represents gradient of expression
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                #row_title = 'Statistically significant genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                column_order = sorder[sorder %in% colnames(heat)],
                show_column_dend = F,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = T,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2',
                
                # specify top and bottom annotations
                top_annotation = colAnn,
                bottom_annotation = boxplotCol)
draw(hmap + genelabels,
     heatmap_legend_side = 'left',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')

#write out
out <- data.frame(gene=colnames(mat), module=moduleColors, gsig_temp=geneTraitSignificance)
write_csv(out, file=paste0('data/',tissue,'_gene_modules.csv') )
