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

#since we're only interested here in the overlap
#we need to filter each of these TOM so they're identical dimensions
#so we start by assembling a list of genes present in all networks

tgenes <- Reduce(intersect, 
                 list(colnames((heart_data$mat)), 
                      colnames(brain_data$mat), 
                      colnames(kidney_data$mat), 
                      colnames(liver_data$mat)))

#then we reduce each TOM and add it to a 3D TOM
nSets <- 4
TOM <- array(0, dim=c(nSets, length(tgenes), length(tgenes)))
TOM[1,,] = heart_data$TOM[which(colnames(heart_data$mat) %in% tgenes),which(colnames(heart_data$mat) %in% tgenes)]
TOM[2,,] = brain_data$TOM[which(colnames(brain_data$mat) %in% tgenes),which(colnames(brain_data$mat) %in% tgenes)]
TOM[3,,] = liver_data$TOM[which(colnames(liver_data$mat) %in% tgenes),which(colnames(liver_data$mat) %in% tgenes)]
TOM[4,,] = kidney_data$TOM[which(colnames(kidney_data$mat) %in% tgenes),which(colnames(kidney_data$mat) %in% tgenes)]



nGenes <- length(tgenes)
setLabels <- c('heart','brain','liver','kidney')
# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets first to calculate scale factor, then to scale
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8)
}
# Scale the TOM based on highest quantile, works better than normalize.quantiles
mset <- which(scaleQuant==max(scaleQuant))
for (set in 1:nSets)
{
  if (set!=mset)
    {
      scalePowers[set] = log(scaleQuant[mset])/log(scaleQuant[set])
      TOM[set, ,] = TOM[set, ,]^scalePowers[set]
    }
}



##QQ plots to ensure even distribution of network edge weights
scaledTOMSamples = list()
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set] 
sizeGrWindow(6,6)

#list of all tissue comparisons
comps <- expand.grid(1:4, 1:4) %>% filter(Var1!=Var2)

#iterate a qq plot for each comparison
for (i in 1:nrow(comps)){
  s1 <- comps[i,'Var1']
  s2 <- comps[i,'Var2']
  pdf(file = paste0("figs/TOMScaling-QQPlot_",setLabels[s1],"_",setLabels[s2],".pdf"), wi = 6, he = 6);
  # qq plot of the unscaled samples
  qqUnscaled = qqplot(TOMScalingSamples[[s1]], TOMScalingSamples[[s2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[s1]), ylab = paste("TOM in", setLabels[s2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
  qqScaled = qqplot(scaledTOMSamples[[s1]], scaledTOMSamples[[s2]], plot.it = FALSE)
  points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
  abline(a=0, b=1, col = "blue")
  legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
  dev.off()
}


#####
# generate final consensus network
#####


consensusTOM <- pmin(TOM[1,,], TOM[2,,], TOM[3,,], TOM[4,,])
# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average")
# At this point, smaller modules are acceptable. Usually mms=30, here we use 10
minModuleSize = 10
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE )
unmergedColors = labels2colors(unmergedLabels)
sizeGrWindow(8,6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Calculate module eigengenes
multiExpr <- vector(mode='list',length=nSets)
multiExpr[[1]] <- list(data=heart_data$mat[,tgenes])
multiExpr[[2]] <- list(data=brain_data$mat[,tgenes])
multiExpr[[3]] <- list(data=liver_data$mat[,tgenes])
multiExpr[[4]] <- list(data=kidney_data$mat[,tgenes])

unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
pdf(file='figs/wgcna_consensus_module_relationships.pdf', width=9, height=6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",xlab = "", sub = "")
abline(h=0.05, col = "red")
dev.off()

merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.05, verbose = 3)
# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs
pdf(file='figs/wgcna_consensus_module_dendrogram.pdf', width=9, height=6)
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

saveRDS(list(consMEs=consMEs, moduleColors=moduleColors, setLabels=setLabels,
             moduleLabels=moduleLabels, consTree=consTree,
             consensusTOM=consensusTOM, multiExpr=multiExpr), file = 'data/wgcna/consensus_wgcna_objs.rds')
cyt <- exportNetworkToCytoscape(consensusTOM,
                                edgeFile = 'data/cytoscape/consensus.edges',
                                nodeFile = 'data/cytoscape/consensus.nodes',
                                weighted = T,
                                threshold = 0.02,
                                nodeNames = tgenes,
                                nodeAttr = moduleColors)
