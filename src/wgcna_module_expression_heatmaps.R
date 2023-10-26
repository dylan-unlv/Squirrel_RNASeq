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

data <- kidney_data

tissue <- 'kidney'
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
                              select(-c(sample, color)),
                            which='col',
                            na_col='white',
                            col=mcolors,
                            simple_anno_size = unit(3, 'cm'), 
                            annotation_label = 'Temperature',
                            annotation_legend_param = list(legend_height=50)) 
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
      gp = gpar(fontsize = 24),
      side = 'top')),
  annotation_width = unit(c(2.0), 'cm'),
  which = 'row')

#in this code snippet, we ‘step through’ the rownames and only retain every (nsamp)th successive label
nsamp <- 8
genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = seq(1, nrow(heat), nsamp),
    labels = rownames(heat)[seq(1, nrow(heat), nsamp)],
    labels_gp = gpar(fontsize = 24, fontface = 'bold'),
    padding = 0.75),
  width = unit(4.0, 'cm') +
    
    max_text_width(
      rownames(heat)[seq(1, nrow(heat), nsamp)],
      gp = gpar(fontsize = 20,  fontface = 'bold')))

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
                  legend_width = unit(20, 'cm'),
                  legend_height = unit(40, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 24, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 24, fontface = 'bold'), plot=T),
                
                # row (gene) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                #row_title = 'Statistically significant genes',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 24,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 20, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(5,'cm'),
                
                # column (sample) parameters
                cluster_columns = FALSE,
                column_order = sorder[sorder %in% colnames(heat)],
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 24, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 20, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2',
                
                # specify top and bottom annotations
                top_annotation = colAnn,
                bottom_annotation = boxplotCol
                
                #legend params
                )

png(filename = paste0('figs/module_heatmaps/',tissue,'_',module,'_expression_heatmap.png'),
    width = 2000, height=3000, units='px')

draw(hmap + genelabels,
     heatmap_legend_side = 'left',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')
dev.off()