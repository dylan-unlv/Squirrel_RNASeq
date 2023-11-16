library(DESeq2)
library(cluster)
library(tximport)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)


##########various clustering things

#import data
dds <- readRDS('data/deseq/all_tiss.dds.rds')
meta <- read.csv('data/deseq/all_tiss_meta.csv')
counts <- counts(dds, normalized=T)
vst <- assay(vst(dds))

#select n genes
topn <- 500
rvars <- rowVars(counts, useNames = T) 
ngenes <- data.frame('gene'=names(rvars), 'vars'=rvars) %>% 
  mutate(idx=row_number()) %>% 
  arrange(vars) %>% head(n=topn) %>% 
  pull(gene)



#select grouping variable (within tissue for instance)
group <- 'tissue' 
groups <- meta %>% pull(get(group)) %>% unique()

for (g in groups){
  
  ###Filter
  samps <- meta %>% filter(get(group)==g) %>% pull(sample) #only samples in group
  
  topn <- 500 #top n by variance (or by DEG eventually) (hypothesis neutral vs hypothesis driven)
  rvars <- rowVars(vst[,which(colnames(vst) %in% samps)], useNames = T) 
  ngenes <- data.frame('gene'=names(rvars), 'vars'=rvars) %>% 
    mutate(idx=row_number()) %>% 
    arrange(-vars) %>% head(n=topn) %>% 
    pull(gene)
  #counts_i <- counts[ngenes,which(colnames(vst) %in% samps)]
  vst_i <- vst[ngenes,which(colnames(vst) %in% samps)]
  vst_i <- vst_i[!is.na(rowVars(vst_i, useNames = T)),]
  vst_i <- vst_i[rowVars(vst_i, useNames = T)>0,]
  
  ###Cluster
  pamClusters <- cluster::pam(vst_i, k = 3) # pre-select k = 9 centers
  pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
  pamClusters$clustering <- factor(pamClusters$clustering,
                                   levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 
                                              'Cluster 5', 'Cluster 6', 'Cluster 7', 'Cluster 8', 'Cluster 9'))

  ###Graph Expr
  
  #annotations
  mcolors<- list(temp=c('SummerActive'="#A50F15",
                        '4C'="#FEE5D9",
                        '12C'="#FCBBA1",
                        '20C'="#FC9272",
                        '25C'="#FB6A4A",
                        '30C'="#DE2D26"))
  colAnn <- HeatmapAnnotation(df= meta %>%  
                                filter(sample %in% colnames(vst_i)) %>% 
                                select(temp) %>% as.data.frame(),
                              which='col',
                              na_col='white',
                              col=mcolors, simple_anno_size = unit(3,'cm'),
                              annotation_name_gp= gpar(fontsize = 35)) 
  boxplotCol <- HeatmapAnnotation(
    boxplot = anno_boxplot(
      vst_i,
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
      vst_i,
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
  
  nsamp <- 2
  genelabels <- rowAnnotation(
    Genes = anno_mark(
      at = seq(1, nrow(vst_i), nsamp),
      labels = rownames(vst_i)[seq(1, nrow(vst_i), nsamp)],
      labels_gp = gpar(fontsize = 6, fontface = 'bold'),
      padding = 0.75),
    width = unit(2.0, 'cm') +
      
      max_text_width(
        rownames(vst_i)[seq(1, nrow(vst_i), nsamp)],
        gp = gpar(fontsize = 6,  fontface = 'bold')))
  
  myCol <- colorRampPalette(c('#ea6bdd', 'black','#50eea7'))(100)
  myBreaks <- seq(-2.5, 2.5, length.out = 100)
  
  #order samples by temp
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
  
  
  
  #gen heatmap
  png(paste0('figs/clusters/',g,'_top_',topn,'_var_pam_heatmap.png'), width = 2000, height = 3000, units = 'px')
  hmap <- Heatmap(t(scale(t(vst_i))),
                  
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
                    title_gp=gpar(fontsize = 22, fontface = 'bold'),
                    labels_gp=gpar(fontsize = 22, fontface = 'bold')),
                  
                  # row (gene) parameters
                  cluster_rows = TRUE,
                  show_row_dend = TRUE,
                  #row_title = 'Statistically significant genes',
                  row_title_side = 'left',
                  row_title_gp = gpar(fontsize = 22,  fontface = 'bold'),
                  row_title_rot = 90,
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 20, fontface = 'bold'),
                  row_names_side = 'left',
                  row_dend_width = unit(75,'mm'),
                  
                  # column (sample) parameters
                  cluster_columns = FALSE,
                  column_order = sorder[sorder %in% colnames(vst_i)],
                  show_column_dend = F,
                  column_title = '',
                  column_title_side = 'bottom',
                  column_title_gp = gpar(fontsize = 22, fontface = 'bold'),
                  column_title_rot = 0,
                  show_column_names = T,
                  column_names_gp = gpar(fontsize = 20, fontface = 'bold'),
                  column_names_max_height = unit(20, 'cm'),
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
    dev.off()
    write.csv(counts_i[!is.na(rowVars(counts_i)),], paste0('data/clusters/',g,'_top_',topn,'_var_pam_counts.csv'), row.names = T)
  }
 
