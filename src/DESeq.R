library(tidyverse)
library(data.table)
library(DESeq2)
library(tximport)
library(RColorBrewer)
library(cluster)
library(circlize)
library(ComplexHeatmap)
library(rtracklayer)


#####
#read in metadata, tximport quants
#####
tissue <- 'heart'
meta <- read_csv('data/temp_metadata.csv') %>% mutate(sample=paste0('CL',sample,str_to_title(tissue),'RNA'))

#tximport quants
files<- paste0('data/quants/',tissue,'/',list.files(path=paste0('data/quants/',tissue,'/'), pattern='*\\.sf', recursive = T))

#sample names from salmon
snames <- tibble(fnames=files) %>% separate(fnames, into=c('a','b','c','sample','e')) %>% pull(sample)


#fix the ones that are wrong
#liver
#snames[which(!(snames %in% meta$sample))] <- c("CL4LiverRNA", "CL278LiverRNA", "CL764CLiverRNA",  "CL76SALiverRNA")
#heart
snames[which(!(snames %in% meta$sample))] <- c("CL278HeartRNA", "CL764CHeartRNA",  "CL76SAHeartRNA" )
#brain
#snames[which(!(snames %in% meta$sample))] <- c("CL4BrainRNA", "CL278BrainRNA", "CL764CBrainRNA", "CL76SABrainRNA"  )
#kidney
#snames[which(!(snames %in% meta$sample))] <- c("CL13KidneyRNA", "CL14KidneyRNA", "CL15KidneyRNA", "CL17KidneyRNA", 
#   "CL19KidneyRNA", "CL20KidneyRNA", "CL21KidneyRNA", "CL22KidneyRNA",
#  "CL23KidneyRNA", "CL278KidneyRNA", "CL27KidneyRNA", "CL29KidneyRNA", 
#  "CL32KidneyRNA", "CL40KidneyRNA", "CL43KidneyRNA", "CL45KidneyRNA",
#  "CL46KidneyRNA", "CL47KidneyRNA", "CL50KidneyRNA", "CL59KidneyRNA", 
#  "CL72KidneyRNA", "CL73KidneyRNA", "CL75KidneyRNA", "CL76SAKidneyRNA")

names(files) <- snames
translate<- read_tsv('data/gid_translation.tsv') %>% drop_na()
setnames(translate, old=c('tnames','gnames'), new=c('TXNAME','GENEID'))
txi.salmon <- tximport(files, type = "salmon", tx2gene = translate)

#make sure meta and counts are in the same order
meta <- meta %>% filter(sample %in% snames)
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample
meta<-meta[colnames(txi.salmon$counts),]



#####
#run DESeq2
#####

dds<-DESeqDataSetFromTximport(txi.salmon,
                              colData=meta,
                              design = ~ temp)

dds <- DESeq(dds)
rld <- vst(dds)

q<-0.01

contrasts <- list(c("SummerActive", "4C"), c("SummerActive", "12C"), c("SummerActive", "20C"),
              c("SummerActive", "25C"), c("SummerActive", "30C"),
              c("30C", "4C"), c("30C", "12C"), c("30C", "20C"),c("30C", "25C"), 
              c("25C", "4C"), c("25C", "12C"), c("25C", "20C"),
              c("20C", "4C"), c("20C", "12C"),
              c("12C", "4C"))


#####
#generate heatmaps
#####
  
sig_genes <- c()
rout <- c()
for (i in contrasts){
  print(i)
  res <- results(dds, contrast=c('temp',i))
  sig_genes <- c(sig_genes, 
                 rownames(subset(res, padj<=q)))
  res$contrast <- paste(i, collapse='_')
  res$gene <- rownames(res)
  rownames(res) <- NULL
  if (length(rout)==0){
    rout <- res
  }
  else{
    rout <- rbind(rout, res)
  }
}
write_csv(as.data.frame(rout), file=paste0('data/', tissue, '_DESeq.csv') )
sig_genes <- unique(sig_genes)
mat <- assay(rld[sig_genes,])
heat <- t(scale(t(mat)))


#fill color gradient
myCol <- colorRampPalette(c('#50eea7', 'black', '#ea6bdd'))(100)
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
                              select(-c(sample, color)),
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
nsamp <- 80
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
                                 levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'))

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
                column_order = sorder,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = FALSE,
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
