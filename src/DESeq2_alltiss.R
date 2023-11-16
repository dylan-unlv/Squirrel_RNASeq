library(DESeq2)
library(tximport)
library(RColorBrewer)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)


#####
#read in metadata, tximport quants
#####

#need this for proper gene names
translate<- read_tsv('data/gid_translation.tsv') %>% drop_na()
setnames(translate, old=c('tnames','gnames'), new=c('TXNAME','GENEID'))

#tissue <- 'kidney'
tissues <- c('heart','brain','liver','kidney')
lsnames <- c()
lsfiles <- c()
lsmeta <- c()
for (tissue in tissues){
  meta <- read_csv('data/temp_metadata.csv') %>% mutate(sample=paste0('CL',sample,str_to_title(tissue),'RNA')) %>% mutate(tissue=tissue)

  #tximport quants
  files<- paste0('data/quants/',tissue,'/',list.files(path=paste0('data/quants/',tissue,'/'), pattern='*\\.sf', recursive = T))

  #sample names from salmon
  snames <- tibble(fnames=files) %>% separate(fnames, into=c('a','b','c','sample','e')) %>% pull(sample)
  #fix the ones that are wrong
  if (tissue=='liver'){snames[which(!(snames %in% meta$sample))] <- c("CL4LiverRNA", "CL278LiverRNA", "CL764CLiverRNA",  "CL76SALiverRNA")
  }else if (tissue=='heart'){snames[which(!(snames %in% meta$sample))] <- c("CL278HeartRNA", "CL764CHeartRNA",  "CL76SAHeartRNA" )
  }else if (tissue=='brain'){snames[which(!(snames %in% meta$sample))] <- c("CL4BrainRNA", "CL278BrainRNA", "CL764CBrainRNA", "CL76SABrainRNA"  )
  }else if (tissue=='kidney'){
      snames[which(!(snames %in% meta$sample))] <- c("CL13KidneyRNA", "CL14KidneyRNA", "CL15KidneyRNA", "CL17KidneyRNA", 
                                               "CL19KidneyRNA", "CL20KidneyRNA", "CL21KidneyRNA", "CL22KidneyRNA",
                                               "CL23KidneyRNA", "CL278KidneyRNA", "CL27KidneyRNA", "CL29KidneyRNA", 
                                               "CL32KidneyRNA", "CL40KidneyRNA", "CL43KidneyRNA", "CL45KidneyRNA",
                                               "CL46KidneyRNA", "CL47KidneyRNA", "CL50KidneyRNA", "CL59KidneyRNA", 
                                               "CL72KidneyRNA", "CL73KidneyRNA", "CL75KidneyRNA", "CL76SAKidneyRNA")
  }
  
  names(files) <- snames

  
  lsfiles <- c(lsfiles, files)
  lsnames <- c(lsnames, snames)
  lsmeta <- rbind(lsmeta, meta)
}


#import all files as one object
txi.salmon <- tximport(lsfiles, type = "salmon", tx2gene = translate)

#####
# Normalized counts from DESeq
#####

#make sure meta and counts are in the same order
lsmeta <- lsmeta %>% filter(sample %in% lsnames)
lsmeta <- as.data.frame(lsmeta)
rownames(lsmeta) <- lsmeta$sample
lsmeta<-lsmeta[colnames(txi.salmon$counts),]

#init deseq
dds<-DESeqDataSetFromTximport(txi.salmon,
                              colData=lsmeta,
                              design = ~ tissue:temp)
dds$group <- factor(paste0(dds$tissue,'.', dds$temp))
design(dds) <-  ~ group
#filter out by genes with low counts and outlier samples
dds <- estimateSizeFactors(dds)
mincount <- 5
gidx <- rowSums( counts(dds, normalized=TRUE) >= mincount ) >= 3

#cluster analysis of samples to identify bad ones
bad_samples <- c('CL59KidneyRNA','CL51HeartRNA','CL278BrainRNA','CL50LiverRNA')
good_samples <- setdiff(dds$sample, bad_samples)

#subset deseq object
dds <- dds[gidx, good_samples]

#run deseq
dds <- DESeq(dds)
vst <- vst(dds)
counts <- counts(dds)

saveRDS(dds, file='data/deseq/all_tiss.dds.rds')
saveRDS(vst, file='data/deseq/all_tiss.vst.rds')
write.csv(data.frame(counts), file = 'data/deseq/all_tiss_counts.csv', row.names = T)
write_csv(lsmeta, file= 'data/deseq/all_tiss_meta.csv')
