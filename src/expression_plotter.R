library(tidyverse)
library(DESeq2)
library(stringr)

gene <- 'Eif4a'


heart_data <- readRDS('data/wgcna/heart_wgcna_objs.rds')$vst
brain_data <- readRDS('data/wgcna/brain_wgcna_objs.rds')$vst
liver_data <- readRDS('data/wgcna/liver_wgcna_objs.rds')$vst
kidney_data <- readRDS('data/wgcna/kidney_wgcna_objs.rds')$vst

hdat <- as.data.frame(t(assay(heart_data[gene,],))) %>% mutate(temp=heart_data$temp) %>% mutate(tissue='heart')
bdat <- as.data.frame(t(assay(brain_data[gene,],))) %>% mutate(temp=brain_data$temp) %>% mutate(tissue='brain')
ldat <- as.data.frame(t(assay(liver_data[gene,],))) %>% mutate(temp=liver_data$temp) %>% mutate(tissue='liver')
kdat <- as.data.frame(t(assay(kidney_data[gene,],))) %>% mutate(temp=kidney_data$temp) %>% mutate(tissue='kidney')

dat <- rbind(hdat, bdat, ldat, kdat) %>% 
  mutate(temp=fct_relevel(temp, c('4C','12C','20C','25C','30C','SummerActive'))) %>% 
  mutate(tissue=str_to_title(tissue))

ggplot(dat)+
  geom_boxplot(mapping=aes(x=temp,y=get(gene), fill=tissue), outlier.alpha = 0)+
  geom_jitter(mapping=aes(x=temp,y=get(gene)), width=0.25, height=0.01, size=0.75)+
  facet_wrap(~tissue)+
  theme_bw()+
  labs(x='Body Temperature', y=paste0(gene,' Expression (VST)'))+
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  guides(fill='none')

ggsave(paste0('figs/gene_expr/', gene, '_expression_4tiss.png' ))
