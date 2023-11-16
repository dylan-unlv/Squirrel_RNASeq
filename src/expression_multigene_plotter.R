library(tidyverse)
library(DESeq2)
library(stringr)


heart_data <- readRDS('data/wgcna/heart_wgcna_objs.rds')$vst
brain_data <- readRDS('data/wgcna/brain_wgcna_objs.rds')$vst
liver_data <- readRDS('data/wgcna/liver_wgcna_objs.rds')$vst
kidney_data <- readRDS('data/wgcna/kidney_wgcna_objs.rds')$vst

#####
# Some function that chooses a set of geness we have quantified

genes <- rownames(heart_data)[grep('Eif4', rownames(heart_data))]




hdat <- as.data.frame(t(assay(heart_data[genes,],))) %>% mutate(temp=heart_data$temp) %>% mutate(tissue='heart')
bdat <- as.data.frame(t(assay(brain_data[genes,],))) %>% mutate(temp=brain_data$temp) %>% mutate(tissue='brain')
ldat <- as.data.frame(t(assay(liver_data[genes,],))) %>% mutate(temp=liver_data$temp) %>% mutate(tissue='liver')
kdat <- as.data.frame(t(assay(kidney_data[genes,],))) %>% mutate(temp=kidney_data$temp) %>% mutate(tissue='kidney')

dat <- rbind(hdat, bdat, ldat, kdat) %>% 
  mutate(temp=fct_relevel(temp, c('4C','12C','20C','25C','30C','SummerActive'))) %>% 
  mutate(tissue=str_to_title(tissue))

#need to reshape data
dat <- dat %>% 
  pivot_longer(-one_of(c('temp','tissue')), names_to = 'gene', values_to = 'expr') %>% 
  mutate(temp=case_when(temp=='SummerActive'~37,
                        temp=='30C'~30, 
                        temp=='25C'~25, 
                        temp=='20C'~20,
                        temp=='12C'~12,
                        temp=='4C'~4))


ggplot(dat, mapping=aes(x=temp, y=expr,group=gene,color=tissue, shape=gene),)+
  stat_summary(geom='line', fun=mean)+
  stat_summary(geom='point', fun=mean, size=1)+
  scale_color_manual(values =  c('#e54ef5','#f0435a', '#b5c95b', '#729fc2'))+
  scale_shape_manual(values=0:20)+
  scale_x_discrete(name='Body Temperature',limits=c(4,12,20,25,30,37))+
  facet_wrap(~tissue)+
  theme_bw()+
  labs(x='Body Temperature', y=paste0(gene,' Expression (VST)'))+
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  guides(color='none')

ggsave(paste0('figs/gene_expr/', paste(genes, collapse='_'), '_expression_4tiss.png' ))
