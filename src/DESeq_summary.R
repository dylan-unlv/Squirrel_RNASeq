library(tidyverse)
library(magrittr)

read_more <- function(x){
  return(read_csv(paste0('data/',x)) %>% mutate(fn=x))
}



brain <- read_csv('data/brain_DESeq.csv')

dat <- list.files(path='data',pattern='_DESeq\\.csv') %>% 
  map_df(., read_more)

dat <- dat  %>% 
  separate(fn, into=c('tissue','trash','trash1')) %>% 
  select(c(-trash, -trash1)) %>% group_by(tissue)


#nonzero counts genes in each tissue
dat  %>% filter(baseMean>0) %>% 
  select(gene, tissue) %>% 
  mutate(gene=as.factor(gene)) %>% 
  summarise(n=n_distinct(gene))

#de genes in each tissue
dat  %>% filter(padj<0.01) %>% 
  select(gene, tissue) %>% 
  mutate(gene=as.factor(gene)) %>% 
  summarise(n=n_distinct(gene))

#maxlfc each tissue
dat  %>% filter(padj<0.01) %>% 
  select(gene, tissue, log2FoldChange) %>%
  mutate(gene=as.factor(gene)) %>% 
  summarise(m=max(log2FoldChange))