library(tidyverse)
library(DESeq2)
library(stringr)
library(scales)
library(car)

#####
#this script outputs a table for each tissue of the most variable genes
#at a given temperature
#meant to be used with expression_plotter.R to viz results

#Load in data
heart_data <- readRDS('data/wgcna/heart_wgcna_objs.rds')$vst
brain_data <- readRDS('data/wgcna/brain_wgcna_objs.rds')$vst
liver_data <- readRDS('data/wgcna/liver_wgcna_objs.rds')$vst
kidney_data <- readRDS('data/wgcna/kidney_wgcna_objs.rds')$vst

hdat <- as.data.frame(t(assay(heart_data))) %>% mutate(temp=heart_data$temp) %>% mutate(tissue='heart')
bdat <- as.data.frame(t(assay(brain_data))) %>% mutate(temp=brain_data$temp) %>% mutate(tissue='brain')
ldat <- as.data.frame(t(assay(liver_data))) %>% mutate(temp=liver_data$temp) %>% mutate(tissue='liver')
kdat <- as.data.frame(t(assay(kidney_data))) %>% mutate(temp=kidney_data$temp) %>% mutate(tissue='kidney')

dat <- rbind(hdat, bdat, ldat, kdat) 

temps <- c('4C','12C','20C','25C','30C','SummerActive')
tissues <- c('heart', 'brain', 'liver', 'kidney')

#add variance within each group
vdat <- dat %>% 
  pivot_longer(!c(temp,tissue), names_to='gene', values_to='expr') %>% 
  group_by(tissue, temp, gene) %>% 
  summarise(variance=var(expr))

#arrange by high variance genes, export entire table
vdat  %>% arrange(-variance) %>% write_csv('data/variance_table.csv')

#select top 0.5% variance genes in each tissue/temp combo
#export table, graph cutoff
for (tiss in tissues){
  for (ttemp in temps){
    tempdat <- vdat %>% filter(tissue==tiss, temp==ttemp, variance>0)
    tmean <- mean(log10(tempdat$variance))
    tsd <- sd(log10(tempdat$variance))
    val <- qnorm(0.995, tmean, tsd)
    tempdat %>% 
      filter(variance > (10**val)) %>% 
      arrange(-variance) %>% 
      write_csv(paste0('data/',tiss,'_',ttemp,'_top_05_variance.csv'))
    
    gg<- ggplot()+
      geom_density(data=tempdat, mapping=aes(x=log10(variance)))+
      geom_vline(xintercept = val, color='red')+
      ggtitle(paste0(str_to_title(tiss), ' ', ttemp, ' Log10(Variance) Density' ))+
      theme_bw()
    ggsave(paste0('figs/variances/',tiss,'_',ttemp,'_top_05_variance.png'))
  }
}


#####
#Levene's Test

genes<-vdat$gene %>% unique()


#I don't feel like doing this concisely, it's going down manual mode
#manual mode did not work, it ran for 2 days and got 21248/133704 results -_-*

gs <- c()
ts <- c()
ps <- c()
fs <- c()
v4 <- c()
v12 <- c()
v20 <- c()
v25 <- c()
v30 <- c()
vSA <- c()

for (tgene in genes){
  for (tiss in tissues){
    tempdat <- dat %>% 
       pivot_longer(!c(temp,tissue), names_to='gene', values_to='expr') %>% 
       group_by(tissue, temp, gene) %>%
       filter(gene==tgene, tissue==tiss)
    res<-leveneTest(expr~temp,data=tempdat)
    gs <- c(gs, tgene)
    ts <- c(ts, tiss)
    ps <- c(ps, res$`Pr(>F)`[[1]])
    fs <- c(fs, res$`F value`[[1]])
    v4 <- c(v4, tempdat %>% filter(temp=='4C') %>% pull(expr) %>% var())
    v12 <- c(v12, tempdat %>% filter(temp=='12C') %>% pull(expr) %>% var())
    v20 <- c(v20, tempdat %>% filter(temp=='20C') %>% pull(expr) %>% var())
    v25 <- c(v25, tempdat %>% filter(temp=='25C') %>% pull(expr) %>% var())
    v30 <- c(v30, tempdat %>% filter(temp=='30C') %>% pull(expr) %>% var())
    vSA <- c(vSA, tempdat %>% filter(temp=='SummerActive') %>% pull(expr) %>% var())
}
}


