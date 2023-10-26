library(tidyverse)
library(fastqcr)

#package sucks and needs to be run from the folder with data
setwd('dylan/data/fastqc/zips/')
qc <- qc_aggregate(qc.dir = getwd())

summary(qc)

#plot total fails warnings etc
gdata <- summary(qc) %>% select(-nb_samples) %>% pivot_longer(!c(module, warned, failed),names_to='message',values_to = 'N') %>% group_by(message)
ggplot(gdata)+
  geom_bar(mapping=aes(x=module, y=N,fill=message),stat='identity', pos='dodge')+
  theme_bw()+
  scale_fill_discrete(labels=c('Fail','Pass','Warn'))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        axis.title.x = element_blank())
ggsave('../../../figs/fastqc_agg.png', width=3000, height=3000, dpi=300, units='px')

#select module, output sample names and over represented seqs
mod <- 'Per sequence GC content'
failed <- summary(qc) %>% filter(module==mod) %>% pull(failed) %>% str_split(', ')
warned <- summary(qc) %>% filter(module==mod) %>% pull(warned) %>% str_split(', ')

#sample names
for (fn in c('fail','warn')){
  fname <- file(paste0('../../../data/fastqc/gc_',fn,'.txt'))
  writeLines(get(eval(paste0(fn,'ed')))[[1]], fname)
  close(fname)
}

#OR Seqs into a fasta
#stupid package can't extract them, have to do myself
samples <- str_split(warned[[1]], '\\.',n = 2) %>% sapply('[[',1)
fasta_content <- c()
for (s in samples){
  lines <- readLines(unz(paste0('../zips/',s,'_fastqc.zip'), paste0(s,'_fastqc/fastqc_data.txt')))
  imin <- FALSE
  ors <- c()
  for (line in lines){
    if ((line=='>>END_MODULE') & (imin)){
      break
    }
    if (imin){
      ors <- c(ors, str_split(line,'\t')[[1]][1])
    }
    if (line=='#Sequence	Count	Percentage	Possible Source'){
      imin <- TRUE
    }
    
    
  }
  if (length(ors)>0){
    for (j in 1:length(ors)){
      fasta_content <- c(fasta_content, paste0('>',s,'_seq',j), ors[j])
    }
  }
  
  }
writeLines(fasta_content,'../gc_warned_orseqs.fasta')
