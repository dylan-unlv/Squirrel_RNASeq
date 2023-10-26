library(tidyverse)

data <- read_delim('data/Result_seal_formatted.csv', delim='\t')

genes <- data %>% pull(SimpleGeneName) %>% unique()
length(genes)

#copy number table
data %>% group_by(SimpleGeneName) %>% mutate(copy=row_number()) %>% filter(copy==max(copy)) %>% ungroup() %>% count(copy)
#data <- data %>% group_by(SimpleGeneName) %>% mutate(copy=row_number()) %>% filter(copy==max(copy)) 
#2 copy number genes
c2_genes <- data %>% group_by(SimpleGeneName) %>% mutate(copy=row_number())  %>% filter(copy==max(copy)) %>% filter(copy==2) %>% pull(SimpleGeneName)

#find proportion of total length assuming no overlap btwn scaffolds
gdata <- data %>% filter(SimpleGeneName %in% c2_genes) %>% 
    group_by(SimpleGeneName) %>% 
    mutate(gene_length=sum(Length)) %>% 
    mutate(proportion = Length/gene_length) %>% 
    arrange(-proportion) %>% 
    mutate(copy=row_number()) %>% 
    filter(copy==min(copy)) 

#max proportion distributions
ggplot(gdata)+
  geom_density(mapping=aes(x=proportion), fill='#f7c081', color='#f7c081', adjust=0.5) +
  theme_bw()

#