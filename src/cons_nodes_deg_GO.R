library(tidyverse)
library(poweRlaw)
library(biomaRt)

data <- read_csv('data/cytoscape/tri_cons_node_analyzed.csv')

thresh <- 0.05
topn <- round(nrow(data)*thresh)
mindeg <- data %>% arrange(-Degree) %>% dplyr::select(name, module, Degree) %>%  head(n=topn) %>% pull(Degree) %>% min()

data <- data %>% mutate(sig=case_when(Degree>=mindeg~'sig', Degree<mindeg~'insig'))

ggplot(data)+
  geom_histogram(mapping=aes(x=Degree, fill=sig),binwidth = 1)+
  scale_fill_manual(values=c('grey','red'))+
  theme_bw()+
  labs(x='Degree',y='Count')+
  guides(fill='none')


###
#GO terms
###
mart = useEnsembl('genes')
ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")
go_terms <- getBM(
  mart=ensembl, attributes=c("hgnc_symbol", "uniprot_gn_id", "uniprot_gn_symbol", "go_id", "namespace_1003", "name_1006"),
  filters="hgnc_symbol", values=data$name)
colnames(go_terms) <- c("name" , "uniprot_gn_id" , "uniprot_gn_symbol", "go_id", "go_cat", "go_desc")
go_terms <- go_terms %>% mutate(name = str_to_title(name))

go_out <- left_join(data %>% filter(sig=='sig'), go_terms %>% 
            dplyr::select(name, go_desc, go_cat, go_id), by='name') %>% 
  dplyr::select(name, go_cat, go_id, go_desc,Degree, module) %>% 
  group_by(name,Degree, module,  go_cat) %>% 
  summarise(go_desc=str_c(go_desc, collapse = ', '), go_id=str_c(go_id, collapse=', ')) %>% 
  filter(go_cat!='') %>% 
  arrange(-Degree)

write_csv(go_out, 'data/tri_consensus_sig_deg_GO.csv')

#data.frame(v=str_split('GO:0006357, GO:0006355, GO:0007283, GO:0006357, GO:0006355, GO:0007283', pattern = ', ')[[1]]) %>% print(row.names=F)












