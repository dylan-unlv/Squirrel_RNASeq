library(tidyverse)
library(DESeq2)
#this script runs simulations to determine which combination of DESeq output
#is the best at recapturing the true distribution of counts in a better scaffolded genome


gene_size <- 100
a <- 0.7
b <- 0.3
nreads <- 250


draw <- function(fc, reads, a, n){
  ret <- c()
  treads <- c()
  for (i in 1:n){
    #transform expected reads to a sample from a normal distribution like with human data
    reads <- as.integer(rnorm(1, mean=reads, sd=sqrt(reads)))
    fcr <- as.integer(fc*reads)
    ret <- c(ret,length(which(runif(fcr, min=0, max=gene_size) <= a*gene_size)))
    treads <- c(treads, fcr)
  }
  return(list(ret, treads, treads-ret))
}

#create simulated count matrix
n<-10
fc_range <- c(0.25, 0.5, 1, 1.5, 2)
read_range <- c(50, 100, 200, 300, 500)

#colnames are samples input to DESeq
scolnames <- paste(rep('seal',n), 1:n, sep='_')
hcolnames <- paste(rep('human',n), 1:n, sep='_')
tcolnames <- c(scolnames, hcolnames)

#rownames are 'genes' we want to test -- combination of fold changes and proportions (a is bigger half, c is control, full counts)
rnames <- expand.grid(paste('fc',as.character(fc_range),sep='_'), c('gene_a','gene_b','gene_c'))
rnames <- paste(rnames$Var1, rnames$Var2,sep='_')

out <- data.frame(lfc=0,condition='test',nread=0)

for (nreads in read_range){

#create input matrix
test <- data.frame(matrix(0, length(rnames), length(tcolnames)), row.names=rnames)
colnames(test) <- tcolnames

#fill human data with control values (samples from known normal dist, mean=nreads, sd=sqrt(nreads))
test[,(n+1):length(test)] <- matrix(as.integer(rnorm((nrow(test)*(length(test) - (n))), mean=nreads, sd=sqrt(nreads))), nrow(test), (length(test) - (n)))


#fill seal data with various fold changes and gene part fractions
for (fc in fc_range){
  ret <- draw(fc, nreads, a, n)
  gene_a <- ret[[1]]
  gene_c <- ret[[2]]
  gene_b <- ret[[3]]
  test[paste0('fc_',fc,'_gene_a'),1:n] <- gene_a
  test[paste0('fc_',fc,'_gene_b'),1:n] <- gene_b
  test[paste0('fc_',fc,'_gene_c'),1:n] <- gene_c
  
}

#count matrix filled with simulated data, now fill metadata for DESeq
meta<-data.frame(row.names = tcolnames, species=c(rep('seal',n),rep('human',n)))

#generate DESeq object
dds <- DESeqDataSetFromMatrix(countData = test,
                              colData = meta, 
                              design = ~ species)
dds <- DESeq(dds)
res <- results(dds)

out <- rbind(out, data.frame(lfc=res$log2FoldChange, condition=rownames(res), nread=nreads))
##lfc doesn't match up very well....? whatever, that's what control is for


}

out <- out %>% separate(condition, c('trash1','fc','trash2','condition'),sep='_')




ggplot(out)+
  geom_line(mapping=aes(x=nread, y=lfc, color=condition))+
  facet_wrap(~fc)





gdata <- expand.grid(reads=read_range, fc=fc_range) #get combo of all read and fc 
gdata <- gdata %>% mutate(hcounts = reads) #make human the baseline counts
gdata <- gdata %>% mutate(scounts = reads*fc) #create true seal counts based on fc
gdata <- gdata %>% mutate(fcreads = as.integer(reads*fc)) 
gdata <- gdata %>% mutate(scountsa = draw(fcreads, a, r)) #use a uniform distribution to draw some percent from a
gdata <- gdata %>% mutate(scountsb = fcreads-scountsa) #rest of counts go to b
gdata <- gdata %>% mutate(true_lfc = log2(scounts/hcounts)) #true log fold change
gdata <- gdata %>% mutate(comb_lfc_a = case_when(scountsa==0~log2(scountsb/hcounts),
                                                scountsb==0~log2(scountsa/hcounts), 
                                                .default = log2(scountsa/hcounts) + log2(scountsb/hcounts)))
gdata <- gdata %>% mutate(comb_lfc_b = case_when(scountsa==0~b*log2(scountsb/hcounts),
                                                 scountsb==0~a*log2(scountsa/hcounts), 
                                                .default = a*log2(scountsa/hcounts) + b*log2(scountsb/hcounts)))
gdata <- gdata %>% mutate(comb_lfc_c = case_when(scountsa==0~log2(b * scountsb/hcounts),
                                                 scountsb==0~log2(a * scountsa/hcounts), 
                                                .default = log2(a * scountsa/hcounts) + log2(b * scountsb/hcounts)))
gdata <- gdata %>% mutate(comb_lfc_d =  case_when(scountsa==0~b*log2(scountsb/hcounts)+ log2(reads*a)/2,
                                                  scountsb==0~a*log2(scountsa/hcounts) + log2(reads*a)/2, 
                                                .default = a*log2(scountsa/hcounts) + b*log2(scountsb/hcounts) + log2(reads*a)/2))
gdata <- gdata %>% pivot_longer(c(true_lfc, comb_lfc_a, comb_lfc_b, comb_lfc_c, comb_lfc_d), names_to = 'strategy', values_to = 'est_fc')
ggplot(gdata)+
  geom_line(mapping=aes(x=reads, y=est_fc, color=strategy))+
  facet_wrap(~fc)+
  labs(caption = 'comb_lfc_a = log2(scountsa/hcounts) + log2(scountsb/hcounts) \n
        comb_lfc_b = a*log2(scountsa/hcounts) + b*log2(scountsb/hcounts) \n
        comb_lfc_c = log2(a * scountsa/hcounts) + log2(b * scountsb/hcounts) \n
       comb_lfc_d = a*log2(scountsa/hcounts) + b*log2(scountsb/hcounts) + log2(reads*a)/2')
  
ggsave('figs/deseq_sims_combd.png', heigh=3000, width=4000, units = 'px', dpi=300)  

