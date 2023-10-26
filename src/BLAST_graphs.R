library(tidyverse)
library(jsonlite)
setwd("/home/dylan/Postdoc/Hindle/dylan")
fjson <- fromJSON('data/fastqc/BLAST_gc_fail.json', flatten=TRUE)
wjson <- read_json('data/fastqc/BLAST_gc_warn.json')

#sample names
print(length(fjson$BlastOutput2$report.results.search.query_title))

