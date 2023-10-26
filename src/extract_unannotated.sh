#!/bin/bash

bedtools getfasta -fi raw_data/ref/ncbi_dataset/data/GCF_016881025.1/GCF_016881025.1_HiC_Itri_2_genomic.fna -bed data/cons_unannotated.bed -nameOnly > data/cons_unannotated.fasta
