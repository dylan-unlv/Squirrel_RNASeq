#!/bin/bash

#####
#first we must create an index, which requires the gff in fasta format. I downloaded a tool to convert, used here
#this step needs to be run only once, if Itri.idx exists, skip to generating quants for the tissue of interest
#####

#fix genome headers to have ONLY scaffold name in fasta file..
#python src/fix_fna_headers.py

##convert gff to fasta
#/data8/han_lab/dbarth/trish_transcriptomics/tools/tophat-2.1.1.Linux_x86_64/gtf_to_fasta raw_data/ref/genomic.filtered.gff raw_data/ref/GCF_016881025.1_HiC_Itri_2_genomic.fxheader.fna data/salmon/salmon_idx/Itri_gff.fasta

##fix effed up headers so that the feature name (gene, rna, transcript, exome, CDS etc ) is the ONLY header in the fasta file
#python src/fix_gff_headers.py

##create decoy file to catch DNA alignments
#cat raw_data/ref/GCF_016881025.1_HiC_Itri_2_genomic.fxheader.fna | cut -d " " -f 1 | sed -n 's/>//p' | tr -d '>' > data/salmon/salmon_idx/decoys.txt

##concatenate the gff fasta with the genome fasta
cat data/salmon/salmon_idx/Itri_salmon_idx.fxheader.fasta raw_data/ref/GCF_016881025.1_HiC_Itri_2_genomic.fxheader.fna > data/salmon/salmon_idx/Itri_salmon_idx.fasta

## once we have the final fasta file, we can create the salmon index
salmon index -t data/salmon/salmon_idx/Itri_salmon_idx.fasta --decoys data/salmon/salmon_idx/decoys.txt -p 70 -k 31 -i data/salmon/salmon_idx/Itri.idx


#####
#now that the index is generated, we generate quants
#####

mkdir -p data/salmon/quants
tissue=kidney #heart brain liver kidney
mkdir -p data/salmon/quants/${tissue}
for fn in raw_data/${tissue}/*1.fq.gz
do 
samp=$(basename ${fn} | sed "s/\_1.fq.gz//")
salmon quant -i data/salmon/salmon_idx/Itri.idx -l A -1 raw_data/${tissue}/${samp}_1.fq.gz -2 raw_data/${tissue}/${samp}_2.fq.gz -p 70 --validateMappings -o data/salmon/quants/${tissue}/${samp}_quant
done
