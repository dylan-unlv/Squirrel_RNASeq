#!/bin/bash

#first get list of unidentified genes with both 'LOC' and 'CUN' prefixes
GNAMES=$(grep -E -i 'LOC|CUN' data/wgcna/networks/consensus_gene_list.csv)

#del output file if exists
rm -f data/cons_unannotated.bed

#next, for each GNAME, search gff for the corresponding rows and take min-max of the ranges
for G in $GNAMES
do
	min=$(fgrep -i $G data/genomic.gff | awk '{print $1 "\t" $4 "\t" $5}' | awk 'NR == 1 || $2 < min {line = $0; min = $2}END{print min}')
	max=$(fgrep -i $G data/genomic.gff | awk '{print $1 "\t" $4 "\t" $5}' | awk 'NR == 1 || $3 > max {line = $0; max = $3}END{print max}')
	chr=$(fgrep -i $G data/genomic.gff -m 1 | awk '{print $1}')
	echo -e $chr'\t'$min'\t'$max'\t'$G >> data/cons_unannotated.bed
done

