#!/bin/bash

mkdir R2_bias/gene_div
SEQ_DIV="/home/projects/cge/people/maloj/src/gene_variance/seqdiv.py"
while read c; do
	ortho_dir="/home/projects/cge/data/projects/5001/Binning_vamb/clusters/$c/ortho_aln"
	echo "geneid variance" > R2_bias/gene_div/$c.gene_div.txt
	ls $ortho_dir | while read f
	do
		div_res=$($SEQ_DIV -g -I fasta $c/ortho_aln/$f)
		gene_diversity=$(echo $div_res | awk '{print $6}')
		geneid=$(echo $f | awk -F".fna" '{print $1}' | sed 's/ortho_group/gene/g')
		echo $geneid $gene_diversity >> R2_bias/gene_div/$c.gene_div.txt
	done 
	echo $c
done < treeclusters.txt