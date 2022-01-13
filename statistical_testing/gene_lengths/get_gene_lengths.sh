#!/bin/bash

while read c; do
	ortho_dir="/home/projects/cge/data/projects/5001/Binning_vamb/clusters/$c/ortho"
	echo "geneid length" > R2_bias/gene_lengths/$c.gene_lengths.txt
	ls $ortho_dir | while read f
	do
		gene_length=$(awk '{/>/&&++a||b+=length()}END{print b/a}' $f)
		geneid=$(echo $f | awk -F".fna" '{print $1}' | sed 's/ortho_group/gene/g')
		echo $geneid $gene_length >> R2_bias/gene_lengths/$c.gene_lengths.txt
	done 
	echo $c
done < treeclusters.txt