#!/bin/bash

CLUSTER=$1
BIN=$2

cd $CLUSTER

# remove stop codons from protein sequences
sed -i "s/\*//g" prodigal/proteomes/${BIN}.predicted_proteins.faa 
#echo $(ls -lh prodigal/proteomes/${BIN}.predicted_proteins.faa)

# load modules
module load perl/5.30.2 anaconda3/4.4.0 java/1.8.0-openjdk interproscan/5.36-75.0

# scan protein sequences from bin agianst InterPro db
interproscan.sh -goterms -pa -f tsv --cpu 20 \
-i prodigal/proteomes/${BIN}.predicted_proteins.faa \
-b interproscan/${BIN}.interproscan

module unload interproscan/5.36-75.0 anaconda3/4.4.0 java/1.8.0-openjdk