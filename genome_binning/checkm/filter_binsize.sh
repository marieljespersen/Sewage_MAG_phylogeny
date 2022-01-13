#!/bin/bash

CLUSTERS="$1"
OUTFILE="$2"

#adding contig lengths to tmp cluster.tsv
awk '{print $2}' $CLUSTERS | awk -F_length_ '{print $2}' | awk -F_ '{print $1}' > contig_lengths.tmp
paste $CLUSTERS contig_lengths.tmp > tmp_clusters.tsv

#counting number of bp in each bin
module load datamash/1.4
datamash groupby 1 sum 3 < tmp_clusters.tsv > binsize.tmp

#creating list of bins containing 10^6 bp or more
awk '$2 >= 1000000 {print $1 }' < binsize.tmp > largebins.tsv

#creating subset of clusters.tsv containing only bins of 10^6 bp or more
touch $OUTFILE
for bin in $(cat largebins.tsv)
do
	grep "^"$bin$'\t' $CLUSTERS >> $OUTFILE
	echo $bin
done

#removing tmp files
rm binsize.tmp tmp_clusters.tsv contig_lengths.tmp
