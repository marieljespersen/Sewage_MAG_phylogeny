#!/bin/bash

module load mafft/7.453 trimal/1.4.1
module unload R gcc perl
module load anaconda3/4.4.0

CLUSTER="C37"
CLUSTER_DIR="/home/projects/cge/data/projects/5001/Binning_vamb/clusters/$CLUSTER/codeml/"

mkdir $CLUSTER_DIR
cd $CLUSTER_DIR
cp /home/projects/cge/data/projects/5001/Binning_vamb/clusters/C5/codeml/codeml.ctl .

START_CODON_SCRIPT="/home/projects/cge/people/maloj/src/dnds/check_start_codon.py"
ALIGNMENTS="/home/projects/cge/data/projects/5001/Binning_vamb/clusters/$CLUSTER/ortho_aln/"
START_CODON_GENES="/home/projects/cge/data/projects/5001/Binning_vamb/clusters/$CLUSTER/codeml/start_codon_genes.txt"

python $START_CODON_SCRIPT -d $ALIGNMENTS -o $START_CODON_GENES
# loop through all genes of cluster and save files for running codeml
touch run_codeml.sh
while read i
do
	# create dir per gene and enter to save one codeml.ctl file per gene
	gene=$(echo $i | sed 's/ortho_group/gene/g')
	mkdir $gene
	cd $gene
	
	# identify path for phylip alignmen
	ortho_phylip=$(echo $i.ortho_phy)
	phy_path=$(echo ../../ortho_phy/$ortho_phylip)
	fasta_file=$(echo $i.aln.fna)
	trimal -in $phy_path -fasta -out $fasta_file

	#remove stop codons from alignments prior to running codeml
	fasta_ns=$(echo $i.nostop.aln.fna)
	python /home/projects/cge/people/maloj/src/dnds/remove_stop_codons.py -f $fasta_file -o $fasta_ns
	
	# identify path for tree file
	treefile=$(ls ../../gene_trees/$gene.treefile)
	
	# name results file
	results_file=$(echo $gene.selection.results)
	
	# change codeml.ctl file to contain input and outputs for the specific gene
	cp ../codeml.ctl .
	sed -i "s|ortho.renamed.aln.fna|$fasta_ns|g" codeml.ctl
	sed -i "s|C5.astral.branch_length.clean.tree|$treefile|g" codeml.ctl
	sed -i "s|selection.results|$results_file|g" codeml.ctl
	cd ..
	
	# if running codeml from within the loop, it stops after one gene, instead save in file
	echo "cd $gene" >> run_codeml.sh
	echo "codeml" >> run_codeml.sh
	echo "cd .." >> run_codeml.sh
	echo $i
done < $START_CODON_GENES