#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 11:13:28 2021

This script removes stop codons from a multi fasta alignment
Takes as inputs: multi fasta file
    
Outputs a fasta file with stop codons removed from the sequences

@author: maloj
"""

import sys, getopt
from Bio import SeqIO

# python /home/projects/cge/people/maloj/src/dnds/remove_stop_codons.py -f /home/projects/cge/data/projects/5001/Binning_vamb/clusters/C5/codeml/genewise/gene_100/gene_100.aln.fna -o /home/projects/cge/data/projects/5001/Binning_vamb/clusters/C5/codeml/genewise/gene_100/gene_100.nostop.aln.fna

# this function was found in: https://www.biostars.org/p/216181/ made by WouterDeCoster 
def identifystops(sequence):
	"""
	Taking a list of codons from a sequence starting with a startcodon
	identifies stopcodons in the same reading frame, converting them to "---"
	"""
	output = [] #An output list to contain the result
	for elem in sequence: #Remember that I split up the sequence in codons, which are now used to loop over
		if elem in ["TAA","TAG","TGA"]:
			output.append("---") #When a stopcodon is found, replace it by a "---" in the output list
		else:
			output.append(elem) #If the codon is no stopcodon, just keep the codon
	return("".join(output)) #Return the output, but join the list to a string so it again is a continuous sequence


def main(argv):
    # get inputs
    try:
        opts, args = getopt.getopt(argv,"hf:o:",["fasta=", "ofile="])
    except getopt.GetoptError:
        print('check_start_codon.py -f <fasta>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('check_start_codon.py -f <fasta>')
            sys.exit()
        elif opt in ("-f", "--fasta"):
            filename = arg
        elif opt in ("-o", "--ofile"):
            outfile = arg

    # read multi fasta
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    
    
    # for each entry in fasta; replace stop codons with ---
    new_fasta=[]
    for record in fasta_sequences:
        codon_list = [str(record.seq)[i:i+3] for i in range(0, len(str(record.seq)), 3)]
        no_stop_seq = identifystops(codon_list)
        #record.seq=Seq(no_stop_seq)
        new_fasta.append('>%s\n%s' % (record.id, no_stop_seq))  

    # write sequences witout stop codons in new file
    with open(outfile, 'w') as f:
        f.write('\n'.join(new_fasta))
    f.close()
 
if __name__ == "__main__":
   main(sys.argv[1:])


  
    











