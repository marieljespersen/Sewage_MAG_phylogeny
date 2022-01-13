#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 14:42:05 2021

This script checks if a multifasta is in frame aligned by checking if any of the sequences has a start codon in the beginning
Takes as inputs: 
    a directory containing multi fasta files
Outputs a file containing the file/gene names of the files that do start with a start codon

@author: maloj
"""

import sys, getopt
from Bio import SeqIO
import os

#python /home/projects/cge/people/maloj/src/dnds/check_start_codon.py -d /home/projects/cge/data/projects/5001/Binning_vamb/clusters/C5/ortho_aln/ -o startcodon_genes.txt
    
def main(argv):
    # get inputs
    DIR = ''
    try:
        opts, args = getopt.getopt(argv,"hd:o:",["directory=", "ofile="])
    except getopt.GetoptError:
        print('check_start_codon.py -d <directory>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('check_start_codon.py -f <directory>')
            sys.exit()
        elif opt in ("-d", "--dir"):
            DIR = arg
        elif opt in ("-o", "--ofile"):
            OUT_FILE = arg
    
    outfile = open(OUT_FILE,'w')            
    for filename in os.listdir(DIR):
        if filename.endswith(".aln"):
            genename = filename.strip('.aln')
            FASTA = os.path.join(DIR, filename)
            print(genename)
            # read multi fasta
            fasta_sequences = SeqIO.parse(open(FASTA),'fasta')
            
            # check if first codon is start codon in any sequence. End if start codon is found
            for fasta in fasta_sequences:
                first_codon = str(fasta.seq)[0:3]
                if(first_codon == "ATG" or first_codon == "TTG" or first_codon == "GTG"):
                    #print(first_codon)
                    outfile.write(genename + "\n")
                    break
        else:
            continue
    
    outfile.close()
    

 
if __name__ == "__main__":
   main(sys.argv[1:])


  
    











