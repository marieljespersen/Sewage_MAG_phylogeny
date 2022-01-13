#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_sonic(i, max, min):
    '''Parse sonic output'''
    
    #min = max - (max * 0.5) # unecessary when min is an input as well?
    # open file
    fh_in = open(i, 'r')
    
    # get header
    h = []
    header = fh_in.readline()
    header = header.rstrip()
    header_ = header.split('\t')
    header_list = header_[4::2]
    for h_element in header_list:
        h.append(h_element.split(".")[0])
    
    # parse data
    ortho = {}
    for line in fh_in:
        line = line.rstrip()
        fields = line.split('\t')
        # if sp_in_group from sonicparanoid is larger than minimum, save samplenames
        if  int(fields[2]) >= int(min):
            ortho[fields[0]] = fields[4::2]
    
    return h, ortho


def write(h, ortho, fasta_folder, fna, o):
    '''Write out ortho genes'''
    
    # read in fastas
    faa_list = []
    for faa in h:
        if fna:
            f = '%s/%s.predicted_proteins.fna' % (fasta_folder, faa)
        else:
            f = '%s/%s.predicted_proteins.faa' % (fasta_folder, faa)
        f_fh = open(f, 'r')
        f_dict = SeqIO.to_dict(SeqIO.parse(f_fh, 'fasta'))
        faa_list.append(f_dict)
    
    # go through ortho
    for key,value in ortho.items():
        n = 1
        if fna:
            f_out = '%s/ortho_group_%s.fna' % (o, key)
        else:
            f_out = '%s/ortho_group_%s.faa' % (o, key)
        fh_out = open(f_out, 'w')
        gene_list = []
        for i,v in enumerate(value):
            id_ = h[i]
            if v == '*':
                if fna:
                    s = SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNNN"), id=id_, description="missing")
                else:
                    s = SeqRecord(Seq("XXXXXXXXXXXXXXXXXXXXX"), id=id_, description="missing")
                n=n+1
            elif v.find(',') > -1:
                if fna:
                    s = SeqRecord(Seq("NNNNNNNNNNNNNNNNNNNNN"), id=id_, description="missing")
                else:
                    s = SeqRecord(Seq("XXXXXXXXXXXXXXXXXXXXX"), id=id_, description="missing")
                n=n+1
            else:
                s = faa_list[i][v]
                s.id = id_
            gene_list.append(s)
        
        # write out
        c=SeqIO.write(gene_list, fh_out, 'fasta')
        

def main(args):
    '''Main program'''
    
    h, ortho = parse_sonic(args.i, args.max, args.min)
    write(h, ortho, args.fasta_folder, args.fna, args.o)

if __name__ == '__main__':
   
   # create the parser
   parser = argparse.ArgumentParser(prog='sonic2fasta.py', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=130), usage='%(prog)s [options]', description='Extract fasta from sonic outputs')
   
   parser.add_argument('--i', help='input tsv', required=True)
   parser.add_argument('--max', help='max number of orthologs', required=True, type=int)
   parser.add_argument('--min', help='min number of orthologs', required=True, type=int)
   parser.add_argument('--fasta_folder', help='folder with genes as fasta', required=True)
   parser.add_argument('--fna', help='use fna instead of faa', default=False, action='store_true')
   parser.add_argument('--o', help='output npz', default=True)
   
   args = parser.parse_args()
   #args = parser.parse_args('--i para/runs/sonic_24320171214_fast_20cpus_ml05/ortholog_groups/ortholog_groups.tsv --max 69 --min 63 --fasta_folder genes --o test'.split())
   
   main(args)


