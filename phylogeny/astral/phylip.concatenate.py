#!/usr/bin/env python

import argparse
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def main(args):
   '''Main program'''
   
   # open out
   fh_out = open(args.o+".phy", 'w')
   fh_out_partition = open(args.o+".partitions.txt", 'w')
   
   # get phylip files to read from file
   files = []
   fh_in = open(args.i, 'r')
   for line in fh_in:
      line = line.rstrip()
      files.append(line)
   
   # parse
   phy_dict = {}
   lengths_list = []
   for i,f in enumerate(files):
      alignment = AlignIO.read(open(f), "phylip-relaxed")
      lengths_list.append(len(alignment[0]))
      if i == 0:
         for record in alignment:
            if record.id in phy_dict:
               raise ValueError('Same key found twice %s\n' % record.id)
            else:
               phy_dict[record.id] = str(record.seq)         
      else:
         for record in alignment:
            if record.id in phy_dict:
               phy_dict[record.id] = phy_dict[record.id]+str(record.seq)
            else:
               raise ValueError('New key found in file %s' % f)
   
   # create seqrecord
   out = []
   for key,value in phy_dict.items():
      phy = SeqRecord(MutableSeq(value), id=key, description="")
      print(phy)
      out.append(phy)
   # write out as phylip
   msa=MultipleSeqAlignment(out)
   c=AlignIO.write(msa, fh_out, "phylip-relaxed")
   fh_out.close()
   
   # write out partition file
   s = 0
   for i,l in enumerate(lengths_list):
      start = s+1
      end = start+l-1
      fh_out_partition.write('%s, genecodon%i = %i-%i\n' % (args.submatrix, i, start, end))
      s=end
   
   fh_out_partition.close()

if __name__ == '__main__':
   parser = argparse.ArgumentParser(prog='phylip.concatenate.py', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=110),
      description='''Concatenate (paste) phylip alignments together ''',
      usage='%(prog)s input output [options]')
   
   parser.add_argument('--i', help='file with list of phylip files to paste', required=True)
   parser.add_argument('--submatrix', help='substitution matrix', default='DNA')
   parser.add_argument('--o', help='output prefix', required=True)
   
   args = parser.parse_args()
   #args = parser.parse_args('--i cds_alignment.phylip.YPTB3958codon1.phy cds_alignment.phylip.YPTB3958codon2.phy cds_alignment.phylip.prfAcodon1.phy --o test'.split())
   
   main(args)


