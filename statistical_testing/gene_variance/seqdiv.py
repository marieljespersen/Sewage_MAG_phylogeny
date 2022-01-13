#!/usr/bin/env python3
# By Anders Gorm Pedersen, gorm@cbs.dtu.dk, CBS, Technical University of Denmark, 2015
# Computes average sequence diversity

import sys, os.path, seqlib
from optparse import OptionParser

# Build commandline parser
parser = OptionParser(usage="usage: %prog [options] [SEQFILE]",
                      version="0.2")

parser.add_option("-I", type="choice", dest="informat",
                  choices=["nexus", "phylip", "fasta", "clustal", "raw", "tab"], metavar="FORM",
                  help="Input format: clustal, fasta, nexus, raw, phylip, tab")

parser.add_option("-g", "--ignoregaps", action="store_true", dest="ignoregaps",
                      help="ignore gappy columns from pairwise comparisons. Note: this means different sets of alignment positions may be used for each pairwise comparison")

parser.add_option("-p", "--printpairs", action="store_true", dest="printpairs",
                      help="print all pairwise distances instead of average and standard error")

parser.add_option("-n", "--noheader", action="store_true", dest="noheader",
                      help="Do not print header (just mean and standard error)")

parser.set_defaults(informat="fasta", ignoregaps=False, printpairs=False, noheader=False)

# Parse commandline, open seqfile and get seqs (if file name is given and file exists)
(options, args) = parser.parse_args()

# If no filename is given: assume stdin
if len(args) < 1:
    filename ="-"
else:
    filename = args[0]
    if not os.path.isfile(filename):
        parser.error("file %s not found." % filename)

# Select appropriate filehandle and read seqs as alignment
seqfile = seqlib.Seqfile(filename, options.informat)
seqs = seqfile.read_alignment()

# Find average sequence diversity and standard deviation
if options.ignoregaps:
    distmatrix = seqs.distmatrix("pdist_ignoregaps")
else:
    distmatrix = seqs.distmatrix("pdist")

# if requested: print all pairwise distances
if options.printpairs:
    print (distmatrix.pairdists())
    
# else: print summary:
else:
    mean, sem = distmatrix.distsummary()
    if options.noheader:
        print ("%f\t%f" % (mean, sem))
    else:
        print ("# Mean pi       Standard error")
        print ("%.6f           %.6f" % (mean, sem))
