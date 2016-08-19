#!/usr/bin/python

import sys

def read_contig(f):
    fout = open(sys.argv[2], "w")
    lenout = open(sys.argv[3], "w")
    line = f.readline()
    number = 0
    while line:
        seq = f.readline().strip()
        #if len(seq) >= 200:
        fout.write(">contig_%d_length_%d\n" % (number,len(seq)))
        fout.write(seq + "\n")
        lenout.write(str(len(seq)) + "\n")
        number += 1
        line = f.readline()
if len( sys.argv) != 4:
    print "Usage: " + sys.argv[0] + " cdbg_copy_number.fa contig.fasta contig_len"
    exit(0)

with file(sys.argv[1]) as f:
    read_contig(f)
#./blasr contig.fasta contig.fasta -out bb -m 0 -nproc 8
