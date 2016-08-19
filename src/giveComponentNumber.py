#!/usr/bin/python

import sys

def read_contig(f):
    fout = open(sys.argv[2], "w")
    lenout = open(sys.argv[3], "w")
    comout = open(sys.argv[4], "w")
    sgout = open(sys.argv[5], "w")
    
    line = f.readline()
    comNum = 0 
    while line:
        cols = line.strip().split("_")
        if len(cols) == 3:
            fout.write(">scaf_%s_length_%s\n" % (cols[1],cols[2]))
            lenout.write(cols[2] + "\n")
            
            if int(cols[2]) >= 150:
                if cols[1] == "0":
                    sgout.write(">seq_%s\t%d\n" % (cols[1], 1))
                else:
                    sgout.write("\n>seq_%s\t%d\n" % (cols[1], 1))
                
                comout.write(">component\t%d\n" % comNum)
                comout.write(cols[1] + "\n\n")
                comNum = comNum + 1
            else:
                if cols[1] == "0":
                    sgout.write(">seq_%s\t%d\n" % (cols[1], 1))
                else:
                    sgout.write("\n>seq_%s\t%d\n" % (cols[1], 2))
        else: 
            sgout.write(line.strip())
            fout.write(line)
        line = f.readline()

if len( sys.argv) != 6:
    print "Usage: " + sys.argv[0] + " 31mer.scaf_seq_with_gaps contigs.fasta contig_len component_5 sg_copy_number.fa"
    exit(0)

with file(sys.argv[1]) as f:
    read_contig(f)
#./blasr contig.fasta contig.fasta -out bb -m 0 -nproc 8
