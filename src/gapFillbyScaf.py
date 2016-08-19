#from tranBlasr2LPAndGraph.py import *
import sys
##############################
#
#import: component_6, sg_copy_number, re_pacbio_read.fasta, map_detail output: gap_filling_after_3rd 
#
##############################
def read_component(f):
    line = f.readline().strip()
    components = {}
    gaps = {}
    while line:
        comId = int(line.split()[1].strip())
        components[comId] = f.readline().strip().split()
        gaps[comId] = f.readline().strip().split()
        line = f.readline().strip()
    return components,gaps

def read_scaf(f):
    line = f.readline().strip()
    scafs = {}
    while line:
        copyNum = int(line.split()[1])
        Id = int(line.split()[0].split('_')[1])
        line = f.readline().strip()
        scafs[Id] = line
        line = f.readline().strip()
    return scafs

def read_read(f):
    line = f.readline().strip()
    reads = {}
    while line:
        Id = int(line.split('_')[1])
        line = f.readline().strip()
        reads[Id] = line
        line = f.readline().strip()
    return reads


def read_map(f):
    line = f.readline().strip()
    line = f.readline().strip()
    mapDetails = {}
    while line:
        con_i = int(line.split()[0])
        con_j = int(line.split()[1])
        readId = int(line.split()[2])
        con_iend = int(line.split()[3])
        con_jbegin = int(line.split()[4])
        key = "%d_%d" % (con_i, con_j)  
        if not mapDetails.has_key(key):
            mapDetails[key] = []
        mapDetails[key].append((readId, con_iend, con_jbegin))
        line = f.readline().strip()
    return mapDetails


if len( sys.argv ) != 6:
    print "Usage: " + sys.argv[0] + "component_6, sg_copy_number, re_pacbio_read.fasta, map_detail, gap_filling_after_3rd"
    exit(0)

#workspace = sys.argv[1]

print "[info] read component from file " + sys.argv[1]
with file( sys.argv[1], 'r' ) as componentFile:
    components, gaps = read_component(componentFile)
print "[info] components size" , len(components) 


print "[info] read sg_copy_number from file " + sys.argv[2]
with file ( sys.argv[2], 'r' ) as seqFile:
    seq = read_scaf(seqFile)
print "[info] seq size(2rd scaf)" , len(seq) 

print "[info] read re_pac from file " + sys.argv[3]
with file ( sys.argv[3], 'r' ) as readFile:
    reads = read_read(readFile)
print "[info] reads size" , len(reads) 

print "[info] read map from file " + sys.argv[4]
with file (sys.argv[4], 'r') as mapFile:
    mapDetails = read_map(mapFile)
print "[info] map details size" , len(mapDetails) 


def suitHole(dis, details):
    tempdis = details[0][2] - details[0][1]
    minD = abs(tempdis-dis)
    readSegment = reads[details[0][0]][details[0][1]+1: details[0][1]]
    for i in range(1,len(details)):
        tempdis = details[i][2] - details[i][1]
        if(tempdis<minD):
            minD = abs(tempdis - dis)
            readSegment = reads[details[0][0]][details[0][1]+1: details[0][1]]
    print minD
    return readSegment 

scaf = {}
#base on seq

fout = open(sys.argv[5],"w")
for com in components:
    print com,components[com]
    temp = seq[int(components[com][0])]
    if len(components[com]) >= 2 :
        for i in range(0, len(components[com])-1 ):
            j=i+1
            con_i = int(components[com][i])
            con_j = int(components[com][j])
            if(int(gaps[com][i]) > 0):
                key = "%d_%d" % (con_i, con_j)
                readsSegment = suitHole(int(gaps[com][i]), mapDetails[key])
                temp = temp + readsSegment
                temp = temp + seq[con_j]
            else:
                temp = temp + seq[con_j][-int(gaps[com][i]):]
    scaf[com] = temp
    fout.write(">scaf_%d\n" % (com))
    fout.write("%s\n" % scaf[com])

