#!/usr/bin/python

import sys
EDGE_LENGTH_CUTOFF = 150
COVRAGE_PERCENT_CUTOFF = 0.8
#COVER_NUM_CUTOFF = 10
COVER_NUM_CUTOFF = 5
SIM_CUTOFF = 85
SCORE_CUTOFF = -600
#SIM_UPPER = 99


###pay attention contigsId from 0 or 1
def read_contigLen(f):
    contigLen = []
    line = f.readline().strip()
    while line:
        contigLen.append( int(line) );
        line = f.readline().strip()
    return contigLen


def read_component(f):
    
    print "read contig len from", sys.argv[1]
    with file(sys.argv[1] ) as contigLenFile:
        contigLen = read_contigLen(contigLenFile)
    print "contig number", len(contigLen)

    fout = open(sys.argv[6], "w")

    contigInComp = {}
    line = f.readline().strip()
    componentNum = 0
    
    print "read component from", sys.argv[2]
    while line:
        componentId = int(line.split()[1].strip())
        componentNum = componentNum + 1
        contigs = f.readline().strip().split()
        gaps    = f.readline().strip().split()
        contigInComp[ int(contigs[0]) ] = (componentId, 0)
        last = 0
        for i in range(1, len(contigs)):
            cur = last + contigLen[ int(contigs[i-1]) ] + int(gaps[i-1])
            contigInComp[ int(contigs[i]) ] = (componentId, cur)
            last = cur
        
        #contigId is from 1 or from 0
        #fout.write(str(last + contigLen[ int(contigs[ len(contigs)-1 ]) -1 ]) + "\n")
        fout.write(str(last + contigLen[ int(contigs[ len(contigs)-1 ])]) + "\n")
        line = f.readline().strip()
    fout = open(sys.argv[7], "w")
    fout.write("EDGE_CLUSTER_NUM=%d\n" % (componentNum))
    return contigLen, contigInComp, componentNum
    
def read_blasr(f):
    with file( sys.argv[2] ) as componentFile:
        contigLen, contigInComp, componentNum = read_component(componentFile)
    graph = {}
    contigGraph = {}
    read2Contigs = {}
    line = f.readline().strip()
    '''
    # format -m 0
    print "readId , contigId, contig in read, readStart, readEnd, readLength, contigStart, contigEnd, contigLength, sim, score"
    '''
    map_detail = {}
    while line:
        # format -m 1
        #print line
        contigId = int(line.split()[1].split("_")[1].strip())
        readId = int(line.split()[0].split("_")[1].strip())
        
        sim = float(line.split()[5].strip())         
        score = float(line.strip().split()[4])
        
        readStart = int(line.split()[9].strip())
        readEnd   = int(line.split()[10].strip())
        readLength = int(line.split()[11].strip())
        
        contigStart = int(line.split()[6].strip())
        contigEnd   = int(line.split()[7].strip())
        contigLength = int(line.split()[1].split("_")[3].strip())
       
        assert contigEnd > contigStart
        mapLength = contigEnd - contigStart

        readDirection = int(line.split()[2].strip()) 
        contigDirection = int(line.split()[3].strip())
       
        #map in single direction
        if (readDirection != 0 or contigDirection != 0):
            line = f.readline();
            continue
        ''' 
        if (contigId == 0):
            print "readId, contigId, readDirection, contigDirection, score, sim, contigStart, contigEnd, contigLength, readStart, readEnd, readLength" ,readId, contigId, readDirection, contigDirection, score, sim, contigStart, contigEnd, contigLength, readStart, readEnd, readLength
        '''
        # rule -600 150 85
        if (score > SCORE_CUTOFF or contigLength < EDGE_LENGTH_CUTOFF or sim < SIM_CUTOFF):
            line = f.readline()
            continue
        
        shouldLength = mapLength + min(contigStart, readStart) + min(contigLength-contigEnd, readLength-readEnd) 
        
        # actual map / should map >  COVRAGE_PERCENT_CUTOFF(0.8)
        if (mapLength < shouldLength * COVRAGE_PERCENT_CUTOFF):
            line = f.readline()
            continue 
        
        #print "after filter"
        if not read2Contigs.has_key(readId):
            read2Contigs[ readId ] = [];
        read2Contigs[ readId ].append(( contigId, readStart - contigStart))
        key = '%d_%d' % (readId, contigId) 
        if map_detail.has_key(key):
            print "error, same read, same contig, twice" , key 
        else:
            map_detail[key] = (readStart-contigStart, readEnd + (contigLength-contigEnd))    
        line = f.readline()
    
    print "read 2 contigs num", len(read2Contigs)

    fout = open(sys.argv[8],"w")
    fout.write("con_i\tcon_j\treadId\tcon_i_end\tcon_j_begin\n")
    for idx in read2Contigs:
        read2Contigs[idx].sort(key=lambda x:x[1])
        for i in range(len( read2Contigs[idx] )):
            contig_i, readPos_i = read2Contigs[idx][i]
            if not contigInComp.has_key( contig_i ):
                continue
            com_i, comPos_i = contigInComp[ contig_i ]
            for j in range(i+1, len( read2Contigs[idx] )):
                contig_j, readPos_j = read2Contigs[idx][j]
                if not contigInComp.has_key( contig_j ):
                    continue
                com_j, comPos_j = contigInComp[ contig_j ]
                if com_j == com_i:
                    continue
                dis = readPos_j - readPos_i + comPos_i - comPos_j
                key = '%d_%d' % (com_i, com_j)
                contigKey = '%d_%d' % (contig_i, contig_j)
                map_left = "%d_%d" % (idx,contig_i)
                map_right = "%d_%d" % (idx,contig_j)
                fout.write("%d\t%d\t%d\t\t%d\t%d\n" % (contig_i,contig_j,idx,map_detail[map_left][1] ,map_detail[map_right][0]))
                if not graph.has_key(key):
                    graph[key] = (dis, 1)
                    contigGraph[contigKey] = (dis, 1)
                else:
                    x, y = graph[key]
                    graph[key] = (x + dis, y + 1)
                    x, y = contigGraph[contigKey]
                    contigGraph[contigKey] = (x + dis, y + 1)
    eleForDelete = []
    print "contig graph"
    print "contig graph size" , len(contigGraph)
    for key in contigGraph:
        contigs = key.split("_")
        con_i = int(contigs[0])
        con_j = int(contigs[1])
        dis, cov = contigGraph[key]
        #print con_i,con_j,dis,cov,dis/cov
    #computer degree 
    inDegree = {}
    outDegree = {}
    print "initial graph"
    print "graph size", len(graph)
    for key in graph: 
        #print key,graph[key]
        components = key.split("_")
        com_i = int(components[0])
        com_j = int(components[1])
        if com_j in inDegree:
            inDegree[com_j] = inDegree[com_j] + 1
        else:
            inDegree[com_j] = 1
    
        if com_i in outDegree:
            outDegree[com_i] = outDegree[com_i] + 1
        else:
            outDegree[com_i] = 1
    '''
    print "indegree"
    print inDegree

    print "outdegree"
    print outDegree
    '''
    for key in graph:
        components = key.split("_")
        com_i = int(components[0])
        com_j = int(components[1])
        dis, cov = graph[key]
        if (cov < COVER_NUM_CUTOFF and (inDegree[com_j] > 1 or outDegree[com_i] > 1)) or (cov < 2):
            #print "cov < 5 delete", com_i, com_j, graph[key]
            eleForDelete.append(key)
    for key in eleForDelete:
        del graph[key]
    #print graph["2_515"]
    print "final graph size" , len(graph)
    return graph, componentNum

if len( sys.argv) != 9:
    print "Usage: " + sys.argv[0] + " contig_len component_n blasr.out contig_arc_graph_after_remove_ambigous_arcs_n position_lp_n.math edge_cluster_len_n scaffold_parameter_n map_detail"
    exit(0)

with file(sys.argv[3]) as f:
    graph, componentNum = read_blasr(f)
    '''
    print "graph"
    for key in graph:
        print key,graph[key]
    '''
    '''
    #filter degree too larger
    #see newARCS3
    '''
#write LP
fout = open(sys.argv[4], "w")
for key in graph:
    components = key.split("_")
    com_i = components[0]
    com_j = components[1]
    dis, cov = graph[key]
    fout.write(com_i + "\t" + com_j + "\t" + str(int(dis/cov)) + "\t" + str(cov) + "\t0.0\n")
fout = open(sys.argv[5], "w")
for i in range(componentNum):
    fout.write("var x_%d;\n" % i)
for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    fout.write("var e_%d_%d;\n" % (com_i, com_j))
    fout.write("var E_%d_%d;\n" % (com_i, com_j))
fout.write("minimize z: \n")
count = 0
for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    fout.write(" E_%d_%d + " % (com_i, com_j))
    count += 1
    if count % 10 == 0:
        fout.write("\n")
fout.write("0;\n")

index = 1
for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    dis, cov = graph[key]
    fout.write("s.t. con%d : x_%d - x_%d + e_%d_%d = %d;\n" % (index , com_j , com_i , com_i , com_j , int(dis/cov)))
    index += 1

for key in graph:
    components = key.split("_")
    com_i = int(components[0])
    com_j = int(components[1])
    fout.write("s.t. con%d : E_%d_%d + e_%d_%d >= 0;\n" % (index , com_i , com_j, com_i , com_j))
    index += 1
    fout.write("s.t. con%d : E_%d_%d - e_%d_%d >= 0;\n" % (index , com_i , com_j, com_i , com_j))
    index += 1
fout.write("end;\n")



