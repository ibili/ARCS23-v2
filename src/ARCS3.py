#!/usr/bin/python

import time
from datetime import datetime
from itertools import ifilter, imap
import string
import os,sys
import getopt

def configfile_parse(f, config):
    transform = {
            'INSERT_SIZE': int, 
            'LINK_QUALITY_PERCENT': float, 
            'PAIR_KMER_CUTOFF': int, 
            'PAIR_READS_CUTOFF': int, 
            'EDGE_LENGTH_CUTOFF': int, 
            'READ_LENGTH_CUTOFF': int,
            'q1': str, 
            'q2': str,
            'pb': str
            }
    library_list = []

    options = {}
    for line in ifilter(lambda x: len(x) > 0, imap(string.strip, f)):
        if line.startswith('#'): continue
        key, val = map(string.strip, line.split('=', 1))
        if transform.has_key(key):
            if options.has_key(key):
                library_list.append(options);
                options = {}
            if key == 'READ_LENGTH_CUTOFF':
                config[key] = transform[key](val) 
            if key == 'pb':
                config[key] = transform[key](val)
            else:
                options[key] = transform[key](val)
    if len(options) > 0:
        library_list.append(options)

    config['library_list'] = library_list

    return config

def command_run(arcs, cmd, args, config):
    proc = '%s %s %s' % (arcs, cmd, args)
    if config.has_key('verbose') and config['verbose']:
        print '-------------------------------------'
        print proc
        print '-------------------------------------'

    start = datetime.now()
    if not config.get('test', False) and os.system(proc) != 0:
        sys.exit(1)
    end = datetime.now()
    if not config.get('test', False) and config.has_key('runtime') and config['runtime']:
        print '%s time is %d seconds' %  (cmd, (end - start).seconds)

def workspace_check(workspace, config):
    if not os.path.exists(workspace) or not os.path.isdir(workspace):
        print 'directory %s  does not exist' % (workspace)
        print "creat this directory? yes/no"
        word = sys.stdin.readline().strip().lower()
        if word == 'y' or word == 'yes':
            os.makedirs(workspace)
        else:
            print "system exit"
            sys.exit(0)

    print 'path : %s' % (workspace)

print ""
print "======================================================================================"
print ""
print "                                 ARCS (Q Version 0.9)                               "
print "        Copyright(c) 2014, Renyu Zhang, Qing Xu, Dongbo Bu. All Rights Reserved.       " 
print ""
print "======================================================================================"
print ""
print datetime.now()
print ""

USAGE = """"USAGE: ARCS.py [options]
 DESCRIPTION:
 OPTIONS:
    -s --configure_file       configure file
    -d --workspace            workspace
    -K --kmer_size            kmer size
    -e --edge_length_cutoff   edge length cutoff
    -O --max_overlap          max overlap to detect conflict
    -h --help                 help information
    -v --version              software version
    -E --kmer_filter          filter low quality kmers
    -V --verbose              verbose model
    -t --test                 test model
    -x --arcs                 arcs path
    -P --scaf_output          after gapfilling output file name
"""
config = {
        'kmer_size': 27, 
        'max_overlap': 200, 
        'runtime': True,  
        'kmer_filter': False, 
        'test': False, 
        'verbose': False,
        'cpu_num': 1
        }

if os.environ.has_key('ARCS_CMD'):
    ARCS_CMD = os.environ['ARCS_CMD']
else:
    ARCS_CMD = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'arcs')

opts, var = getopt.getopt(sys.argv[1:], 
        's:d:K:e:hO:vp:EVtc:x:P:', 
        ['help', 'version', 'configure_file=', 'workspace=', 'kmer_size=', 'edge_length_cutoff=', 'max_overlap', 'CPU','kmer_filter', 'verbose', 'test', '=arcs','scaf_output=']
        )
if len(opts) == 0:
    print USAGE
    sys.exit(0)

for o, a in opts:
    if o == '-h' or o == '--help':
        print USAGE
        sys.exit(0)
    elif o == '-s' or o == '--configure_file':
        config['configure_file'] = a
    elif o == '-d' or o == '--workspace':
        config['workspace'] = a
    elif o == '-P' or o == '--scaf_output':
        config['scaf_output'] = a
    elif o == '-K' or o == '--kmer_size':
        config['kmer_size'] = int(a)
        kmer_size = int(a)
        if kmer_size < 1 or kmer_size > 196:
            print "kmer size must be in [1, 196]"
            print "system exit"
            sys.exit(1)
    elif o == '-e' or o == '--edge_length_cutoff':
        config['edge_length_cutoff'] = int(a)
    elif o == '-p' or o == '--CPU':
        config['cpu_num'] = int(a)
        cpu_num = int(a)
    elif o == '-O' or o == '--max_overlap':
        config['max_overlap'] = int(a)
    elif o == '-E' or o == '--kmer_filter':
        config['kmer_filter'] = True
    elif o == '-x' or o == '--arcs':
        ARCS_CMD = a
    elif o == '-t' or o == '--test':
        config['test'] = True
    elif o == '-V' or o == '--verbose':
        config['verbose'] = True
    elif o == '-v' or o == '--version':
        command_run(ARCS_CMD, '', '-v', {})
        sys.exit(0)

#
# check options
#
if not config.has_key('workspace'):
    print "please specify a workspace!!!"
    sys.exit(1)

if not config.has_key('configure_file'):
    print "please specify a configure file!!!"
    sys.exit(1)

#
# check the workspace
#
workspace_check(config['workspace'], config)

with file(config['configure_file']) as f:
    config = configfile_parse(f, config)

if len(config['library_list']) == 0:
    print "please specify libraries!!!"
    sys.exit(1)

start = datetime.now()

###################################################################
#
# arcs preprocess
#
###################################################################
args = '-K %d -o %s -e %d' % (config['kmer_size'], os.path.join(config['workspace'], 'kmers.arff'), config['READ_LENGTH_CUTOFF'])

#args = '-K %d -o %s -e -1' % (config['kmer_size'], os.path.join(config['workspace'], 'kmers.arff'))
if config['kmer_filter'] and config['kmer_size'] < 33:
    args = '%s -E' % (args)
args = '%s %s %s' % (args, config['library_list'][0]['q1'], config['library_list'][0]['q2'])
command_run(ARCS_CMD, 'preprocess', args, config)

###################################################################
#
# arcs assemble
#
###################################################################
args = '-d %s -K %d -i %s' % (config['workspace'], config['kmer_size'], os.path.join(config['workspace'], 'kmers.arff'))
command_run(ARCS_CMD, 'assemble', args, config)

###################################################################
#
# arcs copy_num_estimate 
#
###################################################################
args = '-s %s -i %s -G %s -C %s' % (os.path.join(config['workspace'], 'contig_parameter'), os.path.join(config['workspace'], 'condensed_de_bruijn_graph_after_trimming.data'), os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_0')) 
command_run(ARCS_CMD, 'copy_num_estimate', args, config)


###############################################
#
#only used generate scaffold_parameter
###############################################
for i, library in enumerate(config['library_list']):
    print "............ iter %d" % (i + 1)

    if not library.has_key('EDGE_LENGTH_CUTOFF'):
        if config.has_key('edge_length_cutoff'):
            library['EDGE_LENGTH_CUTOFF'] = config['edge_length_cutoff']
        else:
            library['EDGE_LENGTH_CUTOFF'] = config['kmer_size']

    ###################################################################
    #
    # arcs scaffold
    #
    ###################################################################
    args = '-d %s -K %d -C %s -f %s -e %d -1 %s -2 %s -L %d -P %f -i %d -r %d -R %d -p %d' % (config['workspace'], config['kmer_size'], os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_%d' % (i)), library['EDGE_LENGTH_CUTOFF'], library['q1'], library['q2'], library['INSERT_SIZE'], library['LINK_QUALITY_PERCENT'], i, library['PAIR_KMER_CUTOFF'], library['PAIR_READS_CUTOFF'], config['cpu_num'])
    command_run(ARCS_CMD, 'scaffold', args, config)

    
'''     
    ###################################################################
    #
    # arcs solveLP
    #
    ###################################################################
    args = '-s %s -i %s -o %s' % (os.path.join(config['workspace'], 'scaffold_parameter_%d' % (i)), os.path.join(config['workspace'], 'position_lp_%d.math' % (i)), os.path.join(config['workspace'], 'edge_cluster_pos_%d' % (i)))
    command_run(ARCS_CMD, 'solveLP', args, config)

    ###################################################################
    #
    # arcs remove_repeats
    #
    ###################################################################
    args = '-d %s -K %d -O %d -i %d' % (config['workspace'], config['kmer_size'], config['max_overlap'], i)
    command_run(ARCS_CMD, 'remove_repeats', args, config)
    
###################################################################
#
# python reverse_filter.py
#
###################################################################
args = '%s %s %s' % (config['workspace'], os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_%d' % len(config['library_list'])))
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'reverse_filter.py'), args, config)

###################################################################
#
# arcs gapfill
#
###################################################################
args = '-s %s -d %s -K %d -O %d -C %s -l %s -I %s -P %s' % (os.path.join(config['workspace'], 'scaffold_parameter_0'), config['workspace'], config['kmer_size'], config['kmer_size'] - 10, os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_last'), os.path.join(config['workspace'], 'condensed_de_bruijn_graph_before_trimming.data'), os.path.join(config['workspace'], 'scaf.after.pair'))
command_run(ARCS_CMD, 'gapfill', args, config)
'''

###################################################################
#
# python giveContigsNumber.py
#
###################################################################
args = '%s %s %s' % (os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'contigs.fasta'), os.path.join(config['workspace'], 'contig_len'))
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'giveContigsNumber.py'), args, config)

###################################################################
#
# run reversePacBioRead.py 
#
###################################################################
args = '%s %s' % (config['pb'], os.path.join(config['workspace'], 're_pacbio_read.fasta'));
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'SPAReversePacBioRead.py'), args, config)

###################################################################
#
# run blasr 
#
###################################################################

args = '%s %s -out %s -m 1 -nproc %d' % (os.path.join(config['workspace'], 're_pacbio_read.fasta'), os.path.join(config['workspace'], 'contigs.fasta'), os.path.join(config['workspace'], 'blasrResult'), 8)
command_run('', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'blasr'), args, config)

###################################################################
#
# python tranBlasr2LPAndGraph.py
#
###################################################################

args = '%s %s %s %s %s %s %s %s' % (os.path.join(config['workspace'], 'contig_len'), os.path.join(config['workspace'], 'component_0'), os.path.join(config['workspace'], 'blasrResult'), os.path.join(config['workspace'], 'contig_arc_graph_after_remove_ambigous_arcs_%d' % i), os.path.join(config['workspace'], 'position_lp_%d.math' % i), os.path.join(config['workspace'], 'edge_cluster_len_%d' % i), os.path.join(config['workspace'], 'parameter_%d' % (i)),os.path.join(config['workspace'], 'map_detail'))
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'tranBlasr2LPAndGraph.py'), args, config)

###################################################################
#
# arcs solveLP
#
###################################################################
args = '-s %s -i %s -o %s' % (os.path.join(config['workspace'], 'parameter_%d' % (i)), os.path.join(config['workspace'], 'position_lp_%d.math' % (i)), os.path.join(config['workspace'], 'edge_cluster_pos_%d' % (i)))
command_run(ARCS_CMD, 'solveLP', args, config)

###################################################################
#
# arcs remove_repeats
#
###################################################################
args = '-d %s -K %d -O %d -i %d' % (config['workspace'], config['kmer_size'], config['max_overlap'], i)
command_run(ARCS_CMD, 'remove_repeats', args, config)

###################################################################
#
# python reverse_filter.py
#
###################################################################
args = '%s %s %s %s' % (config['workspace'], os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_%d' % (i+1)), os.path.join(config['workspace'], 'component_%d' % (i+2)))
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'reverse_filter.py'), args, config)


###################################################################
#
# arcs gapfill by gapFillbyScaf.py
#
###################################################################
args = '%s %s %s %s %s' % (os.path.join(config['workspace'], 'component_%d' % (i+2)), os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 're_pacbio_read.fasta'), os.path.join(config['workspace'], 'map_detail'), os.path.join(config['workspace'], 'gap_filling_after_3rd'))
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'gapFillbyScaf.py'), args, config)


'''
###################################################################
#
# arcs gapfill
#
###################################################################
# notice that component_(i+1) not component_(i+2) resaon see youdao 22,June
args = '-s %s -d %s -K %d -O %d -C %s -l %s -I %s -f %s' % (os.path.join(config['workspace'], 'scaffold_parameter_0'), config['workspace'], config['kmer_size'], config['kmer_size'] - 10, os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_%d' % (i+2)), os.path.join(config['workspace'], 'condensed_de_bruijn_graph_before_trimming.data'), os.path.join(config['workspace'], 'scaf_after_3rd'))
command_run(ARCS_CMD, 'gapfill', args, config)
'''

'''
i = i+1
#######
#unfinish scaffold again
#######
args = '%s %s %s %s %s %s %s' % (os.path.join(config['workspace'], 'contig_len'), os.path.join(config['workspace'], 'component_%d' % i), os.path.join(config['workspace'], 'blasrResult'), os.path.join(config['workspace'], 'contig_arc_graph_after_remove_ambigous_arcs_%d' % i), os.path.join(config['workspace'], 'position_lp_%d.math' % i), os.path.join(config['workspace'], 'edge_cluster_len_%d' % i), os.path.join(config['workspace'], 'scaffold_parameter_%d' % (i)))
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'tranBlasr2LPAndGraph.py'), args, config)

###################################################################
#
# arcs solveLP
#
###################################################################
args = '-s %s -i %s -o %s' % (os.path.join(config['workspace'], 'scaffold_parameter_%d' % (i)), os.path.join(config['workspace'], 'position_lp_%d.math' % (i)), os.path.join(config['workspace'], 'edge_cluster_pos_%d' % (i)))
command_run(ARCS_CMD, 'solveLP', args, config)

#i=0
###################################################################
#
# arcs remove_repeats
#
###################################################################
args = '-d %s -K %d -O %d -i %d' % (config['workspace'], config['kmer_size'], config['max_overlap'], i)
command_run(ARCS_CMD, 'remove_repeats', args, config)

###################################################################
#
# python reverse_filter.py
#
###################################################################
args = '%s %s %s' % (config['workspace'], os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_%d' % (i+1)))
command_run('python', os.path.join(os.path.abspath(os.path.dirname(__file__)), 'reverse_filter.py'), args, config)

###################################################################
#
# arcs gapfill
#
###################################################################

args = '-s %s -d %s -K %d -O %d -C %s -l %s -I %s' % (os.path.join(config['workspace'], 'scaffold_parameter'), config['workspace'], config['kmer_size'], config['kmer_size'] - 10, os.path.join(config['workspace'], 'cdbg_copy_number.fa'), os.path.join(config['workspace'], 'component_last'), os.path.join(config['workspace'], 'condensed_de_bruijn_graph_before_trimming.data'))
command_run(ARCS_CMD, 'gapfill', args, config)
'''
end = datetime.now()
print 'total running time is %d seconds' % ((end - start).seconds)



