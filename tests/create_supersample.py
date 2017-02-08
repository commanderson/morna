#!/usr/bin/env python

import argparse
import gzip
import random
import sqlite3
import sys
import tempfile
from math import log
from itertools import izip

"""
create supersample

Processes an "intropolis-like" format junction-by-junction sample coverage file
and produces a single file with lines in format: 
{chromosome}\t{start_position}\t{end_position}\t{sum_coverage}
where each line represents a unique junction with its first three fields
and the sum total coverage of that junction among all samples
listed in an intput file 

Compatible with pypy for faster speed!
"""


parser = argparse.ArgumentParser()
parser.add_argument('-i','--intropolis', metavar='<str>', type=str,
            required=True,
            help=('path to gzipped file recording junctions across samples '
                  'in intropolis format')
        )
        
parser.add_argument('-s','--sampleids', metavar='<str>', type=str,
            required=True,
            help=('path to file listing sample ids that should be included')
        )
        
parser.add_argument('-o','--output', metavar='<str>', type=str,
            default = 'supersample.qry',
            help=('path to output file')
        )


args = parser.parse_args()
intropolis_like_file = args.intropolis
ids_file = args.sampleids
output_file = args.output

wanted_ids = []

#file of desired ids should be in very simple format,
#one sample index per line with nothing else present
#eg:
#10223
#13224
#14
#...
with open(ids_file) as ids_file_handle:
    for i,line in enumerate(ids_file_handle):
        wanted_ids.append(int(line))

#format for "intropolis-like file" described in construct_tf_idf.py
with gzip.open(intropolis_like_file) as introp_handle,\
     open(output_file,"w") as queryfile_handle:
    for i,line in enumerate(introp_handle):
        if (i % 100 == 0):
                    sys.stdout.write( str(i) + " lines into file\r")
        coverage = 0;
        line_pieces = line.split()
        
        samples_with_junction = (line_pieces[6].split(','))
        samples_with_junction = [int(num) for num in samples_with_junction]
        
        samples_junction_coverages = (line_pieces[7].split(','))
        samples_junction_coverages = [int(num) for num in
                                     samples_junction_coverages]
                                 
        for j, sample in enumerate(samples_with_junction):
            if sample in wanted_ids:
                coverage+=samples_junction_coverages[j]
    
        queryfile_handle.write("\t".join(line_pieces[:3]) + "\t" 
                                + str(coverage) + "\n")
print('')