#!/usr/bin/env python
"""
alignment_comparison

compare a list of .bam alignment files with a master .bam alignment
to determine what percentage of reads in each listed file have an identical
highest-scoring mapping (coordinate + cigar string) in the master file

Note the concept of fractional mappings - if master file has 4 equally high scoring mappings, and 3 of these are present in a compared file, then .75 of an identical mapping is present for that read in that file.

"""
import argparse
import subprocess
import sys
from itertools import groupby

_help_intro = \
"""alignment_comparison determines what portion of reads in each of a list of bam files have an identical mapping (coordinate + cigar string) present in a
master bam file. 
Requires that reads be sorted by read name aka the QNAME field-
use samtools sort -n to get these
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_help_intro)
    
    parser.add_argument('-m', '--master', metavar='<file>', type=str,
            required=True,
            help=('path to .bam alignment file which serves as master')
        )
    parser.add_argument('-s', '--subfiles', metavar='<m1>',
            type=str, nargs='+',
            required=False,
            help='paths to .bam alignment files to evaluate'
        )
    parser.add_argument('-v', '--verbose', action='store_const',
            const=True,
            default=False,
            help='be talkative'
        )
        
    args = parser.parse_args()
    
    read_dict = {}
    process1 = subprocess.Popen(["samtools","view","-F","4",
                                 args.master], stdout=subprocess.PIPE)
                                 
    sub_iterators = []
    current_read_groups = []
    for file in args.subfiles:
        process_sub = subprocess.Popen(["samtools","view","-F","4",
                                     file], stdout=subprocess.PIPE)
        sub_iterators.append(groupby(iter(process_sub.stdout.readline, ''),
                                             lambda x: x.split()[0]))
    for iterator in sub_iterators:
        current_read_groups.append(next(iterator))
        
    for name, group in groupby(iter(process1.stdout.readline, ''),
                                             lambda x: x.split()[0]):
        max_alignment_score=-99999999999
        alignments_list = []
        for alignment in group:
            fields = alignment.split()
            score = int(fields[11].strip("AS:i"))
            #Is it a max scoring alignment?
            if score > max_alignment_score:
                max_alignment_score = score
                alignments_list = [(fields[0],fields[3],fields[5],)]
            #if equals current, add to list
            elif score == max_alignment_score:
                alignments_list.append((fields[0],fields[3],fields[5],))
            #else do nothing
        
        if args.verbose:
            print "--------------------------------------------------"
            sys.stdout.write("Alignments in this group:\n")
            for entry in alignments_list:
                sys.stdout.write(str(entry) + "\n")
            sys.stdout.write("Max score: " + str(max_alignment_score) + "\n")
        for i,file in enumerate(args.subfiles):
            #if the current group from that file has the same read name
            if current_read_groups[i][0] == name:
                #then we need to look at the alignments in that group in subfile
                #to see how many from our list are matched
                #TODO: FIND HIGHEST SCORING LIST HERE TOO?
                matching = 0.0
                for sub_alignment in current_read_groups[i][1]:
                    sub_fields = sub_alignment.split()
                    if (sub_fields[0],sub_fields[3],
                                        sub_fields[5],) in alignments_list:
                        matching += 1
                matching /= 
                    
                
        #for name, group in groupby(iter(process2.stdout.readline, ''),
        #                                     lambda x: x.split()[0]):
        