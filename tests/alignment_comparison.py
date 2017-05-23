#!/usr/bin/env python
"""
alignment_comparison

compare a list of .bam alignment files with a master .bam alignment
to determine what percentage of reads in each listed file have an identical
mapping (coordinate + cigar string) in the master file
"""
import argparse
import subprocess
import sys
from itertools import groupby

_help_intro = \
"""alignment_comparison determines what portion of reads in each of a list of bam files have an identical mapping (coordinate + cigar string) present in a
master bam file.
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
    
    args = parser.parse_args()
    
    read_dict = {}
    process1 = subprocess.Popen(["samtools","view","-F","4",
                                 args.master], stdout=subprocess.PIPE)
    for name, group in groupby(iter(process1.stdout.readline, ''),
                                             lambda x: x.split()[0]):
        max_alignment_score=0
        for read in group:
            #Is it a max scoring alignment? 
            #If replaces current, new list
            #if equals current, add to list
            #else do nothing
            print read.split()[:3]
        print "-------------------------"

    for i,line in enumerate(iter(process1.stdout.readline, '')):
        if (i%100 == 0):
            sys.stdout.write(str(i) + " lines into master file processing\r")
        fields = line.split()
        #name = (fields[0],fields[1])
        #rep = (fields[2],fields[3],fields[5])
        #The "key" which uniquely represents a read is the SAM fields 
        #QNAME (query template name)
        #FLAG (bitwise flag about the mapping)
        #RNAME (reference sequence name)
        #A given combination of these uniquely represents a specific read
        #The mapping of that read is represented in the SAM fields
        #POS (1-based leftmost mapping position of read)
        #CIGAR (details of exact base mapping)
        #If two reads have identical POS and CIGAR fields, 
        #they have mapped identically.
        read_dict[(fields[0],fields[1],fields[2])] = (fields[3],fields[5])
        
    for file in args.subfiles:
        sys.stdout.write("\n")
        process2 = subprocess.Popen(["samtools","view","-F","4",
                                     args.master], stdout=subprocess.PIPE)
        num_mapped_reads = 0
        num_identical_mappings = 0
        for i,line in enumerate(iter(process2.stdout.readline, '')):
            if (i%100 == 0):
                sys.stdout.write(str(i) + " lines into " + file 
                                                        + " processing\r")
            num_mapped_reads += 1
            fields = line.split()
            if read_dict[(fields[0],fields[1],fields[2])] == (fields[3],
                                                                fields[5]):
                num_identical_mappings += 1
        print i
        print("For file " + file + " there are " + str(num_mapped_reads) + " mapped reads, " + str(num_identical_mappings) + " identical (" + str(float(num_identical_mappings)/num_mapped_reads) + ")")
    