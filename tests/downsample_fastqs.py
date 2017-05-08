#!/usr/bin/env python
"""
downsample_fastqs

beginning with a set of 2 paired-end fastq files, downsample randomly to a given number and output new fastq files which are thus downsampled

"""
import argparse
import gzip
#import os
#import re
import sys
from itertools import takewhile,repeat
from random import sample,seed,shuffle
#from collections import defaultdict
#from math import log, sqrt, ceil

_help_intro = \
"""downsample_fastqs takes a pair of paired-end fastq files and 
downsamples the reads to a requested number
"""

def gzipped_file_linecount(filename):
    f = gzip.open(filename, 'r')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen if buf )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_help_intro)
    parser.add_argument('-1', '--left', metavar='<file>', type=str,
            required=True,
            help=('path to gzipped fastq file for left reads')
        )
    parser.add_argument('-2', '--right', metavar='<file>', type=str,
            required=True,
            help=('path to gzipped fastq file for right reads')
        )
    parser.add_argument('-o', '--original-reads', metavar='<int>',
            type=int,
            required=False,
            default=None,
            help=('number of paired-end reads in the initial fastq files; '
                  ' will be calculated if necessary')
        )
    parser.add_argument('-d', '--downsamples', metavar='<int>',
            type=int,
            required=True,
            help=('number of reads to include in downsampled files; must be '
                  'less than the number of original reads')
        )
    parser.add_argument('-s', '--seed', metavar='<str>', type=str,
            required=False,
            default="867-5309",
            help=('Seed the random module for reproducibility')
        )
    
    args = parser.parse_args()
    seed(args.seed)
    if args.original_reads is None:
        print("Counting reads in " + args.left + " since it was not specified.")
        file_lines = gzipped_file_linecount(args.left)
        if not file_lines % 4 == 0:
            print("It's concerning that file lines isn't divisble by 4, " 
                                    + "since each read should be 4 lines")
        args.original_reads = file_lines/4
        print("found file linecount of " + str(file_lines) 
                + " which means " + str(args.original_reads) + " reads.")
    #To construct a list of reads to retain in our downsampled files,
    #We randomly sample downsides indexes 
    #from the range 0 to original_reads-1,
    #Then sort these so we can use them as we loop
    #through the input file
    print("Generating subsample index list")
    keeplist = sorted(sample(range(args.original_reads),
                                    args.downsamples))
    print("List generated, beginning downsampling")
    #print keeplist
                    
    with gzip.open(args.left) as leftin,\
            open(args.left + ".downsamp", "w") as leftout:
        index = 0
        for i,line in enumerate(leftin):
            if i%10000 == 0:
                sys.stdout.write("Made it to line " + str(i) 
                                    + " in left file.\r")
                #sys.stdout.flush()
            if i/4 == keeplist[index]:
                sys.stdout.write("\nkeeping line " + str(i) 
                    + " since keeplist[index] is "+ str(keeplist[index]) +"\r")
                leftout.write(line)
                if (i-3)%4 == 0:
                    try:
                        index += 1
                        keeplist[index[
                    except IndexError:
                        #if the list of junction to match is empty, 
                        #stop looping
                        break
            
    sys.stdout.write("\n")
    print("Finished left file, begin right file")

    with gzip.open(args.right) as rightin,\
         open(args.right + ".downsamp", "w") as rightout:
        index = 0
        for i,line in enumerate(rightin):
            if i%10000 == 0:
                sys.stdout.write("Made it to line " + str(i) 
                                    + " in right file.\r")
                #sys.stdout.flush()
            if i/4 == keeplist[index]:
                rightout.write(line)
                if (i-3)%4 == 0:
                    try:
                        index += 1
                        keeplist[index]
                    except IndexError:
                        #if the list of junction to match is empty, 
                        #stop looping
                        break