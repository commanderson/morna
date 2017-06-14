#!/usr/bin/env python
"""
create_masterfile

given a sorted list of read names to retain, and a sorted sam or bam file to 
filter, output a still-sorted sam file with all reads not appearing in the list 
removed to stdout

"""
import argparse
import gzip
import sys
import subprocess

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
    parser.add_argument('-r', '--readnames', metavar='<file>', type=str,
            required=True,
            help=('path to list of sorted readnames')
        )
    parser.add_argument('-b', '--bamfile', metavar='<file>', type=str,
            required=True,
            help=('path to SORTED .bam file you want filtered')
        )
    args = parser.parse_args()

    process1 = subprocess.Popen(["samtools","view",
                                 args.bamfile], stdout=subprocess.PIPE)
                                 
    header_process = subprocess.Popen(["samtools","view","-H",
                                 args.bamfile], stdout=subprocess.PIPE)
    #must reproduce header
    for header_line in iter(header_process.stdout.readline, ''):
        sys.stdout.write(header_line)
            
    with open(args.readnames) as keep_readnames:
        last_line_retained = False
        current_readname = keep_readnames.readline().strip()
        lines_to_keep = []
        bam_iterator = iter(process1.stdout.readline, '')
        current_bamline = next(bam_iterator)
        
        while current_readname and current_bamline:
            #print current_readname + "||" + current_bamline.split()[0]
            #if this is a keeper
            if current_bamline.split()[0] == current_readname:
                #add current bamline to kept
                lines_to_keep.append(current_bamline)
                #mark that we retained a line
                last_line_retained = True
                #and move forward in bam
                try:
                    current_bamline = next(bam_iterator)
                except StopIteration:
                    current_bamline = ""
            #Otherwise,
            #if we are just finishing a run of keepers
            elif last_line_retained == True:
                #output our kept lines and empty them
                for keeper in lines_to_keep:
                    sys.stdout.write(keeper)
                lines_to_keep = []
                #mark that we had a mismatch
                last_line_retained = False
                #and move forward in keep list
                current_readname = keep_readnames.readline().strip()
            #if we're a miss and last comparison was a miss
            else:
            #keep moving forward in bamline
                try:
                    current_bamline = next(bam_iterator)
                except StopIteration:
                    current_bamline = ""
                    
        #one final flush of retained lines
        for keeper in lines_to_keep:
                    sys.stdout.write(keeper)
