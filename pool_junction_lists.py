#!/usr/bin/env python
"""
pool_junction_lists

Combines junction data from input list of intropolis-like files to 
create a single pooled junctions file from which 
a master morna index can be made
"""


import argparse
import gzip
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-f','--files', nargs='+', metavar='<files>', type=str,
            required=True,
            help=('paths to gzipped files recording junctions across samples '
                  'in intropolis format, space separated')
        )
parser.add_argument('-o', '--outfile', metavar='<file>', type=str,
            required=True,
            help='path to output file'
        )

junction_line_offsets = {}

if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.outfile, "w+") as output_handle:
        current_offset = 0
        for numfile,file in enumerate(args.files):
            output_handle.seek(0)
            with gzip.open(file) as input_handle:
                for i, line in enumerate(input_handle):
                    if (i % 100 == 0):
                        sys.stdout.write(str(i) 
                            + " lines into file " + file + "\r")
                    tokens = line.strip().split("\t")
                    if numfile == 0: #always just spit first input into output
                        output_handle.write(line)
                        
                        #out.write(" ".join(
                        #                    [" ".join(tokens[:3]),
                        #                    tokens[-2],
                        #                    tokens[-1]]))
                        #out.write"\n"
                        junction_line_offsets[
                                    ' '.join(tokens[:3])] = current_offset
                        current_offset += len(line)
                        
                    else: #
                        try:
                            output_handle.seek(
                                junction_line_offsets[' '.join(tokens[:3])])
                            print("already have line for \n" 
                                        + ' '.join(tokens[:3]))
                            print("it's: \n" + output_handle.readline())
                        except KeyError:
                            output_handle.seek(current_offset)
                            output_handle.write(line)
                            junction_line_offsets[
                                    ' '.join(tokens[:3])] = current_offset
                            current_offset += len(line)
                
                print "-----------------------------"
