#!/usr/bin/env python
"""
pool_junction_lists

Parses a SORTED combined junction list with final column 
Combines junction data from input list of intropolis-like files to 
create a single pooled junctions file from which 
a master morna index can be made
"""


import argparse
import gzip
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='<file>', type=str,
            required=True,
            help='path to gzipped input file of combined junctions'
            'eg: sorted_cjs.tsv.gz'
        )
parser.add_argument('-s','--sources', nargs='+', metavar='<files>', type=str,
            required=True,
            help=('list of sources which must exactly match the last column '
                  'values of the input file; this helps align the sample ids'
                  'eg: tcga sra gtex')
        )
parser.add_argument('-m','--metafiles', nargs='+', metavar='<files>', type=str,
            required=False,
            default=None,
            help=('list of metadata files matching up with sources in order;'
                  'must have tab separated first field with sample id, then'
                  'remainder of line is any and all metadata that should be'
                  'associated with that sample id, such as accession number'
                  'and other useful details to include with search results.'
                  'If included, output will include a combined metadata file'
                  'named source1_source2..._metadata.tsv ')
        )
parser.add_argument('-n', '--increment', metavar='<int>', type=int,
            required=True,
            help='since each input file has its own integer sample ids '
            'which start from 0, increment the ids in each subsequent file '
            'by this much for each file it is after the first, which must '
            'be made enough to separate the sample id ranges entirely'
            'eg: since TCGA has 11349 samples, GTEX 9662 and SRA 20854,'
            '25000 is an appropriate and readable increment'
        )
parser.add_argument('-o', '--outfile', metavar='<file>', type=str,
            required=True,
            help='path to output file, which will be gzipped'
            'eg: pooled_junctions.tsv.gz'
        )


#junction_line_offsets = {}
junction_presence = {}
stored_tokens=""

if __name__ == '__main__':
    args = parser.parse_args()
    
    #this source offset dictionary is used to alter sample ids
    source_offset = {}
    for i,source in enumerate(args.sources):
        source_offset[source] = i * args.increment
        
    with gzip.open(args.input) as input_handle,\
     gzip.open(args.outfile, "w") as out_handle:
        for i, line in enumerate(input_handle):
            if (i % 1000 == 0):
                sys.stdout.write(str(i) 
                    + " lines into combined file " + args.input + "\r")
                sys.stdout.flush()
            tokens = line.strip().split("\t")
            #we'll have to adjust the sample ids by an appropriate offset
            offset = source_offset[tokens[-1]]
            #so break them out into numbers, add offset, and back to string!
            tokens[-3] = ",".join(str(int(x) + offset) 
                    for x in tokens[-3].split(","))
                    
            #if this line isn't the same junction as the stored, write stored
            if stored_tokens[:3] != tokens[:3]:
                #write accumulated junction line, leaving out source indicator
                if i>0:
                    out_handle.write("\t".join(stored_tokens[:-1])+"\n") 
                #and move current tokens into stored
                stored_tokens = tokens
                
            #otherwise, this is a continuation of a junction from another source
            else:
                #add newly offset sample ids and coverages
                stored_tokens[-3] = stored_tokens[-3] + "," + tokens[-3]
                stored_tokens[-2] = stored_tokens[-2] + "," + tokens[-2]
        #finish off with final outfile write
        out_handle.write("\t".join(stored_tokens[:-1])+"\n")

if len(args.metafiles) > 0:
    with open("_".join(args.sources) + "_metadata.tsv","w") as meta_out:
        for i,file in enumerate(args.metafiles):
            with open(file) as metafile:
                for line in metafile:
                    tokens = line.strip().split("\t")
                    tokens[0]= str(int(tokens[0]) + i * args.increment)
                    meta_out.write("\t".join(tokens) + "\n")