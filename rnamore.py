#!/usr/bin/env python
"""
rnamore

Uses mmh3 to implement feature hashing for dimensionality reduction of 
intropolis junction vectors by sample and uses Spotify's annoy library for
obtaining nearest neighbors.

Requires mmh3 (pip install mmh3) and Spotify's annoy (pip install annoy)
"""
import argparse
import mmh3
from annoy import AnnoyIndex

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--intropolis', type=str, required=False,
            default=None,
            help='path to gzipped intropolis; specify only when indexing'
        )
    parser.add_argument('-x', '--annoy-idx', type=str, required=False,
            default='jx.idx',
            help='path to junction index for either writing or reading'
        )
    parser.add_argument('--features', type=int, required=False,
            default=3000,
            help='dimension of feature space'
        )
    parser.add_argument('--n-trees', type=int, required=False,
            default=200,
            help='number of annoy trees'
        )
    parser.add_argument('--search-k', type=int, required=False,
            default=None,
            help='larger = more accurage search'
        )
    args = parser.parse_args()
    if args.intropolis:
        # Index
        
        with open(agrs.intropolis) as introp_file_handle:
        for i, introp_line in enumerate(junc_file_handle):
            if (i % 100 == 0):
                print( str(i) + " lines into 2nd pass (tf-idf) writing")
            introp_line_pieces = introp_line.split()
            new_line_pieces = introp_line_pieces[:6]
        
            #This time we will write the non tf-idf score line pieces first
            tf_idf_file_handle.write("\t".join(new_line_pieces))
        
            samples_junction_coverages = (introp_line_pieces[7].split(','))
            samples_junction_coverages = [int(num) for num in
                                     samples_junction_coverages]
                                     
            num_samples_with_junction = len(samples_junction_coverages)
            idf_value = norm_log(1.0 + (num_samples/num_samples_with_junction))
        
    else:
        # Search