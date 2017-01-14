#!/usr/bin/env python

import re
import argparse
import random
#import linecache

def tf_idf_score(term, document, corpus):
    """ Produces a tf-idf score for the given term in 
        question: string with question to be printed to console
        yes: automatically answers yes to question
        Return value: float
    """
    return (term_frequency(term, document) 
    * inverse_document_frequency(term, corpus))

#######################
#Main body begins here#
#######################

#Set seed for reproducible results!
random.seed(8675309)
#argparse expects argument f for filepath of junctions file
#junctions file in tab separated format with fields:
#1. chromosome
#2. intron start position (1-based; inclusive)
#3. intron end position (1-based; inclusive)
#4. strand (+ or -)
#5. donor dinucleotide (e.g., GT)
#6. acceptor dinucleotide (e.g., AG)
#7. comma-separated list of indices of samples in which junction was found
#8. comma-separated list of coverages for spanning reads
#    in the corresponding samples from field 7

#First set up parser and parse in input file's name.
parser = argparse.ArgumentParser()
parser.add_argument("-f", nargs="?")


args = parser.parse_args()
junctions_file = args.f

#Next we need a list of all the samples; we'll go by indices now
#and look up in the idmap file later.

#Initialize some storage spaces
samples = []
junctions = []

with open(junctions_file) as file_handle:
    for i, line in enumerate(file_handle):
        line_pieces = line.split()
        id_string = '_'.join(line_pieces[0:4])
        junctions.append(id_string)
        for index in line_pieces[6].split(','):
            if not index in samples:
                samples.append(index)
#print(sorted(samples))
#print i

tf_idf_matrix = ([[0 for x in range(len(samples))] 
                for x in range(i)])

#print tf_idf_matrix
with open(junctions_file) as file_handle:
    for i, line in enumerate(file_handle):
        line_pieces = line.split()
        for j, sample in enumerate(samples):##TODO: REAL WORK STARTS HERE
            tf_idf_matrix[i,j] = tf_idf_score(junctions[i], samples[j], 
                                            sum_of_frequencies)
        