#!/usr/bin/env python

import re
import argparse
import random
#import linecache

def tf_idf_score(term, document, corpus):
    """ Produces a tf-idf score for the given term in the document
        given the corpus
        term: in this case, an exon-exon junction
        document: a sample with coverage information for each junction
        corpus: all the samples taken together as a collection of coverages
        Return value: float
    """
    return (term_frequency(term, document) 
    * inverse_document_frequency(term, corpus))
    
def match(lst, value):
    """ Produces the index in a UNIQUE list whose contents matches a value
        lst: list to search within
        value: the value to match
        Return value: list of ints, error string if not found or multiple found
    """
    matches = [i for i, x in enumerate(lst) if x==value]
    if (len(matches) == 0):
        return "Err: " + str(value) + " not found in " + str(lst)
    elif (len(matches) > 1):
        return ("Err: " + str(value) + " occurs " + str(len(matches)) 
                + " times in" + str(lst) + "\nFound at indices: " 
                + str(matches))
    else:
        return matches[0]
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
all_samples = []
max_sample_index = []
#junctions = []

with open(junctions_file) as file_handle:
    for i, line in enumerate(file_handle):
        line_pieces = line.split()
        #id_string = '_'.join(line_pieces[0:4])
        #junctions.append(id_string)
        samples_with_junction = (line_pieces[6].split(','))
        samples_with_junction = [int(num) for num in samples_with_junction]
        for sample in samples_with_junction:
            if not sample in all_samples:
                all_samples.append(sample)

all_samples = sorted(all_samples)                
print(all_samples)
print i

#keeping the matrix at minimal size is tricky. 
#Because all sample indices might not be present in a continuous range, we treat #them more like ids to be looked up later. The rows of the matrix correspond 1:1 #with the lines (each representing a junction) in the input file, and there is #one column for each unique sample index, but there is no guarantee that e.g. #sample index 12 corresponds to column 12 of the matrix; instead, the columns #correspond to the matching sample index from the set of all indices sorted in #increasing order - column 0 being the lowest sample index, column 9 the #10th-lowest, and so on.

tf_idf_matrix = ([[0 for x in range(len(all_samples)+1)] 
                for x in range(i+1)])

#print tf_idf_matrix

print( "tf_idf matrix is "+str(len(tf_idf_matrix))+" rows by "
        + str(len(tf_idf_matrix[0]))+" columns")
        
with open(junctions_file) as file_handle:
    for i, line in enumerate(file_handle):
        line_pieces = line.split()
        samples_with_junction = (line_pieces[6].split(','))
        samples_with_junction = [int(num) for num in samples_with_junction]
        samples_junction_coverages = (line_pieces[7].split(','))
        samples_junction_coverages = [int(num) for num in
                                     samples_junction_coverages]
        
        for j, sample in enumerate(samples_with_junction):
            print ("junction "+str(i) + " in sample" 
                    + str(samples_with_junction[j]))
            sample_index = match(all_samples, sample)
            print("store that in column " + str(sample_index))
            tf_idf_matrix[i][sample_index] = samples_junction_coverages[j]
#            tf_idf_matrix[i,j] = tf_idf_score(junctions[i], samples[j], 
#                                            sum_of_frequencies)
print tf_idf_matrix