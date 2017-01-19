#!/usr/bin/env python

import re
import argparse
import random
from math import log
from itertools import izip
#import linecache

def uniquify(seq):
    """Returns a version of list seq with duplicate elements reduced to unique
        seq: a list
        Return value: a list
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]
    
def norm_log(input):
    """ Returns natural log of input unless input is 0;
        then it returns 0.0
        input: integer or float input
        Return value: a float
    """
    if (input == 0):
        return 0.0
    else:
        return log(float(input))
    
def tf_idf_score(frequency_matrix):
    """ Produces a sample-by-junction tf-idf matrix,
        given a sample-by-junction raw frequency matrix
        frequency_matrix: a matrix where each row represents a unique exon-exon
        junction, and each column corresponds to a different sample; cells 
        contain raw coverage number of reads spanning that junction 
        in the corresponding sample
        Return value: a matrix of the same dimensions contining tf-idf scores
        calculated using junctions as terms and samples as documents.
    """
    #As of now, the term frequency measure is 1+log(junction coverage),
    #and the inverse document frequency is:
    # log(1+(number of samples/number of samples where junction is present)
    num_samples = float(len(frequency_matrix[0]))
    junction_present_counts = [0 for x in range(len(frequency_matrix)+1)]
    for i,row in enumerate(frequency_matrix):
        samples_with_junction = len([x for x in row if x > 0])
#        print samples_with_junction
#        print row
        for j,frequency in enumerate(row):
            frequency_matrix[i][j] = ((1+ norm_log(frequency_matrix[i][j])) *
                            norm_log(1.0 + (num_samples/samples_with_junction)))        
    return frequency_matrix
    
def match(lst, value):
    """ Produces the index in a UNIQUE list whose contents matches a value
        lst: list to search within
        value: the value to match
        Return value: list of ints, error string if not found or multiple found
    """
    matches = [i for i, x in enumerate(lst) if x==value]
    if (len(matches) < 1):
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

#We must calculate the number of samples to do IDF calculations
all_samples = []

with open(junctions_file) as junc_file_handle:
    for i, line in enumerate(junc_file_handle):
        if (i % 100 == 0):
            print( str(i) + " lines into sample count")
        line_pieces = line.split()
        #id_string = '_'.join(line_pieces[0:4])
        #junctions.append(id_string)
        samples_with_junction = (line_pieces[6].split(','))
        samples_with_junction = [int(num) for num in samples_with_junction]
        for sample in samples_with_junction:
            all_samples.append(sample)

unique_samples = uniquify(all_samples)
num_samples = len(unique_samples)

##################################
#pass1: create intropolis of IDFs#
##################################

idf_file = junctions_file + ".idf"

with open(junctions_file) as junc_file_handle, \
     open(idf_file, 'w') as idf_file_handle:
    for i, line in enumerate(junc_file_handle):
        if (i % 100 == 0):
            print( str(i) + " lines into 1st pass (idf) writing")
        line_pieces = line.split()
        samples_with_junction = (line_pieces[6].split(','))
        num_samples_with_junction = len(samples_with_junction)
        new_line_pieces = line_pieces[:6]
        idf = norm_log(1.0 + (num_samples/num_samples_with_junction))
        new_line_pieces.append(str(idf)+"\n")
        idf_file_handle.write("\t".join(new_line_pieces))
            
###########################################
#pass2: create intropolis of tf-idf scores#
###########################################

tf_idf_file = junctions_file + ".tf_idf"

with open(junctions_file) as junc_file_handle, \
     open(idf_file) as idf_file_handle, \
     open(tf_idf_file, 'w') as tf_idf_file_handle:
    i=0;    #since izip is already playing with iterators I'll track iterations
    for introp_line, idf_line in izip(junc_file_handle,
                                                idf_file_handle):
        if (i % 100 == 0):
            print( str(i) + " lines into 2nd pass (tf-idf) writing")
        introp_line_pieces = introp_line.split()
        idf_line_pieces = introp_line.split()
        new_line_pieces = introp_line_pieces[:6]
        idf_value = float(idf_line_pieces[6])
        #This time we will write the non tf-idf score line pieces first
        tf_idf_file_handle.write("\t".join(new_line_pieces))
        
        samples_with_junction = (introp_line_pieces[6].split(','))
        samples_with_junction = [int(num) for num in samples_with_junction]
        
        samples_junction_coverages = (introp_line_pieces[7].split(','))
        samples_junction_coverages = [int(num) for num in
                                     samples_junction_coverages]
        for j, coverage in enumerate(samples_junction_coverages):
            tf_idf_score = (1+ norm_log(coverage) 
                            * idf_value)
            tf_idf_file_handle.write("\t" + tf_idf_score)
        
        tf_idf_file_handle.write("\n")
        i+=1    #keep this at the end of the loop to track iteration!
 
quit()
#Next, we'll write a file to hold Inverse Document Frequency info
#Keeping the matrix at minimal size is tricky. 
#Because all sample indices might not be present in a continuous range, we treat
#them more like ids to be looked up later. The rows of the matrix correspond 1:1 
#with the lines (each representing a junction) in the input file, and there is 
#one column for each unique sample index, but there is no guarantee that e.g. 
#sample index 12 corresponds to column 12 of the matrix; instead, the columns 
#correspond to the matching sample index from the set of all indices sorted in 
#increasing order - column 0 being the lowest sample index, column 9 the 
#10th-lowest, and so on.

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
        print("junction "+ str(i) + " in " + str(len(samples_with_junction)) 
                + "samples.")
        for j, sample in enumerate(samples_with_junction):
#            print ("junction "+str(i) + " in sample" 
#                    + str(samples_with_junction[j]))
            sample_index = match(all_samples, sample)
#            print("store that in column " + str(sample_index))
            tf_idf_matrix[i][sample_index] = samples_junction_coverages[j]
#            tf_idf_matrix[i,j] = tf_idf_score(junctions[i], samples[j], 
#                                            sum_of_frequencies)
#print tf_idf_matrix

tf_idf_matrix = tf_idf_score(tf_idf_matrix)
#for row in tf_idf_matrix:
#    print row

outfile = workfile+".deduped"
fho=open(outfile, "w")
fho.write(orig_header)