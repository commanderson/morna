#!/usr/bin/env python

import argparse
import gzip
import random
import sqlite3
import sys
import tempfile
from math import log
from itertools import izip

def magic_sql_write(cursor, sample_id, tf_idf_score):
    """ Wrapper function for writing to our temp sql db
            Note: the sample id should have a table already!
        cursor: cursor object for writing to the db
        sample_id: the sample id we are writing for
        tf_idf_score: the score we are recording. Note: 0.0 is empty
        Return value: null
    """
    #since everything's been converted to numeric types we're assuming
    #no problems with injection;
    #also, standard float precision warnings apply.
    cursor.execute("INSERT INTO sample_%d VALUES (%f)"
                    % (sample_id,tf_idf_score))

def dot_product(v1,v2):
    """returns the dot product of two lists
        v1: first vector, of equal length to v2
        v2: second vector, of equal length to v1
        Return value: float equal to
        v1[0]*v2[0] + v1[1]*v2[1] + ... v1[n]*v2[n]
    """
    return sum([i*j for (i, j) in zip(v1, v2)])
    
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
parser.add_argument('-f', '--file', metavar='<str>', type=str,
            required=True,
            help=('path to gzipped file recording junctions across samples '
                  'in intropolis format')
        )
parser.add_argument('-v', '--verbose', action='store_const',
            const=True,
            default=False,
            help='be talkative'
        )
parser.add_argument('-o', '--output', metavar='<str>', type=str,
            required=False,
            default="tf_idf_temp",
            help='output file basename'
        )

args = parser.parse_args()
junctions_file = args.file

#Next we need a list of all the samples; we'll go by indices now
#and look up in the idmap file later.

#We must determine the number of samples to do IDF calculations
all_samples = set()

#########################
#pass1: count of samples#
#########################

with gzip.open(args.file) as introp_file_handle:
    if args.verbose:
        for i, line in enumerate(introp_file_handle):
            if i % 100 == 0:
                sys.stdout.write(
                    str(i) + " lines into sample count, " 
                    + str(len(all_samples)) + " samples so far.\r"
                )
            all_samples.update([int(num) for num 
                                in(line.split('\t')[-2].split(','))])
    else:
        for i, line in enumerate(introp_file_handle):
            all_samples.update([int(num) for num 
                                in(line.split('\t')[-2].split(','))])
num_samples = len(all_samples)
if args.verbose:
    print '\nThere are {} samples.'.format(num_samples)
    print(all_samples)
###########################################
#pass2: create intropolis of tf-idf scores#
###########################################

tf_idf_file = args.file + ".tf_idf"

with gzip.open(args.file) as junc_file_handle, open(
        tf_idf_file, 'w') as tf_idf_file_handle:
    if args.verbose:
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
            idf_value = log((float(num_samples)/num_samples_with_junction))
        
            first_one = True    #must format commas appropriately!
            for j, coverage in enumerate(samples_junction_coverages):
                tf_idf_score = ((coverage) * idf_value)
                if (first_one):   
                    tf_idf_file_handle.write("\t" + str(tf_idf_score))
                    first_one = False
                else:
                    tf_idf_file_handle.write("," + str(tf_idf_score))
            tf_idf_file_handle.write("\n")
    else:
        for i, introp_line in enumerate(junc_file_handle):
            introp_line_pieces = introp_line.split()
            new_line_pieces = introp_line_pieces[:6]
        
            #This time we will write the non tf-idf score line pieces first
            tf_idf_file_handle.write("\t".join(new_line_pieces))
        
            samples_junction_coverages = (introp_line_pieces[7].split(','))
            samples_junction_coverages = [int(num) for num in
                                         samples_junction_coverages]
                                     
            num_samples_with_junction = len(samples_junction_coverages)
            idf_value = log((float(num_samples + 1)/num_samples_with_junction))
        
            first_one = True    #must format commas appropriately!
            for j, coverage in enumerate(samples_junction_coverages):
                tf_idf_score = ((coverage) * idf_value)
                if (first_one):   
                    tf_idf_file_handle.write("\t" + str(tf_idf_score))
                    first_one = False
                else:
                    tf_idf_file_handle.write("," + str(tf_idf_score))
            tf_idf_file_handle.write("\n")
        
chosen_sample_ids = random.sample(all_samples,3)
print("randomly chosen sample ids are: " + str(chosen_sample_ids))

# Each chosen sample id will get its own table in a sqlite db
#We will treat these like lists onto which we can 
#add a tf-idf score for each junction.

if args.verbose:
    print("Now writing " + str(len(all_samples)) + " sample tables.")
##TODO: Make the db files temporary files
conn = sqlite3.connect(args.output + '.db')
cursor = conn.cursor()
if args.verbose:
    for i,sample_id in enumerate(all_samples):
        print("Creating table number " + str(i) + ", sample " + str(sample_id))
        #I don't think we worry about injection, but if we do I'm
        #currently assuming our int conversion prevents it 
        cursor.execute("CREATE TABLE sample_%d (tf_ifd float)" 
                        % sample_id)
else:
    for i,sample_id in enumerate(all_samples):
        cursor.execute("CREATE TABLE sample_%d (tf_ifd float)" 
                        % sample_id)

with gzip.open(args.file) as junc_file_handle, \
     open(tf_idf_file) as tf_idf_file_handle:
    i = 0;
    if args.verbose:
        for introp_line, score_line in izip(junc_file_handle, 
                                            tf_idf_file_handle):
            if (i % 100 == 0):
                print( str(i) + " lines into 3rd pass (inversion into sql)")
            samples_with_junction = ((introp_line.split())[6]).split(',')
            samples_with_junction = [int(num) for num in samples_with_junction]
        
            tf_idfs_of_samples = ((score_line.split())[6]).split(',')
            tf_idfs_of_samples = [float(num) for num in tf_idfs_of_samples]
        
            for j, sample_id in enumerate(all_samples):
                if sample_id in samples_with_junction:
                    #print("SPECIAL")
                    idx = match(samples_with_junction, sample_id)
                    magic_sql_write(cursor, sample_id,
                     tf_idfs_of_samples[
                                        match(samples_with_junction, sample_id)
                                        ]
                                    )
                else:
                    #print("DEFAULT")
                    magic_sql_write(cursor, sample_id, 0.0)
            i+=1
    else:
        for introp_line, score_line in izip(junc_file_handle, 
                                            tf_idf_file_handle):
            samples_with_junction = ((introp_line.split())[6]).split(',')
            samples_with_junction = [int(num) for num in samples_with_junction]
        
            tf_idfs_of_samples = ((score_line.split())[6]).split(',')
            tf_idfs_of_samples = [float(num) for num in tf_idfs_of_samples]
        
                
            for j, sample_id in enumerate(all_samples):
                if sample_id in samples_with_junction:
                    #print("SPECIAL")
                    idx = match(samples_with_junction, sample_id)
                    magic_sql_write(cursor, sample_id,
                     tf_idfs_of_samples[
                                        match(samples_with_junction, sample_id)
                                        ]
                                    )
                else:
                    #print("DEFAULT")
                    magic_sql_write(cursor, sample_id, 0.0)
            i+=1
conn.commit()

for sample_id in all_samples:
    print("Contents of sample " + str(sample_id))
    for row in cursor.execute(
    #"SELECT * FROM sample_%d ORDER BY ROWID ASC LIMIT 1" %sample_id):
    "SELECT * FROM sample_%d ORDER BY ROWID ASC" %sample_id):
        print row
conn.close()        
#for sample_id in chosen_sample_ids:
#    cursor.execute("DROP TABLE sample_%d" %sample_id)
    
