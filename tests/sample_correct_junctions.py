#!/usr/bin/env python
"""
sample_correct_junctions.py

manually calculates the exact junctions and coverages for a given sample, 
and verifies that the information retrieved from the provided .junc.mor 
junctions-by-sample morna db matches up.
"""

import argparse
import gzip
import re
import sqlite3
import sys
from BitVector import BitVector

parser = argparse.ArgumentParser()
parser.add_argument('-d','--database', metavar='<file>', type=str,
            required=True,
            help=('path to junction database file for which you wish'
                  'to test the sample')
        )
parser.add_argument('-s','--sample', metavar='<int>', type=int,
            required=True,
            help=('sample id of sample to compare')
        )
parser.add_argument('-i','--intropolis', metavar='<file>', type=str,
            required=True,
            help=('path to gzipped file recording junctions across samples '
                  'in intropolis format')
        )

args = parser.parse_args()
print("starting analysis of " + args.intropolis + " for sample " 
            + str(args.sample))


if True:
    introp_bv = BitVector(size = 0)
    introp_coverages = ""
    with gzip.open(args.intropolis) as introp_handle:
        for i,line in enumerate(introp_handle):
            if (i % 100 == 0):
                print "On junction " + str(i)
            line_pieces = line.split()
    
            samples_with_junction = (line_pieces[6].split(','))
            samples_with_junction = [int(num) for num in samples_with_junction]
    
            samples_junction_coverages = (line_pieces[7].split(','))
            samples_junction_coverages = [int(num) for num in
                                         samples_junction_coverages]
            if args.sample in samples_with_junction:
                introp_bv = introp_bv + BitVector(intVal = 1)
                idx = samples_with_junction.index(args.sample)
                introp_coverages += (str(samples_junction_coverages[idx]) + ',')
            else:
                introp_bv = introp_bv + BitVector(intVal = 0)

    introp_coverages = introp_coverages[:-1] #remove the trailing comma
    print "Intropolis Junctions"
    #print introp_bv           
    print "Length: " + str(len(introp_bv))
    print "Num. Present Junctions: " + str(sum(introp_bv))
    #print "Coverages: " + introp_coverages
    print("Length of coverages (should agree): " 
            + str(len(introp_coverages.split(','))))

conn=sqlite3.connect(args.database)
c=conn.cursor()

for line in c.execute(("SELECT * FROM sample_%d") % args.sample ):
    db_rle_juncs = line[0]
    db_coverages = line[1]

shr_str = db_rle_juncs
db_bv=BitVector(size=0)
while len(shr_str)>0:
    #print "starting with shr-str of:" +shr_str
    m = re.search("(.*)([io])(\d+)$",shr_str)
    if m.groups()[1] == 'i':
            bit='1'
    elif m.groups()[1] == 'o':
            bit='0'
    else: 
        raise ValueError("Found non i-o bit in parsing " + shr_str)
        
    digit = m.groups()[2]
    if not (m.groups()[2].isdigit()):
        raise ValueError("Found non digit runlength in parsing " + shr_str)
    db_bv = BitVector(bitstring=(int(m.groups()[2])*bit)) + db_bv
    shr_str = m.groups()[0]
    
print "DB Junctions"
#print db_bv
print "Length: " + str(len(db_bv))
print "Num. Present Junctions: " + str(sum(db_bv))
#print "Coverages: " + introp_coverages
print("Length of coverages (should agree): " 
        + str(len(db_coverages.split(','))))

if (db_bv == introp_bv):
    print("The junction-presenece bitvectors derived from the original file"
        +" and the DB match.")
else:
    print("ERR:The junction-presenece bitvectors derived from the original file"
        +" and the DB DO NOT match.")
