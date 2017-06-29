#!/usr/bin/env python
"""
sample_correct_junctions.py

manually calculates the exact junctions and coverages for a given sample, 
and verifies that the information retrieved from the provided .junc.mor 
junctions-by-sample morna db matches up.
"""

import argparse
import cPickle
import gzip
import mmh3
import re
import sqlite3
import sys
from BitVector import BitVector

def isbase64(s):
    """ Checks whether a string is a valid base64 number 
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63
        
        s: string of base 64 number in format above to increment
        
        return value: boolean True (if valid by format) or False
    """
    for i in s:
        if ord(i) < 48 or ord(i) > 111:
            return False
    return True
    
def decode_64(s):
    """ Decodes an input base 64 number formatted
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63 into decimal int

        s: string of base 64 number in format above to convert
        
        return value: int decimal number form of input
    """
    return sum([(ord(s[idx])-48)
                    *(64**(len(s)-idx-1)) for idx in (xrange(len(s)))])


parser = argparse.ArgumentParser()
parser.add_argument('-d','--database', metavar='<file>', type=str,
            required=True,
            help=('path to junction database file for which you wish'
                  'to test the sample; if sharded provide only basename!')
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
parser.add_argument('-u', '--unsharded', action='store_const',
            const=True,
            default=False,
            help='not checking sharded db; -d argument is exactly literal'
        )
parser.add_argument('-x', '--exclude-introp', action='store_const',
            const=True,
            default=False,
            help='skip checking intropolis file (to quickly get db info)'
        )
parser.add_argument('-pv', '--print-vectors', action='store_const',
            const=True,
            default=False,
            help='print the entire vector of junction presence/absence '
                 'and coverages (only recommended on short test data)'
        )


args = parser.parse_args()
print("starting analysis of " + args.intropolis + " for sample " 
            + str(args.sample))

introp_bv = BitVector(size = 0)
introp_coverages = ""
full_length = 0
if not args.exclude_introp:
    with gzip.open(args.intropolis) as introp_handle:
        for i,line in enumerate(introp_handle):
            if (i % 100 == 0):
                sys.stdout.write("On junction " + str(i) + "\r")
                sys.stdout.flush()
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
    print "Intropolis Junctions:"
    if args.print_vectors:
        print introp_bv           
    full_length = len(introp_bv)
    print "Length: " + str(full_length)
    print "Num. Present Junctions: " + str(sum(introp_bv))
    if args.print_vectors:
        print "Coverages:\n" + introp_coverages
    print("Length of coverages (should agree): " 
            + str(len(introp_coverages.split(','))))
    if not (sum(introp_bv) == len(introp_coverages.split(','))):
        print("ERR: The number of present junctions and the length "
        + "of the coverages list don't match in INTROPOLOLIS-LIKE FILE!.")

                
if not args.unsharded:# figure out which shard has our sample in sharded case
    print("Sample id is " + str(args.sample))
    shard_id = mmh3.hash(str(args.sample)) % 100
    print "shard_id is " + format(shard_id, '02d')
    args.database = args.database + ".sh" + format(shard_id, '02d') + ".junc.mor"
    
conn=sqlite3.connect(args.database)
c=conn.cursor()
print("Connected to " + args.database)
db_rle_juncs=[]
db_coverages=[]
for line in c.execute(("SELECT * FROM sample_%d") % args.sample ):
    db_rle_juncs.append(line[0])
    db_coverages.append(line[1])

shr_str = "".join(db_rle_juncs)
db_coverages = "".join(db_coverages)
db_bv=BitVector(size=0)
i = 0
while len(shr_str)>0:
    #print "starting with shr-str of:" +shr_str
    m = re.search("(.*)!([0-o]+)$",shr_str)
    if i % 2 == 0:
            bit='1'
    else:
            bit='0'
      
    if not isbase64(m.group(2)):
        raise ValueError("Found invalid runlength '" +  m.group(2) + "' in parsing " + shr_str)
    db_bv =  BitVector(bitstring=(decode_64(m.group(2))*bit)) + db_bv
    shr_str = m.group(1)
    i+=1
    
print "DB Junctions"
if args.print_vectors:
    print db_bv
db_length = len(db_bv)
print "Pre-padding Length: " + str(len(db_bv))
db_bv += BitVector(size = (full_length - db_length))
print "Post-padding Length: " + str(len(db_bv))
print "Num. Present Junctions: " + str(sum(db_bv))
if args.print_vectors:
    print "Coverages: " + db_coverages
db_coverages = db_coverages[:-1] #remove the trailing comma 

print("Length of coverages (should agree): " 
        + str(len(db_coverages.split(','))))
if not (sum(db_bv) == len(db_coverages.split(','))):
    print("ERR: The number of present junctions and the length "
    + "of the coverages list don't match in the DB.")
if (db_coverages == introp_coverages):
    print("The comma-separated coverage lists match")
else:
    print("ERR:The comma-separated coverage lists DO NOT match")
    #print "INTROP COVS: " + introp_coverages
    #print "DB COVS: " + db_coverages
if (db_bv == introp_bv):
    print("The junction-presenece bitvectors derived from the original file"
        +" and the DB match.")
else:
    print("ERR:The junction-presenece bitvectors derived from the original file"
        +" and the DB DO NOT match.")
    print "introp_bv: "
    #print introp_bv
    print "db_bv: "
    #print db_bv