#!/usr/bin/env python
"""
junction_integrity.py

tests each table in the provided .junc.mor junctions-by-sample morna db
for logical consistency; the number of junctions indicated in the 
run-length-encoded junctions field must agree with the number of coverages int
the comma-separated coverages field, and all run-length encoded junctions lists
must be of a length equal to the number of junctions and must be in a valid
format, with strictly alternating runs of 'o' and 'i'
"""

import argparse
import re
import sys
import sqlite3

parser = argparse.ArgumentParser()
parser.add_argument('-d','--database', metavar='<str>', type=str,
            required=True,
            help=('path to junction database file for which you wish'
                  'to test integrity')
        )
parser.add_argument('-j','--junctions', metavar='<int>', type=int,
            required=True,
            help=('number of junctions originally provided as input')
        )

args = parser.parse_args()
print "parsed args"
conn=sqlite3.connect(args.database)
c=conn.cursor()
print "connection established"
tables=[]
for line in c.execute("SELECT name FROM sqlite_master WHERE type='table'"):
    tables.append(str(line[0]))
print "tables populated"
problem = False
for table in sorted(tables):
    sys.stdout.write("Checking table " + table + "\n")
    last_one = 2
    min_num_juncs = 0
    num_pos_juncs = 0
    num_covs = 0
    for juncs,covs in c.execute(("SELECT * FROM %s")%table):
        num_covs += len(covs.strip(',').split(','))
        shr_str=str(juncs)
        while len(shr_str)>0:
            m = re.search("^([io])(\d+)(.*)$",shr_str)
            min_num_juncs += int(m.group(2))
            if m.group(1) == 'i':
                bit=1
                num_pos_juncs += int(m.group(2))
            else:
                bit=0
            if last_one==bit:
                problem = True
                print("\nERR:With len(shr_str) " + str(len(shr_str)) 
                        + " we found 2 " + str(bit) + "s in a row")
            else:
                last_one=bit
            shr_str = m.group(3)
    if not num_pos_juncs==num_covs:
        print("\nERR:In table " + table + " found num_pos_juncs " 
                + str(num_pos_juncs) + " but num_covs " 
                + str(num_covs))
        problem = True
    if min_num_juncs>args.junctions:
        print("\nERR:In table " + table + " found min_num_juncs "
                + str(min_num_juncs) + " which is greater than arg junctions"
                + " (" + str(args.junctions)+")")
        problem = True
		    
if problem:
    print("\nERR(s) occurred; see above messages")
else:
    print("\nIntegrity checks out!")
