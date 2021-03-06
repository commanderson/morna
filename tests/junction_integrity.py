#!/usr/bin/env python
"""
junction_integrity.py

tests each table in the provided .junc.mor junctions-by-sample morna db
for logical consistency; the number of junctions indicated in the 
run-length-encoded junctions field must agree with the number of coverages in
the comma-separated coverages field, and all run-length encoded junctions lists
must be of a length equal to the number of junctions and must be in a valid
format, run lengths separated by ! characters
"""

import argparse
import re
import sys
import sqlite3

def decode_64(s):
    """ Decodes an input base 64 number formatted
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63 into decimal int

        s: string of base 64 number in format above to convert
        
        return value: int decimal number form of input
    """
    return sum([(ord(s[idx])-48)
                    *(64**(len(s)-idx-1)) for idx in (xrange(len(s)))])
    
if __name__ == '__main__':

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
    #print "parsed args"
    conn=sqlite3.connect(args.database)
    c=conn.cursor()
    print "connection established to " + args.database
    tables=[]
    for line in c.execute("SELECT name FROM sqlite_master WHERE type='table'"):
        tables.append(str(line[0]))
    #print "tables populated"
    problem = False
    for table in sorted(tables):
        sys.stdout.write("Checking table " + table + "\r")
        sys.stdout.flush()
        #last_one = 2
        min_num_juncs = 0
        num_pos_juncs = 0
        num_covs = 0
        for juncs,covs in c.execute(("SELECT * FROM %s")%table):
            num_covs += len(covs.strip(',').split(','))
            shr_str=str(juncs)
            i=0
            while len(shr_str)>0:
                m = re.search("(.*)!([0-o]+)$",shr_str)
                min_num_juncs += decode_64(m.group(2))
                if i%2==0:
                    num_pos_juncs += decode_64(m.group(2))
                i += 1
                shr_str = m.group(1)
        if not num_pos_juncs==num_covs:
            print("\nERR:In table " + table + " found num_pos_juncs " 
                    + str(num_pos_juncs) + " but num_covs " 
                    + str(num_covs))
            problem = True
        if min_num_juncs>args.junctions:
            print("\nERR:In table " + table + " found min_num_juncs "
                    + str(min_num_juncs) 
                    + " which is greater than arg junctions"
                    + " (" + str(args.junctions)+")")
            problem = True
            
    if problem:
        print("\nERR(s) occurred; see above messages")
    else:
        print("\nIntegrity checks out!")
