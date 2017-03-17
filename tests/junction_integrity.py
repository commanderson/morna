#!/usr/bin/env python

import argparse
import re
import sqlite3
from BitVector import BitVector

parser = argparse.ArgumentParser()
parser.add_argument('-d','--database', metavar='<str>', type=str,
            required=True,
            help=('path to junction database file for which you wish'
                  'to test integrity')
        )
args = parser.parse_args()

conn=sqlite3.connect(args.database)
c=conn.cursor()
tables=[]
for line in c.execute("SELECT name FROM sqlite_master WHERE type='table'"):
    tables.append(str(line[0]))

for table in sorted(tables):
    print "Checking table " + table
    for juncs,covs in c.execute(("SELECT * FROM %s")%table):
        num_juncs=0
        shr_str=str(juncs)
        last_one = 'x'
        while len(shr_str)>0:
            m = re.search("(.*)([io])(\d+)$",shr_str)
            if m.groups()[1] == 'i':
                bit='1'
                num_juncs+=int(m.groups()[2])
            else:
                bit='0'
            if last_one==bit:
                print "With shr_str " + shr_str + " we found 2 " + bit + "s in a row"
            else:
                last_one=bit
            shr_str = m.groups()[0]
        if not num_juncs==len(covs.split(',')):
		print "In table " + table + "found numjuncs " +str(num_juncs) + " but len(covs.split(',')) " +str(len(covs.split(',')))
