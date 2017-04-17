#!/usr/bin/env python
"""
shard_check.py

checks that each sample id in a list is present in its appropriate shard in a given database
OR
checks that each sample id in a given sharded db file belings in that shard
"""

import argparse
#import gzip
import mmh3
import re
import sqlite3
import sys
#from BitVector import BitVector


parser = argparse.ArgumentParser()
parser.add_argument('-d','--database', metavar='<file>', type=str,
            required=True,
            help=('path to basename for junction database file')
        )
parser.add_argument('-f','--shard', metavar='<int>', type=int,
            required=False,
            default=None, 
            help=('specific shard to check ( will check all if left blank)')
        )
parser.add_argument('-s','--samples', metavar='<int>', type=int,
            nargs='+',
            required=False,
            default=None, 
            help=('sample ids of sample to compare')
        )
parser.add_argument('-m','--missing', metavar='<int>', type=int,
            nargs='+',
            required=False,
            default=None, 
            help=('If you cannot find a sample table where you expect it, '
                  'specify only database and add it here to search for it')
        )

        
args = parser.parse_args()

if args.samples is not None:
    print("Specific sample(s) mode")
    for sample in args.samples:
        shard_id = mmh3.hash(str(sample)) % 100
        print "shard_id of " + str(sample) + " is " + format(shard_id, '02d')
        db = args.database + ".sh" + format(shard_id, '02d') + ".junc.mor"
        conn=sqlite3.connect(db)
        c=conn.cursor()
        print("Connected to " + db)
        c.execute(("select count(*) from sqlite_master where"
                +" type='table' and name='sample_%d'") % sample )
        if not c.fetchone()[0] == 1:
            print("ERR: table sample_" + str(sample) + " not found in db " 
                    + db + " (Shard_id: " + format(shard_id, '02d'))
            ##TODO: ERROR MESSAGE IF THIS IS NOT 1

else:
    if args.shard is not None:
        print("Specific shard mode")
        shard_id = args.shard
        db = args.database + ".sh" + format(shard_id, '02d') + ".junc.mor"
        conn=sqlite3.connect(db)
        c=conn.cursor()
        print("Connected to " + db)
        tables=[]
        for line in c.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"):
            tables.append(str(line[0]))
        for table in sorted(tables):
            #print table
            #TODO: use re to strip out sample name 
            #and hash it to see if it agrees with shard_id
            m = re.search("sample_([0-9]*)$", table)
            sample_id = m.group(1)
            hashval = mmh3.hash(sample_id) % 100
            if hashval != shard_id:
                print("ERR: table " + table + " found in shard " + db 
                    + " but the shard id of " + sample_id 
                    + " is " + str(hashval))
        
    else:
        print("Whole database mode")                                            
        print("Checking EVERY " + args.database + ".shXX.junc.mor file...")
        all_samples = []
        if args.missing is not None:
            print("Searching for missing samples: " 
                        + ",".join(str(x) for x in args.missing))
        for shard_id in xrange(0,100):
            db = args.database + ".sh" + format(
                                                shard_id, '02d') + ".junc.mor"
            conn=sqlite3.connect(db)
            c=conn.cursor()
            print("Connected to " + db)
            tables=[]
            for line in c.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"):
                tables.append(str(line[0]))
            for table in sorted(tables):
                #print format(shard_id, '02d') + ": table " + table
                all_samples.append(table)
                m = re.search("sample_([0-9]*)$", table)
                sample_id = m.group(1)
                print sample_id
                if int(sample_id) in args.missing:
                    print("ALERT: Found missing table " + table 
                            + " in shard " + db)
                hashval = mmh3.hash(sample_id) % 100
                if hashval != shard_id:
                    print("ERR: table " + table + " found in shard " 
                                    + db + " but the shard id of " 
                                    + sample_id + " is " + str(hashval))
                    
        #print((sorted(all_samples)))
            
conn.close()