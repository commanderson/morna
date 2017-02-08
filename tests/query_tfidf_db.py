import argparse
import gzip
import random
import sqlite3
import sys
import tempfile
from math import log, sqrt
from itertools import izip

def dot_product(v1, v2):
    """returns the dot product of two lists
        v1: first vector, of equal length to v2
        v2: second vector, of equal length to v1
        Return value: float equal to
        v1[0]*v2[0] + v1[1]*v2[1] + ... v1[n]*v2[n]
    """
    return sum([float(i)*j for (i, j) in zip(v1, v2)])
    
def euclidean_norm(v1):
    """returns the euclidean norm of a vector
        v1: first vector, of equal length to v2
        Return value: float equal to
        sqrt(v1[0]^2 + v1[1]^2 + ... + v1[n]^2)
    """
    return sqrt(sum([float(i)**2 for i in v1]))
    
    
#######################
#Main body begins here#
#######################

#Set seed for reproducible results!
random.seed(8675309)

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--database', metavar='<str>', type=str,
            required=True,
            help=('path to .db file created by construc_tf_idf.py')
        )
        
parser.add_argument('-q', '--query', metavar='<str>', type=str,
            required=True,
            help=('path to junctions file for querying the db')
        )

parser.add_argument('-f', '--format', metavar='<choice>', type=str,
            required=False,
            default='sam',
            help=('format of query; one of {sam, bed, raw}')
        )
parser.add_argument('-v', '--verbose', action='store_const',
            const=True,
            default=False,
            help='be talkative'
        )
parser.add_argument('-o', '--output', metavar='<str>', type=str,
            required=False,
            default="query_results",
            help='output file name'
        )

args = parser.parse_args()

conn = sqlite3.connect(args.database)
c = conn.cursor()
c2 = conn.cursor()


query_vector = []#TODO: actually take in query vector

top_20_sims = [[-2] for _ in xrange(20)]
top_20_vectors = [[] for _ in xrange(20)]
for table_line in c.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY ROWID ASC;"
        ):
    table_id = table_line[0].strip("()',")
    print("Contents of table " + table_id)
    table_vector=[]
    for row in c2.execute(
    #"SELECT * FROM sample_%d ORDER BY ROWID ASC LIMIT 1" %sample_id):
        "SELECT * FROM %s ORDER BY ROWID ASC" %table_id
        ):
        table_vector.append(float(row[0].strip("()',")))
    
    #TODO: work this logic out
    cosine_sim = (dot_product(table_vector, query_vector) /                  
        (euclidean_norm(table_vector) * euclidean_norm(vector))
    if cosine_sim "belongs on the list!!":
        #insert it into top_20_sims
        #insert vector into corresponding place in top 20 vectors
conn.close()      