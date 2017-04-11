#!/usr/bin/env python
"""
morna

Uses mmh3 to implement feature hashing for dimensionality reduction of 
intropolis junction vectors by sample and uses Spotify's annoy library for
obtaining nearest neighbors.

Requires mmh3 (pip install mmh3) and Spotify's annoy (pip install annoy)


morna index


"""
import argparse
import bisect
import gzip
import cPickle
import mmh3
import os
import re
import sqlite3
import subprocess
import sys
from annoy import AnnoyIndex
from BitVector import BitVector
from collections import defaultdict
from math import log, sqrt
from utils import *

_help_intro = \
"""morna searches for known RNA-seq samples with exon-exon junction expression
patterns similar to those in a query sample.
"""

            
def encode_64(num):
    """ Encodes an input number as a base 64 number,
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63

        num: integer number to change to base 64
        
        return value: string of base 64 number in format above
    """
    s = [chr(48 + num % 64)]
    num /= 64
    while num > 0:
        s.append(chr(48 + num % 64))
        num /= 64
    return ''.join(s[::-1])

def decode_64(s):
    """ Decodes an input base 64 number formatted
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63 into decimal int

        s: string of base 64 number in format above to convert
        
        return value: int decimal number form of input
    """
    return sum([(ord(s[idx])-48)
                    *(64**(len(s)-idx-1)) for idx in (xrange(len(s)))])
    
def increment_64(s):
    """ Increments an input base 64 number formatted
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63 into decimal int

        s: string of base 64 number in format above to increment
        
        return value: string of base 64 number s + 1 in format above
    """
    return encode_64(decode_64(s) + 1)
    
def magnitude(v1):
    """ Calculate the magnitude of a vector
    
        v1 vectors/list of numerical values
    """ 
    mag = [(a)**2 for a in v1]
    mag = sqrt(sum(mag))
    return mag
    
def dot_prod(v1,v2):
    """Calculate dot product for 2 vectors
    
        v1 and v2: vectors of matching length
    """
    return sum([i*j for (i, j) in zip(v1, v2)])
    
def cosine_distance(v1,v2):
    """Calculate dot product for 2 vectors
    
        v1 and v2: vectors of matching length
    """
    cosine_similarity = dot_prod(v1, v2)/ (magnitude(v1) * magnitude(v2))
    return 1-cosine_similarity

class MornaIndex(AnnoyIndex):
    """ Augments AnnoyIndex to accommodate parameters needed by morna

        Four indexes are written when invoking save(): an Annoy index
        (basename.annoy.mor), a dictionary that maps each junction to the 
        number of samples in which it's found (basename.freq.mor), a stats 
        file including the total number of samples and the size of the 
        hashed feature space (basename.stats.mor), and a 100-sharded 
        sqlite database mapping each sample to the junctions it contained 
        and their repective read coverages (basename.shXX.junc.mor).
        An additional, optional index in the form of a sqlite database 
        containing metadata info associated with sample ids is included if 
        the -m (--metadata) option is specified with an appropriately-
        formatted metadata file.

        See https://github.com/spotify/annoy/blob/master/annoy/__init__.py
        for original class. 
    """
    def __init__(self, sample_count, basename, dim=3000, 
                 sample_threshold=100, metafile=None, buffer_size=1024):
        super(MornaIndex, self).__init__(dim, metric='angular')
        # Store the total number of samples represented in the index
        self.sample_count = sample_count
        #Store basename for index files
        self.basename = basename
        #Store the map between provided sample ids and internal ids,
        #which are an integer range from 0-number_of_samples-1
        self.internal_id_map = {}
        self.new_internal_id = 0
        # Store numbers of samples in which each junction is found
        self.sample_frequencies = defaultdict(int)
        # Store low-dimensional representations of samples
        self.sample_feature_matrix = defaultdict(lambda:
                                                [0.0 for _ in xrange(dim)]
                                        )
        #Store filename of metadata database
        self.metafile = metafile
        # Dimension of low-dimensional feature feature vector
        self.dim = self.dimension_count = dim
        # Minimum number of samples in which junction should be found
        self.sample_threshold = sample_threshold
        # list of junction ids per sample database connection;
        #Will remove old DB if it already exists!
        for shard_id in range(0,100):
            try:
                os.remove(self.basename  + ".sh" 
                + format(shard_id, '02d') + ".junc.mor")
            except OSError:
                pass
        self.buffer_size = buffer_size
        self.junction_conns={}
        self.cursors={}
        #self.junction_conn = sqlite3.connect(self.basename + '.junc.mor')
        #self.junction_conn.isolation_level = None
        #self.junc_cursor = self.junction_conn.cursor()
        #self.junc_cursor.execute("BEGIN")
        self.junc_id = -1
        self.last_present_junction = defaultdict(lambda:-1)
        self.jns_write_buffer = defaultdict(list)
        self.cov_write_buffer = defaultdict(list)                
    
    #@profile 
    def add_junction(self, junction, samples, coverages):
        """ Adds contributions of a single junction to feature vectors
            and also constructs junctions-by-sample database.
            
            junction: string representation of junction
                (typically chromosome, 1-based inclusive start position,
                 1-based inclusive end position)
            samples: list of integer indices of samples
            coverages: list of integer coverages corresponding to samples

            No return value.
        """
        
        self.junc_id += 1
        
        
        for i,sample_id in enumerate(samples):
            try:
                internal_id = self.internal_id_map[sample_id]
            except KeyError:
                self.internal_id_map[sample_id] = self.new_internal_id
                self.new_internal_id+=1
                internal_id = self.internal_id_map[sample_id]
            try:
                shard_id = mmh3.hash(str(internal_id)) % 100
                
                #junc_cursor = self.junction_conns[shard_id].cursor()
                junc_cursor = self.cursors[shard_id]
            except KeyError:
                #Make the db
                self.junction_conns[shard_id] = sqlite3.connect(self.basename 
                                           + ".sh" + format(shard_id, '02d') 
                                           + ".junc.mor")
                self.junction_conns[shard_id].isolation_level = None
                self.cursors[shard_id] = self.junction_conns[shard_id].cursor()
                junc_cursor = self.cursors[shard_id]
                junc_cursor.execute("BEGIN")
        
            if self.last_present_junction[internal_id] == -1: 
                #if this is the first time we've seen this sample,
                #we'll have to make a table for it.
                junc_cursor.execute(("CREATE TABLE sample_%d "         
                                        +"(junctions TEXT,"
                                         +"coverages TEXT)") 
                                            % internal_id)                
                   
                #All start with buffers now
                brand_new_junctions = "!1"
                if self.junc_id >0:
                    brand_new_junctions = ("." + encode_64(self.junc_id) 
                                            + "!1")
                
                #self.junc_cursor.execute(sql, 
                #            [brand_new_junctions,str(coverages[i])])
                self.jns_write_buffer[internal_id].append(brand_new_junctions)
                self.cov_write_buffer[internal_id].append(coverages[i])
                
            else: #if there is already a table for this sample
                #this will imply self.last_present_junction[internal_id] isn't -1
                #first we check if we're on a run of 1s for this sample
                if self.last_present_junction[internal_id] == self.junc_id - 1:
                    #If so, we need to get the current run length
                    #If there's stuff in the buffer, we can check there
                    if (self.jns_write_buffer[internal_id]):
                        m = re.search("!([0-o]+)$",
                                  self.jns_write_buffer[internal_id][-1])
                        new_run_length = increment_64(m.group(1))
                        new_terminus = ('!' + new_run_length)
                        #and even update the run itself in buffer
                        self.jns_write_buffer[internal_id][-1]=(new_terminus)
                        self.cov_write_buffer[internal_id].append(coverages[i])
                    #Otherwise, we actually have to update the run in db :/    
                    #This is an edge case that rarely applies, 
                    #but it's definitely the worst case performance-wise.
                    else:
                        current_junctions=""
                        junc_cursor.execute(("SELECT ROWID FROM sample_%d"
                                 + " ORDER BY ROWID DESC LIMIT 1") % internal_id)
                        last_rowid = junc_cursor.fetchone()[0]
                        junc_cursor.execute(
                        ("SELECT junctions FROM sample_%d "
                                + "ORDER BY ROWID DESC LIMIT 1") % internal_id)
                        current_junctions = junc_cursor.fetchone()[0]
                        m = re.search("(.*?)!([0-o]+)$", current_junctions)
                        new_run_length = increment_64(m.group(2))
                        new_junctions = (m.group(1) + '!' +
                                         new_run_length)
                        
                        sql = (("UPDATE sample_%d SET junctions = ?," 
                                + " coverages = coverages||?||','"
                                + "WHERE ROWID = %d") 
                                                    % (internal_id,last_rowid))
                        junc_cursor.execute(sql,
                            [new_junctions,coverages[i]])
                    
                #Otherwise, provided we're not on a run of 1s
                else: 
                    self.jns_write_buffer[internal_id].append(
                            "." + encode_64((self.junc_id 
                                - self.last_present_junction[internal_id]) - 1))
                    self.jns_write_buffer[internal_id].append("!1")
                    self.cov_write_buffer[internal_id].append(coverages[i])
                
                #We will write buffer contents when buffer gets too big,
                #and then empty the sample's buffer
                if (sys.getsizeof(self.jns_write_buffer[internal_id]) 
                                                > self.buffer_size):
                    #print "let's write buffer for sample " + str(internal_id)
                    
                    sql = ((
                    "INSERT INTO sample_%d VALUES (?, ?)")
                    % internal_id)
                    #TODO: format juncs buffer entries as always having an .X!X 
                    #format, so each has a corresponding single cov entry? Right 
                    #now we are content with variable length rows in db
                    junc_cursor.execute(sql, 
                            ["".join(self.jns_write_buffer[internal_id]),
                            (",".join(str(c) for c in         
                                    self.cov_write_buffer[internal_id]) + ",")
                            ]
                    
                    )
                    #self.junction_conn.commit()
                    self.jns_write_buffer[internal_id] = []
                    self.cov_write_buffer[internal_id] = []
            #No matter what case we hit, each sample must update its 
            #last_present_junction to the current.
            self.last_present_junction[internal_id] = self.junc_id
        
        #if (self.junc_id % 1000 == 0):
        #    self.junction_conn.commit()
        #    print >>sys.stderr, (
         #                   'Wrote junc_id {} into tables'
          #              ).format(self.junc_id)
        
        if len(samples) < self.sample_threshold:
            return
        self.sample_frequencies[junction] += len(samples)
        #right now we hash on 'chromosome start stop'
        #maybe strand someday but shouldn't matter
        hashed_value = mmh3.hash(junction)
        multiplier = (-1 if hashed_value < 0 else 1)
        hashed_value = hashed_value % self.dim
        idf_value = log(float(self.sample_count) 
                        / self.sample_frequencies[junction]
                        )
        for sample, coverage in zip(samples, coverages):
            tf_idf_score = (coverage * idf_value)
            #previously used 1+ norm_log(int(coverage))
            self.sample_feature_matrix[
                    int(sample)][
                    hashed_value] += (multiplier * tf_idf_score)

    def build(self, n_trees, verbose=False):
        """ Adds sample-feature matrix to Annoy index before building.

            n_trees: number of trees in Annoy index
            verbose: write status updates to stderr

            No return value.
        """
        if verbose:
            for i, sample in enumerate(self.sample_feature_matrix):
                self.add_item(sample, self.sample_feature_matrix[sample])
                if not (i % 100):
                    print >>sys.stderr, (
                            'Added {} samples to Annoy index so far.'
                        ).format(i+1)
            print >>sys.stderr, (
                    'Added a total of {} samples to Annoy index.'
                ).format(i+1)
        else:
            for sample in self.sample_feature_matrix:
                self.add_item(sample, self.sample_feature_matrix[sample])
        super(MornaIndex, self).build(n_trees)

    def save(self, basename):
        """ Saves index information to hard disk in three files;
            The annoy index as basename.annoy.mor,
            the sample count and stats information as basename.stats.mor,
            the sample frequency dicitonary as basename.freq.mor,
            and the metadata from metafile as basename.meta.mor

            basename: path used as base of the filename for each saved file
            
            No return value.
        """
        # Save Annoy index first
        super(MornaIndex, self).save(basename + '.annoy.mor')
        # Write total number of samples
        with open(basename + ".stats.mor", 'w') as stats_stream:
            print >>stats_stream, str(self.sample_count)
            print >>stats_stream, str(self.dim)
        # Pickle sample frequency dictionary
        with open(basename + ".freq.mor", 'w') as pickle_stream:
            cPickle.dump(self.sample_frequencies, pickle_stream,
                         cPickle.HIGHEST_PROTOCOL)
        # Pickle sample id to interal id map
        with open(basename + ".map.mor", 'w') as pickle_stream:
            cPickle.dump(self.internal_id_map, pickle_stream,
                         cPickle.HIGHEST_PROTOCOL)
        #Now, deal with junction db:
        #We don't need to pad out all junction lists with runs of 0s,
        #we just assume all junctions after the end of the coverage are 0s
        #We know all shards have been created now,
        # so getting shard cursor is easier
        
        for internal_id in self.last_present_junction.keys():
            ##First, add terminal 0s to each buffer (if needed)
            #if self.last_present_junction[internal_id]<self.junc_id:
            #    self.jns_write_buffer[internal_id].append('o'+str(self.junc_id
            #                    - self.last_present_junction[internal_id]))
                #print("gonna add " + str(additional_junctions))
            
            shard_id = mmh3.hash(str(internal_id)) % 100
            junc_cursor = self.cursors[shard_id]
            
            #Write each buffer out, skipping only empty buffers
            #(NOTE: Not adding terminal 0s just infeerring them. If we did want 
            #terminal 0s, must use different sql if there are no coverages to 
            #add, aka the only-terminal-0s case)
            
            if (self.jns_write_buffer[internal_id]):
                sql = ((
                    "INSERT INTO sample_%d VALUES (?, ?)")
                    % internal_id)
                junc_cursor.execute(sql, 
                        ["".join(self.jns_write_buffer[internal_id]),
                        (",".join(str(c) for c in         
                                self.cov_write_buffer[internal_id]) + ",")
                        ]
                )
                
        #with that, we have finished writing our run-length-encoded junctions!
        for conn in self.junction_conns.values():
            conn.commit()
            conn.execute("VACUUM")
            conn.close()
        
        if self.metafile:
            #Parse file at metafaile and save it into a metadata db
            #Metafile format (whitespace separated)
            #First column: sample id.
            #Any number of additional columns: keywords describing sample,
            #separated by whitespace.
            conn = sqlite3.connect(basename + '.meta.mor')
            cursor = conn.cursor()
        
            #this next part catches if you are overwriting an old index db;
            #in that case it drops the old table before proceeding
            cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
            +" AND name='metadata'")
            if cursor.fetchone():
                cursor.execute("DROP TABLE metadata")
        
            cursor.execute(
                       "CREATE TABLE metadata (sample_id real, keywords text)"
                          )
        
            with open(self.metafile) as metafile_handle:
                for i,line in enumerate(metafile_handle):
                    cursor.execute("INSERT INTO metadata VALUES ('%s','%s')"
                        % (line.split(None,1)[0],line.split(None,1)[1]))
            conn.commit()
            conn.close()

class MornaSearch(object):
    """A class for searching our morna indexes, 
        
        Three index files are needed when initializing: an Annoy index
        (basename.annoy.mor), a dictionary that maps each junction to the 
        number of samples in which it's found (basename.freq.mor), and
        a stats file including the total number of samples (basename.stats.mor).
    """
    def __init__(self, basename):
        # Load AnnoyIndex with AnnoyIndex class, not MornaIndex
        #self.dim = self.dimension_count = dim
        self.basename = basename
        
        with open(basename + ".stats.mor") as stats_stream:
            self.sample_count = int(stats_stream.readline())
            self.dim = int(stats_stream.readline())
        
        self.query_sample = [0.0 for _ in xrange(self.dim)]
        
        self.annoy_index = AnnoyIndex(self.dim)
        self.annoy_index.load(basename + '.annoy.mor')
        
        with open(basename + ".freq.mor") as pickle_stream:
            self.sample_frequencies = cPickle.load(pickle_stream)
            
        with open(basename + ".map.mor") as pickle_stream:
            self.internal_id_map = cPickle.load(pickle_stream)
            
    def inverse_lookup(self, internal_id):
        """ Finds the sample_id key in self.internal_id_map corresponding 
            to the provided internal id; this mapping should always be 1:1
            
            internal_id: an integer in the range 0-(numberOfSamples-1)

            return value: the input file sample_id corresponding to this 
            internal id
        """
        match = None
        found_one_already = False
        for sample_id in self.internal_id_map.keys():
            if self.internal_id_map[sample_id] == internal_id:
                match = sample_id
                if found_one_already == True:
                    raise RuntimeError(str(internal_id) 
                     + " does not have unique mapping in self.internal_id_map.")
                found_one_already = True
        return match
                
                
    def update_query(self, junction):
        """ Updates the query with a junction in an additive fashion

            junction: A junction with which to update the query,
            in format 'chromosome start end coverage'
            
            No return value.
        """
        hashable_junction = ' '.join(map(str, junction[:3]))
        
        if (self.sample_frequencies[hashable_junction] == 0.0):
            idf_value = 0
        else:
            idf_value = log(float(self.sample_count)
                            / self.sample_frequencies[hashable_junction])

        hash_value = mmh3.hash(hashable_junction)
        multiplier = (-1 if hash_value < 0 else 1)
        self.query_sample[hash_value % self.dim] += (
                    multiplier * (junction[3] * idf_value)
                )
        
        
            
    def search_nn(self, num_neighbors, search_k, 
                    include_distances=True, meta_db=False):
        """ Uses the current query_sample to query the annoy index

            num_neighbors: an integer number of nearest neighbors to return
            search_k: an integer number of nodes to check while searching 
            the trees in the index; larger values lead to slower, more accurate 
            searches.
            include_distances: bolean indicator - should the distance to each 
            neighbor be included along with the sample id?
            meta_db: bolean indicator - should an associated metadata db be 
            used to add metadata keywords to results?
            
            Return value: a list of the sample ids of the neareast 
            neighbors of the query, optionally joined in a tuple by 
            a corresponding list of distances and/or a list of metadata keywords 
            found in a supplied metadata db.
        """
        if include_distances:
            results = self.annoy_index.get_nns_by_vector(
                                                 [feature for feature in
                                                 self.query_sample], 
                                                 num_neighbors,
                                                 search_k,
                                                 include_distances
                                                )
        else:
            results = (self.annoy_index.get_nns_by_vector(
                                                 [feature for feature in
                                                 self.query_sample], 
                                                 num_neighbors,
                                                 search_k,
                                                 include_distances
                                                ),)
        if meta_db:
            meta_results = ['' for _ in results[0]]
            conn = sqlite3.connect(self.basename + ".meta.mor")
            cursor = conn.cursor()
            for i,internal_id in enumerate(results[0]):
                sample_id = self.inverse_lookup(internal_id)
                cursor.execute(
                    'SELECT keywords FROM metadata WHERE sample_id=?',
                                                 (str(sample_id),))
                meta_results[i] = cursor.fetchone()
            return results + (meta_results,)
        else:
            return results

            
    def exact_search_nn(self, num_neighbors, include_distances=True,
                         meta_db=False):
        """ Manually search the annoy index for the nearest neighbors to 
            the current query_sample.

            num_neighbors: an integer number of nearest neighbors to return
            include_distances: bolean indicator - should the distance to each 
            neighbor be included along with the sample id?
            meta_db: bolean indicator - should an associated metadata db be 
            used to add metadata keywords to results?
            
            Return value: a list of the sample ids of the neareast 
            neighbors of the query, optionally joined in a tuple by 
            a corresponding list of distances and/or a list of metadata keywords 
            found in a supplied metadata db.
        """ 
        neighbor_indexes=[]
        neighbor_distances=[]

        for i in range(0,self.sample_count):
            current_distance = cosine_distance(
                self.annoy_index.get_item_vector(i),
                self.query_sample)
            insert_point = bisect.bisect_left(neighbor_distances,
                                          current_distance)
            if insert_point < num_neighbors:
                neighbor_distances.insert(insert_point,current_distance)
                neighbor_indexes.insert(insert_point,i)
            if len(neighbor_distances) > num_neighbors:
                neighbor_distances = neighbor_distances[0:num_neighbors]
                neighbor_indexes = neighbor_indexes[0:num_neighbors]

        results = (neighbor_indexes,)
        if include_distances:
            results += (neighbor_distances,)
        
        if meta_db:
            meta_results = ['' for _ in results[0]]
            conn = sqlite3.connect(self.basename + ".meta.mor")
            cursor = conn.cursor()
            for i,internal_id in enumerate(results[0]):
                sample_id = self.inverse_lookup(internal_id)
                cursor.execute(
                    'SELECT keywords FROM metadata WHERE sample_id=?',
                                                 (str(sample_id),))
                meta_results[i] = cursor.fetchone()
            return results + (meta_results,)
        else:
            return results
            
            
    def search_member_n(self, id, num_neighbors, search_k, 
                    include_distances=True, meta_db=False):
        """ Uses a given element of the annoy index to query it

            id: the id in the annoy index (corresponds to sample id) of the 
            element to use as the query
            num_neighbors: an integer number of nearest neighbors to return
            search_k: an integer number of nodes to check while searching 
            the trees in the index; larger values lead to slower, more accurate 
            searches.
            include_distances: bolean indicator - should the distance to each 
            neighbor be included along with the sample id?
            meta_db: bolean indicator - should an associated metadata db be 
            used to add metadata keywords to results?
            
            Return value: a list of the sample ids of the neareast 
            neighbors of the query, optionally joined in a tuple by 
            a corresponding list of distances and/or a list of metadata keywords 
            found in a supplied metadata db.
        """
        if include_distances:
            results = self.annoy_index.get_nns_by_item(
                                                 id, 
                                                 num_neighbors,
                                                 search_k,
                                                 include_distances
                                                )
        else:
            results = (self.annoy_index.get_nns_by_item(
                                                 id, 
                                                 num_neighbors,
                                                 search_k,
                                                 include_distances
                                                ),)
        if meta_db:
            meta_results = ['' for _ in results[0]]
            conn = sqlite3.connect(self.basename + ".meta.mor")
            cursor = conn.cursor()
            for i,internal_id in enumerate(results[0]):
                sample_id = self.inverse_lookup(internal_id)
                cursor.execute(
                    'SELECT keywords FROM metadata WHERE sample_id=?',
                                                 (str(sample_id),))
                meta_results[i] = cursor.fetchone()
            return results + (meta_results,)
        else:
            return results

def add_search_parameters(subparser):
    """ Adds command-line parameters associated with search subcommand

        subparser: subparser of argparse.ArgumentParser

        No return value.
    """
    subparser.add_argument('-x', '--basename', metavar='<idx>', type=str,
            required=True,
            help='path to junction index basename for search'
        )
    subparser.add_argument('-v', '--verbose', action='store_const',
            const=True,
            default=False,
            help='be talkative'
        )
    subparser.add_argument('--search-k', metavar='<int>', type=int,
            required=False,
            default=100,
            help='a larger value makes for more accurate search'
        )
    subparser.add_argument('-f', '--format', metavar='<choice>', type=str,
            required=False,
            default='sam',
            help=('one of {sam, bed, raw}')
        )
    subparser.add_argument('-d', '--distances', action='store_const',
            const=True,
            default=False,
            help='include distances to nearest neighbors'
        )
    subparser.add_argument('-m', '--metadata', action='store_const',
            const=True,
            default=False,
            help='display results mapped to metadata'
        )
    subparser.add_argument('-c','--converge', metavar='<int>', type=int,
            required=False,
            default=None,
            help=('attempt to converge on a solution with concordance equal'
                  'to the given argument as a percentage (truncated to 0-100)')
        )
    subparser.add_argument('-q','--query-id', metavar='<int>', type=int,
            required=False,
            default=None,
            help=('if provided, search is for nearest neighbors to the sample'
                  'already in the index with this id')
        )
    subparser.add_argument('-e', '--exact', action='store_const',
            const=True,
            default=False,
            help='search for exact nearest neighbor to query within morna index'
                 'rather than using annoy hyperplane division algorithm'
        )
    subparser.add_argument('-r','--results', metavar='<int>', type=int,
            required=False,
            default=20,
            help=('the number of nearest neighbor results to return')
        )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_help_intro, 
                formatter_class=help_formatter)
    subparsers = parser.add_subparsers(help=(
                'subcommands; add "-h" or "--help" '
                'after a subcommand for its parameters'),
                dest='subparser_name'
            )
    index_parser = subparsers.add_parser('index', help='creates a morna index')
    search_parser = subparsers.add_parser('search',
                                            help='searches a morna index')
    
    align_parser = subparsers.add_parser('align',
                                            help='performs alignment')
    # Add command-line arguments
    index_parser.add_argument('--intropolis', metavar='<file>', type=str,
            required=True,
            help=('path to gzipped file recording junctions across samples '
                  'in intropolis format')
        )
    index_parser.add_argument('-x', '--basename', metavar='<str>', type=str,
            required=False,
            default='morna',
            help='basename path of junction index files to create'
        )
    index_parser.add_argument('--features', metavar='<int>', type=int,
            required=False,
            default=3000,
            help='dimension of feature space'
        )
    index_parser.add_argument('--n-trees', metavar='<int>', type=int,
            required=False,
            default=200,
            help='number of annoy trees'
        )
    index_parser.add_argument('-s','--sample-count', metavar='<int>', type=int,
            required=False,
            default=None,
            help=('optionally specify number of unique samples to speed '
                  'indexing')
        )
    index_parser.add_argument('-t', '--sample-threshold', metavar='<int>',
            type=int,
            required=False,
            default=100,
            help=('minimum number of samples in which a junction should '
                  'appear for it to be included in morna index')
        )
    index_parser.add_argument('-b', '--buffer-size', metavar='<int>',
            type=int,
            required=False,
            default=1024,
            help=('size in bytes of each entry (one per sample) in the '
             'buffer for writing per-sample junctions to the' 
             '.juncs.mor database')
        )
    index_parser.add_argument('-v', '--verbose', action='store_const',
            const=True,
            default=False,
            help='be talkative'
        )
    index_parser.add_argument('-m', '--metafile',  metavar='<file>', type=str,
            required=False,
            default=None,
            help=('path to metadata file with sample index in first column'
                  'and other keywords in other columns, whitespace delimited')
        )
    add_search_parameters(search_parser)
    add_search_parameters(align_parser)
    align_parser.add_argument('--aligner', '-a', metavar='<exe>', type=str,
            required=False,
            default='STAR',
            help='STAR or HISAT(2) binary'
        )
    align_parser.add_argument('-1', '--mates-1', metavar='<m1>',
            type=str, nargs='+',
            required=False,
            help='files with #1 mates'
        )
    align_parser.add_argument('-2', '--mates-2', metavar='<m2>',
            type=str, nargs='+',
            required=False,
            default=None,
            help='files with #2 mates'
        )
    align_parser.add_argument('-U', '--unpaired', metavar='<r>',
            type=str, nargs='+',
            required=False,
            default=None,
            help='files with unpaired reads'
        )
    align_parser.add_argument('-i', '--index', metavar='<idx>', type=str,
            required=True,
            help=('index basename or directory; use same assembly as was used '
                  'for junctions indexed by morna\'s similarity search')
        )
    align_parser.add_argument('--aligner-args',
            type=str,
            required=False,
            default=None,
            help='additional arguments to pass to aligner; use quote string'
        )
    align_parser.add_argument('-p1', '--pass1-sam', metavar='<sam>',
            type=str,
            required=False,
            default="pass1.sam",
            help='filename for pass 1 alignment file output by aligner'
        )
    align_parser.add_argument('--junction-filter', type=str, required=False,
        default='.05,5',
        help='Two part junction filter settings separated by comma. Only retain' 
             'junctions for alignment found in at least {first part} proportion'
             'of result samples, or junctions with at least {second part}'
             'coverage in any one result sample.')
    args = parser.parse_args()
    if args.subparser_name == 'index':
        # Index
        if not args.sample_count:
            # Count number of samples
            samples = set()
            with gzip.open(args.intropolis) as introp_file_handle:
                if args.verbose:
                    for i, line in enumerate(introp_file_handle):
                        if i % 100 == 0:
                            sys.stdout.write(
                                str(i) + " lines into sample count, " 
                                + str(len(samples)) + " samples so far.\r"
                            )
                            sys.stdout.flush()
                        samples.update(line.split('\t')[-2].split(','))
                else:
                    for i, line in enumerate(introp_file_handle):
                        samples.update(line.split('\t')[-2].split(','))
            args.sample_count = len(samples)
        if args.verbose:
            print '\nThere are {} samples.'.format(args.sample_count)
            
        morna_index = MornaIndex(args.sample_count, args.basename, 
                                    dim=args.features,
                                    sample_threshold=args.sample_threshold, 
                                    metafile=args.metafile,
                                    buffer_size=args.buffer_size)
        with gzip.open(args.intropolis) as introp_file_handle:
            if args.verbose:
                for i, line in enumerate(introp_file_handle):
                    if i % 1000 == 0:
                        sys.stdout.write(str(i) + " lines into index making\r")
                        sys.stdout.flush()
                    tokens = line.strip().split('\t')
                    morna_index.add_junction(
                            ' '.join(tokens[:3]),
                            map(int, tokens[-2].split(',')), # samples
                            map(int, tokens[-1].split(',')), # coverages
                        )
            else:
                for line in introp_file_handle:
                    tokens = line.strip().split('\t')
                    morna_index.add_junction(
                            ' '.join(tokens[:3]),
                            map(int, tokens[-2].split(',')), # samples
                            map(int, tokens[-1].split(',')) # coverages
                        )
        if args.verbose: 
            print 'Finished making index; now building'
        morna_index.build(args.n_trees, verbose=args.verbose)
        morna_index.save(args.basename)
    else:
        junction_stream = sys.stdin 
        #this will be changed to the output sam file 
        #from first pass alignment in the morna align case
        
        if args.subparser_name == 'align':
            # Check command-line parameters
            if args.unpaired is not None:
                if args.mates_1 is not None or args.mates_2 is not None:
                    raise RuntimeError(
                            'Cannot align both paired and unpaired reads at '
                            'once.'
                        )
            elif (args.mates_1 is not None and args.mates_2 is None 
                    or args.mates_2 is not None and args.mates_1 is None):
                raise RuntimeError(
                        'If analyzing paired-end reads, must specify paths '
                        'to both mates files.'
                    )
            elif len(args.mates_1) != len(args.mates_2):
                raise RuntimeError(
                        'The numbers of #1-mate and #2-mate files must be '
                        'equal.'
                    )
                    
            command = []
                
            if args.aligner == "STAR":
                command.append("STAR")
                command += ["--genomeDir", args.index]
                command += ["--outFileNamePrefix", args.pass1_sam]
                command.append("--readFilesIn")
                if args.mates_1 is not None and args.mates_2 is not None:
                    command += args.mates_1
                    command += args.mates_2
                else:
                    command += args.unpaired
                    
                    
            elif args.aligner == "hisat2":
                command.append("hisat2")
                command += ["-x", args.index]
                command += ["-S",args.pass1_sam]
                if args.mates_1 is not None and args.mates_2 is not None:
                    command += ["-1", ",".join(args.mates_1)]
                    command += ["-2", ",".join(args.mates_2)]
                else:
                    command += ["-U", ",".join(args.unpaired)]            
            else:
                raise RuntimeError(
                        'Currently supported aligner options are STAR and'
                        'hisat2.'
                    )
            if args.aligner_args is not None:
                for arg in args.aligner_args.split(" "):
                    command.append(arg)
            if args.verbose:
                print("Calling: " 
                        + " ".join(command))
            ret = subprocess.check_call(command)
            junction_stream = open(args.pass1_sam)
            args.format = "sam"
            


        searcher = MornaSearch(basename=args.basename)
        # Search
        # The search by query id case cuts out early
        if args.query_id is not None:
            results = searcher.search_member_n(args.query_id,
                                    args.results, args.search_k,
                                    include_distances=args.distances,
                                     meta_db = args.metadata)
            print results
            quit()
        # Read BED or BAM
        if args.format == "sam":
            junction_generator = junctions_from_sam_stream(junction_stream)
        elif args.format == "bed":
            junction_generator = junctions_from_bed_stream(junction_stream)
        else:
            assert args.format == "raw"
            junction_generator = junctions_from_raw_stream(junction_stream)
        if (args.converge) and (args.format == 'sam'):
            threshold = args.converge
            backoff = 100
            old_results = []
            if args.verbose:
                for i, junction in enumerate(junction_generator):
                    if (i % 100 == 0):
                        sys.stderr.write( str(i) 
                                         + " junctions into query sample\r")
                        sys.stderr.flush()
                    searcher.update_query(junction)
                    if (i == backoff):
                        backoff += backoff
                        sys.stderr.write("\n")
                        results = (searcher.search_nn(args.results, 
                                    args.search_k,
                                    include_distances=args.distances,
                                     meta_db = args.metadata))
                        shared = 0
                        for result in results[0]:
                            if (result in old_results):
                                shared+=1
                        sys.stderr.write(str(shared) + 
                                    "Shared junctions with old results")
                        if (float(shared)/len(results[0]) >= threshold/100.0):
                            sys.stderr.write("Converged after " + str(i) 
                                                + "junctions.\n")
                            print results
                            
                            quit()
                        else:
                            sys.stderr.write("Not converged after " + str(i) 
                                                + "junctions.\n")
                            old_results = results[0]
                sys.stderr.write("No convergence, but here's results:")
                print results
            else:
                for i, junction in enumerate(junction_generator):
                    searcher.update_query(junction)
                    if (i == backoff):
                        backoff += backoff
                        results = (searcher.search_nn(args.results, 
                                    args.search_k,
                                    include_distances=args.distances,
                                     meta_db = args.metadata))
                        shared = 0
                        for result in results[0]:
                            if (result in old_results):
                                shared+=1
                        if (float(shared)/len(results[0]) >= threshold/100.0):
                            print results
                            quit()
                        else:
                           old_results = results[0]
                print results
            
        else:
            if args.verbose:
                for i, junction in enumerate(junction_generator):
                    if (i % 100 == 0):
                        sys.stderr.write( str(i) 
                                         + " junctions into query sample\r")
                        sys.stdout.flush()
                    searcher.update_query(junction)
            else:
                for i, junction in enumerate(junction_generator):
                    searcher.update_query(junction)
            
            if args.verbose:
                sys.stderr.write("\n")
            if args.exact:
                results = searcher.exact_search_nn(args.results,
                 include_distances=args.distances, meta_db = args.metadata)
            else:
                results = (searcher.search_nn(args.results, args.search_k,
                            include_distances=args.distances,
                            meta_db = args.metadata))
            print results
            
        if args.subparser_name == 'align':
            for result in results[0]:
                shard_id = mmh3.hash(str(result)) % 100
                print "shard_id is " + str(shard_id)
                database = args.basename + ".sh" + str(shard_id) + ".junc.mor"
                conn=sqlite3.connect(database)
                c=conn.cursor()
                print("Connected to " + database)
                db_rle_juncs=[]
                db_coverages=[]
                for line in c.execute(
                    ("SELECT * FROM sample_%d") % result):
                        db_rle_juncs.append(line[0])
                        db_coverages.append(line[1])
                print("-------------------------------------------------")
                print("Internal id " + str(result) + " junctions/converages:")
                print db_rle_juncs
                print db_coverages
                conn.close()
            