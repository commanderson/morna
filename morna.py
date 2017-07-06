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
from math import log, sqrt, ceil
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
    """ Calculate dot product for 2 vectors
    
        v1 and v2: vectors of matching length
    """
    return sum([i*j for (i, j) in zip(v1, v2)])
    
def cosine_distance(v1,v2):
    """ Calculate cosine distance for 2 vectors
    
        v1 and v2: vectors of matching length
    """
    cosine_similarity = dot_prod(v1, v2)/ (magnitude(v1) * magnitude(v2))
    return 1-cosine_similarity

def results_output(results):
    """ output morna search results in human-readable fashion
    
        results: tuple of lists, definitely including sample ids and
        also optionally adding distances and metadata
    """
    
    for i in range(len(results[0])):
        sys.stdout.write(str(i+1) + ".")
        for list in results:
            sys.stdout.write("\t" + str(list[i]))
        sys.stdout.write("\n")
    
def running_sum(rls):
    """ Generate list of present junction ids from our run-length_list
    
        rls: list of run length strings encoded in our base64 format
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63
    """
    tot = 0
    for i, item in enumerate(rls):
        length = decode_64(item)
        if i % 2:
            tot += length
        else:
            for i in xrange(length):
                yield tot + i
            tot += length

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
        #Store the total number of samples represented in the input data
        self.sample_count = sample_count
        # The number actually included in the index may differ (if every
        # junction for a sample fails to meet the sample threshold, a sample
        # would be excluded from the final index entirely)
        # the final value of self.new_internal_id (+1) is the actual number 
        # of samples included in the index
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
        #before being included in the index
        self.sample_threshold = sample_threshold
        
        #count the number of junctions left out due to sample threshold
        self.skipped = 0
        
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
    def update_junction_dbs(self,junction,samples,coverages):
        """ Updates junction dbs with information from this junction 
            as appropriate
            
            junction: string representation of junction
                (typically chromosome, 1-based inclusive start position,
                 1-based inclusive end position)
            samples: list of integer indices of samples
            coverages: list of integer coverages corresponding to samples

            No return value.
        """ 
        
        for i,sample_id in enumerate(samples):
            try:
                shard_id = mmh3.hash(str(sample_id)) % 100
                
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
        
            if self.last_present_junction[sample_id] == -1: 
                #if this is the first time we've seen this sample,
                #we'll have to make a table for it.
                junc_cursor.execute(("CREATE TABLE sample_%d "         
                                        +"(junctions TEXT,"
                                         +"coverages TEXT)") 
                                            % sample_id)                
                   
                #All start with buffers now
                brand_new_junctions = "!1"
                if self.junc_id >0:
                    brand_new_junctions = ("!" + encode_64(self.junc_id) 
                                            + "!1")
                
                #self.junc_cursor.execute(sql, 
                #            [brand_new_junctions,str(coverages[i])])
                self.jns_write_buffer[sample_id].append(brand_new_junctions)
                self.cov_write_buffer[sample_id].append(coverages[i])
                
            else: #if there is already a table for this sample
                #this implies self.last_present_junction[sample_id] isn't -1
                #first we check if we're on a run of 1s for this sample
                #TODO: fix assumption that self.junc_id -1 is always previous one; it might be that we skipped last one due to threshold but are still on a run!
                if self.last_present_junction[sample_id] == self.junc_id - 1:
                    #If so, we need to get the current run length
                    #If there's stuff in the buffer, we can check there
                    if (self.jns_write_buffer[sample_id]):
                        m = re.search("!([0-o]+)$",
                                  self.jns_write_buffer[sample_id][-1])
                        new_run_length = increment_64(m.group(1))
                        new_terminus = ('!' + new_run_length)
                        #and even update the run itself in buffer
                        self.jns_write_buffer[sample_id][-1]=(new_terminus)
                        self.cov_write_buffer[sample_id].append(coverages[i])
                    #Otherwise, we actually have to update the run in db :/    
                    #This is an edge case that rarely applies, 
                    #but it's definitely the worst case performance-wise.
                    else:
                        current_junctions=""
                        junc_cursor.execute(("SELECT ROWID FROM sample_%d"
                                 + " ORDER BY ROWID DESC LIMIT 1") % sample_id)
                        last_rowid = junc_cursor.fetchone()[0]
                        junc_cursor.execute(
                        ("SELECT junctions FROM sample_%d "
                                + "ORDER BY ROWID DESC LIMIT 1") % sample_id)
                        current_junctions = junc_cursor.fetchone()[0]
                        m = re.search("(.*?)!([0-o]+)$", current_junctions)
                        new_run_length = increment_64(m.group(2))
                        new_junctions = (m.group(1) + '!' +
                                         new_run_length)
                        
                        sql = (("UPDATE sample_%d SET junctions = ?," 
                                + " coverages = coverages||?||','"
                                + "WHERE ROWID = %d") 
                                                    % (sample_id,last_rowid))
                        junc_cursor.execute(sql,
                            [new_junctions,coverages[i]])
                    
                #Otherwise, provided we're not on a run of 1s
                else: 
                    self.jns_write_buffer[sample_id].append(
                            "!" + encode_64((self.junc_id 
                                - self.last_present_junction[sample_id]) - 1))
                    self.jns_write_buffer[sample_id].append("!1")
                    self.cov_write_buffer[sample_id].append(coverages[i])
                
                #We will write buffer contents when buffer gets too big,
                #and then empty the sample's buffer
                if (sys.getsizeof(self.jns_write_buffer[sample_id]) 
                                                > self.buffer_size):
                    #print "let's write buffer for sample " + str(sample_id)
                    
                    sql = ((
                    "INSERT INTO sample_%d VALUES (?, ?)")
                    % sample_id)
                    #TODO
                    #NOT: format juncs buffer entries as always having an .X!X 
                    #format, so each has a corresponding single cov entry? Right 
                    #now we are content with variable length rows in db
                    junc_cursor.execute(sql, 
                            ["".join(self.jns_write_buffer[sample_id]),
                            (",".join(str(c) for c in         
                                    self.cov_write_buffer[sample_id]) + ",")
                            ]
                    
                    )
                    #self.junction_conn.commit()
                    self.jns_write_buffer[sample_id] = []
                    self.cov_write_buffer[sample_id] = []
            #No matter what case we hit, each sample must update its 
            #last_present_junction to the current.
            self.last_present_junction[sample_id] = self.junc_id
            
            
    def add_junction(self,junction,samples,coverages):
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
        
        self.sample_frequencies[junction] += len(samples)
        
        self.update_junction_dbs(junction,samples,coverages)
        
        if len(samples) < self.sample_threshold:
            self.skipped += 1
            return
        
        #right now we hash on 'chromosome start stop'
        #maybe strand someday but shouldn't matter
        hashed_value = mmh3.hash(junction)
        multiplier = (-1 if hashed_value < 0 else 1)
        hashed_value = hashed_value % self.dim
        idf_value = log(float(self.sample_count) 
                        / self.sample_frequencies[junction]
                        )

        for sample_id, coverage in zip(samples, coverages):
            try:
                internal_id = self.internal_id_map[sample_id]
            except KeyError:
                self.internal_id_map[sample_id] = self.new_internal_id
                self.new_internal_id+=1
                internal_id = self.internal_id_map[sample_id]

            tf_idf_score = (coverage * idf_value)
            #previously used 1+ norm_log(int(coverage))
            self.sample_feature_matrix[
                    int(internal_id)][
                    hashed_value] += (multiplier * tf_idf_score)

    def build(self, n_trees, verbose=False):
        """ Adds sample-feature matrix to Annoy index before building.

            n_trees: number of trees in Annoy index
            verbose: write status updates to stderr

            No return value.
        """
        #TODO: add handling for silly case when NO samples actually made it
        if self.new_internal_id == 0:
            raise ValueError("No internal ids were assigned, indicating that "
                           + "no samples were added to the index. Likely "
                           + "caused when no junctions pass the sample "
                           + "threshold.")
        if verbose:
            for i, internal_id in enumerate(self.sample_feature_matrix):
                self.add_item(internal_id, 
                        self.sample_feature_matrix[internal_id])
                if not (i % 100):
                    sys.stderr.write(
                        'Added {} samples to Annoy index so far.\r'.format(i+1))
            print >>sys.stderr, (
                    '\nAdded a total of {} samples to Annoy index.'
                ).format(i+1)
            #print >>sys.stderr, (
            #        'Final new internal id value: {}'
            #    ).format(self.new_internal_id)
            print >>sys.stderr, (
                    '{} junctions skipped for not meeting sample threshold'
                ).format(self.skipped)
                    
        else:
            for internal_id in self.sample_feature_matrix:
                self.add_item(internal_id, 
                        self.sample_feature_matrix[internal_id])
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
        # Write total number of samples in input, 
        # total number of samples RETAINED, and dimensionality to
        # the stats file so search can recover them
        with open(basename + ".stats.mor", 'w') as stats_stream:
            print >>stats_stream, str(self.sample_count)
            print >>stats_stream, str(self.new_internal_id)
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
        
        for sample_id in self.last_present_junction.keys():
            ##First, add terminal 0s to each buffer (if needed)
            #if self.last_present_junction[internal_id]<self.junc_id:
            #    self.jns_write_buffer[internal_id].append('o'+str(self.junc_id
            #                    - self.last_present_junction[internal_id]))
                #print("gonna add " + str(additional_junctions))
            
            shard_id = mmh3.hash(str(sample_id)) % 100
            junc_cursor = self.cursors[shard_id]
            
            #Write each buffer out, skipping only empty buffers
            #(NOTE: Not adding terminal 0s just inferring them. If we did want 
            #terminal 0s, must use different sql if there are no coverages to 
            #add, aka the only-terminal-0s case)
            
            if (self.jns_write_buffer[sample_id]):
                sql = ((
                    "INSERT INTO sample_%d VALUES (?, ?)")
                    % sample_id)
                junc_cursor.execute(sql, 
                        ["".join(self.jns_write_buffer[sample_id]),
                        (",".join(str(c) for c in         
                                self.cov_write_buffer[sample_id]) + ",")
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
    """A class for searching our morna indexes for approximate nearest neighbors 
        
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
            self.index_size = int(stats_stream.readline())
            self.dim = int(stats_stream.readline())
        
        self.query = defaultdict(int)
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
                
                
    def old_update_query(self, junction):
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
        
    def update_query(self, junction):
        """ Updates the query dictionary with a junction's coverage info

            junction: A junction in list format with which to update the query,
            in format 'chromosome start end coverage [other stuff]'
            
            No return value.
        """
        #TODO: check if junction is in dictionary of junctions in sample
        #if it isn't, don't include
        self.query[tuple(junction[:3])] += int(junction[3])
                
    def finalize_query(self):
        """Uses the 'self.query' dictionary to create a final query sample 
           stored in self.query_sample; doing it this way avoids hashing at 
           every step
           
           No return value.
        """
        self.query_sample = [0.0 for _ in xrange(self.dim)]
        for junction in self.query.keys():
            hashable_junction = ' '.join(str(_) for _ in junction)        
            if (self.sample_frequencies[hashable_junction] == 0.0):
                idf_value = 0
            else:
                idf_value = log(float(self.sample_count)
                                / self.sample_frequencies[hashable_junction])

            hash_value = mmh3.hash(hashable_junction)
            multiplier = (-1 if hash_value < 0 else 1)
            self.query_sample[hash_value % self.dim] += (
                        multiplier * (self.query[junction] * idf_value)
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
            
            Return value: a list of the internal ids of the neareast 
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
            
            Return value: a list of the internal ids of the neareast 
            neighbors of the query, optionally joined in a tuple by 
            a corresponding list of distances and/or a list of metadata keywords 
            found in a supplied metadata db.
        """ 
        neighbor_indexes=[]
        neighbor_distances=[]

        for i in range(0,self.index_size):
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
            
            
    def search_member_n(self, query_id, num_neighbors, search_k, 
                    include_distances=True, meta_db=False):
        """ Uses a given element of the annoy index to query it

            id: the sample_id  of the element to use as the query
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
        print("querying by sample id " + str(query_id))
        try:
            internal_id = self.internal_id_map[query_id]
        except KeyError:
            raise ValueError("Querying sample id " + str(query_id) 
            + " is not possible because no internal id is mapped to that " 
            + "sample id. Likely no sample with that id was included " 
            + "in the index.")
        print("this is internal id " + str(internal_id))
        if include_distances:
            results = self.annoy_index.get_nns_by_item(
                                                 internal_id, 
                                                 num_neighbors,
                                                 search_k,
                                                 include_distances
                                                )
        else:
            results = (self.annoy_index.get_nns_by_item(
                                                 internal_id, 
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

def count_samples(introp_file_handle,verbose):
    """ Count distinct sample ids in samples-by-junction file 
        at 'intropolis' file path
        
        introp_file_handle: file handle for tab delimited samples-by-junction 
         (intropolis format) file.
            tab-delimited columns:
            1. chromosome where junction is located e.g. chr10
            2. 1-based start position e.g. 101039923
            3. 1-based end position e.g. 101040322
            4. forward or reverse strand, one of '+' or '-'
            5. donor dinucleotide (e.g., GT)
            6. acceptor dinucleotide (e.g., AG)
            7. comma-separated list of indexes of samples in which junction was 
                found
            8. comma-separated list of corresponding numbers of reads mapping 
                across junction in samples from field 7
        verbose: boolean True or False indicating whether progress updates 
         should be written to stdout
    """
    samples = set()
    if verbose:
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
    return len(samples)
    
def go_index(intropolis,basename,features,n_trees,sample_count,
                sample_threshold,buffer_size,verbose,metafile):
    """ Runs Morna index to create an index at 'basename' (all directories 
        must already exist) based on the samples-by-junction file 
        at 'intropolis'.
    """
    if not sample_count:
        with gzip.open(intropolis) as introp_file_handle:
            sample_count = count_samples(introp_file_handle,verbose)
    if verbose:
        print '\nThere are {} samples.'.format(sample_count)
        
    morna_index = MornaIndex(sample_count, basename, 
                                dim=features,
                                sample_threshold=sample_threshold, 
                                metafile=metafile,
                                buffer_size=buffer_size)
    with gzip.open(intropolis) as introp_file_handle:
        if verbose:
            for i, line in enumerate(introp_file_handle):
                if i % 1000 == 0:
                    sys.stdout.write(str(i) 
                            + " lines into index making\r")
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
    if verbose: 
        print 'Finished making index; now building'
    morna_index.build(n_trees, verbose=verbose)
    morna_index.save(basename)

def add_test_parameter(subparser):
    """ Adds test parameter for running tests instead of executing
    """ 
    subparser.add_argument('--test', action='store_const',
            const=True,
            default=False,
            help='run tests'
        )

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
    subparser.add_argument('-c','--convergence-backoff', metavar='<int>', 
            type=int,
            required=False,
            default=None,
            help=('attempt to converge on a solution with this backoff interval'
                  '(check after c, then 2c, then 4c, and stop when no change)')
        )
    subparser.add_argument('-ch','--checkpoint', metavar='<int>', type=int,
            required=False,
            default=0,
            help=('Do not start backoff until this point (check at checkpoint,'
                  ' then checkpoint + c, then continue exponential backoff')
        )
    subparser.add_argument('-q','--query-id', metavar='<int>', type=int,
            required=False,
            default=None,
            help=('if provided, search is for nearest neighbors to the sample'
                  'already in the index with this sample id')
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
    add_test_parameter(index_parser)
    add_test_parameter(search_parser)
    add_test_parameter(align_parser)
    
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
            default='hisat2',
            help='One of STAR, hisat, or hisat2'
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
            help='filename for first pass alignment file output by aligner'
        )
    align_parser.add_argument('-p2', '--pass2-sam', metavar='<sam>',
            type=str,
            required=False,
            default="pass2.sam",
            help='filename for second pass alignment file output by aligner'
        )
    align_parser.add_argument('--junction-filter', type=str, required=False,
        default=".05,5",
        help='Two part junction filter settings separated by comma. Only retain' 
             'junctions for alignment found in at least {first part} proportion'
             'of result samples, or junctions with at least {second part}'
             'coverage in any one result sample.'
        )
    align_parser.add_argument('--junction-file', type=str, metavar='<gz>',
        required=True,
        help='path to gzipped file with junction name elements in same order '
             'as file used to create index (that file will do in a pinch)'
        )
    align_parser.add_argument('-pa', '--previous-alignment',
            action='store_const',
            const=True,
            default=False,
            help='use file at pass1_sam arg as a previous alignment '
                 '(skip first alignment pass and go straight to search)'
        )
    align_parser.add_argument('-sf','--splicefile', type=str, metavar='<gz>',
        required=False,
        default=None,
        help='path to output intropolis-like file with junctions retained from '
             'search results (only change to intropolis-like format is addition'
             'of ranking of samples which provided junction'
        )
    args = parser.parse_args()
    
    
    if args.test:
        print("Tests are a work in progress")
        if args.subparser_name == 'index':
            del sys.argv[1:]
            import unittest
            import tempfile
            import shutil

            class IndexTestGo(unittest.TestCase):
                """Index Tests go() """
                
                def setUp(self):
                    # Set up temporary directory
                    self.temp_dir_path = tempfile.mkdtemp()
                    self.input_file = os.path.join(self.temp_dir_path,
                                                               'junctions.temp')
                    self.meta_file = os.path.join(self.temp_dir_path,
                                                               'meta.temp')
                    self.output_file = os.path.join(self.temp_dir_path,
                                                                 'index.temp')
                    self.generic_input_string = (
                        'chr10\t100138610\t100139562\t-\tGC\tAG\t1\t1\n'
                        'chr10\t100347256\t100362229\t-\tGC\tAG\t2\t3\n'
                        'chr10\t100526502\t100529392\t-\tAT\tAC\t3\t1\n'
                        'chr10\t100526535\t100527418\t+\tAT\tAC\t4\t1\n'
                        'chr10\t100548141\t100548202\t+\tGC\tAG\t5\t1\n'
                        'chr10\t100947776\t101440795\t-\tGC\tAG\t6\t1\n'
                        'chr10\t100962117\t100963108\t+\tGC\tAG\t7\t1\n'
                        'chr10\t100979320\t100982343\t+\tGC\tAG\t8\t1\n'
                        'chr10\t101003908\t101007637\t-\tGC\tAG\t9\t1\n'
                        'chr10\t101035534\t101036515\t-\tGT\tAG\t10\t1\n'
                        'chr10\t101039923\t101040322'
                        '\t-\tGC\tAG\t1,2,3,4\t1,1,1,1\n'
                        'chr10\t101354229\t101407426\t-\tAT\tAC\t5\t1\n'
                        'chr10\t101579666\t101583591\t-\tAT\tAC\t6\t1\n'
                        'chr10\t101669187\t102046680\t+\tAT\tAC\t7,8\t1,1\n'
                        'chr10\t101885932\t101918495\t-\tGC\tAG\t9\t1\n'
                        'chr10\t102078553\t102078686\t-\tGC\tAG\t10\t1\n'
                        'chr10\t102161920\t102162275\t+\tGC\tAG'
                        '\t1,2,3,4,5,6,7,8,9,10\t2,1,1,2,2,1,2,3,8,2\n'
                        'chr10\t102424521\t102425341\t-\tGC\tAG\t1\t1\n'
                        'chr10\t102479334\t102813545\t+\tGT\tAG\t2\t1\n'
                        'chr10\t102594066\t102615267\t+\tGT\tAG'
                        '\t3,4,5,6,7,8,9,10\t1,4,1,1,1,7,1,1\n'
                    )
                    
                    self.lossy_input_string = (
                        'chr10\t100138610\t100139562\t-\tGC\tAG\t10\t1\n'
                        'chr10\t100347256\t100362229\t-\tGC\tAG\t9\t3\n'
                        'chr10\t100526502\t100529392\t-\tAT\tAC\t8\t1\n'
                        'chr10\t100526535\t100527418\t+\tAT\tAC\t7\t1\n'
                        'chr10\t100548141\t100548202\t+\tGC\tAG\t6\t1\n'
                        'chr10\t100947776\t101440795\t-\tGC\tAG\t5\t1\n'
                        'chr10\t100962117\t100963108\t+\tGC\tAG\t4\t1\n'
                        'chr10\t100979320\t100982343\t+\tGC\tAG\t3\t1\n'
                        'chr10\t101003908\t101007637\t-\tGC\tAG\t2\t1\n'
                        'chr10\t101035534\t101036515\t-\tGT\tAG\t1\t1\n'
                        'chr10\t101039923\t101040322\t-\tGC\tAG\t'
                        '10,9,8,7\t1,1,1,1\n'
                        'chr10\t101354229\t101407426\t-\tAT\tAC\t6\t1\n'
                        'chr10\t101579666\t101583591\t-\tAT\tAC\t5\t1\n'
                        'chr10\t101669187\t102046680\t+\tAT\tAC\t4,3\t1,1\n'
                        'chr10\t101885932\t101918495\t-\tGC\tAG\t2\t1\n'
                        'chr10\t102078553\t102078686\t-\tGC\tAG\t1\t1\n'
                        'chr10\t102161920\t102162275\t+\tGC\tAG\t'
                        '10,9,8,7,6,5,4\t2,1,1,2,2,1,2\n'
                        'chr10\t102424521\t102425341\t-\tGC\tAG\t10\t1\n'
                        'chr10\t102479334\t102813545\t+\tGT\tAG\t9\t1\n'
                        'chr10\t102594066\t102615267\t+\tGT\tAG\t'
                        '4,3,2,1\t1,7,1,1\n'
                    )
                    self.meta_string = (
                        '10\tSample 10\tThe 10th sample\n'
                        '1\tSample 1\tThe 1st sample\n'
                        '2\tSample 2\tThe 2nd sample\n'
                        '3\tSample 3\tThe 3rd sample\n'
                        '4\tSample 4\tThe 4th sample\n'
                        '5\tSample 5\tThe 5th sample\n'
                        '6\tSample 6\tThe 6th sample\n'
                        '7\tSample 7\tThe 7th sample\n'
                        '8\tSample 8\tThe 8th sample\n'
                        '9\tSample 9\tThe 9th sample\n'
                    )
                def test_sample_count(self):
                    """ Fails if sample count is inaccurate
                    """
                    with open(self.input_file, 'w') as input_stream:
                        input_stream.write(
                            self.generic_input_string
                        )
                    with open(self.input_file) as input:
                        sample_count = count_samples(input,False)
                    self.assertEqual(sample_count, 10)
                    
                #def test_add_junction(self):
                #    """ Fails if junction does not add correctly
                #    """
                
                def test_simple_indexing(self):
                    """ Test a runthrough on simple sample data with 
                        all junctions included
                    """
                    with gzip.open(self.input_file, 'w') as input_stream:
                        input_stream.write(
                            self.generic_input_string
                        )
                    with open(self.meta_file, 'w') as meta_stream:
                        meta_stream.write(
                            self.meta_string
                        )
                    
                    index_path = os.path.join(self.temp_dir_path ,"tempIndex")
                    go_index(intropolis=self.input_file, basename=index_path,
                        features=3000,n_trees=20,
                        sample_count=10,sample_threshold=1,
                        buffer_size=1024,verbose=False,
                        metafile=self.meta_file)
                    
                    test_index=AnnoyIndex(3000)
                    test_index.load(index_path+".annoy.mor")
                    
                    self.assertEqual(test_index.get_n_items(), 10)
                    
                    expected_results = (
                        [[0, 2, 3, 1, 4, 5, 6, 7, 8, 9],
                        [1, 2, 3, 0, 4, 5, 6, 7, 8, 9],
                        [2, 3, 0, 1, 7, 6, 4, 5, 8, 9],
                        [3, 7, 2, 0, 1, 6, 4, 5, 8, 9],
                        [4, 7, 3, 2, 6, 5, 8, 9, 0, 1],
                        [5, 7, 3, 2, 6, 4, 8, 9, 0, 1],
                        [6, 7, 3, 2, 4, 5, 8, 9, 0, 1],
                        [7, 6, 3, 2, 4, 5, 8, 9, 0, 1],
                        [8, 7, 3, 2, 6, 4, 5, 9, 0, 1],
                        [9, 7, 3, 2, 6, 4, 5, 8, 0, 1]]
                    )

                    for i in range(10):
                        self.assertEqual(
                            test_index.get_nns_by_item(i,10,search_k=100,
                                include_distances=False),
                            expected_results[i])
                
                def test_metadata_indexing(self):
                    """ Test a runthrough on simple sample data with 
                        all junctions included for proper metadata mapping
                    """
                    with gzip.open(self.input_file, 'w') as input_stream:
                        input_stream.write(
                            self.generic_input_string
                        )
                    with open(self.meta_file, 'w') as meta_stream:
                        meta_stream.write(
                            self.meta_string
                        )
                    
                    index_path = os.path.join(self.temp_dir_path ,"tempIndex")
                    go_index(intropolis=self.input_file, basename=index_path,
                        features=3000,n_trees=20,
                        sample_count=10,sample_threshold=1,
                        buffer_size=1024,verbose=False,
                        metafile=self.meta_file)
                    
                    test_index=AnnoyIndex(3000)
                    test_index.load(index_path+".annoy.mor")
                    
                    self.assertEqual(test_index.get_n_items(), 10)
                    
                    expected_results = (
                        ([0, 2, 3, 1, 4, 5, 6, 7, 8, 9],
                        [(u'Sample 1\tThe 1st sample\n',),
                         (u'Sample 3\tThe 3rd sample\n',),
                         (u'Sample 4\tThe 4th sample\n',),
                         (u'Sample 2\tThe 2nd sample\n',), 
                         (u'Sample 5\tThe 5th sample\n',), 
                         (u'Sample 6\tThe 6th sample\n',), 
                         (u'Sample 7\tThe 7th sample\n',), 
                         (u'Sample 8\tThe 8th sample\n',), 
                         (u'Sample 9\tThe 9th sample\n',), 
                         (u'Sample 10\tThe 10th sample\n',)],)
                    )
                    
                    searcher = MornaSearch(basename=index_path)
                    
                    results = searcher.search_member_n(1,
                                        10, 100,
                                        include_distances=False,
                                         meta_db = True)
                    self.assertEqual(results, expected_results)
                    
                def test_lossy_indexing(self):
                    """ Test a runthrough on simple sample data with 
                        many junctions excluded by sample threshold
                    """
                    with gzip.open(self.input_file, 'w') as input_stream:
                        input_stream.write(
                            self.lossy_input_string
                        )
                    with open(self.meta_file, 'w') as meta_stream:
                        meta_stream.write(
                            self.meta_string
                        )
                    
                    index_path = os.path.join(self.temp_dir_path ,"tempIndex")
                    go_index(intropolis=self.input_file, basename=index_path,
                        features=3000,n_trees=20,
                        sample_count=10,sample_threshold=4,
                        buffer_size=1024,verbose=False,
                        metafile=self.meta_file)
                    
                    test_index=AnnoyIndex(3000)
                    test_index.load(index_path+".annoy.mor")
                    
                    self.assertEqual(test_index.get_n_items(), 10)
                    
                    expected_results = (
                        [[0, 3, 1, 2, 4, 5, 6, 7, 8, 9],
                        [1, 2, 0, 3, 4, 5, 6, 7, 8, 9],
                        [1, 2, 0, 3, 4, 5, 6, 7, 8, 9],
                        [0, 3, 1, 2, 4, 5, 6, 7, 8, 9],
                        [4, 5, 0, 3, 6, 1, 2, 7, 8, 9],
                        [4, 5, 0, 3, 6, 1, 2, 7, 8, 9],
                        [6, 8, 9, 7, 4, 5, 0, 3, 1, 2],
                        [7, 8, 9, 6, 0, 1, 2, 3, 4, 5],
                        [8, 9, 7, 6, 0, 1, 2, 3, 4, 5],
                        [8, 9, 7, 6, 0, 1, 2, 3, 4, 5]]
                    )

                def test_lose_sample_indexing(self):
                    """ Test a runthrough on simple sample data with 
                        enough junctions excluded by sample threshold
                        that some samples are entirely lost
                    """
                    with gzip.open(self.input_file, 'w') as input_stream:
                        input_stream.write(
                            self.lossy_input_string
                        )
                    with open(self.meta_file, 'w') as meta_stream:
                        meta_stream.write(
                            self.meta_string
                        )
                    
                    index_path = os.path.join(self.temp_dir_path ,"tempIndex")
                    go_index(intropolis=self.input_file, basename=index_path,
                        features=3000,n_trees=20,
                        sample_count=10,sample_threshold=6,
                        buffer_size=1024,verbose=False,
                        metafile=self.meta_file)
                    
                    test_index=AnnoyIndex(3000)
                    test_index.load(index_path+".annoy.mor")
                    
                    self.assertEqual(test_index.get_n_items(), 7)
                    
                    expected_results = (
                        [[0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6],
                        [0, 1, 2, 3, 4, 5, 6]]
                    )


                    for i in range(10):
                        self.assertEqual(
                            test_index.get_nns_by_item(i,10,search_k=100,
                                include_distances=False),
                            expected_results[i])
                    
                def tearDown(self):
                    # Kill temporary directory
                    shutil.rmtree(self.temp_dir_path)

            unittest.main()
            
    else:
    
        if args.subparser_name == 'index':
            #index case
            go_index(args.intropolis,args.basename,args.features,args.n_trees,
                    args.sample_count,args.sample_threshold,args.buffer_size,
                    args.verbose,args.metafile)
        else:
            #either align or search case
            junction_stream = sys.stdin 
            #this will be changed to the output sam file 
            #from first pass alignment in the morna align case
        
            if args.subparser_name == 'align':
                if args.previous_alignment:
                    args.format = "sam"
                    junction_stream = open(args.pass1_sam)
                else:
                    # Check command-line parameters
                    if args.unpaired is not None:
                        if args.mates_1 is not None or args.mates_2 is not None:
                            raise RuntimeError(
                                    'Cannot align both paired & unpaired reads'
                                    'at once.'
                                )
                    elif (args.mates_1 is not None and args.mates_2 is None 
                            or args.mates_2 is not None 
                                                and args.mates_1 is None):
                        raise RuntimeError(
                                'If analyzing paired-end reads, must specify '
                                'paths to both mates files.'
                            )
                    elif len(args.mates_1) != len(args.mates_2):
                        raise RuntimeError(
                                'The numbers of #1-mate and #2-mate files must'
                                'be equal.'
                            )
                    
                    command = []
                
                    if args.aligner == "STAR":
                        command.append("STAR")
                        command += ["--genomeDir", args.index]
                        command += ["--outFileNamePrefix", args.pass1_sam]
                        command.append("--readFilesIn")
                        if (args.mates_1 is not None 
                                                and args.mates_2 is not None):
                            command += args.mates_1
                            command += args.mates_2
                        else:
                            command += args.unpaired
                    
                    
                    elif (args.aligner == "hisat2") or (args.aligner 
                                                            == "hisat"):
                        command.append(args.aligner)
                        command += ["-x", args.index]
                        command += ["-S",args.pass1_sam]
                        if (args.mates_1 is not None 
                                                and args.mates_2 is not None):
                            command += ["-1", ",".join(args.mates_1)]
                            command += ["-2", ",".join(args.mates_2)]
                        else:
                            command += ["-U", ",".join(args.unpaired)]            
                    else:
                        raise RuntimeError(
                                'Currently supported aligner options are STAR,'
                                'hisat, and hisat2.'
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
                #MAYBE BABY
                results = searcher.search_member_n(args.query_id,
                                        args.results, args.search_k,
                                        include_distances=args.distances,
                                         meta_db = args.metadata)
                results_output(results)
                quit()
            # Read BED or BAM
            if args.format == "sam":
                junction_generator = junctions_from_sam_stream(junction_stream)
            elif args.format == "bed":
                junction_generator = junctions_from_bed_stream(junction_stream)
            else:
                assert args.format == "raw"
                junction_generator = junctions_from_raw_stream(junction_stream)
            if (args.convergence_backoff) and (args.format == 'sam'):
                backoff = args.convergence_backoff
                checkpoint = args.checkpoint
                old_results = [-1 for _ in range(args.results)]
                if args.verbose:
                    for i, junction in enumerate(junction_generator):
                        if (i % 1000 == 0):
                            sys.stderr.write( str(i) 
                                             + " junctions into query sample\r")
                            sys.stderr.flush()
                        if (" ".join(junction[:3]) in 
                                searcher.sample_frequencies):
                            searcher.update_query(junction)
                        if (i == checkpoint):
                            checkpoint += (backoff)
                            backoff += backoff
                            sys.stderr.write("\n")
                            searcher.finalize_query()
                            results = (searcher.search_nn(args.results, 
                                        args.search_k,
                                        include_distances=args.distances,
                                         meta_db = args.metadata))
                            same = True
                            for j,result in enumerate(results[0]):
                                if not (result == old_results[j]):
                                    same = False

                            if same:
                                sys.stderr.write("Converged after "+ str(i) 
                                                    + " junctions.\n")
                                results_output(results)
                            
                                quit()
                            else:
                                sys.stderr.write("Not converged after " + str(i) 
                                                    + " junctions.\n")
                                sys.stderr.write("Old results:\n" 
                                    + str(old_results) + "\n")
                                sys.stderr.write("New results:\n" 
                                    + str(results[0]) + "\n")
                                old_results = results[0]
                    sys.stderr.write("No convergence after " + str(i) 
                                            + " junctions, but here's results:\n")
                    searcher.finalize_query()
                    results = (searcher.search_nn(args.results, 
                                    args.search_k,
                                    include_distances=args.distances,
                                     meta_db = args.metadata))
                    results_output(results)
                else:
                    for i, junction in enumerate(junction_generator):
                        #TODO: propagate this to see that we get dead on results
                        #with annoy
                        if (" ".join(junction[:3]) in 
                                searcher.sample_frequencies):
                            searcher.update_query(junction)
                        if (i == checkpoint):
                            checkpoint += backoff
                            backoff += backoff
                            sys.stderr.write("\n")
                            searcher.finalize_query()
                            results = (searcher.search_nn(args.results, 
                                        args.search_k,
                                        include_distances=args.distances,
                                         meta_db = args.metadata))
                            same = True
                            for j,result in enumerate(results[0]):
                                if not (result == old_results[j]):
                                    same = False

                            if same:
                                results_output(results)
                                quit()
                            else:
                                old_results = results[0]
                    results_output(results)
            
            else:
                if args.verbose:
                    for i, junction in enumerate(junction_generator):
                        if (i % 1000 == 0):
                            sys.stderr.write( str(i) 
                                             + " junctions into query sample\r")
                            sys.stderr.flush()
                        if (" ".join(junction[:3]) in 
                                searcher.sample_frequencies):
                            searcher.update_query(junction)
                    searcher.finalize_query()
                else:
                    for i, junction in enumerate(junction_generator):
                        if (" ".join(junction[:3]) in 
                                searcher.sample_frequencies):
                            searcher.update_query(junction)
                    searcher.finalize_query()
                
                if args.verbose:
                    sys.stderr.write("\n")
                if args.exact:
                    results = searcher.exact_search_nn(args.results,
                     include_distances=args.distances, meta_db = args.metadata)
                else:
                    results = (searcher.search_nn(args.results, args.search_k,
                                include_distances=args.distances,
                                meta_db = args.metadata))
                results_output(results)
            
            if args.subparser_name == 'align':
                filter = args.junction_filter.split(",")
                frequency_filter = float(filter[0])
                coverage_filter = int(filter[1])
                max_junctions = -1
                #print presence_filter
                #print coverage_filter
                result_sample_ids=[]
                result_juncs = []
                result_covrs = []
    #            Two part junction filter settings separated by comma. 
    #            Only retain junctions for alignment found in at least 
    #            {first part} proportion of result samples, or 
    #            junctions with at least {second part} coverage 
    #            in any one result sample.
                for result in results[0]:
                    #the results are internal ids; we need sample ids
                    sample_id = searcher.inverse_lookup(result)
                    result_sample_ids.append(sample_id)
                    shard_id = mmh3.hash(str(sample_id)) % 100
                    print "shard_id is " + format(shard_id, '02d')
                    database = (args.basename + ".sh" + format(shard_id, '02d') 
                                                                + ".junc.mor")
                    conn=sqlite3.connect(database)
                    c=conn.cursor()
                    #print("Connected to " + database)
                    #The following makes this_one_juncs a list 
                    #of '!'-separated lists of based64 encoded run lengths
                    #and it makes this_one_covrs a list of ','-separated lists
                    #of coverage integers
                    this_one_juncs, this_one_covrs = zip(*list(c.execute(
                            ("SELECT * FROM sample_%d") % sample_id
                        )))
                
                    #Then we join them, making this_one_juncs a string which is
                    #a single !-separated list of based64 encoded run lengths
                    #and this_one_covrs a string representing a ,-separated list
                    #of coverage integers
                    this_one_juncs="".join(this_one_juncs)
                    this_one_covrs="".join(this_one_covrs)
                    #then we split both on their separator and extract junction 
                    #indexes from this_one_juncs, adding a list of integer 
                    #junction indexes to result_juncs and a 
                    #list of integer coverages to result_covrs
                    result_juncs.append([j for j in
                                     running_sum(this_one_juncs.split("!"))])
                    result_covrs.append(this_one_covrs.strip(",").split(","))
                    conn.close()
                print "result_juncs lengths: "
                print [len(ls) for ls in result_juncs]
                print "result_covrs lengths: "
                print [len(ls) for ls in result_covrs]
            
                retain_junctions = set()
            
                frequency_counts = defaultdict(int)
                if args.splicefile is not None:
                    found_in_map = defaultdict(list) # maps junction index to list
                                                     # of results it was found in
                                             
                    #we need results -> sample ids
                    #result_sample_ids[result_number] = result sample id
                    #this next gets the minimum number of samples with a junction
                    #to pass the frequency filter; we can stop doing anything
                    #for that junction once we have passed this threshold
                    min_count = int(ceil(frequency_filter * len(result_juncs)))
                    for i, junction_list in enumerate(result_juncs):
                        for index in junction_list:
                            frequency_counts[index]+=1
                            found_in_map[index].append(i)
                    for i, junction_list in enumerate(result_juncs):
                        for index in frequency_counts:
                            if frequency_counts[index] >= min_count:
                                retain_junctions.add(index)
                else:
                    #this next gets the minimum number of samples 
                    #with a given junction to pass the frequency filter;
                    # we can stop doing anything for that junction 
                    #once we have passed this threshold
                    min_count = int(ceil(frequency_filter * len(result_juncs)))
                    for i, junction_list in enumerate(result_juncs):
                        for index in junction_list:
                            if frequency_counts[index]<min_count:
                                frequency_counts[index]+=1
                            elif frequency_counts[index]==min_count:
                                frequency_counts[index]+=1
                                #print("Retaining junction " + str(index) 
                                #        + " for passing frequency filter")
                                retain_junctions.add(index)
                            
                for i, coverages_list in enumerate(result_covrs):
                    for j, coverage in enumerate(coverages_list):
                        if int(coverage) >= coverage_filter:
                            #Mark the jth (0 based) index 
                            #from the list of junctions 
                            #indexes present in this sample for retention since 
                            #its corresponding coverage passes the filter.
                            retain_junctions.add(result_juncs[i][j])
                            #print("Retaining junction " 
                            #           + str(result_juncs[i][j]) 
                            #        + " for passing coverage filter")
                print("Number of retained junctions: " 
                                            + str(len(retain_junctions)))
                #~/Downloads/hisat2-2.0.5/hisat2 
                #-x /Users/andechri/Downloads/alignment/hg38/genome 
                #-S ${i}hisat_out.sam -1 ${i}62.filtered.fastq 
                #-2 ${i}63.filtered.fastq 
                #--novel-splice-infile ~/Downloads/splicesites.txt            
            
                if args.splicefile is not None:
                    with open(args.splicefile, "w") as splices,\
                                  gzip.open(
                                  args.junction_file) as junction_names:
                        ordered_junctions = sorted(retain_junctions)
                        sys.stderr.write(str(len(ordered_junctions)) 
                                        + " junctions to begin with\n")
                        junction_index = ordered_junctions.pop(0)
                
                        sys.stderr.flush()                
                        for i, line in enumerate(junction_names):
                            if i == junction_index:
                            
                                #result_sample_ids maps result index 
                                #(0-19 for 20 results) to the sample ids 
                                #those results correspond to
                                #retain sample ids will be, crucially, 
                                #in the same order as the indexes in
                                # found_in_map[i], so they will correspond 
                                #1:1 with the result numbers
                                retain_sample_ids = []
                                for result_index in found_in_map[i]:
                                    retain_sample_ids.append(result_sample_ids[
                                                                  result_index])
               #chr1 438529 544346 - GC AG 2407,3250,5952,6791 2,1,2,3	[3]
                                tokens = line.strip().split("\t")
                                tokens[1] = str(int(tokens[1])-2)
                                ####TODO: use just found in map to reconstruct!
                                #new_samples = retain_sample_ids
                                #new_coverages =old_covs[
                                #        old_samples.index(sample_id)]
                                ####
                                old_samples = [int(x) 
                                                for x in tokens[6].split(",")]
                                old_covs = [int(x) 
                                                for x in tokens[7].split(",")]
                                new_samples = []
                                new_covs = []
                                for sample_id in retain_sample_ids:
                                    if sample_id in old_samples:
                                        new_samples.append(sample_id)
                                        new_covs.append(old_covs[
                                                old_samples.index(sample_id)])
                                tokens[6] = ",".join([str(_) 
                                                        for _ in new_samples])
                                tokens[7] = ",".join([str(_) for _ in new_covs])

                                splices.write("\t".join(tokens) + "\t" 
                                        + str(found_in_map[i]) + "\n")
                                try:
                                    junction_index = ordered_junctions.pop(0)
                                except IndexError:
                                    #if the list of junction to match is empty, 
                                    #stop looping
                                    break 
                                #sys.stderr.write(str(len(ordered_junctions)) 
                                #        + " junctions remain\r")
                else:
                    with open(args.pass1_sam 
                            + "_splicesites.txt", "w") as splices,\
                                  gzip.open(
                                        args.junction_file) as junction_names:
                        ordered_junctions = sorted(retain_junctions)
                        sys.stderr.write(str(len(ordered_junctions)) 
                                        + " junctions to begin with\n")
                        junction_index = ordered_junctions.pop(0)
                
                        sys.stderr.flush()                
                        for i, line in enumerate(junction_names):
                            if i == junction_index:
                                tokens = line.strip().split("\t")[:4]
                                tokens[1] = str(int(tokens[1])-2)
                                splices.write("\t".join(tokens) + "\n")
                                try:
                                    junction_index = ordered_junctions.pop(0)
                                except IndexError:
                                    #if the list of junction to match is empty, 
                                    #stop looping
                                    break 
                                #sys.stderr.write(str(len(ordered_junctions)) 
                                #        + " junctions remain\r")
                        
                    print ""
            
                    command = []
                
                    if args.aligner == "STAR":
                        command.append("STAR")
                        command += ["--genomeDir", args.index]
                        command += ["--outFileNamePrefix", args.pass2_sam]
                        command.append("--readFilesIn")
                        if (args.mates_1 is not None 
                                                and args.mates_2 is not None):
                            command += args.mates_1
                            command += args.mates_2
                        else:
                            command += args.unpaired
                    
                    
                    elif (args.aligner == "hisat2") or (args.aligner 
                                                                == "hisat"):
                        command.append(args.aligner)
                        command += ["-x", args.index]
                        command += ["-S",args.pass2_sam]
                        if (args.mates_1 is not None 
                                                and args.mates_2 is not None):
                            command += ["-1", ",".join(args.mates_1)]
                            command += ["-2", ",".join(args.mates_2)]
                        else:
                            command += ["-U", ",".join(args.unpaired)]            
                    else:
                        raise RuntimeError(
                                'Currently supported aligner options are '
                                'STAR, hisat, and hisat2.'
                            )
                    if args.aligner_args is not None:
                        for arg in args.aligner_args.split(" "):
                            command.append(arg)

                    command += ["--novel-splicesite-infile", args.pass1_sam 
                                                        + "_splicesites.txt"]
                    if args.verbose:
                        print("Calling: " 
                                + " ".join(command))
                    ret = subprocess.check_call(command)
                    junction_stream.close()
