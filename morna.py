#!/usr/bin/env python
"""
morna

Uses mmh3 to implement feature hashing for dimensionality reduction of 
intropolis junction vectors by sample and uses Spotify's annoy library for
obtaining nearest neighbors.

Requires mmh3 (pip install mmh3) and Spotify's annoy (pip install annoy)
"""
import argparse
import bisect
import gzip
import cPickle
import mmh3
import sqlite3
import sys
from annoy import AnnoyIndex
from collections import defaultdict
from math import log, sqrt
from utils import *

_help_intro = \
"""morna searches for known RNA-seq samples with exon-exon junction expression
patterns similar to those in a query sample.
"""

            
def euclidean_distance(v1, v2):
    """ Calculate the eculcidean distance between 2 vectors/points
    
        v1 and v2: vectors/points of matching length
    """ 
    distance = [(a - b)**2 for a, b in zip(v1, v2)]
    distance = sqrt(sum(distance))
    return distance
    
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

        Four index files are written when invoking save(): an Annoy index
        (basename.annoy.mor), a dictionary that maps each junction to the 
        number of samples in which it's found (basename.freq.mor), a stats 
        file including the total number of samples (basename.stats.mor),
        and a sqlite database containing metadata info associated with 
        sample ids.

        See https://github.com/spotify/annoy/blob/master/annoy/__init__.py
        for original class. 
    """
    def __init__(self, sample_count, basename, dim=3000, 
                 sample_threshold=100, metafile=None):
        super(MornaIndex, self).__init__(dim, metric='angular')
        # Store the total number of samples represented in the index
        self.sample_count = sample_count
        #Store basename for index files
        self.basename = basename
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
        # list of junction ids per sample database connection
        self.junction_conn = sqlite3.connect(self.basename + '.junc.mor')
        self.junc_cursor = self.junction_conn.cursor()
        self.junc_id = -1
        
    def add_junction(self, junction, samples, coverages):
        """ Adds contributions of a single junction to feature vectors

            junction: string representation of junction
                (typically chromosome, 1-based inclusive start position,
                 1-based inclusive end position)
            samples: list of integer indices of samples
            coverages: list of integer coverages corresponding to samples

            No return value.
        """
        
        self.junc_id += 1
        
        for i,sample_id in enumerate(samples):
            self.junc_cursor.execute(("CREATE TABLE IF NOT EXISTS sample_%d "         
                                        +"(junction_id real, coverage real)") 
                            % sample_id)
            self.junc_cursor.execute("INSERT INTO sample_%d VALUES (%d, %d)"
                    % (sample_id,self.junc_id,coverages[i]))
        
        if (self.junc_id % 1000 == 0):
            self.junction_conn.commit()
            print >>sys.stderr, (
                            'Wrote junc_id {} into tables'
                        ).format(self.junc_id)
        
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
        self.junction_conn.commit()
        self.junction_conn.close()
        
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
            for i,sample_id in enumerate(results[0]):
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
            for i,sample_id in enumerate(results[0]):
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
            for i,sample_id in enumerate(results[0]):
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
            help=('attempt to converge on a solution with concordance equal'
                  'to the given argument as a percentage (truncated to 0-100)')
        )
    subparser.add_argument('-e', '--exact', action='store_const',
            const=True,
            default=False,
            help='search for exact nearest neighbor to query within morna index'
                 'rather than using annoy hyperplane division algorithm'
        )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_help_intro, 
                formatter_class=help_formatter)
    subparsers = parser.add_subparsers(help=(
                'subcommands; add "-h" or "--help" '
                'after a subcommand for its parameters')
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
            default='STAR'
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
    align_parser.add_argument('-x', '--index', metavar='<idx>', type=str,
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
                                    metafile=args.metafile)
        with gzip.open(args.intropolis) as introp_file_handle:
            if args.verbose:
                for i, line in enumerate(introp_file_handle):
                    if i % 100 == 0:
                        sys.stdout.write(str(i) + " lines into index making\r")
                    tokens = line.strip().split('\t')
                    morna_index.add_junction(
                            ' '.join(tokens[:3]),
                            map(int, tokens[-2].split(',')), # samples
                            map(int, tokens[-1].split(',')) # coverages
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
        searcher = MornaSearch(basename=args.basename) 
        # Search
        # The search by query id case cuts out early
        if args.query_id:
            results = searcher.search_member_n(args.query_id,20, args.search_k,
                                    include_distances=args.distances,
                                     meta_db = args.metadata)
            print results
            quit()
        # Read BED or BAM
        if args.format == 'sam':
            junction_generator = junctions_from_sam_stream(sys.stdin)
        elif args.format == 'bed':
            junction_generator = junctions_from_bed_stream(sys.stdin)
        else:
            assert args.format == 'raw'
            junction_generator = junctions_from_raw_stream(sys.stdin)
        if (args.converge) and (args.format == 'sam'):
            threshold = args.converge
            backoff = 100
            old_results = []
            if args.verbose:
                for i, junction in enumerate(junction_generator):
                    if (i % 100 == 0):
                        sys.stderr.write( str(i) 
                                         + " junctions into query sample\r")
                    searcher.update_query(junction)
                    if (i == backoff):
                        backoff += backoff
                        sys.stderr.write("\n")
                        results = (searcher.search_nn(20, args.search_k,
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
                        results = (searcher.search_nn(20, args.search_k,
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
                    searcher.update_query(junction)
            else:
                for i, junction in enumerate(junction_generator):
                    searcher.update_query(junction)
            
            if args.verbose:
                sys.stderr.write("\n")
            if args.exact:
                results = searcher.exact_search_nn(20,
                 include_distances=args.distances, meta_db = args.metadata)
            else:
                results = (searcher.search_nn(20, args.search_k,
                            include_distances=args.distances,
                            meta_db = args.metadata))
            print results
            
            