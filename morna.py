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
from math import log
from utils import *

_help_intro = \
"""morna searches for known RNA-seq samples with exon-exon junction expression
patterns similar to those in a query sample.
"""

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
    def __init__(self, sample_count, dim=3000, 
                 sample_threshold=100, metafile = None):
        super(MornaIndex, self).__init__(dim, metric='angular')
        # Store the total number of samples represented in the index
        self.sample_count = sample_count
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
    
    def add_junction(self, junction, samples, coverages):
        """ Adds contributions of a single junction to feature vectors

            junction: string representation of junction
                (typically chromosome, 1-based inclusive start position,
                 1-based inclusive end position)
            samples: list of integer indices of samples
            coverages: list of integer coverages corresponding to samples

            No return value.
        """
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
        #Parse file at metafaile and save it into a metadata db
        #Metafile format (whitespace separated)
        #First column: sample id.
        #Any number of additional columns: keywords describing sample,
        #separated by whitespace.
        conn = sqlite3.connect(basename + '.meta.mor')
        cursor = conn.cursor()
        cursor.execute("CREATE TABLE metadata (sample_id real, keywords text)" )
        print self.metafile
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
        
        
            
    def search_nn(self, num_neighbors, search_k, include_distances=True):
        """ Uses the current query_sample to query the annoy index

            num_neighbors: an integer number of nearest neighbors to return
            search_k: an integer number of nodes to check while searching 
            the trees in the index; larger values lead to slower, more accurate 
            searches.
            include_distances: bolean indicator - should the distance to each 
            neighbor be included along with the sample id?
            
            Return value: a list of the sample ids of the neareast 
            neighbors of the query, optionally joined in a tuple by 
            a corresponding list of distances.
        """
        if include_distances:
            return self.annoy_index.get_nns_by_vector(
                                                 [feature for feature in
                                                 self.query_sample], 
                                                 num_neighbors,
                                                 search_k,
                                                 include_distances
                                                )
        else:
            return [self.annoy_index.get_nns_by_vector(
                                                 [feature for feature in
                                                 self.query_sample], 
                                                 num_neighbors,
                                                 search_k,
                                                 include_distances
                                                )]

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
    # Add command-line arguments
    index_parser.add_argument('--intropolis', metavar='<str>', type=str,
            required=True,
            help=('path to gzipped file recording junctions across samples '
                  'in intropolis format')
        )
    index_parser.add_argument('-x', '--basename', metavar='<str>', type=str,
            required=False,
            default='morna',
            help='basename path of junction index files to create'
        )
    search_parser.add_argument('-x', '--basename', metavar='<str>', type=str,
            required=True,
            help='basename path of junction index to search'
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
    search_parser.add_argument('--search-k', metavar='<int>', type=int,
            required=False,
            default=100,
            help='a larger value makes for more accurate search'
        )
    search_parser.add_argument('-f', '--format', metavar='<choice>', type=str,
            required=False,
            default='sam',
            help=('one of {sam, bed, raw}')
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
    search_parser.add_argument('-v', '--verbose', action='store_const',
            const=True,
            default=False,
            help='be talkative'
        )
    search_parser.add_argument('-c','--converge', metavar='<int>', type=int,
            required=False,
            default=None,
            help=('attempt to converge on a solution with concordance equal'
                  'to the given argument as a percentage (truncated to 0-100)')
        )
    index_parser.add_argument('-m', '--metafile',  metavar='<str>', type=str,
            required=False,
            default=None,
            help=('path to metadata file with sample index in first column'
                  'and other keywords in other columns, whitespace delimited')
        )
    args = parser.parse_args()

    if hasattr(args, 'intropolis'):
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
        
        if hasattr(args, 'metafile'):
            with open(args.metafile) as metafile_handle:
                    for i,line in enumerate(metafile_handle):
                        if (i == 0):
                            continue
            
        morna_index = MornaIndex(args.sample_count, dim=args.features,
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
        # Search
        # Read BED or BAM
        if args.format == 'sam':
            junction_generator = junctions_from_sam_stream(sys.stdin)
        elif args.format == 'bed':
            junction_generator = junctions_from_bed_stream(sys.stdin)
        else:
            assert args.format == 'raw'
            junction_generator = junctions_from_raw_stream(sys.stdin)
        
        
        searcher = MornaSearch(basename=args.basename) 
        
        if hasattr(args, 'converge'):
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
                                    include_distances=True))
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
            if args.verbose:
                sys.stderr.write("\n")
            results = (searcher.search_nn(20, args.search_k,
                       include_distances=True))
            
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
            results = (searcher.search_nn(20, args.search_k,
                       include_distances=True))
            if args.metafile:
                with open(args.metafile) as metafile_handle:
                    for i,line in enumerate(metafile_handle):
                        if (i == 0):
                            continue
                        print(line.partition(' '))
                        #print(results[0])
                        if int(line.partition(' ')[0]) in results[0]:
                            print line 
            else:
                print results