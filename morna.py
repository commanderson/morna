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
import sys
from annoy import AnnoyIndex
from collections import defaultdict
from math import log
from utils import *

class MornaIndex(AnnoyIndex):
    """ Augments AnnoyIndex to accommodate parameters needed by morna

        Three index files are written when invoking save(): an Annoy index
        (basename.annoy.mor), a dictionary that maps each junction to the 
        number of samples in which it's found (basename.freq.mor), and
        a stats file including the total number of samples.

        See https://github.com/spotify/annoy/blob/master/annoy/__init__.py
        for original class. 
    """
    def __init__(self, sample_count, dim=3000, sample_threshold=100):
        super(MornaIndex, self).__init__(dim, metric='angular')
        # Store numbers of samples in which each junction is found
        self.sample_frequencies = defaultdict(int)
        # Store low-dimensional representations of samples
        self.sample_feature_matrix = defaultdict(lambda:
                                                [0.0 for _ in xrange(dim)]
                                        )
        # Dimension of low-dimensional feature feature vector
        self.dim = self.dimension_count = 3000
        # Minimum number of samples in which junction should be found
        self.sample_threshold = sample_threshold

    def add_junction(junction, samples, coverages):
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
        value = mmh3.hash(junction)
        multiplier = (-1 if value < 0 else 1)
        hashed_value = hashed_value % args.features
        idf_value = log(float(num_samples) / self.sample_frequencies[junction])                                           
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
                self.add_item(sample, sample_feature_matrix[sample])
                if not (i % 100):
                    print >>sys.stderr, (
                            'Added {} samples to Annoy index so far.'
                        ).format(i+1)
            print >>sys.stderr, (
                    'Added a total of {} samples to Annoy index.'
                ).format(i+1)
        else:
            for sample in self.sample_feature_matrix:
                self.add_item(sample, sample_feature_matrix)
        super(MornaIndex, self).build(n_trees)

    def save(self, basename):
        # Save Annoy index first
        super(MornaIndex, self).save(basename + '.annoy.mor')
        # Write total number of samples
        with open(basename + ".stats.mor", 'w') as stats_stream:
            print >>stats_stream, str(self.sample_count)
        # Pickle sample frequency dictionary
        with open(basename + ".freq.mor", 'w') as pickle_stream:
            cPickle.dump(self.sample_frequencies, pickle_stream,
                         cPickle.HIGHEST_PROTOCOL)

class MornaSearch(object):
    def __init__(self, dim=3000, basename):
        # Load AnnoyIndex with AnnoyIndex class, not MornaIndex

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--intropolis', type=str, required=False,
            default=None,
            help='path to gzipped intropolis; specify only when indexing'
        )
    parser.add_argument('-x', '--basename', type=str, required=False,
            default='morna',
            help='path to junction index for either writing or reading'
        )
    parser.add_argument('--features', type=int, required=False,
            default=3000,
            help='dimension of feature space'
        )
    parser.add_argument('--n-trees', type=int, required=False,
            default=200,
            help='number of annoy trees'
        )
    parser.add_argument('--search-k', type=int, required=False,
            default=None,
            help='larger = more accurage search'
        )
    parser.add_argument('--format', type=str, required=False,
            default='sam',
            help=('one of {sam, bed, raw}')
        )
    parser.add_argument('-s','--sample-count', type=int, required=False,
            default=None,
            help=('optionally specify number of unique sample ids')
        )
    parser.add_argument('-t', '--sample-threshold', type=int, required=False,
            default=100,
            help='minimum number of samples in which a junction should appear '
                 'for it to be included in morna index')
    parser.add_argument('-v', '--verbose', action='store_const', const=True,
            default=False,
            help='be talkative'
        )
        
    args = parser.parse_args()

    if args.intropolis:
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
                                + str(len(all_samples)) + " samples so far.\r"
                            )
                        samples.update(line.rpartition('\t')[-1].split(','))
                else:
                    for line in enumerate(introp_file_handle):
                        samples.update(line.rpartition('\t')[-1].split(','))
            args.sample_count = len(samples)
            print 'There are {} samples.'.format(args.sample_count)
        morna_index = MornaIndex(args.sample_count, dim=args.features,
                                    sample_threshold=args.sample_threshold)
        with gzip.open(args.intropolis) as introp_file_handle:
            if args.verbose:
                for i, line in enumerate(introp_file_handle):
                    if i % 100 == 0:
                        sys.stdout.write(str(i) + " lines into index making\r")
                    tokens = line.strip().split('\t')
                    morna_index.add_junction(
                            ' ' .join(tokens[:3]),
                            map(int, tokens[-2].split(',')), # samples
                            map(int, tokens[-1].split('',)) # coverages
                        )
            else:
                for line in introp_file_handle:
                    tokens = line.strip().split('\t')
                    morna_index.add_junction(
                            ' ' .join(tokens[:3]),
                            map(int, tokens[-2].split(',')), # samples
                            map(int, tokens[-1].split('',)) # coverages
                        )
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
                
        num_samples = 0
        with open(args.basename + ".stats.mor") as num_fh:
            for introp_line in enumerate(introp_file_handle):
                num_samples = int(line)
        
        junction_sketch = madoka.CroquisUint16()
        junction_sketch.load(arg.basename + ".sketch.mor")
        
        query_sample = [0 for _ in xrange(args.features)]
        for i, junction in enumerate(junction_generator):
            if (i % 100 == 0):
                sys.stderr.write( str(i) + " junctions into query sample\r")
            hashable_junction = ' '.join(map(str, junction[:3]))
            
            
            num_samples_with_junction = junction_sketch[
                                                hashable_junction]
            idf_value = log(float(num_samples)/num_samples_with_junction)
            hash_value = mmh3.hash(hashable_junction)
            multiplier = (-1 if hash_value < 0 else 1)
            query_sample[hash_value % args.features] += (
                        multiplier * (junction[3] * idf_value)
                    )
            
        print("")
        annoy_index = AnnoyIndex(args.features)
        annoy_index.load(args.annoy_idx)
        print annoy_index.get_nns_by_vector(
                                                [feature for feature in
                                                 query_sample], 
                                                 500,
                                                 args.search_k,
                                                 include_distances=True
                                            )
