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
import cPickle
import gzip
import madoka
import mmh3
import sys
from annoy import AnnoyIndex
from collections import defaultdict
from math import log


def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers
        md: an MD string (example: 33A^CC).
        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
        if i == 0: continue
        if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
                md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
            md_group.append(char)
    if md_group:
        md_to_parse.append(''.join(md_group))
    return [char for char in md_to_parse if char != '0']

def indels_junctions_exons_mismatches(
            cigar, md, pos, seq, drop_deletions=False, junctions_only=False
        ):
    """ Finds indels, junctions, exons, mismatches from CIGAR, MD string, POS
    
        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence
        drop_deletions: drops deletions from coverage vectors iff True
        junctions_only: does not populate mismatch list
        Return value: tuple
            (insertions, deletions, junctions, exons, mismatches).
        Insertions is a list of tuples (last genomic position before insertion, 
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). Junctions is a list
            of tuples (junction start position (inclusive),
                       junction end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)). Mismatches is a list
            of tuples (genomic position of mismatch, read base)
    """
    insertions, deletions, junctions, exons, mismatches = [], [], [], [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    seq_size = len(seq)
    cigar_chars, cigar_sizes = [], []
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            aligned_base_cap = int(cigar[cigar_index])
            aligned_bases = 0
            while True:
                try:
                    aligned_bases += int(md[md_index])
                    if aligned_bases <= aligned_base_cap:
                        md_index += 1
                except ValueError:
                    # Not an int, but should not have reached a deletion
                    assert md[md_index] != '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
                    if not junctions_only:
                        mismatches.append(
                                (pos + aligned_bases,
                                    seq[seq_index + aligned_bases])
                            )
                    correction_length = len(md[md_index])
                    m_length = aligned_base_cap - aligned_bases
                    if correction_length > m_length:
                        md[md_index] = md[md_index][:m_length]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += correction_length
                        md_index += 1
                if aligned_bases > aligned_base_cap:
                    md[md_index] = aligned_bases - aligned_base_cap
                    break
                elif aligned_bases == aligned_base_cap:
                    break
            # Add exon
            exons.append((pos, pos + aligned_base_cap))
            pos += aligned_base_cap
            seq_index += aligned_base_cap
        elif cigar[cigar_index+1] == 'N':
            skip_increment = int(cigar[cigar_index])
            # Add junction
            junctions.append((pos, pos + skip_increment,
                            seq_index, seq_size - seq_index))
            # Skip region of reference
            pos += skip_increment
        elif cigar[cigar_index+1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append(
                    (pos - 1, seq[seq_index:seq_index+insert_size])
                )
            seq_index += insert_size
        elif cigar[cigar_index+1] == 'D':
            assert md[md_index] == '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
            # Deletion
            delete_size = int(cigar[cigar_index])
            md_delete_size = len(md[md_index+1])
            assert md_delete_size >= delete_size
            deletions.append((pos, md[md_index+1][:delete_size]))
            if not drop_deletions: exons.append((pos, pos + delete_size))
            if md_delete_size > delete_size:
                # Deletion contains a junction
                md[md_index+1] = md[md_index+1][delete_size:]
            else:
                md_index += 2
            # Skip deleted part of reference
            pos += delete_size
        else:
            # Soft clip
            assert cigar[cigar_index+1] == 'S'
            # Advance seq_index
            seq_index += int(cigar[cigar_index])
        cigar_index += 2
    '''Merge exonic chunks/deletions; insertions/junctions could have chopped
    them up.'''
    new_exons = []
    last_exon = exons[0]
    for exon in exons[1:]:
        if exon[0] == last_exon[1]:
            # Merge ECs
            last_exon = (last_exon[0], exon[1])
        else:
            # Push last exon to new exon list
            new_exons.append(last_exon)
            last_exon = exon
    new_exons.append(last_exon)
    return insertions, deletions, junctions, new_exons, mismatches

def dummy_md_index(cigar):
    """ Creates dummy MD string from CIGAR in case of missing MD.
        cigar: cigar string
        Return value: dummy MD string
    """
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    cigar_index = 0
    max_cigar_index = len(cigar)
    md = []
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            try:
                if type(md[-1]) is int:
                    md[-1] += int(cigar[cigar_index])
                else:
                    md.append(int(cigar[cigar_index]))
            except IndexError:
                md.append(int(cigar[cigar_index]))
            cigar_index += 2
        elif cigar[cigar_index+1] in 'SIN':
            cigar_index += 2
        elif cigar[cigar_index+1] == 'D':
            md.extend(['^', 'A'*int(cigar[cigar_index])])
            cigar_index += 2
        else:
            raise RuntimeError(
                        'Accepted CIGAR characters are only in [MINDS].'
                    )
    return ''.join(str(el) for el in md)

def junctions_from_raw_stream(raw_stream):
    """ Generates junctions from raw stream
    
        Raw format has four fields: junction's (chromosome, start (1-based,
        inclusive), end (1-based, inclusive)), coverage
        
        Yield value: tuple (chrom, start position, end position, coverage)
    """
    for line in raw_stream:
        tokens = line.strip().split('\t')
        yield (tokens[0], tokens[1], tokens[2], int(tokens[3]))

def junctions_from_bed_stream(bed_stream):
    """ Generates junctions from BED stream

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions

        Yield value: tuple (chrom, start position, end position, coverage)
            representing a junction. Start position is 1-based and inclusive;
            end position is 1-based and inclusive. 
    """
    for line in bed_stream:
        tokens = line.rstrip().split('\t')
        if len(tokens) < 12:
            continue
        chrom = tokens[0]
        chrom_start = int(tokens[1])
        chrom_end = int(tokens[2])
        coverage = int(tokens[4])
        block_sizes = tokens[10].split(',')
        block_starts = tokens[11].split(',')
        # Handle trailing commas
        try:
            int(block_sizes[-1])
        except ValueError:
            block_sizes = block_sizes[:-1]
        try:
            int(block_starts[-1])
        except ValueError:
            block_starts = block_starts[:-1]
        block_count = len(block_sizes)
        if block_count < 2:
            # No junctions
            continue
        assert block_count == len(block_starts)
        junctions = []
        # First block characterizes junction on left side of junction
        junctions.append(chrom_start + int(block_starts[0]) 
                                + int(block_sizes[0]))
        for i in xrange(1, block_count - 1):
            # Any intervening blocks characterize two junctions
            junction_start = chrom_start + int(block_starts[i])
            junctions.append(junction_start)
            junctions.append(junction_start + int(block_sizes[i]))
        # Final block characterizes junction on right side of junction
        junctions.append(chrom_start + int(block_starts[-1]))
        for i in xrange(len(junctions)/2):
            yield (chrom, junctions[2*i]+1, junctions[2*i+1]+1, coverage)

def junctions_from_sam_stream(sam_stream):
    """ Generates junctions from BED stream

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions

        Yield value: tuple (chrom, start position, end position, 1)
            representing a junction. Start position is 1-based and inclusive;
            end position is 1-based and inclusive.
    """
    for line in sam_stream:
        if line[0] == '@': continue
        try:
            tokens = line.strip().split('\t')
            flag = int(tokens[1])
            if flag & 4:
                continue
            name = tokens[0]
            rname = tokens[2]
            cigar = tokens[5]
            pos = int(tokens[3])
            seq = tokens[9]
            flag = int(tokens[1])
            if 'N' not in cigar or flag & 256:
                continue
            #md = [token[5:] for token in tokens if token[:5] == 'MD:Z:'][0]
            _, _, junctions_to_add, _ = indels_junctions_and_exons(
                        cigar, dummy_md_index(cigar), pos, seq
                    )
            for junction in junctions_to_add:
                yield (rname,) + junction + (1,)
        except IndexError:
            print >>sys.stderr, ('Error found on line: ' + line)
            raise

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
    parser.add_argument('-s','--num-samples', type=int, required=False,
            default=None,
            help=('optionally specify number of unique sample ids')
        )
        
    args = parser.parse_args()

    if args.intropolis:
        #First, sample count to facilitate IDF calculation
        #list is initialized with a value such that new elements will
        #ALWAYS sort to the left of it; this eliminates a check for empty list 
        #from sorted_unique_insert which I'm trying to optimize
        all_samples = set()
        
        num_samples = args.num_samples
        
        if not args.num_samples:
            with gzip.open(args.intropolis) as introp_file_handle:
                for i, line in enumerate(introp_file_handle):
                    if (i % 100 == 0):
                        sys.stdout.write(str(i) + " lines into sample count, " 
                                + str(len(all_samples)) + " samples so far.\r")
                    line_pieces = line.split()
                    samples_with_junction = (line_pieces[6].split(','))
                    for sample in samples_with_junction:
                        all_samples.add(sample)

            print("")
            num_samples = len(all_samples)
            print("Samples counted:")
            
        print("There are " + str(num_samples) + " samples.")
        # Index
        sample_feature_matrix = defaultdict(
                                lambda: [0.0 for _ in xrange(args.features)])
        
        with open(args.basename + ".stats.mor", 'w') as num_fh:
            num_fh.write(str(num_samples))
        
        junction_sketch = madoka.CroquisUint16()
        
        with gzip.open(args.intropolis) as introp_file_handle:
            for i, introp_line in enumerate(introp_file_handle):
                if (i % 100 == 0):
                    sys.stdout.write(str(i) + " lines into index making\r")
                introp_line_pieces = introp_line.split()
                #right now we hash on 'chromosome start stop'
                #maybe strand someday but shouldn't matter
                hashable_junction = " ".join(introp_line_pieces[:3])
                hashed_value = mmh3.hash(hashable_junction)
                multiplier = (-1 if hashed_value < 0 else 1)
                hashed_value = hashed_value % args.features
                samples_junction_coverages = (introp_line_pieces[7].split(','))
                samples_junction_coverages = [int(num) for num in
                                             samples_junction_coverages]
                                     
                num_samples_with_junction = len(samples_junction_coverages)
                junction_sketch[hashable_junction
                                        ] += num_samples_with_junction
                idf_value = log(float(num_samples)/num_samples_with_junction)
                                                            
                for sample_index, coverage in zip(
                                    introp_line_pieces[6].split(','),
                                    samples_junction_coverages
                                ):
                    tf_idf_score = (coverage * idf_value)
                    #previously used 1+ norm_log(int(coverage))
                    sample_feature_matrix[
                            int(sample_index)][
                            hashed_value] += (multiplier * tf_idf_score)
        
        junction_sketch.save(arg.basename + ".sketch.mor")
        
        print("")
        
        annoy_index = AnnoyIndex(args.features)
        i=0;  
        for sample_index in sample_feature_matrix:
            if (i % 100 == 0):
                    sys.stdout.write( str(i) + " samples into annoy building\r")
            annoy_index.add_item(i,sample_feature_matrix[sample_index])
            i+=1
        print("")
        annoy_index.build(args.n_trees) # n trees specified by args
        annoy_index.save(args.basename + '.annoy.mor')
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
