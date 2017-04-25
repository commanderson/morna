#!/usr/bin/env python
"""
junction_recovery_performance.py

Outputs precision, recall, and related performance measurements of a spliced
aligner given a SAM file T with true junctions (typically from a simulator like
Flux) and SAM lines Y with spliced alignments retrieved by the aligner.

A line in an input BED file characterizes an junction by decking it with blocks.
A single junction defines two junctions, one on either side of it. A true
positive is an junction recovered by the aligner, a false positive is an junction
that appears in Y but not in T, and a false negative is an junction that appears
in T but not in Y.

SAM lines are read from stdin.

Output (written to stdout)
----------------------------
Two columns delimited by tabs, where the first column characterizes the numbers
in the second column. The first column is given below.

true junction count
retrieved junction count
true positive count
false positive count
false negative count
precision
recall
"""

import sys
from collections import defaultdict
import re

def junctions_from_bed_stream(bed_stream):
    """ Converts BED to dictionary that maps RNAMES to sets of junctions.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions.

        Return value: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an junction on RNAME. Each tuple is of the form
            (start position, end position).
    """
    to_return = defaultdict(int)
    for line in bed_stream:
        tokens = line.rstrip().split('\t')
        if len(tokens) < 12:
            continue
        chrom = tokens[0]
        chrom_start = int(tokens[1])
        chrom_end = int(tokens[2])
        block_sizes = tokens[10].split(',')
        block_starts = tokens[11].split(',')
        coverage = tokens[5]
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
            to_return[(chrom, junctions[2*i]+1, junctions[2*i+1]+1)] += 1
    return to_return

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

def indels_junctions_and_exons(cigar, md, pos, seq):
    """ Computes indels, junctions, and exons from CIGAR, MD string,
        and POS of a given alignment.

        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence

        Return value: tuple (insertions, deletions, junctions, exons). Insertions
            is a list of tuples (last genomic position before insertion, 
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). junctions is a list
            of tuples (junction start position (inclusive),
                       junction end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)).
    """
    insertions, deletions, junctions, exons = [], [], [], []
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
                    if aligned_bases + len(md[md_index]) > aligned_base_cap:
                        md[md_index] = md[md_index][
                                            :aligned_base_cap-aligned_bases
                                        ]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += len(md[md_index])
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
            if md_delete_size > delete_size:
                # Deletion contains an junction
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
    return insertions, deletions, junctions, new_exons

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

def junctions_from_sam_stream(sam_stream):
    """ Writes output that maps QNAMES to exon-exon junctions overlapped.

        sam_stream: where to find retrieved alignments in SAM form

        Return value: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a set of tuples, each
            denoting an junction on RNAME. Each tuple is of the form
            (start position, end position).
    """
    junctions = defaultdict(int)
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
            _, _, junctions_to_add, _ = indels_junctions_and_exons(cigar,
                                        dummy_md_index(cigar), pos, seq)
            for junction in junctions_to_add:
                junctions[(rname,) + junction[:2]] += 1
        except IndexError:
            print >>sys.stderr, ('Error found on line: ' + line)
            raise
    return junctions

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--true-junctions-bed', type=str, required=True, 
        help='Full path of BED file containing true junctions')
    args = parser.parse_args(sys.argv[1:])
    with open(args.true_junctions_bed) as true_junctions_bed_stream:
        true_junctions = junctions_from_bed_stream(true_junctions_bed_stream)
    retrieved_junctions = junctions_from_sam_stream(sys.stdin)
    (false_positive_count, false_negative_count, true_junction_count, 
        retrieved_junction_count) = [0]*4
    true_ones = set(
                [junction_and_coverage[0] for junction_and_coverage
                    in true_junctions.items()
                    if junction_and_coverage[1] == 1]
            )
    retrieved_ones = set(
                [junction_and_coverage[0] for junction_and_coverage
                    in retrieved_junctions.items()
                    if junction_and_coverage[1] == 1]
            )
    true_one_count = len(true_ones)
    retrieved_one_count = len(retrieved_ones)
    retrieved_true_one_count = len(true_ones.intersection(set(retrieved_junctions.keys())))
    recall_ones = float(retrieved_true_one_count) / true_one_count
    true_retrieved_one_count = len(retrieved_ones.intersection(set(true_junctions.keys())))
    precision_ones = float(true_retrieved_one_count) / len(retrieved_ones)
    print '# junctions in sim each covered by one read\t%d' % true_one_count
    print '# junctions in SAM each covered by one read\t%d' % retrieved_one_count
    print '# coverage-1 junctions recovered in SAM\t%d' % retrieved_true_one_count
    print '# coverage-1 junctions in SAM that were true positives\t%d' % true_retrieved_one_count
    print 'precision\t%.12f' % precision_ones
    print 'recall\t%.12f' % recall_ones
    print 'fscore\t%.12f' % (2 * precision_ones * recall_ones / (precision_ones + recall_ones))
