#!/usr/bin/env python
"""
process_junctionlist

Takes a master junction list created by morna align -sf and turns 
it into an appropriately formatted junction list with given filtering restrictions
"""


import argparse
import sys
import ast

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='<file>', type=str,
        required=True,
        help='path to master junction list input file from '
                 'morna align -sf; format is tab separated as follows:'
                 'chr, start, end, strand, donor, acceptor, ...'
                 'sample_ids, coverages, [ranks of results with junction]'
        )
parser.add_argument('-n', '--minnum', metavar='<int>', type=int,
        required=False,
        default=-1,
        help='minimum number of results a junction must have appeared in'
                 'junctions not meeting this are filtered from output'
                 
        )
parser.add_argument('-r', '--requiredrank', metavar='<int>', type=int,
        required=False,
        default=-1,
        help='A junction must have appeared in a rank at least this low; '
                 'junctions not meeting this are filtered from output'
        )
parser.add_argument('-c', '--mincoverage', metavar='<int>', type=int,
        required=False,
        default=-1,
        help='minimum coverage a junction must have in any sample; '
                 'junctions not meeting this are filtered from output'
        )
parser.add_argument('--junction-filter', type=str, required=False,
        default=None,
        help='Two part junction filter settings separated by comma. Only retain' 
             'junctions for alignment found in at least {first part} number'
             'of result samples, OR junctions with at least {second part}'
             'coverage in any one result sample. Mimics morna align argument'
             'except first part is quantity not proportion;'
             'example value is 1,5'
        )
args = parser.parse_args()

with open(args.input) as input_handle:
    for i, line in enumerate(input_handle):
        #if i % 1000 == 0:
        #    sys.stderr.write(str(i) 
        #            + " lines into index making\r")
        tokens = line.strip().split('\t')
        in_results = ast.literal_eval(tokens[-1])
        if args.junction_filter is not None:
            filter = args.junction_filter.split(",")
            frequency_filter = float(filter[0])
            coverage_filter = int(filter[1])
            if len(in_results) < 2:
                covs = [ast.literal_eval(tokens[-2])]
            else:
                covs = list(ast.literal_eval(tokens[-2]))
            if ((len(in_results) >= float(filter[0])) 
                    or (min(covs) >= int(filter[1]))):
               print "\t".join(tokens[:4])
               continue 
        if args.minnum>-1:
            if len(in_results) < args.minnum:
                continue
        if args.requiredrank>-1:
            if min(in_results) > args.requiredrank:
                continue
        if args.mincoverage>-1:
            if len(in_results) < 2:
                covs = [ast.literal_eval(tokens[-2])]
            else:
                covs = list(ast.literal_eval(tokens[-2]))
            if min(covs) < args.mincoverage:
                continue
        print "\t".join(tokens[:4])