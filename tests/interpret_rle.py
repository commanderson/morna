#!/usr/bin/env python
""" interpret_rle.py
    interpret an input string as a !-separated list of run lengths 
    encoded in our base64 format using 0-9 as digits 0-9 
    and characters ':' thru 'o'as digits 10-63
"""

import argparse
import re
from BitVector import BitVector

def decode_64(s):
    """ Decodes an input base 64 number formatted
        using 0-9 as digits 0-9 and characters ':' thru 'o'
        as digits 10-63 into decimal int

        s: string of base 64 number in format above to convert
        
        return value: int decimal number form of input
    """
    return sum([(ord(s[idx])-48)
                    *(64**(len(s)-idx-1)) for idx in (xrange(len(s)))])


            
parser = argparse.ArgumentParser()
parser.add_argument('-s','--string', metavar='<str>', type=str,
            required=True,
            help=('String to interpret as run-length-encoded'
                  'junction list')
        )

args = parser.parse_args()

rls = args.string.split("!")
tot = 0
jlist = []
for i, item in enumerate(rls):
    length = decode_64(item)
    if i % 2 == len(jlist) % 2:
        tot += length
    else:
        for i in xrange(length):
            jlist.append(tot + i)
        tot += length
print "Input string: " + args.string
print "List of junctions: " + jlist
