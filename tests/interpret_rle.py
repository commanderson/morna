#!/usr/bin/env python

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
                  'bit vector')
        )
args = parser.parse_args()

shr_str=args.string
print shr_str
bv=BitVector(size=0)
while len(shr_str)>0:
    print "starting with shr-str of:" +shr_str
    m = re.search("(.*)([!.])([0-o]+)$",shr_str)
    shr_str = m.groups()[0]
    if m.groups()[1] == '!':
            bit='1'
    else:
            bit='0'
    bv = BitVector(bitstring=(decode_64(m.groups()[2])*bit)) + bv
print "Original String:"
print args.string
print "Bit vector:"
print bv
