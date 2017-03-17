#!/usr/bin/env python

import argparse
import re
from BitVector import BitVector

parser = argparse.ArgumentParser()
parser.add_argument('-s','--string', metavar='<str>', type=str,
            required=True,
            help=('String to interpret as run-length-encoded'
                  'bit vector')
        )
args = parser.parse_args()

shr_str=args.string
bv=BitVector(size=0)
while len(shr_str)>0:
    print "starting with shr-str of:" +shr_str
    m = re.search("(.*)([io])(\d+)$",shr_str)
    shr_str = m.groups()[0]
    if m.groups()[1] == 'i':
            bit='1'
    else:
            bit='0'
    bv = BitVector(bitstring=(int(m.groups()[2])*bit)) + bv
print "Original String:"
print args.string
print "Bit vector:"
print bv
