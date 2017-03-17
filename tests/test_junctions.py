#!/usr/bin/env python
import re
from BitVector import BitVector
orig_str="o3i4o1i1o1i1"
shr_str=orig_str
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
print orig_str
print bv
