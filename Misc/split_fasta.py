#!/usr/bin/env python

'Split a multi-FASTA file into multiple FASTA files.'

import sys
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('rU'),
                        help='FASTA file to split')
args = parser.parse_args()

header_re = re.compile('^>(\w+)\s.*')
first_line = True
outfile = None

for line in args.infile:
    if first_line:
        first_line = False
        if not header_re.match(line):
            raise Exception("This doesn't look like a FASTA file!")
    
    match = header_re.match(line)
    if match:
        if outfile is not None:
            outfile.close()
        outfile = open(match.group(1)+".fa", 'w')
        
    outfile.write(line)
    
outfile.close()