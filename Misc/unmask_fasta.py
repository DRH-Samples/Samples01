#!/usr/bin/env python

''' Convert all soft-masked (lower case) sequence letters in a FASTA file to upper case.
Actually, just skips lines starting with '>' and converts everything else to upper case.

'''

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('r'),
                        help='FASTA file to convert to upper-case')
args = parser.parse_args()

for line in args.infile:
    if line[0] == '>':
        sys.stdout.write(line)
    else:
        sys.stdout.write(line.upper())
