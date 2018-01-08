#!/usr/bin/env python

'''
Convert a Geo SOFT file into a text file for input into Bowtie2.
Parses out each 24-nt sequence from the input file and writes it
to stdout "Count" number of times, which is the raw abundance indicated
in the input file on the same line as the sequence. Ignores sequences
that are not 24-nt in length.
'''

import argparse
import sys
import re

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType('r'), help='SOFT file to process')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_cl()
    for line in args.infile:
        #print(line)
        match = re.match(r'[ACGT]', line , re.M)
        if match:
            seq, n = line.split()
            if len(seq) == 24:
                for x in range(int(n)):
                    print seq

