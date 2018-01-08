#!/usr/bin/env python

'''
Based on an input FASTA file, calculates a chromosome sizes file of the type
needed by e.g. wigToBigWig, i.e. a tab-delimited list.
'''

from __future__ import print_function
import sys
import argparse
from Bio import SeqIO

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin,
                        help='FASTA input file to convert [default: STDIN]')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_cl()
    for record in SeqIO.parse(args.infile, "fasta"):
        print('\t'.join((record.id, str(len(record.seq)))))
    