#!/usr/bin/env python

''' Scans for poly-(A) (or poly-(T)) and outputs a BED file. '''

from __future__ import print_function
import sys
import argparse
from Bio import SeqIO

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='BedGraph input file [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='BedGraph output file [default: STDOUT]')
    parser.add_argument('--min_length', type=int, default=6,
                        help='minimum number of consecutive As [default: 6]')
    return parser.parse_args()

def scan(infile, outfile, min_length):
    for seq in SeqIO.parse(infile, "fasta"):
        countA = 0
        countT = 0
        for i in range(len(seq)):
            if seq[i] == 'A':
                countA = countA + 1
            else:
                if countA >= min_length:
                    recordCount(countA, i, '+', seq, outfile)
                countA = 0    
            if seq[i] == 'T':
                countT = countT + 1
            else:
                if countT >= min_length:
                    recordCount(countT, i, '-', seq, outfile)
                countT = 0

def recordCount(count, i, strand, seq, outfile):
    print('%s\t%i\t%i\t%s\t%i\t%s' % (seq.name, i-count, i, 'polyA', count, strand), file=outfile)

if __name__ == '__main__':
    args = parse_cl()
    scan(args.infile, args.outfile, args.min_length)