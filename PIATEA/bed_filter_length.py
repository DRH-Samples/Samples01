#!/usr/bin/env python

'''
Default usage (only usage available for now): Specify minimum length.
'''

from __future__ import print_function
import sys
import argparse

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='BedGraph input file [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='BedGraph output file [default: STDOUT]')
    parser.add_argument('--min_length', type=int, default=0,
                        help='minimum allowed feature length (bp) [default: 0]')
    parser.add_argument('--min_score', type=float, default=0.0,
                        help='minimum allowed feature length (bp) [default: 0.0]')
    return parser.parse_args()

def filter_length(infile, outfile, min_length, min_score):
    for line in infile:
        if '#' == line[0]:              # keep comment lines
            print(line, file=outfile, end='')
        (chrom, start, end, name, score) = (line.split())[0:5]
        start = int(start); end = int(end); score = float(score)
        if (end - start > min_length and score >= min_score):
            print(line, file=outfile, end='')

if __name__ == '__main__':
    args = parse_cl()
    filter_length(args.infile, args.outfile, args.min_length, args.min_score)