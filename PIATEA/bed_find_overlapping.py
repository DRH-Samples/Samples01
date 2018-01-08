#!/usr/bin/env python

''' Report any overlapping BED records. Uses 0-based, half-open indexing (the same as UCSC Browser).'''

from __future__ import print_function
import sys
import argparse
import collections

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='BedGraph input file [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='BedGraph output file [default: STDOUT]')
    return parser.parse_args()

def process(args):
    start0 = -1
    end0 = -1
    chrom0 = None
    completedChroms = []
    found = False
    for line in args.infile:
        d = line.strip().split('\t')
        chrom = d[0]
        start = int(d[1])
        end   = int(d[2])
        if chrom in completedChroms:
            exit("ERROR! Chroms not sorted. Sort first with (e.g.): sort -k1,1 -k2,2n bedfile >sorted_bedfile)\n" + line)
        if chrom != chrom0:
            completedChroms.append(chrom0)
            chrom0 = chrom
        else:
            if start < start0:
                exit("ERROR! Starts not sorted. Sort first with (e.g.): sort -k1,1 -k2,2n bedfile >sorted_bedfile)\n" + line)
            if start < end0:
                found = True
                print(line, file=args.outfile, end='')
                None
        start0 = start
        end0 = end
    if not found:
        print('No overlapping records found.', file=args.outfile)

if __name__ == '__main__':
    args = parse_cl()
    process(args)