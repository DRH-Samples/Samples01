#!/usr/bin/env python

'''
Convert BED to BEDGRAPH by counting the number of BED segments overlapping each base
(then compressing equal-valued segments to one line).

Uses 0-based, half-open coordinates for consistency with the UCSC browser
(see http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms).

Note: This may use a lot of memory.
'''

from __future__ import print_function
import sys
import array
import argparse

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='LASTZ output file to convert [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='LASTZ output file to convert [default: STDOUT]')
    parser.add_argument('--subtract_1', action='store_true', default=False,
                        help='subtract 1 from all positions (to eliminate self-alignment) [default: False]')
    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_cl()
    counts = {}
    
    for line in args.infile:
        if '#' == line[0]:              # ignore comment lines
            continue
        (chrom, start, end) = (line.split())[0:3]
        start = int(start); end = int(end)
        if chrom not in counts:
            counts[chrom] = array.array('I')  # unsigned int (long)
        if len(counts[chrom]) < end:
            counts[chrom].extend([0]*(end-len(counts[chrom])))
        for i in range(start,end):
            counts[chrom][i] = counts[chrom][i]+1
    
    for chrom in counts:
        a = counts[chrom]
        i = 0
        x = a[i]
        for j in range(len(a)):
            if x != a[j] or len(a)-1 == j:
                if args.subtract_1: x=x-1
                if x > 0:
                    print('%s\t%i\t%i\t%i' % (chrom,i,j,x), file=args.outfile)
                i = j
                x = a[i]
    