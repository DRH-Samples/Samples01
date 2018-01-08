#!/usr/bin/env python

''' Certain tracks in the VEGI browser, such as ChainD (orthology chain depth) and sncRNAS (sRNA),
have a strange stutter where every other base has height 0. This script fills in these single base
gaps with the value from the previous base. '''

from __future__ import print_function
import sys
import array
import argparse


def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='BedGraph input file [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='BedGraph output file [default: STDOUT]')
    return parser.parse_args()

    
def read_values(infile):
    values = {}
    for line in infile:
        if '#' == line[0] or 'track' == line[0:5]:              # ignore comment lines
            print('WARNING: Ignoring line: '+line, sys.stderr)
            continue
        (chrom, start, end, val) = (line.split())[0:4]
        start = int(start); end = int(end); val = float(val)
        if chrom not in values:
            values[chrom] = array.array('f')  # float
        if len(values[chrom]) < end:
            values[chrom].extend([0]*(end-len(values[chrom])))
        for i in range(start,end):
            values[chrom][i] = values[chrom][i]+val
    return values

def scan(values):
    LENGTH = 100
    results = {}
    for chrom in values:
        results[chrom] = array.array('f')
        results[chrom].extend([0]*(len(values[chrom])))
        for i in range(1,len(values[chrom])):
            if values[chrom][i] == 0:
                results[chrom][i] = values[chrom][i-1]
            else:
                results[chrom][i] = values[chrom][i]
    return results


def compress_and_write(results, outfile):
    for chrom in results:
        a = results[chrom]
        i = 0
        x = a[i]
        for j in range(len(a)):
            if x != a[j] or len(a)-1 == j:
                #if args.subtract_1: x=x-1
                if x > 0:
                    print('%s\t%i\t%i\t%.1f' % (chrom,i,j,x), file=outfile)
                i = j
                x = a[i]





if __name__ == '__main__':
    args = parse_cl()
    values = read_values(args.infile)
    results = scan(values)
    compress_and_write(results, args.outfile)