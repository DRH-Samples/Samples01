#!/usr/bin/env python

'''
Attempt to measure consistency of input BedGraph data in 2 ways:
1. Compare the mean over segments to the left and the right of each location.
2. Compare the variation of a segment, e.g. using standard or median absolute variation.
'''

from __future__ import print_function
import sys
import array
import argparse
import scipy.stats
import math

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
    #START = False
    for chrom in values:
        results[chrom] = array.array('f')
        results[chrom].extend([0]*(len(values[chrom])))
        for i in range(len(values[chrom])):
            # t-test
            #if i-LENGTH>=0 and i+LENGTH<=len(values[chrom]):
            #    a = values[chrom][i-LENGTH:i]
            #    b = values[chrom][i:i+LENGTH]
            #    if math.fsum(a) + math.fsum(b) == 0:
            #        t_result = -99
            #        v_result = -99
            #    else:
            #        START = True
            #        if max(a) == min(a) and max(b) == min(b):
            #            if a[0] == b[0]:
            #                t_result = -1
            #                v_result = -1
            #            else:
            #                t_result = 99
            #                v_result = 99
            #        else:
            #            t_result = math.floor(-math.log10((scipy.stats.ttest_ind(a,b))[1]))
            #    if START:
            #        print('%s\t%i\t%i\t%i' % (chrom,i,i+1,t_result))
            # variation
            if i-LENGTH/2>0 and i-LENGTH/2<len(values[chrom]):
                c = values[chrom][i-LENGTH/2:i+LENGTH/2]
                if math.fsum(c) == 0:
                    v = -1
                else:
                    #START = True
                    v = scipy.stats.variation(c)
                results[chrom][i] = v
                #if START:
                    #print('%s\t%i\t%i\t%-.2f' % (chrom,i,i+1,v_result))
    return results

def write_results(results, outfile):
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
    write_results(results, args.outfile)
    #raw_input("Press any key to exit.")
