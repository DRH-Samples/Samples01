#!/usr/bin/env python

'''
Convert LASTZ 'general format' output file to BED (chrom, start, end).
Note: Default coordinates (zstart1, zstart2) are origin-zero half-open, which is
the same as the UCSC browser, so don't convert (see
http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/README.lastz-1.02.00a.html#fmt_general).
'''

from __future__ import print_function
import sys
import argparse

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='LASTZ output file to convert [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='LASTZ output file to convert [default: STDOUT]')
    parser.add_argument('--sequence', choices=['query','target'], default='query',
                        help='which sequence to use, query or target [default: query]')
    parser.add_argument('--min_identity', type=float, default=0,
                        help='minimum percent identity [default: 0]')
    parser.add_argument('--max_identity', type=float, default=100,
                        help='minimum percent identity [default: 100]')
    return parser.parse_args()

if __name__ == '__main__':

    args = parse_cl()
    
    line1 = args.infile.readline()
    HLABELS = "#score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,idPct,coverage,covPct"
    if HLABELS != ','.join(line1.split()):
        exit("This doesn't look like LASTZ general output; use LASTZ option --format=general.")
    
    for line in args.infile:
        t = line.split()
        ident = float(t[12].rstrip('%'))
        if args.min_identity < ident and ident <= args.max_identity:
            # query
            q_seq = t[6]
            q_strand = t[7]
            if '-' == q_strand:
                q_start = int(t[8]) - int(t[10])
                q_end   = int(t[8]) - int(t[9])
            else:
                assert '+' == q_strand
                (q_start, q_end) = t[9], t[10]
            q_name = '%s:%s-%s_%s%%' % (q_seq, q_start, q_end, ident)
                
            # target
            t_seq = t[1]
            t_strand = t[2]
            if '-' == t_strand:
                t_start = int(t[3]) - int(t[5])
                t_end   = int(t[3]) - int(t[4])
            else:
                assert '+' == t_strand
                (t_start, t_end) = t[4], t[5]
            t_name = '%s:%s-%s_%s%%' % (t_seq, t_start, t_end, ident)
            
            if 'query' == args.sequence:
                print('%s '*4 % (q_seq, q_start, q_end, t_name), file=args.outfile)
            else:
                assert 'target' == args.sequence
                print('%s '*4 % (t_seq, t_start, t_end, q_name), file=args.outfile)
        None