#!/usr/bin/env python


'''
Converts BigBed format to BED format. Excludes bases with value less than the given
threshold and combines each retained contiguous segments into a single line.
'''

from __future__ import print_function
import argparse
import sys

def convert(infile, outfile, threshold):
    seqid0 = None
    start0 = None
    end0 = None     # needed if a segment runs to end of a sequence
    bed_id = 1
    for line in infile:
        if line[0] == '#':
            None
            #sys.stderr.write('DEBUG: Discarding: ' + line)
        else:
            tup = line.strip().split('\t')
            if len(tup) != 4:
                sys.stderr.write('ERROR: Invalid line: ' + line)
                exit(1)
            (seqid, start, end, value) = tup
            value = float(value)
            if seqid0 is None:
                seqid0 = seqid
            if seqid0 == seqid:
                if value >= threshold:
                    end0 = end
                    if start0 is None:
                        start0 = start
                else:
                    if start0 is not None:
                        print('\t'.join([seqid0, start0, end0]), file=outfile)
                        start0 = None
            else:
                if start0 is not None:  # print last segment of sequence if needed
                    print('\t'.join([seqid0, start0, end0]), file=outfile)
                    start0 = None
                seqid0 = seqid
                if value >= threshold:
                    start0 = start
                    end0 = end
    if start0 is not None:  # print last segment of last sequence if needed
        print('\t'.join([seqid0, start0, end0]), file=outfile)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType('r'),
                    help='BigBed input file to be converted')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                    help='Bed output file')
    parser.add_argument('threshold', type=float,
                    help='minimum value to be included in output')
    args = parser.parse_args()
    convert(args.infile, args.outfile, args.threshold)
    