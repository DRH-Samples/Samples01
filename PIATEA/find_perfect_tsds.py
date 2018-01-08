#!/usr/bin/env python

'''
Looks for perfect TSDs.
Input: BED files output from termini_search.py along with the corresponding 
Output: A new BED file with all perfect non-overlapping TSDs within the specified region
'''

from __future__ import print_function
import argparse
#import os.path
from Bio import SeqIO

class Termini():
    def __init__(self):
        self.start = None
        self.end = None
        self.id = None

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('sequence_names', nargs='+',
                        help='names of sequences to search; for each name, there must be ' +
                        'a corresponding FASTA file and a [name].termini_search.bed file')
    parser.add_argument('--sequence_dir', default='.',
                        help='directory containing sequence files [default: cwd]')
    parser.add_argument('--bed_dir', default='.',
                        help='directory containing bed files [default: cwd]')
    parser.add_argument('--inner_margin', type=int, default='5',
                        help='# bp inside termini to search [default: 10]')
    parser.add_argument('--outer_margin', type=int, default='20',
                        help='# bp inside termini to search [default: 10]')
    parser.add_argument('--max_tsd_size', type=int, default='20',
                        help='# bp inside termini to search [default: 10]')
    return parser.parse_args()

def find_tsds(seqfilepath, bedfilepath):
    seq = readseq(seqfilepath)
    for termini in parsepairs(bedfilepath):
        find_tsds(seq, termini)


def read_seq(infile):
    result = list(SeqIO.parse(infile, "fasta"))
    if len(result) == 0:
        raise Exception("No sequence found in file %s. Are you sure this is a FASTA file?" % (infile.name))
    if len(result) > 1:
        raise Exception("More than one sequence found in file %s." % (infile.name))
    return result[0]

if __name__ == '__main__':
    args = parse_cl()
    for seqname in args.sequence_names:
        seqfilepath = os.path.join(args.sequence_dir, seqname)
        if not os.path.exists(seqfilepath):
            print('Sequence file does not exist: ' + seqfilepath, sys.stderr)
            exit(1)
        bedfilepath = os.path.join(args.bed_dir, seqname + '.termini_search.bed')
        if not os.path.exists(bedfilepath):
            print('BED file does not exist: ' + bedfilepath, sys.stderr)
            exit(1)
        find_tsds(seqfilepath, bedfilepath)
