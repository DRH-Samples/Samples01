#!/usr/bin/env python

'''
- Filters a multi-fasta file based on the sequence ID
- Filter can be a string (substring or exact match) or a file
- By default, prints sequences that do match the filter string
    - Optionally, does the opposite
'''

import sys
import re
import argparse
from Bio import SeqIO


def parse_cl():
    parser = argparse.ArgumentParser(description='Filter a multi-FASTA file. Either --filterString or --filterFile are required.')
    parser.add_argument('fasta', type=str, 
                            help='multi-FASTA file to filter')
    parser.add_argument('--filterString', type=str, 
                            help='substring (or string) to search for')
    parser.add_argument('--filterFile', type=argparse.FileType('rU'),
                            help='rather than matching a substring, the ID must exactly match the specified string [default: False]')
    parser.add_argument('--exact', action="store_true",
                            help='rather than matching a substring, the ID must exactly match the specified string  [default: False]')
    parser.add_argument('--exclude', action="store_true",
                            help='exclude the matching sequences instead of including them [default: False]')
    args = parser.parse_args()
    if args.filterString is None and args.filterFile is None:
        parser.print_help()
        exit()
    return args

def run(args):
    patterns = set()
    if args.filterFile is not None:
        for line in args.filterFile:
            patterns.add(get_pattern(args, line.rstrip()))
    else:
        patterns.add(get_pattern(args, args.filterString))
    
    filt(args, patterns)
    
def get_pattern(args, s):
    if args.exact:
        return re.compile('^{}$'.format(s))
    else:
        return re.compile('.*{}.*'.format(s))

def filt(args, patterns):
    handle = open(args.fasta, 'rU')
    for record in SeqIO.parse(handle, 'fasta'):
        if args.exclude:
            found = False
            for pat in patterns:
                if pat.match(record.id):
                    found = True
                    break
            if not found:
                SeqIO.write(record, sys.stdout, 'fasta')
        else:
            for pat in patterns:
                if pat.match(record.id):
                    SeqIO.write(record, sys.stdout, 'fasta')
                
    handle.close()


if __name__ == '__main__':
    args = parse_cl()
    run(args)