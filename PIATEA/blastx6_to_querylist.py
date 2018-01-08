#!/usr/bin/env python

'''
Given the output of BLASTX (-outfmt 6), outputs a unique list of query IDs 

'''

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('rU'), 
                        help='BLASTX -outfmt 6 output')
args = parser.parse_args()

queries = set()
for line in args.infile:
    queries.add(line.split()[0])
for q in sorted(queries):
    print(q)