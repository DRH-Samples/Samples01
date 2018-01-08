#!/usr/bin/env python

# Make a BED track of the position of all N's in the assembly.

import sys
import argparse
from Bio import SeqIO

desc='Make a BED file with the positions of all Ns in the supplied FASTA file.'
argParser = argparse.ArgumentParser(description=desc)
argParser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),              
                    default=sys.stdin, help='FASTA input file [default: STDIN]')
argParser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),             
                    default=sys.stdout, help='BED output file [default: STDOUT]')
args = argParser.parse_args()

out = args.outfile
for record in SeqIO.parse(args.infile, "fasta"):
    locations = []
    start = -1
    for index, letter in enumerate(record):
        if 'N' == letter and start < 0:
            start = index
        if start >= 0:
            if 'N' != letter:
                locations.append([start, index-1])
                start = -1
    out.write("track name=Ns description=\"Ns\" \n");
    for location in locations:
        out.write("{0}\t{1}\t{2}\n".format(record.id, location[0], location[1]))

