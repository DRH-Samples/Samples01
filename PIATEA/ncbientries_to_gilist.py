#!/usr/bin/env python

'''
Given a list of entries genereated by a blastdbcmd query, generate a unique
list of GIs that can be used as input into blast using -gilist or -negative_gilist.
The format is simply one GI per line.

E.g.,
$ blastdbcmd -db refseq_protein -entry all | grep -f te_strings.txt >te_entries.txt
$ entries_to_gilist.py te_entries.txt >te_gilist.txt
$ blastx ... -negative_gilist te_gilist.txt

'''

from __future__ import print_function
import sys
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('rU'), nargs='?', 
                        help='output of blastdbcmd, e.g. "blastdbcmd -db refseq_protein -entry all | grep -f te_strings.txt" [default: STDIN]')
args = parser.parse_args()

rexp = re.compile('gi\|(\d+)\|')
gis = set()
for line in args.infile:
    for match in rexp.finditer(line):
        gis.add(match.group(1))
for gi in sorted(gis, key=int):
    print(gi)
