#!/usr/bin/env python

'''
Given a classified RepeatModeler library ("consensi.fa.classified") and
a library of non-TE proteins (e.g., RefSeq protein), generates a new file
in which unclassified ("Unknown") TEs that align to a RefSeq protein are
filtered out.
'''

import sys
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=argparse.FileType('rU'),
                        help='consensi.fa.classified file to filter')
parser.add_argument('db', type=str,
                        help='name of preformatted blastx database of non-TE proteins')
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='output file [default: STDOUT]')
args = parser.parse_args()