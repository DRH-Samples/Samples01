#!/usr/bin/env python

'''
Convert the output of HelitronScanner to BED. E.g.,

>scaffold_1
23054:22088 [12:16] 121244:102969 [15:13] 1195407:1194529 [13:11] 2076048:2075181 [13:12] ...

'''

from __future__ import print_function
import sys
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile',
                        type=argparse.FileType('rU'), 
                        help='output of HelitronScanner pairends command')
args = parser.parse_args()

chrRe = re.compile('^>(\S*)\s*')
coordsRe = re.compile('(\d+)\:(\d+) \[(\d+)\:(\d+)\]')
chrom = ''
for line in args.infile:
    line = line.strip()
    if line == '':
            continue
    match = chrRe.match(line)
    if match:
            chrom = match.group(1)
    elif coordsRe.match(line):
        if chrom == '':
            exit('Chrom not set at line:\n' + line)
        for match in coordsRe.finditer(line):
            # Note: I verified that the positions are 1-based using the HelitronScanner 'draw' program.
            left = int(match.group(1))
            right = int(match.group(2))
            score = min(int(match.group(3)), int(match.group(4)))
            if left < right:
                strand = '+'
                start = left-1  # 0-based for BED
                end = right
            else:
                strand = '-'
                start = right-1 # 0-based for BED
                end = left  
            print("\t".join(map(str,(chrom, start, end, 'DNA/Helitron', score, strand ))))
            None
    else:
        exit('Unexpected line:\n' + line)
    
