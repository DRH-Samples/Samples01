#!/usr/bin/env python

# From Chaux et al. 2012, use data in supplemental file "Mob DNA 2012 de la Chaux.xlsx"
# (opened in Excel and then saved as a CSV)
# to relabel the FASTA library in "Mob DNA 2012 de la Chaux.txt" 
# from: >repeatname
# to:   >repeatname#class/subclass

import sys
import argparse
from Bio import SeqIO

argParser = argparse.ArgumentParser()
argParser.add_argument('--fasta', type = argparse.FileType('r'),  help = 'FASTA input file')
argParser.add_argument('--csv',   type = argparse.FileType('rU'), help = 'CSV input file') #rU needed here: Excel f-s up line endings
args = argParser.parse_args()

names = {}
for line in args.csv:
    rname, rsubclass, rclass, etc = line.split(',', 3)
    names[rname] = rname
    if rclass:
        names[rname] = names[rname] + '#' + rclass
        if rsubclass:
            names[rname] = names[rname] + '/' + rsubclass
    
for record in SeqIO.parse(args.fasta, "fasta"):
    record.id = names[record.id]
    record.description = ''
    SeqIO.write(record, sys.stdout, "fasta")
