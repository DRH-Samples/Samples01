#!/usr/bin/env python

'''
Convert LASTZ output file to BED for self-LASTZ.
LASTZ output file must be in default general format, which is a tab-separated file
with the following columns:
#score  name1   strand1 size1   zstart1 end1    name2   strand2 size2   zstart2 end2    identity        idPct   coverage        covPct
See file:///usr/local/lastz-distrib-1.02.00/README.lastz.html#fmt_general.
'''

import sys

HLABELS = "#score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,idPct,coverage,covPct"

with open(sys.argv[1], 'r') as infile:
    line1 = infile.readline()
    if HLABELS != ','.join(line1.split()):
        exit("This doesn't look like LASTZ default general output; use --format=general.")
    for line in infile:
        t = line.split()
        score = min(int(t[0])/100, 1000)                                        # fit score into range 1-1000
        strand = t[7]                                                           # if strand is '+' it is a direct repeat, if '-' it is an inverted repeat
        bSizes  = '%d,%d' % (int(t[5])-int(t[4]), int(t[10])-int(t[9]))
        bStarts = '0,%d'  % (int(t[9])-int(t[4]))
        label = '%s_(%s)' % (t[11], t[12])
        bed = '%s '*12 % (t[1], t[4], t[10], label, score, strand, 0, 0, 0, 2, bSizes, bStarts)
        print bed
        None
        
None