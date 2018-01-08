#!/usr/bin/env python

'''
Convert a TRF dat output file to BED format. Also, can split results based on
specifications, e.g. copies-#, in which case 2 BED files are output, one with
results that match the specifications, one with results that do not match.


Example head of TRF dat file:
================================================================================

Tandem Repeats Finder Program written by:

Gary Benson
Program in Bioinformatics
Boston University
Version 4.07b


Sequence: Chr4 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02

Mirror:Test doughoen$ head -30 chr4.fas.2.7.7.80.10.50.500.dat
Tandem Repeats Finder Program written by:

Gary Benson
Program in Bioinformatics
Boston University
Version 4.07b


Sequence: Chr4 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02



Parameters: 2 7 7 80 10 50 500


6319 6395 39 2.0 39 84 5 102 41 15 12 29 1.85 TAAACTGATACATGCCGAATTGCACATGATAAATATAAC TAAACTGATACATGCACAATTGCACATGATAAATATAACTAAATTGATACATGCCGATTTGCGCATGATAATTATAA
38405 38502 26 3.8 26 98 0 187 45 8 30 15 1.75 AATGAACTGAGAGACAAGAGATTAGG AATGAACTGAGAGACAAGAGATTAGGAATGAACTGAGAGACAAGAGATTAGGAATGAACTGAGAGACAAGAGATTTGGAATGAACTGAGAGACAAGAG
...

================================================================================


Fields per line (above) are:
start end period_size copies pattern_size match_% indel_% score

'''

import sys
import re

def main():
    with open(sys.argv[1], 'rU') as datfile:
        parse(datfile)
    
def parse(datfile):
    seq_re = re.compile('^Sequence: (\w+) .*')
    dat_re = re.compile('^' + '\d+ '*3 + '\d+\.\d ' + '\d+ '*8 + '\d+\.\d+ \w+ \w+$')
    seqname = None
    data_mode = False
    for line in datfile:
        if data_mode:
            if dat_re.match(line) is None:
                None    #debug
            to_bed(line, seqname)
        else:
            if seqname is None:
                m = seq_re.match(line)
                if m is not None:
                    seqname = m.group(1)
                    continue
            if dat_re.match(line):
                if seqname is None:
                    raise Exception('Could not parse sequence name.')
                to_bed(line, seqname)
                data_mode = True

def to_bed(line, seqname):
    #(start, end, period, n, pattern, match, indel, score) = line.split()
    s = line.split()
    bed = '%s '*3 % (seqname, s[0], s[1])
    print bed
    

# DO IT!
main()
    