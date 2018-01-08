#!/usr/bin/env python

'''
Convert a TRF dat output file to BED format. Also, can split results based on
the 'copies' found by TRF, in which case 2 BED files are output, one with
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

from __future__ import print_function
#import sys
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('datfile',
                        help='TRF data output file (use -d option when running TRF)')
    parser.add_argument('outfile',
                        help='used to generate output file name(s)')
    parser.add_argument('--copies', type=float, default=0,
                        help='exclude tandem repeats with less than the ' +
                        'specified (float) number of copies')
    parser.add_argument('--allow_no_data', default=False, action="store_true", 
                        help='Do not throw an exception if there is no data ' +
                        '(still throws an exception if there is no header).')
    return parser.parse_args()

def parse(datfile, copies, outfile_lte, outfile_gt, allow_no_data):
    seq_re = re.compile('^Sequence: (\S+)')
    dat_re = re.compile('^' + '\d+ '*3 + '\d+\.\d ' + '\d+ '*8 + '\d+\.\d+ \w+ \w+$')
    seqname = None
    data_mode = False
    data_found = False
    header_found = False
    for line in datfile:
        if data_mode:
            if dat_re.match(line) is None:
                data_mode = False
                seqname = None
            else:
                to_bed(line, seqname, copies, outfile_lte, outfile_gt)
        if not data_mode:
            if seqname is None:
                m = seq_re.match(line)
                if m is not None:
                    seqname = m.group(1)
                    header_found = True
                    continue
            if dat_re.match(line):
                if seqname is None:
                    raise Exception('Could not parse sequence name.')
                to_bed(line, seqname, copies, outfile_lte, outfile_gt)
                data_mode = True
                data_found = True
    if not header_found or not (allow_no_data or data_found):
        raise Exception("Unable to parse file %s. Are you sure this is a TRF DAT file?" % (datfile.name))

def to_bed(line, seqname, copies, outfile_lte, outfile_gt):
    #(start, end, period, n, pattern, match, indel, score) = line.split()
    s = line.split()
    if len(s) < 15:
        raise Exception('Unexpected line: ', line)
    if float(s[3]) <= copies:
        outfile = outfile_lte
    else:
        outfile = outfile_gt
    bed = '%s\t'*3 % (seqname, int(s[0])-1, s[1])    # subtract 1 from start b/c BED is 0-based
    print(bed, file=outfile)

def run(datfilename, outfilename, copies, allow_no_data):
    with open(datfilename, 'r') as datfile:
        if copies < 1:
            with open(outfilename, 'w') as outfile:
                parse(datfile, copies, outfile, outfile, allow_no_data)
        else:
            with open(outfilename+'_excluded', 'w') as outfile_lte:
                with open(outfilename, 'w') as outfile_gt:
                    parse(datfile, args.copies, outfile_lte, outfile_gt, allow_no_data)
    

if __name__ == '__main__':
    args = get_args()
    run(args.datfile, args.outfile, args.copies, args.allow_no_data)
    