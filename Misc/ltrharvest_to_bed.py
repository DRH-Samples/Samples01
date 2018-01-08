#!/usr/bin/env python

'''
Convert an LTRharvest GFF to a BED file with PIATEA labels: TSD, LTR, inside.
(Actually, I got the GFF from Glenn; not sure if it's the direct output of LTRharvest...)
'''


from __future__ import print_function
import argparse
import sys


def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='GFF file from LTRharvest [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='Bed output file [default: STDOUT]')
    return parser.parse_args()


def process(infile, outfile):
    state = 'outside'
    idd = 1 
    insideStart = None
    for line in infile:
        if '#' == line[0] or 'track' == line[0:5]: continue                     # ignore comment lines
        d = line.split()
        start = int(d[3]) - 1
        if d[2] == 'target_site_duplication':
            if state == 'outside':  
                print('\t'.join((d[0], str(start), d[4], 'TSD|left|LTR_TE|'+str(idd), '0', '.' )), file=outfile)
                state = 'left_tsd'
            elif state == 'right_ltr':  
                print('\t'.join((d[0], str(start), d[4], 'TSD|right|LTR_TE|'+str(idd), '0', '.' )), file=outfile)
                idd = idd+1
                state = 'outside'
            else:
                exit('Error: Unexpected TSD from state ' + state)
        if d[2] == 'long_terminal_repeat':
            if state == 'left_tsd':    
                print('\t'.join((d[0], str(start), d[4], 'LTR|left|LTR_TE|'+str(idd), '0', '.' )), file=outfile)
                insideStart = d[4]
                state = 'inside'
            elif state == 'inside':
                print('\t'.join((d[0], insideStart, str(start), 'inside|-|LTR_TE|'+str(idd), '0', '.' )), file=outfile)
                print('\t'.join((d[0], str(start), d[4], 'LTR|right|LTR_TE|'+str(idd), '0', '.' )), file=outfile)
                insideStart = None
                state = 'right_ltr'
            else:
                exit('Error: Unexpected left LTR from state ' + state)        


if __name__ == '__main__':
    args = parse_cl()
    process(args.infile, args.outfile)