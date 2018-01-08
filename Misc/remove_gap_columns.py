#!/usr/bin/env python

'''
Removes columns from an amino acid multiple alignment that have gaps in at least the specified percentage of sequences.
'''

from __future__ import print_function
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import single_letter_alphabet

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType('r'),
                        help='Original FASTA alignment file')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help='Edited FASTA alignment file')
    parser.add_argument('--pct', type=int, default=50,
                        help='Remove columns with gaps in at least this percentage of sequences [default: 50]. \
                             Note: Gaps must be represented by "-".')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_cl()
    align = AlignIO.read(args.infile, 'fasta')
    #print(align)
    nSeqs = len(align)
    filtered = None
    for i in range(align.get_alignment_length()):
        nGaps = (align[:, i]).count('-')
        pctGaps = 100*nGaps/nSeqs
        if pctGaps < args.pct:
            if filtered == None:
                filtered = align[:, i:i+1]
            else:
                filtered = filtered + align[:, i:i+1]
    #print(filtered)
    AlignIO.write(filtered, args.outfile, 'fasta')
    
    
    