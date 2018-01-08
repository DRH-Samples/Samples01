#!/usr/bin/env python

'''
Converts LTR_FINDER output (tabular format, option "-w2") to BED format.
'''

import argparse
import re
import sys

# E.g.)
# index   SeqID         Location        LTR len Inserted element len    TSR     PBS     PPT             RT      IN (core)   IN (c-term)     RH      Strand  Score   Sharpness       Similarity
# [ 1]    scaffold_2    31789-43543     975,975 11755                   ATACT   N-N     42489-42503     N-N     N-N         N-N             N-N     +       6       0.529,0.571     0.997
global pattern
pattern = "^\[\s*(\d+)\]\s+(\S+)\s+(\d+)-(\d+)\s+(\d+),(\d+)\s+\d+\s+(\w+)" + \
          "\s+\S+"*6 + "\s+(\S+)\s+(\S+)\s+\S+\s+\S+$";

# Note: BED format is 0-indexed, includes the start position, excludes the end position.
def print_bed(match):
    te_id = "LTR_TE|" + match.group(1)
    seq_id = match.group(2)
    te_start        = int(match.group(3)) - 1
    te_end          = int(match.group(4)) - 1
    left_ltr_len    = int(match.group(5))
    right_ltr_len   = int(match.group(6))
    tsd             = match.group(7)
    strand          = match.group(8)
    score           = match.group(9)
    
    left_ltr_id = '|'.join(['LTR|left', te_id])
    left_ltr_end = te_start + left_ltr_len
    print '%s\t'*6 % (seq_id, te_start, left_ltr_end, left_ltr_id, score, strand)
    
    right_ltr_id = '|'.join(['LTR|right', te_id])
    right_ltr_start = te_end - right_ltr_len + 1
    print '%s\t'*6 % (seq_id, right_ltr_start, te_end+1, right_ltr_id, score, strand)
    
    inside_id = '|'.join(['inside|-', te_id])
    print '%s\t'*6 % (seq_id, left_ltr_end, right_ltr_start, inside_id, score, strand)
    
    if tsd != 'N':      # 'N' indicates TSD not found
        left_tsd_id = '|'.join(['TSD|left', te_id])
        print '%s\t'*6 % (seq_id, te_start-len(tsd), te_start, left_tsd_id, score, strand)
        right_tsd_id = '|'.join(['TSD|right', te_id])
        print '%s\t'*6 % (seq_id, te_end+1, te_end+len(tsd)+1, right_tsd_id, score, strand)
        

if __name__ == '__main__':
    cl_parser = argparse.ArgumentParser()
    cl_parser.add_argument('ltr_file', type=argparse.FileType('r'),
                    help='LTR_FINDER tabular output file to be converted (use -w2 option when running LTR_FINDER)')
    cl_args = cl_parser.parse_args()
    infile = cl_args.ltr_file
    expr = re.compile(pattern)
    for line in infile:
        match = expr.match(line)
        if match:
            print_bed(match)
        else:
            #None
            sys.stderr.write("Discarding: " + line)

