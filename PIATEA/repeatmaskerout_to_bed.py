#!/usr/bin/env python

'''
Converts a RepeatMasker ".out" file to BED format.
'''

from __future__ import print_function
import sys
import argparse

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('rmfile', type=argparse.FileType('r'),
                        help='RepeatMasker ".out" file')
    parser.add_argument('--minscore', type=int, default='0',
                        help='minimum RepeatMasker "SW score" for result to be included in BED file [default: 0]')
    parser.add_argument('--maxdivergence', type=int, default='100',
                        help='minimum RepeatMasker "SW score" for result to be included in BED file [default: 100]')
    parser.add_argument('--maxtruncation', type=int, default='100',
                        help='minimum RepeatMasker "SW score" for result to be included in BED file [default: 100]')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_cl()
    for lineno in range(3):
        print('INFO: Discarding header line #%s: %s' % (lineno, args.rmfile.readline().strip()), file=sys.stderr)
    print('# The "name" field has the following format: "class/family|repeat|%divergence_%deletions_%insertions|repstart-end_left_%truncation"')
    discarded = 0
    for line in args.rmfile:
        f = line.split()
        score0 = int(f[0])
        score = min(score0/10, 1000)
        strand = f[8]
        if strand == 'C':
            strand = '-'
            (repleft, repend, repstart) = f[11:14]
        else:
            (repstart, repend, repleft) = f[11:14]
        repleft = repleft[1:-1]
        (repleft, repend, repstart) = (int(repleft), int(repend), int(repstart)) 
        truncation = 100*repleft/(repleft+repend-repstart+1)
        if score0 < args.minscore or truncation > args.maxtruncation or float(f[3]) > args.maxdivergence:
            discarded = discarded+1
            continue    
        desc = 'class/family|repeat|%divergence_%deletions_%insertions|repeat:start-end_left_%coverage'
        name = '%s|%s|%s_%s_%s|%s-%s_%s_%s' % (f[10], f[9], f[1], f[2], f[3], repstart, repend, repleft, truncation)
        print('\t'.join([f[4], str(int(f[5])-1), f[6], name, str(score), strand]))
    print('Discarded %s records. Done.' % (discarded), file=sys.stderr)
