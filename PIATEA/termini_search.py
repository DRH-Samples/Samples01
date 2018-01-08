#!/usr/bin/env python

"""
Run a sliding window self-LASTZ scan to identify potential terminal TE repeats, e.g. LTRs and TIRs.
Input:  A FASTA sequence file to be scanned, e.g. a chromosome assembly.
Output: A BED file, written to STDOUT.

Revisions:
1. Originally specified the windows using the subrange action built into lastz, "[start..end]". However,
    it turns out that this adds ~5 sec. overhead to each call to lastz, increasing the running time by
    about an order of magnitude. (Not sure why as it's not clearly stated in the docs what steps are skipped
    for this action.) So instead, revised to manually chunk the input sequence into files and feed them to
    LASTZ, then adjust the coordinates.

"""

from __future__ import print_function
import sys
import tempfile
import re
import time
import argparse
import os
from subprocess import call
from Bio import SeqIO

#START_POSITION = 9400001           # reverse strand error (chr1)
#START_POSITION = 1450001           # reverse strand error (scaffold_1)

hit_id_plus = 1                     # TO DO: make not global
hit_id_minus = 1                    # TO DO: make not global

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infilename',
                        help='FASTA file to search')
    parser.add_argument('--outfilename',
                        help='BED file for output')
    parser.add_argument('--max_te_length', type=int, default=50000,
                        help='maximum allowed TE length (bp) [default 50000]')
    parser.add_argument('--window_length', type=int, default=300000,
                        help='window length (bp); should be at least 2x max_te_length [default 300000]')
    parser.add_argument('--format_type', choices=['simple', 'blocks'], default='simple',
                        help='output each terminus as a separate BED line ("simple") or as "blocks" [default "simple"]')
    parser.add_argument('--gap', default='400,30',
                        help='LASTZ gap penalties [default "400,30"]')
    parser.add_argument('--ydrop',
                        help='LASTZ y-drop threshold [default: not set => open+300*extend = 9400 for default gap settings]')   # see http://galaxylocal.genenetwork.org:8080/tool_runner?tool_id=lastz_paired_reads_wrapper
    parser.add_argument('--start', type=int, default=0,
                        help='start at this position in the sequence (bp); zero-based indexing [default 0]')
    parser.add_argument('--alignments', action='store_true',
                        help='output alignment files (requires an extra scan per window)')
    parser.add_argument('--no_bed', action='store_true',
                        help='do not output bed files (requires an extra scan per window)')
    args = parser.parse_args()
    if not (args.alignments) and args.no_bed:
        print("Cannot have --alignments=False and --no_bed=True because nothing will be output!", file=sys.stderr)
        exit(1)
    cfg = Config()
    cfg.infilename = args.infilename
    cfg.outfilename = args.outfilename
    cfg.max_te_length = args.max_te_length
    cfg.window_length = args.window_length
    cfg.format_type = args.format_type
    cfg.gap = args.gap
    cfg.ydrop = args.ydrop
    cfg.start = args.start
    cfg.alignments = args.alignments
    cfg.no_bed = args.no_bed
    return cfg

class Config():
    def __init__(self):
        self.infilename = None
        self.outfilename = None
        self.max_te_length = 50000
        self.window_length = 300000
        self.format_type = 'simple'
        self.gap = '400,30'
        self.ydrop = None
        self.start = 0
        self.alignments = False
        self.no_bed = False

def read_seq(infile):
    result = list(SeqIO.parse(infile, "fasta"))
    if len(result) == 0:
        raise Exception("No sequence found in file %s. Are you sure this is a FASTA file?" % (infile.name))
    if len(result) > 1:
        raise Exception("More than one sequence found in file %s." % (infile.name))
    return result[0]

def write_subseq(outfile, seq, start, end):
    subseq = seq[start:end+1]                                                                       # adjust end b/c biopython appears to use non-inclusive end (right?)
    SeqIO.write(subseq, outfile, "fasta")
    outfile.flush()                                                                                 # needed!

def calculate_windows(seq_length, cfg):
    windows = []
    for start in range(cfg.start, seq_length, cfg.window_length - cfg.max_te_length):
        windows.append((start, min(start + cfg.window_length-1, seq_length)))
    return windows

def to_bed(outfile, lastzfile, cut_position, window, cfg):
    """Convert LASTZ output (default general format) to BED."""
    global hit_id_plus
    global hit_id_minus
    line1 = lastzfile.readline()
    if line1 != ("#score\tname1\tstrand1\tsize1\tzstart1\tend1\tname2\tstrand2"+
                 "\tsize2\tzstart2\tend2\tidentity\tidPct\tcoverage\tcovPct\n"):
        raise Exception("Unexpected first line of LASTZ output file:\n"+line1)
    empty = True
    for line in lastzfile:
        empty = False
        (score, name1, strand1, size1, zstart1, end1, name2, strand2, size2,
         zstart2, end2, identity, idPct, coverage, covPct) = line.split()
        zstart1 = int(zstart1)
        zstart2 = int(zstart2)
        size2 = int(size2)
        end1 = int(end1)
        end2 = int(end2)
        zstart1 = zstart1 + window[0]
        end1 = end1 + window[0]
        if (strand2 == '+'):
            zstart2 = zstart2 + window[0]
            end2 = end2 + window[0]
        elif (strand2 == '-'):
            right2 = window[1] - zstart2 + 1
            left2 = window[1] - end2 + 1
            zstart2 = left2
            end2 = right2
        else:
            raise Exception("Unexpected strand: '%s' " % strand2)
        real_end =   max(zstart1, zstart2, end1, end2)
        real_start = min(zstart1, zstart2, end1, end2)
        if (real_end > cut_position) and (real_end - real_start <= cfg.max_te_length):
            if (strand2 == '+'):
                hit_id = "%s%s" % (hit_id_plus, '+')
                hit_id_plus = hit_id_plus + 1
            else:
                hit_id = "%s%s" % (hit_id_minus, '-')
                hit_id_minus = hit_id_minus + 1
            if end2 != real_end:
                print("\nDEBUG: end2 not max: ID %s %d-%d, %d-%d\n\t\t\t\t\t"
                                 % (hit_id, zstart1, end1, zstart2, end2), file=sys.stderr)
            score = min(int(score)/100, 1000)                                                           # cram LASTZ score into range 1-1000
            if cfg.format_type == 'blocks':
                bSizes  = '%d,%d' % (end1-zstart1, end2-zstart2)
                bStarts = '0,%d'  % (zstart2-zstart1)
                label = '%s%s(%s)' % (identity, strand2, idPct)
                print('%s\t'*12 % (name1, zstart1, end2, label, score, strand2, 0, 0, 0, 2, bSizes, bStarts), file=outfile)
            else: # cfg.format == 'simple'
                print('%s\t'*6 % (name1, zstart1, end1, hit_id, score, strand1), file=outfile)
                print('%s\t'*6 % (name1, zstart2, end2, hit_id, score, strand2), file=outfile)
    if empty:
        print('WARNING!! No LASTZ results for %s in window [%d..%d].'
                         % (cfg.infilename, window[0]+1, window[1]+1), file=sys.stderr)                              # convert to 1-based indexing: easier to read

def run(outfilename, cfg):
    bedfile = None
    alnfile = None
    if not cfg.no_bed:
        bedfile = open(outfilename+'.bed', 'w')
    if cfg.alignments:
        alnfile = open(outfilename+'.aln', 'w')
    overlap_position = cfg.start-1
    whole_seq = read_seq(cfg.infilename)                                                        # read full sequence
    for window in calculate_windows(len(whole_seq.seq), cfg):
        print('DEBUG: Scanning window [%d..%d]...' % (window[0]+1, window[1]+1), file=sys.stderr)     # convert to 1-based indexing: easier to read
        with tempfile.NamedTemporaryFile() as subseq_file:                                      # automatic deletion on 'with' context exit
            write_subseq(subseq_file, whole_seq, window[0], window[1])
            time0 = time.time()
            scan(subseq_file.name, bedfile, alnfile, overlap_position, window, cfg)
            print('\tCompleted in %0.2f seconds.\n' % (time.time()-time0), file=sys.stderr)
        overlap_position = window[1]                                                            # do not use hits that terminate before the end of the previous scan window
    if not cfg.no_bed:
        bedfile.close()
    if cfg.alignments:
        alnfile.close()

def scan(seqfilename, bedfile, alnfile, overlap_position, window, cfg):
    common = ['lastz', seqfilename, '--self', '--nomirror', '--ambiguous=iupac', '--gap='+cfg.gap]      # or, --ambiguous=n
    if cfg.ydrop is not None:
        common.append('--ydrop='+cfg.ydrop)
    if not cfg.no_bed:
        with tempfile.NamedTemporaryFile() as lastz_general:  
            call(common + ['--format=general', '--output='+lastz_general.name])
            to_bed(bedfile, lastz_general, overlap_position, window, cfg)
    if cfg.alignments:
        with open(os.devnull, 'w') as FNULL:
            with tempfile.NamedTemporaryFile() as lastz_text:  
                call(common + ['--format=text', '--output='+lastz_text.name], stdout=FNULL)     # this format prints tons of blank lines; redirect them to /dev/null
                alnfile.write(lastz_text.read())                                                        # concat new results to alignment file
                alnfile.flush()

'''
TO DO:
Filter results. Look for termini vs. overlapping vs. potential palindromes, etc.
Put results into different files, starting with 1 file for clear potential termini,
one file for everything else (for now).
'''

if __name__ == '__main__':
    cfg = parse_cl()
    if cfg.outfilename is not None:
        outfilename = cfg.outfilename
    else:
        outfilename = cfg.infilename + '.self-lastz'
    run(outfilename, cfg)



