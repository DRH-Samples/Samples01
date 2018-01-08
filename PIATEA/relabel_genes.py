#!/usr/bin/env python

'''
Input:
(1) *Sorted* GTF file listing the exons of a gene and "gene_id".
    Note: If necessary, first sort the GTF file using:
    'sort -k1,1 -k4,4n -k5,5n fgh_scaffolds_1-8.gtf <file>'.

Output: BED file labelled for use in PIATEA as: intron, exon, or intergenic.
'''

from __future__ import print_function
import argparse
import sys

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genes_file', type=argparse.FileType('r'),
                        help='Sorted GTF input file with exon locations')
    parser.add_argument('--chrom_file', type=argparse.FileType('r'),
                        help='Tab-separated value file with chromosome lengths')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='Bed output file [default: STDOUT]')
    return parser.parse_args()

def read_chrom_lengths(chromFile):
    chromLengths = {}
    for line in chromFile:
        chrom, length = line.split()
        chromLengths[chrom] = length
    return chromLengths

def process(genesFile, chromLengths, outfile):
    completed = {}
    prev = None
    for line in genesFile:
        if '#' == line[0] or 'track' == line[0:5]: continue                     # ignore comment lines
        d = line.split()
        if d[2] != 'exon':
            continue
        chrom = d[0]
        if chrom in completed:
            exit("Is the genes file sorted? Non-contiguous chromosome at: " + line)
        start = int(d[3])-1                                                          # -1 b/c BED is 0-indexed
        if prev != None and chrom == prev[1] and start < prev[3]:
            exit("Is the genes file sorted? Start before end at: " + line)
        end = int(d[4])
        strand = d[6]
        geneId = d[9]
        if prev == None:
            print('\t'.join(map(str, (chrom, 0, start, 'intergenic', 0, '.'))), file=outfile)
        elif (prev[0] == geneId):
            print('\t'.join(map(str, (chrom, prev[3], start, 'intron', 0, strand))), file=outfile)
        elif (prev[1] == chrom):
            print('\t'.join(map(str, (chrom, prev[3], start, 'intergenic', 0, '.'))), file=outfile)
        else:
            print('\t'.join(map(str, (prev[1], prev[3], chromLengths[prev[1]], 'intergenic', 0, '.'))), file=outfile)
            print('\t'.join(map(str, (chrom, 0, start, 'intergenic', 0, '.'))), file=outfile)
            completed[prev[1]] = 1
        print('\t'.join(map(str, (chrom, start, end, 'exon', 0, strand))), file=outfile)
        prev = (geneId, chrom, start, end)
    print('\t'.join(map(str, (prev[1], prev[3], chromLengths[chrom], 'intergenic', 0, '.'))), file=outfile)


if __name__ == '__main__':
    args = parse_cl()
    chromLengths = read_chrom_lengths(args.chrom_file)
    process(args.genes_file, chromLengths, args.outfile)
    
    