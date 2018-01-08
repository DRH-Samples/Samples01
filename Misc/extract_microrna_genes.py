#!/usr/bin/env python

# Extract microRNA genes from a the set of NCBI GFFs found here:
# /Users/doughoen/Work/Data/Arabidopsis/TAIR9/NCBI/Genes/gff3 .
# Renames chromosomes in accordance with VEGI standards.

#import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

allMirL = []
overlapMirL = []
allAllMirSRs = []
allOverlapMirSRs = []
for n in range(1, 6):
    #if n == 1: iterator = SeqIO.parse(sys.argv[1], 'gb')
    #else: break
    iterator = SeqIO.parse('chr{0}.gb'.format(n), 'gb')
    original = iterator.next()
    try:
        iterator.next()
    except StopIteration:
        pass # do nothing, this is expected
    else:
        raise Error('more than one sequence in ' + gbFile)
    
    seq = Seq('')
    allMirSR = SeqRecord(seq, 'Chr{0}'.format(n))
    allAllMirSRs.append(allMirSR);
    overlapMirSR = SeqRecord(seq, 'Chr{0}'.format(n))
    allOverlapMirSRs.append(overlapMirSR);
    previousEnd = None
    previousIsMir = False
    previousId = None
    for feature in original.features:
        if 'gene' == feature.type:
            if previousIsMir:
                if feature.location.start <= previousEnd:
                    print("5' overlap: " + previousId)
            previousIsMir = False
            previousId = None
            if 'gene' in feature.qualifiers:
                if feature.qualifiers['gene'][0].startswith('MIR'):
                    tag = '{0}={1}'.format(
                                feature.qualifiers['locus_tag'][0],
                                feature.qualifiers['gene'][0])
                    allMirL.append(tag)
                    feature.qualifiers = { feature.qualifiers['locus_tag'][0] :
                                           feature.qualifiers['gene'][0] }
                    allMirSR.features.append(feature)
                    if feature.location.start <= previousEnd:
                        overlapMirSR.features.append(feature)
                        overlapMirL.append(tag)
                    previousIsMir = True
                    previousId = tag
            previousEnd = feature.location.end
            
print 'All:'
for mir in allMirL:
    print mir
print '\nOverlapping:'
for mir in overlapMirL:
    print mir

GFF.write(allAllMirSRs, open('mir_all.gff', 'w'))
GFF.write(allOverlapMirSRs, open('mir_overlap.gff', 'w'))