#!/usr/bin/env python

"""
Input: FASTA file containing a single DNA sequence.
Performs BLASTN of the sequence against itself.
Output: GFF(?) of local non-self alignments.
"""

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

class SelfBlastn:
    def search(self, infile):
        OUTFILE = "self-blastn-out.xml"
        cline = NcbiblastnCommandline(query=infile, subject=infile, outfmt=5, out=OUTFILE)
        stdout, stderr = cline()
        resultHandle = open(OUTFILE)
        blastRecord = NCBIXML.read(resultHandle)
        if len(blastRecord.alignments) != 1:
            raise Exception("Not exactly 1 alignment")
        alignment = blastRecord.alignments[0]
        hsps = iter(alignment.hsps)
        selfHsp = hsps.next()
        if (selfHsp.query_start != 1 or                                         # sanity check
            selfHsp.sbjct_start != 1 or
            selfHsp.query_end != blastRecord.query_length or
            selfHsp.sbjct_end != blastRecord.query_length or
            selfHsp.positives != blastRecord.query_length or
            selfHsp.strand != (None, None) or
            selfHsp.gaps != 0):
                raise Exception("Self-hit exception.")
        selfAlignments = []
        for hspA in hsps:
            hspB = hsps.next()      # each alignment should result in 2 symmetric HSPs (b/c it's self-BLAST)
    # TODO: deal with strands!!
            if (hspA.strand != (None, None) or hspB.strand != (None, None)):    # sanity check
                    raise Exception("Unexpected strand")
            if (hspA.query_start != hspB.sbjct_start or                         # sanity check
                hspA.query_end != hspB.sbjct_end or
                hspA.score != hspB.score):
                    raise Exception("Unexpected HSP pair")
            termini = ((hspA.sbjct_start, hspA.sbjct_end), (hspA.query_start, hspA.query_end))
            print termini
            selfAlignments.append(termini)
        return selfAlignments


    
# Main
if __name__ == "__main__":
    import sys      # for command-line arguments
    if len(sys.argv) != 2:
        print "Usage: <program> fasta_file_name"
        exit(0)
    infile = sys.argv[1]
    SelfBlastn().search(infile)