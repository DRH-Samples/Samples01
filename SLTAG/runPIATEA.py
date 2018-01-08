#!/usr/bin/env python

import sys      # for command-line arguments
import time     # for performance metrics
from Bio import SeqIO
from selfblastn01 import SelfBlastn
from stag14 import Grammar

if __name__ != "__main__":
    exit("ERROR: Module run-PIATEA is only used to run PIATEA from the command-line.")
    
else:    
    start = time.clock()
    if len(sys.argv) != 2:
        exit("Usage: run-PIATEA-<xx> infile.fa")
        
    infile = sys.argv[1]
    selfBlast = SelfBlastn(infile)
    dna = str(SeqIO.parse(infile, "fasta").next().seq)
    grammar = Grammar(selfBlast, dna)
    grammar.parse()
    
    end = time.clock()
    print; print "Elapsed time: %f s" % (end - start)
    #raw_input("Press <return> to exit: ")
    print