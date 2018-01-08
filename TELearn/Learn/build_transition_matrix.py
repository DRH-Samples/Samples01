#!/usr/bin/env python

'''
Input: BED file with N states (BED names)
Output: N*N transisiton probability matrix (from*to)
'''

import sys

f = open(sys.argv[1])

background = "background"

counts = {} # state -> state -> count
counts[background] = {}
counts[background][background] = 0

prevState = background  # not always correct, but vanishingly small error and simpler
prevScaffold = None
prevEnd = 0
for line in f:
    scaffold, start, end, name = line.split()[0:4]
    state = name.split('|')[0]
    state = state.split('/')[0]
    start = int(start) ; end = int(end)
    length = int(end)-int(start)
    
    # if state previously unseen, add to counts matrix
    if state not in counts:
        a={}
        a[state] = 0
        for state0 in counts:
            counts[state0][state] = 0
            a[state0] = 0
        counts[state] = a
    
    if scaffold == prevScaffold:
        if start < prevEnd:
            exit("Not sorted by start")
        if start > prevEnd: # count gap
            counts[prevState][background] += 1
            counts[background][background] += start - prevEnd
            counts[background][state] += 1
        else: # contiguous 
            counts[prevState][state] += 1
    else:   # new scaffold
        # refactor duplicate code
        if start > 0: # count gap
            counts[background][background] += start
            counts[background][state] += 1
    counts[state][state] += length
        
    prevState = state
    prevScaffold = scaffold
    
for k in counts.keys():
    print k, counts[k]
