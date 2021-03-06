#!/usr/bin/env python

''' Filter sliding-window LASTZ results and remove any overlapping hits.
    Input is the output generated by termini_search.py, namely a BED file with
    pairs of lines representing (subject, query) hit pairs. First, pairs are 
    filtered based on several criteria (see arguments help for details). Then,
    ovelaps between hits that pass all criteria are resolved, as follows:
    1. Clusters of pairs with any overlap between any hit (single-link) are formed.
    2. In each cluster, the best pair in the cluster is determined, as follows:
        2.1. The highest scoring pair and any other pair with score not less than
            --score_buffer_pct % lower score are determined.
        2.2. If there is only one such pair, it is the best pair. Otherwise,
            the subset of these pairs with total hit length (sum of lengths of the 2 hits)
            not less than --hit_length_buffer_pct % of the longest total hit length
            are determined.
        2.3. If there is only one such pair, it is the best pair. Otherwise,
            the longest pair (total hit length plus separation) is the best pair.
    3. The best pair is kept, and any pair overlapping either of its hits is
        truncated so that it no longer overlaps. If the truncation of any hit is more than
        --max_truncation, the corresponding pair is discarded.
    4. Steps #1 to #3 are repeated until there are no pairs remaining in the cluster.
    
    Overlap resolution is required for HMM input. However, the filter is also capable of
    identifying direct overlapping repeats and palindromes (inverted overlapping repeats),
    which may be useful in identifying certain types of elements (e.g., Helitrons). For these
    types, overlaps are not currently resolved.
'''

from __future__ import print_function
import sys
import argparse
import collections
import math

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='BedGraph input file [default: STDIN]')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help='BedGraph output file [default: STDOUT]')
    parser.add_argument('--type', choices=['LTR', 'TIR', 'direct_overlapping', 'palindrome'], default='LTR',
                        help='type of local similarity to look for [default: LTR]')
    
    # length
    parser.add_argument('--min_pair_length', type=int, default=0,
                        help='minimum allowed length from start of first hit to end of second hit (bp) [default: 0]')
    parser.add_argument('--max_pair_length', type=int, default=20000,
                        help='maximum allowed length from start of first hit to end of second hit (bp) [default: 20000]')
    parser.add_argument('--min_hit_separation', type=int, default=50,
                        help='minimum amount of separation between end of first and start of second hit (bp); ' +
                        ' only applies to LTR and TIR [default: 50]')
    parser.add_argument('--max_hit_length', type=int, default=2500,
                        help='maximum length of any hit (bp) [default: 2500]')
    
    # score
    parser.add_argument('--score_filter', choices=['filter_1', 'filter_2'], default='filter_2',
                        help='score filter funtion. [default: filter_2]' +
                        'filter_1: score > min_score ' +
                        'filter_2: score > accept_score or (score>min_score and score*score_multiplier>te_length)')
    parser.add_argument('--min_score', type=int, default=70,
                        help='(for score filter_1 and filter_2) no hits below this score are used [default: 70]')
    parser.add_argument('--accept_score', type=int, default=200,
                        help='(for score filter_2) no hits above this score are rejected based on score [default: 200]')
    parser.add_argument('--score_multiplier', type=int, default=50,
                        help='(for score filter_2) pairs between min_score and accept score are rejected or accepted iff: ' +
                             'score*score_multiplier > te_length [default: 50]')
    
    # overlap resolution
    parser.add_argument('--hit_overlap_pct', type=int, default=50,
                        help='For LTR & TIR, this is the maximum allowed amount of overlap (% of total length) ' +
                        ' between hits of same pair (but note that it is effectively 0 if min_hit_separation != 0)' +
                        '; for overlapping_direct & palindrome, it is the minimum overlap [default: 50]')
    parser.add_argument('--score_buffer_pct', type=int, default=25,
                        help='overlapping pairs with score not more than this % less than max score are ' +
                        ' evaluated based on pair length rather than score [default: 25]')
    parser.add_argument('--hit_length_buffer_pct', type=int, default=15,
                        help='overlapping pairs with scores within score_buffer_pct (above) and total hit (not pair) length ' +
                        ' within this % of the longest total hit length in the group are decided based on longest pair [default: 15]')
    parser.add_argument('--max_truncation', type=int, default=30,
                        help='if any hit (not pair) length is trimmed by more than this % (due to overlaps), the corresponding pair is discarded [default: 30]')

    return parser.parse_args()

class Hit:
    def __init__(self, line):
        d = line.strip().split('\t')
        self.chrom  = d[0]
        self.start  = int(d[1])
        self.end    = int(d[2])
        self.id     = d[3]
        self.score  = float(d[4])
        self.strand = d[5]
        self.orientation = self.id[-1]
        self.originalLength = self.length()
    
    def __repr__(self):
        return '\t'.join(str(x) for x in 
            [self.chrom, self.start, self.end, self.id, self.score, self.strand])
    
    def trimStart(self, x):
        self.start = self.start + x
        self.originalLength = self.length()
    
    def trimEnd(self, y):
        self.end = self.end - y
        self.originalLength = self.length()
    
    def length(self):
        return self.end - self.start
    
    def separation(self, other):
        s = sorted([self, other], key=lambda hit: hit.start)
        return s[1].start - s[0].end
    
    def overlaps(self, other):
        return self.separation(other) < 0
    
    def overlapLength(self, other):
        return -1 * self.separation(other)
    
    # Truncate any parts overlapping either part of hitPair and return the total percent that was truncated.
    def truncateBoth(self, hitPair):
        self.truncate(hitPair.a)
        self.truncate(hitPair.b)
        return 100 * (self.originalLength - self.length()) / self.originalLength
    
    # Truncate any parts overlapping other.
    def truncate(self, other):
        if other.start <= self.start and other.end > self.start:
            if other.end < self.end:
                self.start = other.end
            else:                                                               # set both start & end to -1 if overlap is complete (should not be in result, but just in case...)
                self.start = -1
                self.end = -1
        elif other.start > self.start and other.start < self.end:
            self.end = other.start

class HitPair:
    def __init__(self, a, b):       # sort a,b by start,end
        if a.id != b.id: exit('ERROR: IDs do not match!%s\n%s' % (a, b))        # validate that IDs are the same
        s = sorted([a,b], key=lambda hit: hit.start)
        self.a = s[0]
        self.b = s[1]
        self.start = min(a.start, b.start)
        self.end = max(a.end, b.end)
        self.orientation = a.orientation
        self.score = a.score
        self.trimmedSelfOverlap = 0
    
    def __repr__(self):
        return '\n'.join([self.a.__repr__(), self.b.__repr__()])

    def trimSelfOverlap(self):
        overlap = self.selfOverlap()
        if overlap > 0:
            x = int(math.ceil(overlap/2))
            y = overlap - x
            self.a.trimEnd(x)
            self.b.trimStart(y)
            self.trimmedSelfOverlap = overlap
        
    def selfOverlap(self):
        return self.a.overlapLength(self.b)
    
    def selfOverlapPct(self):
        return 100*self.selfOverlap()*2/self.totalHitLength()
    
    def trimmedSelfOverlapPct(self):
       return 100.0*self.trimmedSelfOverlap/self.originalHitLength()

    def length(self):
        return self.end - self.start
    
    def totalHitLength(self):
        return self.a.length() + self.b.length()
    
    def originalHitLength(self):
        return self.a.originalLength + self.b.originalLength
    
    def hitSeparation(self):
        return self.a.separation(self.b)
    
    def overlaps(self, other):
        return (self.a.overlaps(other.a) or self.a.overlaps(other.b) or
                self.b.overlaps(other.a) or self.b.overlaps(other.b))
    
    # Truncate any parts overlapping other and return the max percent of a or b that were truncated.
    def truncateAll(self, other):
        x = self.a.truncateBoth(other)
        y = self.b.truncateBoth(other)
        return max(x,y)

class HitPairCollection:
    def __init__(self):
        self.hitPairs = []
        self.end = -1
    
    def overlaps(self, hitPair):
        for p in self.hitPairs:
            if p.overlaps(hitPair):
                return True
        return False
    
    def add(self, hitPair):
        self.hitPairs.append(hitPair)
        self.hitPairs.sort(key=lambda x:x.start)
        if hitPair.end > self.end:
            self.end = hitPair.end
        None
    
    def addAll(self, hitPairCollections):
        for coll in hitPairCollections:
            for pair in coll.hitPairs:
                self.add(pair)

def readPairs(args):
    allPairs = []
    for line1 in args.infile:
        a = Hit(line1)
        try:
            b = Hit(args.infile.next())
        except StopIteration:
            exit('ERROR: Unexpected end of file. Last line has no matching pair.')
        pair = HitPair(a,b)
        if args.type == 'LTR' or args.type == 'TIR':
            pair.trimSelfOverlap()                                              # if hits overlap, truncate
        if passesFilters(pair, args):
            allPairs.append(pair)
        if not selfOverlapIsAcceptable(pair, args):
            None
    return allPairs

def process(args):
    sortedPairs = sorted(readPairs(args), key=lambda pair:pair.start)
    remaining = []
    for pair in sortedPairs:
        overlapping = []
        remove = []
        for coll in remaining:
            if coll.end <= pair.start:                                          # '<=' not just '<' b/c of 0-based indexing
                remove.append(coll)
                output(coll, args)
            elif coll.overlaps(pair):
                remove.append(coll)
                overlapping.append(coll)
        for coll in remove:                                                     # can't do it in-place in the loop b/c it messes up the iteration!!
            remaining.remove(coll)
        
        if len(overlapping) > 0:
            merged = HitPairCollection()
            merged.addAll(overlapping)
            merged.add(pair)
            remaining.append(merged)
        else:
            newColl = HitPairCollection()
            newColl.add(pair)
            remaining.append(newColl)
    for coll in remaining:
        output(coll, args)

def passesFilters(pair, args):
    if args.type == 'LTR' or args.type == 'TIR':
        return (
            typeIsCorrect(pair, args) and
            selfOverlapIsAcceptable(pair, args) and
            scoreIsAcceptable(pair, args) and
            pairLengthIsAcceptable(pair, args) and
            hitLengthIsAcceptable(pair, args) and
            separationIsAcceptable(pair, args))
    else:
        assert args.type == 'direct_overlapping' or args.type == 'palindrome'
        return (
            typeIsCorrect(pair, args) and
            selfOverlapIsAcceptable(pair, args) and
            scoreIsAcceptable(pair, args) and
            pairLengthIsAcceptable(pair, args) and
            hitLengthIsAcceptable(pair, args))

def typeIsCorrect(pair, args):
    if pair.orientation == '+':
        return args.type == 'LTR' or args.type == 'direct_overlapping'
    else:
        assert pair.orientation == '-'
        return args.type == 'TIR' or args.type == 'palindrome'

def selfOverlapIsAcceptable(pair, args):
    if args.type == 'LTR' or args.type == 'TIR':
        return pair.trimmedSelfOverlapPct() < args.hit_overlap_pct
    else:
        assert args.type == 'direct_overlapping' or args.type == 'palindrome'
        return pair.selfOverlapPct() > args.hit_overlap_pct

def scoreIsAcceptable(pair, args):
    if args.score_filter == 'filter_1':
        return pair.score >= args.min_score
    else:
        assert args.score_filter == 'filter_2'
        return (
            pair.score >= args.accept_score or (
            pair.score > args.min_score and pair.score*args.score_multiplier > pair.length()))               

def pairLengthIsAcceptable(pair, args):
    return pair.length() >= args.min_pair_length and pair.length() <= args.max_pair_length

def hitLengthIsAcceptable(pair, args):
    return max(pair.a.length(), pair.b.length()) < args.max_hit_length

def separationIsAcceptable(pair, args):
    return pair.hitSeparation() >= args.min_hit_separation

def output(hitPairCollection, args):
    pairs = sorted(hitPairCollection.hitPairs, key=lambda pair: pair.score, reverse=True)
    while len(pairs) > 0:
        # Find all pairs with scores 'close' to the top score
        topPair = None
        topScoring = []
        for p in pairs:
            if p.a.id == '8-':
                None
            if (len(topScoring) == 0 or
                topScoring[0].score - p.score < (args.score_buffer_pct/100.0)*topScoring[0].score):     # must be '100.0' not '100' or else integer division is used, so the test always fails
                    topScoring.append(p)
            else:
                break
        if len(topScoring) == 1:
            topPair = topScoring[0]
        else:
            maxTotalHitLength = 0
            for p in topScoring:
                maxTotalHitLength = max(p.totalHitLength(), maxTotalHitLength)
            for p in topScoring:
                if topPair is None or (
                    100*(maxTotalHitLength-p.totalHitLength())/maxTotalHitLength < args.hit_length_buffer_pct
                    and p.length() > topPair.length()):
                        topPair = p
        # Remove all pairs overlapping by >33% or truncate if overlapping by less.
        pairs.remove(topPair)                                                   # always remove at least 1 pair, guarantee loop completion
        remove = []
        for p in pairs:
            t = p.truncateAll(topPair)
            if t > args.max_truncation:                                         # careful: if truncation is 100%, positions are set to (-1,-1), which is bound to cause downstream errors
                remove.append(p)
        for coll in remove:                                                     # can't do it in-place in the loop b/c it messes up the iteration!!
            pairs.remove(coll)
        print(topPair, file=args.outfile)



if __name__ == '__main__':
    args = parse_cl()
    process(args)
    