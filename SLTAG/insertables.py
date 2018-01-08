#!/usr/bin/env python

'''
TAG tuples (i,j,k,l) stored in a convenient way to create combinations for inserting
into 'middle' sequence.
'''

import math
import itertools
import time
from sortedcollection import SortedCollection   # efficient sorted lists
from operator import itemgetter

# Like TupleList, but indexed by (yield, i), rather than just yield
# (where i is the first index of the tuple).
# TODO: TupleList and InsertablesList subclasses of the same class.
# Note: Expect the number of allowable insertions to be small (once the real grammar is implemented),
# since we only allow complete TEs to be inserted, which should be relatively rare. Therefore,
# insertions are implemented by iterating through the InsertablesList rather than the Tree tuples.
# (See insertion methods in Tree classes.)
# Note: Current implementation (v6) does not use the i index, since when doing insertions it iterates
# through all insertables, so indexing is done for the B1 List instead.
class InsertablesCollection():

    def __init__(self, dna):
        self.dna = dna
        self.byYield    = []
        self.D          = {}                                                 # (yield, i) => set(Tuple)
        self.L          = SortedCollection(key=lambda tup: tup.getStart())   # sorted by i
        
        # Combinations of tuples, stored in a collection sorted by start, of
        # collections sorted by end, of non-overlapping combinations of Tuples.
        # Can be accessed with like: groups = C[start][end], which will return
        # an collection of all combinations of non-overlapping insertables between (start, end).
        #self.C = CombinationCollection() 
        
    def __repr__(self):
        return self.D.values().__repr__()
    
    def beginStage(self):
        self.tuplesForStage = []

    def addToStage(self, tup):
        self.tuplesForStage.append(tup)
    
    def getTuplesForStage(self):
        return self.tuplesForStage

    def filterStage(self, aFilter):
        self.tuplesForStage = aFilter.filter(self.tuplesForStage)
    
    def completeStage(self):
        self.addAll(self.tuplesForStage)
        self.tuplesForStage = None
    
    def addAll(self, tuples):
        for tup in tuples:
            self.add(tup)
        #if len(self.D.values()) > 12:
        #    None
            print "\t\tAdded insertable:" + tup.derivation.toString(self.dna)
    
    def add(self, tup):
        (i,j,k,l) = tup.coordinates
        
        # add to byYield
        try:
            setByYield = self.byYield[tup.getYield()]
        except IndexError:
            setByYield = set()
            self.byYield.extend([set() for x in range(len(self.byYield), tup.getYield())])  # fill in any missing entries up to current tuple yield (which may never be used)
            self.byYield.append(setByYield)
        setByYield.add(tup)
        
        # add to D
        tupleSet = set()
        try:
            tupleSet = self.D[tup.getYield(), i]
        except KeyError:
            self.D[tup.getYield(), i] = tupleSet
        tupleSet.add(tup)
        
        # add to sorted list
        self.L.insert_right(tup)
 
    # return all tuples (i,j,k,l) such that i>=start and l<=end, *** sorted by i ***
    def getTuplesInSegment(self, start, end):
        tuples = []
        if (len(self.L) == 0) or (start > self.L[-1].getStart()):
            return tuples
        left = self.L.find_ge(start)  # tuple with lowest i >= start
        for tup in self.L[self.L.index(left):]:
            (i,j,k,l) = tup.coordinates
            if (tup.getStart() > end):
                return tuples
            if (tup.getEnd() > end):
                continue
            tuples.append(tup)
        return tuples
    
    # Return a list of all tuples
    def getAllTuples(self):
        tuples = []
        for tupleSet in self.D.itervalues():
            tuples.extend(tupleSet)
        return tuples
    
    # Return a list of tuples for the given yield
    def getTuplesOfYield(self, tYield):
        try:
            return list(self.byYield[tYield])
        except IndexError:
            return []                           # no parse entries for stage
    
    # Return a list of tuples for the given yield
    def getTuplesOfYieldAndI(self, tYield, i):
        if (tYield, i) in self.D:
            return self.D[tYield, i]
        else:
            return set()
        
    # Return a list of all combinations of non-overlapping insertables in region (start, end)
    def combinations(self, start, end):
        time0 = time.clock()
        segment = self.getTuplesInSegment(start, end)
        nonoverlapping = []
        limit = min(stag13.Probabilities.MAX_INSERTIONS_PER_MIDDLE+1, len(segment)+1)
        for k in range(1, limit):
            for combination in itertools.combinations(segment, k):
                try:
                    tup0 = None
                    for tup in combination:
                        if tup0 is None:
                            tup0 = tup
                        else:
                            if tup.overlaps(tup0):
                                raise ValueError()
                            tup0 = tup
                    nonoverlapping.append(combination)
                except ValueError: None
        #print "\t\t---> insertion combinations (%s..%s): %-3s %0.1e s" % (
        #    start, end, len(nonoverlapping), time.clock()-time0)
        return nonoverlapping

if __name__ == "__main__":

    #test = [(1,2),(1,3),(2,3),(2,4),(3,5),(4,7),(6,7),(6,8),(7,8),(7,9)]
    test = [(0,1),(0,3),(1,2),(3,4),(4,5),(5,6)]
    
    #for k in range(1, len(test)+1):
    #    combinations = itertools.combinations(test, k)
    #    for o in combinations: print o
    #    None
    
    insertables = InsertablesCollection()
    for t in test:
        insertables.add(stag13.Tuple((t[0],t[1],t[1],t[1]), None, None, '', 1))
    
    for c in insertables.combinations(3,5):
        print c
