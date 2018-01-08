#!/usr/bin/env python

"""
Modified from Rajasekeran 1996 and Rajasekeran et al. 2010.

Modifications:
* Use one list per tree rather than per node.
* Trees are not in normal form, but instead we know from the grammar how many terminals
    each tree produces (e.g., direct repeat trees produce 2 terminals at a time).
* For now there is only one tree, so no need to use a matrix to track the trees.

Grammar:

   X            X   Y       Y   Y
 / |            |   | \     |   |
b  X*           Y   b  Y    X   0
   | \
   X  b
   
   B1          <-- A1 -->   A2

Where b is any nucleotide (A,C,T, or G) and adjunction takes place only at X*.

This grammar recognizes any direct repeat (min length 1) surrounding any
sequence (min length 0). E.g. "gtag" is recognized, "tag" is not recognized.

(New feature:) Tree A2 allows the insertion of a B1 tree anywhere. (Correct??)

v08: Allow mismatches in LTRs; and tuple probabilities / scores.

v10: Add filtering for LTRs:
        For each i,j keep k,l only if P(i,j,k,l) >= a * max P(i,j,k',l'); l'-k' = l-k.
                   Similarly, for each k,l keep i,j only if P(i,j,k,l) >= a * max P(i',j',k,l); j'-i' = j-i.
        Note: The coefficient 'a' will determine the maximum number of (e.g.) consecutive mismatches.
        E.g., Roughly speaking, if a = {(mismatch P)/(match P)}^n, then at most n consecutiive mismatches will survive filtering.
        I.e., For each tuple length and each left substring (i,j), keep right substrings only if they are not too much
            less probable than the best matching right substring. Ditto for left substrings.
        
v11: Add filters for exact matches.
        Either (i,j,k,l) == (i',j',k',l'), or (i,l) == (i',l') and j == k and j' == k'.
        This is primarily to prevent unlikely insertions from propogating.
    Also, add probabilities for operations (adjoin, substitute, insert).

v12: Self-BLAST for LTRs:
        1. Do self-BLASTS (see self-blastn.py).
        2. Simple probability scheme for adjoining (B1 tree):
            i.   In alignment => v. high probability
            ii.  Not in alignment => v. low probability
        3. Null Model Filters
            There is currently no way to eliminate v. low probability alignments, except if there is a better alignment to the same (i,j),
            which doesn't eliminate low probability alignments from regions with no alignment at all. This is what null model filters
            will do. They provides a baseline probability against which to measure any tuple. Eventually, there should be a different null
            models for each tree (or tuple list?). Initially, just implement one for the B1 tree adjoinables list: RandomNullFilter.
            
            *** Unfortunately, this approach proved too slow (under the current, unoptimized implementation), to allow parsing of realistic length TEs
                (e.g. 5 kbp). So instead, we will limit searches for LTRs to the areas identified by BLASTN.
            
        Note: Moved sequence input from command-line to a file because from now on will be working
            with longer, more realistic sequences (which won't easily fit on the command-line).
        
    Tree A1 no longer parses every possible middle sequence.
            
    Fixed: match/mismatch probability table.
    
v13: Changed how "middle" sequence is handled. Motivation: middle is essentially context-insensitive, so there's no need to continually parse every
        possibility. Instead, use an HMM-like approach and simply 'scan' the sequence when needed. For now, the scan simply consists of generating
        a string of whatever length is needed to join the LTRs (only when the final length of the contiguous sequence matches the stage).
            - As a result, need to change how "evaluate" is done. Rather than direclty accessing the A1 list from B1, have B1 ask A1 for the requisite
                sequence, which is then generated on the fly.
    
TO DO next: - Allow gaps (use BLASTN results?). Also allow N's (allow them to match or not?)
            - Find TEs embedded in genomic sequence, not just full TEs.

TO DO: Note: Was going to implement this in v12, but realized that the following approach will likely not work after implementing
        self-BLAST requirement for LTRs (since not all DNA will be parsed using the LTRs). So do self-BLAST first...
    Add TSDs. Possible TSDs are already being parsed by tree B1, so instead of recalculating them, just give a 'bonus' (high probability)
        for inserting a B1 tree (and later, other TE termini) into a B1 tree (of the right length(s)).

------------------------------------------------------------------------

Features (to do):

* handle 'N' (and other ambiguity codes?)
    * will probably be easy to handle with score-based base pairing (to allow for mismatches & indels)

Ideas (to do):

* Include an additional (5th) value in the tuple, specifying the probability of the tuple.
    E.g., At first, record the match likelihood, with something like:
        1 = match, 0.3 = transition, 0.2 = transversion, 0.1 = gap.
        * Can we incorporate affine gap scoring?
    * Only 'grow' the string where probability exceeds some threshold, either static
        (e.g., >0.5) or relative (e.g., among the top Q scores at that length).
        

Optimizations (to do):

* Merge contiguous segments of single, self-adjoining trees (esp. LTR) in a single
    step. Track such contiguous 'homogeneous' segments separately from those
    containing merges of multiple trees.
        * Akin to extension for BLAST or scanning for HMMs.
        * Similarly, don't store all possible sub-parses. We know that any subparse
            is allowable.
                * Complication: What about probabilities? If we ever need to
                    use a substring, won't know probability... Not sure how to
                    deal with this.
                    * Solution: If this is a rare event, just compute when & if needed.
                        Won't cost more than computing before needed.
        * Without this, we have many sets of trees for repeats at each size, because
            every substring is a repeat! This is why memory & processing are high.
        * Complication: how to proceed in steps?
            * Perhaps use stages instead of gradually progressing in stepwise increasing
                length increments. Take advantage of 
        * Complication: how to deal with non-exact matches?
            * At first, just do exact matches? At least this will cover most of the
            LTRs (and TIRs?), which should be the easiest to parse, rather than
            taking the longest.


"""

import abc                                      # abstract base classes
import re                                       # regular expressions
import math
import time                                     # to measure performance
from operator import attrgetter                 # for sorting
import insertables



class Probabilities:
    
    MIN_INSERTABLE_LENGTH       = 100
    MAX_INSERTIONS_PER_MIDDLE   = 5
    ALPHA_MATCH                 = 1e-9
    ALPHA_INSERTABLES           = 0.5
    ALPHA_IDENTICAL_COORDS      = 0.999     # don't set to 1 b/c floating point imprecision may causes issues??
    
    A1_PROBABILITY = 0.25     # "Middle" sequence. Should be > mismatch P but << than match P. ??? Needs to be 0.25 (to sum to 1) ???
    
    # Blastn alignment adjustments
    #class Blastn:
    #    MISSING     = 1e-2
    #    FOUND       = 1 - MISSING
    
    # Match/mismatch probabilities
    class B1Matrix:
        # Should the sum add to one, or the sum of the square roots (since there are 2 terminals being generated)?
        M = math.log(0.90/4  # match         
        I = 0.12/4  # transition      # Needs to be < a1Probability^2 ??
        V = 0.08/8  # transversion    # Needs to be < a1Probability^2 ??
        X = {
            ('a','a'): M, ('c','c'): M, ('g','g'): M, ('t','t'): M,     # match
            ('a','g'): I, ('g','a'): I, ('c','t'): I, ('t','c'): I,     # transitions
            ('a','t'): V, ('t','a'): V, ('g','t'): V, ('t','g'): V,     # transversions
            ('a','c'): V, ('c','a'): V, ('c','g'): V, ('g','c'): V,     # transversions
        }
    
    # Operations
    # TO DO: How exactly should operation probabilities add to 1? Seems like it should be
    # normalized somehow, but how??? E.g., Is B1_ADJOIN an equivalent operation to A1_SUBSTUTUTE?
    # Perhaps simply adjoin, which extends both substrings (i,j) and (k,l), should be
    # equivalent to 2 substitutions, which extends only one substring (think in terms of
    # number of new "links" formed between sequences)?
    class Operations:
        INSERT          = 1e-2
        B1_ADJOIN       = 0.95
        B1_SUBSTITUTE   = 1 - B1_ADJOIN - INSERT     # square-root of B1_ADJOIN to normalize
        A1_SUBSTITUTE   = 1 - INSERT



    # Class methods
    @staticmethod
    def calcSubstitutionProb(length):
        assert length > 0
        prob = math.pow(Probabilities.A1_PROBABILITY, length)                     # P (emissions)
        return prob * math.pow(Probabilities.Operations.A1_SUBSTITUTE, length)    # P (substitutions)



class Filter:
    @abc.abstractmethod
    def filter(self, tuples):
        """Filter tuples."""
        return



# Filter tuples that have identical coordinates, but probability' < alpha*probability.
# Identical coordinates means that (i,j,k,l) == (i',j',k',l'), or
# (i,l) == (i',l') and j == k and j' == k' (but j might not == j' nor k == k').
# (The latter case is for contiguous substrings.)
class IdenticalCoordinatesFilter(Filter):
    
    # Generate set of tuples used for filtering.
    # Input: Iterator over tuples, all to be used for filtering.
    def __init__(self, alpha, tuples):
        self.alpha = alpha

        # construct data structures:
        # 1. j == k: (i,l) => maxP
        # 2. j != k: (i,j,k,l) => maxP
        self.ilMap = {}
        self.ijklMap = {}
        for tup in tuples:
            (i,j,k,l) = tup.coordinates
            if j == k:
                try:
                    self.ilMap[(i,l)] = max(self.ilMap[(i,l)], tup.getProbability())
                except KeyError:
                    self.ilMap[(i,l)] = tup.getProbability()
            else:
                try:
                    self.ijklMap[(i,j,k,l)] = max(self.ijklMap[(i,j,k,l)], tup.getProbability())
                except KeyError:
                    self.ijklMap[(i,j,k,l)] = tup.getProbability()
    
    # Perform filtering based on the tuples that must be previously passed to setTuples().
    # Returns the tuples that pass through the filter.
    def filter(self, tuples):
        filtered = []
        for tup in tuples:
            (i,j,k,l) = tup.coordinates
            if j == k:
                if tup.getProbability() >= self.alpha * self.ilMap[(i,l)]:
                    filtered.append(tup)
            else:
                if tup.getProbability() >= self.alpha * self.ijklMap[(i,j,k,l)]:
                    filtered.append(tup) 
        return filtered



# Filter tuples out that have the same i,j, but probability' < alpha*probability
class MatchFilter(Filter):
    
    def __init__(self, alpha):
        self.alpha = alpha
    
    def filter(self, tuples):
        # construct data structure: (i,j) => maxP
        ijMap = {}
        for tup in tuples:
            (i,j,k,l) = tup.coordinates
            try:
                ijMap[(i,j)] = max(ijMap[(i,j)], tup.getProbability())
            except KeyError:
                ijMap[(i,j)] = tup.getProbability()
        
        # filter
        filtered = []
        for tup in tuples:
            (i,j,k,l) = tup.coordinates
            if tup.getProbability() >= self.alpha * ijMap[(i,j)]:
                filtered.append(tup)

        return filtered
    
        # TO DO: don't just throw away these data structures and return a list... they may be useful to keep
    
    
    
# Filter out tuples that have P' < alpha*P, where
# P' is the probability of the tuple and
# P is the average probability of a random DNA sequence of the same yield as the tuple,
# i.e. 0.25^yield
class InsertablesFilter(Filter):
    
    def __init__(self, alpha):
        self.alpha = alpha

    def filter(self, tuples):
        filtered = []
        for tup in tuples:
            if (tup.getProbability() * Probabilities.Operations.INSERT >=
                self.alpha * Probabilities.calcSubstitutionProb(tup.getEnd() - tup.getStart())):
                    filtered.append(tup)
            else:
                None
        return filtered



# Container for initial and auxiliary trees
class Grammar:
    
    def __init__(self, selfBlastn, dna):
        
        # BLASTN-defined regions in which to look for termini, of form [(start0,end0),(start1,end1)],
        # which define the two blocks that align to one another.
        self.selfBlastn = selfBlastn
        self.dna = dna
        
        # Initial Trees
        self.treeA1 = TreeA1(self)
        
        # Auxilliary Tres
        self.treeB1 = TreeB1(self)
        
        # Insertables, which are common to all trees
        self.insertables = insertables.InsertablesCollection(dna)
        
        # Probabilities & Filters
        #self.probabilities = Probabilities()
        self.matchFilter = MatchFilter(Probabilities.ALPHA_MATCH)
        self.insertablesFilter = InsertablesFilter(Probabilities.ALPHA_INSERTABLES)
    
    # Get all insertables with i>=start, l<=end, sorted by i
    def getInsertablesCombinations(self, start, end):
        return self.insertables.combinations(start, end)
    
    def parse(self):
        dna = self.dna
        lcDna = dna.lower()
        if not self.validateDNA(lcDna):
            return
        
        print "\nPartial contiguous parse(s):"
        for stage in range(1, len(dna)+1):
            if (stage == 12):
                None
            print "Stage %s:" % stage
            self.treeB1.beginStage()
            #self.treeA1.beginStage()
            self.insertables.beginStage()
            self.treeB1.evaluate(lcDna, stage)
            #self.treeA1.evaluate(lcDna, stage)
            self.filter()
            self.treeB1.completeStage()
            #self.treeA1.completeStage()
            self.insertables.completeStage()
            #self.printParses(stage)
            #self.printInsertables(stage)
            None
        
        self.printInsertables(len(dna))
        self.recognizeInput(dna)
    
    def filter(self):        
        start = time.clock()
        
        # Filter insertables: remove insertables that are not much more probable than random
        self.insertables.filterStage(self.insertablesFilter)
        
        # Filter LTRs: remove low probability alignments where better ones exist
        self.treeB1.filterStage(self.matchFilter)
            
        # Filter all tuples for low probability derivations
        allTuples = []
        allTuples.extend(self.insertables.getTuplesForStage())
        #allTuples.extend(self.treeA1.nonadjoinables.getTuplesForStage())
        allTuples.extend(self.treeB1.adjoinables.getTuplesForStage())
        allTuples.extend(self.treeB1.nonadjoinables.getTuplesForStage())
        identicalCoordFilter = IdenticalCoordinatesFilter(
            Probabilities.ALPHA_IDENTICAL_COORDS, allTuples)
        self.treeB1.filterStage(identicalCoordFilter)
        #self.treeA1.filterStage(identicalCoordFilter)
        self.insertables.filterStage(identicalCoordFilter)
        
        print "\t%-18s\t %0.1e s" % ("Filter:", time.clock()-start)
        
    def validateDNA(self, lcDna):
        if re.compile('[acgt]').match(lcDna):
            return True
        else:
            print "ERROR: Illegal characters found in input sequence."
            return False
        
    def recognizeInput(self, dna):
        #completeParses = self.insertables.getTuplesOfYield(len(dna))
        completeParses = self.treeB1.adjoinables.getTuplesOfYield(len(dna))     # Note: Always of form (0,p,p,len(dna))
        for tup in completeParses:
            (i,j,k,l) = tup.coordinates
            assert j == k
        nComplete = len(completeParses)
        if nComplete == 0:
            print "\n---> Input is NOT recognized by the grammar.\n"
            #self.dump()
            return False
        else:
            print "\n---> Input IS recognized by the grammar.\n"
            print "\n%s Complete Parses:" % nComplete
            self.printDerivations(completeParses)
            return True
    
    def printParses(self, stage):
        parses = self.treeB1.getTuplesOfYield(stage)     # Note: Always of form (0,p,p,len(dna))
        print "\n%s Parses:" % len(parses)
        self.printDerivations(parses)
    
    def printInsertables(self, stage):
        insertables = self.insertables.getTuplesOfYield(stage)     # Note: Always of form (0,p,p,len(dna))
        print "\n%s Insertables:" % len(insertables)
        self.printDerivations(insertables)
    
    def printDerivations(self, tuples):
        tuples.sort(key=attrgetter('logP'), reverse=False)
        for count in range(0, len(tuples)):
            tup = tuples[count]
            (i,j,k,l) = tup.coordinates
            (s1, s2) = (dna[i:j], dna[k:l])
            print "\nDerivation #%s: (%d,%d,%d,%d) %s+%s %s" % (
                    len(tuples)-count,
                    i,j,k,l,s1,s2, math.pow(10, tup.logP))
            print "=" * (len(dna)+48)
            #print; print tup.derivationString
            print tup.derivation.toString(dna)
    
    def dump(self):
        print "\n\n\n" + "*" * 30
        print "DUMP:"
        print "\nInsertables:\n%s" % self.insertables
        print "\nTree B1:\n%s" % self.treeB1.__repr__()
        print "\nTree A1:\n%s" % self.treeA1.__repr__()

    def isInsertable(self, tup, dna):
        (i,j,k,l) = tup.coordinates
        if j == k and l-i >= min(Probabilities.MIN_INSERTABLE_LENGTH, len(dna)):
            #(s1, s2) = (dna[i:j], dna[k:l])
            #print "%s+%s \t length %d \t (%d,%d,%d,%d)" % (s1,s2,l-i,i,j,k,l)           
            return True
        return False



# Tuple coordinates plus other data, such as the sub-tuples that formed this tuple for back-tracing.
class Tuple:
    
    def __init__(self, (i,j,k,l), leftTuple, rightTuple, label, probability=None):
        '''If proabability is specified, it is used. Otherwise, leftTuple and rightTuple must not
        be none and probability is the product of their probabilities'''
        self.coordinates = (i,j,k,l)
        self.label = label
        #self.derivationString = label + ":"
        
        if probability is not None:
            self.setProbability(probability)
        else:
            self.setProbability(leftTuple.getProbability() * rightTuple.getProbability())
            
        leftDerivation = None
        if leftTuple is not None:
            leftDerivation = leftTuple.derivation
        #    self.derivationString = self.derivationString + leftTuple.derivationString + ","
        #else:
        #    self.derivationString = self.derivationString + "E,"
        rightDerivation = None
        if rightTuple is not None:
            rightDerivation = rightTuple.derivation
        #    self.derivationString = self.derivationString + rightTuple.derivationString
        #else:
        #    self.derivationString = self.derivationString + "E"
        self.derivation = DerivationTree((i,j,k,l), label, leftDerivation, rightDerivation, self.getProbability())  # TO DO: Ugly to store this in both Tuple & DerivationTree
  
    def __repr__(self):
        return "%s %s" % (self.coordinates.__repr__(), self.getProbability())
    
    def __eq__(self, other):
        if isinstance(other, Tuple):
            return (self.coordinates == other.coordinates and
                  self.label == other.label and
                  self.logP == other.logP)                                    # TO DO: Also test for derivation equality (otherwise some derivation trees may be missed)
        else:
            return False
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __hash__(self):
        return hash((self.coordinates, self.logP))
    
    def getStart(self):
        (i,j,k,l) = self.coordinates
        return i
    
    def getEnd(self):
        (i,j,k,l) = self.coordinates
        return l
    
    def overlaps(self, other):
        if self.getStart() <= other.getStart():
            return self.getEnd() > other.getStart()     # Not if end == start because i specifies the index *before* the start)
        else:
            return self.getStart() < other.getEnd()
    
    def scaleProbability(self, coeff):
        self.setProbability(self.getProbability() * coeff)
    
    def getYield(self):
        (i,j,k,l) = self.coordinates
        return j-i+l-k
    
    def getProbability(self):
        return math.pow(10, self.logP)
    
    def setProbability(self, prob):
        self.logP = math.log10(prob)



# To track derivation trees
class DerivationTree:
    
    def __init__(self, (i,j,k,l), label, left, right, probability):
        self.left = left
        self.right = right
        self.coordinates = (i,j,k,l)
        self.label = label
        self.probability = probability
    
    # TO DO: __repr__(), not toString()
    def toString(self, dna):
        (i,j,k,l) = self.coordinates
        (s1, s2) = (dna[i:j], dna[k:l])
        if self.left is None or self.right is None:
            if s2 == "":
                return "\n%s:   %s" % (s1, self.label)
            else:
                return "\n%s+%s: %s" % (s1, s2, self.label)
        else:
            fill = " " * (len(dna) - (j-i+l-k))
            return "\n\t%s+%s %s %s %s %s %s %s" % (
                s1, s2, fill, self.label, self.coordinates, self.probability,
                self.left.toString(dna), self.right.toString(dna))


# Lists (actually, sets) of tuples, indexed by yield
class TupleList:
    
    def __init__(self):
        # TO DO: integrate these data structures into one (to save memory and a bit of time?)
        self.allTuples      = []
        self.byYield        = {}     
        self.byYieldAndJ    = {}
        self.byYieldAndL    = {}
        self.byIJKL         = {}    # for adjoin
        
    def __repr__(self):
        return self.byYield.__repr__()
    
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
    
    # TO DO: shouldn't allow both add & addToStage -- confusing. Make all TupleLists filtered and disallow direct usage of add().
    def add(self, tup):
        (i,j,k,l) = tup.coordinates
        
        # add to allTuples
        self.allTuples.append(tup)
        
        # add to byYield
        setByYield = set()
        try:
            setByYield = self.byYield[tup.getYield()]
        except KeyError:
            #self.byYield.extend([set() for x in range(len(self.byYield), tup.getYield())])  # fill in any missing entries up to current tuple yield (which may never be used)
            #self.byYield.append(setByYield)
            self.byYield[tup.getYield()] = setByYield
        setByYield.add(tup)
        
        # add to byYieldAndJ
        setByJ = set()
        try:
            setByJ = self.byYieldAndJ[tup.getYield(), j]
        except KeyError:
            self.byYieldAndJ[tup.getYield(), j] = setByJ
        setByJ.add(tup)
        
        # add to byYieldAndL
        setByL = set()
        try:
            setByL = self.byYieldAndL[tup.getYield(), l]
        except KeyError:
            self.byYieldAndL[tup.getYield(), l] = setByL
        setByL.add(tup)
        
        # add to byIJKL
        setByIJKL = set()
        try:
            setByIJKL = self.byIJKL[(i,j,k,l)]
        except KeyError:
            self.byIJKL[(i,j,k,l)] = setByIJKL
        setByIJKL.add(tup)
 
    # Return a list of all tuples
    def getAllTuples(self):
        return self.allTuples
        #tuples = []
        #for y in range(len(self.byYield)):
        #    tuples.extend(self.getTuplesOfYield(y))
        #return tuples
    
    # Return a list of tuples for the given yield
    def getTuplesOfYield(self, tYield):
        try:
            return self.byYield[tYield]
        except KeyError:
            return []                           # no parse entries for stage
        
    def getTuplesOfYieldAndJ(self, tYield, j):
        if (tYield, j) in self.byYieldAndJ:
            return self.byYieldAndJ[tYield, j]
        else:
            return set()
    
    def getTuplesOfYieldAndL(self, tYield, l):
        if (tYield, l) in self.byYieldAndL:
            return self.byYieldAndL[tYield, l]
        else:
            return set()
    
    def getTuplesOfIJKL(self, (i,j,k,l)):
        if ((i,j,k,l)) in self.byIJKL:
            return self.byIJKL[(i,j,k,l)]
        else:
            return set()



# Abstract Tree class
class Tree(object):
    __metaclass__ = abc.ABCMeta
    
    #name        = ""
    #nTerminals  = 0             # number of terminals 
    #list                        # list of substrings that can be parsed (i,j,k,l)       
    
    #def __init__(self, name):
    #    self.name = name
    
    def beginStage(self):
        for tuples in self.getTupleLists():
            tuples.beginStage()
    
    def filterStage(self, aFilter):
        for tuples in self.getTupleLists():
            tuples.filterStage(aFilter)
    
    def completeStage(self):
        start = time.clock()
        for tuples in self.getTupleLists():
            tuples.completeStage()
        print "\t%s%-18s\t %0.1e s" % (self.name, " completeStage:", time.clock()-start)
    
    @abc.abstractmethod
    def getTupleLists(self):
        """Return all tuple lists for this tree."""
        return
    
    @abc.abstractmethod
    def evaluate(self, dna, stage):
        """Update list according to terminals (dna) and stage (l)."""
        return



# Tree for LTRs.
class TreeB1(Tree):
    
    def __init__(self, grammar):
        self.grammar = grammar
        self.adjoinables = TupleList()
        self.nonadjoinables = TupleList()
        self.name = "B1"
        
    def __repr__(self):
        return "Adjoinables:\n" + self.adjoinables.__repr__() + "\nNon-adjoinables:\n" + self.nonadjoinables.__repr__()
    
    def __str__(self):
        return self.__repr__()
    
    def getTupleLists(self):
        return (self.adjoinables, self.nonadjoinables)
    
    def getAllTuples(self):
        tuples = []
        tuples.extend(self.adjoinables.getAllTuples())
        tuples.extend(self.nonadjoinables.getAllTuples())
        return tuples
        
    def getTuplesOfYield(self, tYield):
        tuples = []
        tuples.extend(self.adjoinables.getTuplesOfYield(tYield))
        tuples.extend(self.nonadjoinables.getTuplesOfYield(tYield))
        return tuples

    def evaluate(self, dna, stage):
        if (stage < 2):
            return
        elif (stage == 2):
            self.initialize(dna)
        else:
            self.adjoin(stage, dna)
            self.substitute(stage, dna)
            self.insert(stage, dna)
    
    def initialize(self, dna):
        start = time.clock()
        nPairs = 0
        for i in range(len(dna)):
            for j in range(i+1, len(dna)):
                if not self.grammar.selfBlastn.containsPair(i+1, j+1):
                    continue
                #if (dna[i:i+1] == dna[j:j+1]):
                prob = Probabilities.B1Matrix.X[dna[i:i+1], dna[j:j+1]]
                #if self.grammar.selfBlastn.containsPair(i+1, j+1):             # this approach proved to be too slow
                #    prob = prob * Probabilities.Blastn.FOUND
                #else:
                #    prob = prob * Probabilities.Blastn.MISSING
                tup = Tuple((i,i+1,j,j+1), None, None, "B  ", prob)
                self.adjoinables.addToStage(tup)
                if self.grammar.isInsertable(tup, dna):
                    self.grammar.insertables.addToStage(tup)
                nPairs = nPairs + 1
        print "\tB1 initialize:\t%-8s %0.1e s" % (nPairs, time.clock()-start)
    
    # Tree B1: Adjoin B1 to itself. Result is still adjoinable.
    # For this tree, self adjoin always results in extending by 2 (identical) terminals,
    # i.e. (i,j,k,l) + (j,j+1,l,l+1) => (i,j+1,k,l+1)
    def adjoin(self, stage, dna):
        #if (stage == 12):
        #    None
        start = time.clock()
        nAdjoins = 0
        self.adjoinables.beginStage()
        for t1 in self.adjoinables.getTuplesOfYield(stage-2):
            (i,j,k,l) = t1.coordinates
            if (j == k):    # don't extend j beyond k (i<=j<=k<=l)
                continue
            for t2 in self.adjoinables.getTuplesOfIJKL((j,j+1,l,l+1)):
                newTuple = Tuple((i,j+1,k,l+1), t1, t2, "Ba ")
                newTuple.scaleProbability(Probabilities.Operations.B1_ADJOIN)
                self.adjoinables.addToStage(newTuple)
                if self.grammar.isInsertable(newTuple, dna):
                    self.grammar.insertables.addToStage(newTuple)
                nAdjoins = nAdjoins + 1
        print "\tB1 adjoin:\t%-8s %0.1e s" % (nAdjoins, time.clock()-start)
    
    # Tree A1: Substitute into B1. Result is no longer adjoinable. This is the "middle" part of the TE, between TIRs.
    # For now, substitute every allowed A1 (even though we know most will not be used).
    # Could iterate through a1Tree.nonadjoinables instead, but this one should be shorter and therefore faster
    # (In general, may want to keep track of the list lengths and iterate through the shorter one -- TODO.)
    def substitute(self, stage, dna):
        start = time.clock()
        nSubstitutions = 0
        for termini in self.adjoinables.getAllTuples():
            (i,j,k,l) = termini.coordinates
            if stage == 26:# and (i,j,k,l) == (0,3,9,12):
                None
            if l-i == stage and j<k:
                for middle in self.grammar.treeA1.getTuplesForSegment(j,k):        
                    newTuple = Tuple((i,k,k,l), termini, middle, "Bs ")
                    newTuple.scaleProbability(Probabilities.Operations.B1_SUBSTITUTE)
                    self.nonadjoinables.addToStage(newTuple)
                    if self.grammar.isInsertable(newTuple, dna):
                        self.grammar.insertables.addToStage(newTuple)
                    nSubstitutions = nSubstitutions + 1
        print "\tB1 substitute:\t%-8s %0.1e s" % (nSubstitutions, time.clock()-start)
    
    # Tree A2: Insert permitted tuples inside either or both TIRs. Permitted tuples are currently
    # contiguous parses that
    # Iterate through allowed insertions, which should be much fewer in number than the tuples in this tree.
    # Note: Current implementation inserts to the right of the current position in the growing left or right TIR,
    # which has the effect of allowing insertions to the right of the right TIR, but not to the left of the left TIR.
    # This may eventually need to be changed.
    def insert(self, stage, dna):
        start = time.clock()
        nInsertions = 0
        for insertable in self.grammar.insertables.getAllTuples():
            (s,t,u,v) = insertable.coordinates
            if (t != u):
                raise Exception     # sanity check
            remainingYield = stage - insertable.getYield()
            for left in self.adjoinables.getTuplesOfYieldAndJ(remainingYield, s):    # left TIR
                (i,j,k,l) = left.coordinates
                if v > k:
                    continue    # not enough room between TIRs to insert (TODO?: optimize lookup to avoid these cases)
                newTuple = Tuple((i,v,k,l), insertable, left, "Bil")
                newTuple.scaleProbability(Probabilities.Operations.INSERT)
                self.adjoinables.addToStage(newTuple)
                if self.grammar.isInsertable(newTuple, dna):
                        self.grammar.insertables.addToStage(newTuple)
                nInsertions = nInsertions + 1
            for right in self.adjoinables.getTuplesOfYieldAndL(remainingYield, s):    # right TIR
                (i,j,k,l) = right.coordinates
                newTuple = Tuple((i,j,k,v), insertable, right, "Bir")
                newTuple.scaleProbability(Probabilities.Operations.INSERT)
                self.adjoinables.addToStage(newTuple)
                if self.grammar.isInsertable(newTuple, dna):
                        self.grammar.insertables.addToStage(newTuple)
                nInsertions = nInsertions + 1
        print "\tB1 insert:\t%-8s %0.1e s" % (nInsertions, time.clock()-start)



# Tree for random internal sequence (between termini).
# Note: As currently implemented, this tree does very little and could easily be implemented
# more efficiently. It's true purpose is to have a second tree in the grammar to see if the
# parser will work with >1 trees, especially with different numbers of terminals in each tree.
class TreeA1():
    
    def __init__(self, grammar):
        self.grammar = grammar
        #self.nonadjoinables = TupleList()
        self.name = "A1"
        
    def __repr__(self):
        return self.nonadjoinables.__repr__()
    
    def getTuplesForSegment(self, start, end):
        # insertables in each combination must be sorted by i for the following to work
        time0 = time.clock()
        combinations = self.grammar.getInsertablesCombinations(start, end)
        tuples = []
        
        # with no insertions
        prob = Probabilities.calcSubstitutionProb(end-start)
        tuples.append(Tuple((start, end, end, end), None, None, "A  ", prob))
        
        # with insertions
        for combn in combinations:    # each group forms one tuple
            join = None
            pJoin = 1
            startScan = start
            for insertable in combn:
                # 1. Scan to insert (if needed)
                if insertable.getStart() != startScan:
                    endScan = insertable.getStart()
                    pScan = Probabilities.calcSubstitutionProb(endScan-startScan)
                    scan = Tuple((startScan, endScan, endScan, endScan), None, None, "As ", pScan)
                    pJoin = pJoin * pScan
                    join = Tuple((start, endScan, endScan, endScan), join, scan, "Aj ", pJoin)
                
                # 2. Insert
                pJoin = pJoin * insertable.getProbability() * Probabilities.Operations.INSERT
                endInsert = insertable.getEnd()
                join = Tuple((start, endInsert, endInsert, endInsert), join, insertable, "Ai ", pJoin)
                startScan = endInsert
                
            # 3. Scan to end
            if join.getEnd() < end:
                endScan = end
                pScan = Probabilities.calcSubstitutionProb(endScan-startScan)
                scan = Tuple((startScan, endScan, endScan, endScan), None, None, "As ", pScan)
                pJoin = pJoin * pScan
                join = Tuple((start, endScan, endScan, endScan), join, scan, "Aj ", pJoin)
            tuples.append(join)
            #print join.derivation.toString(dna)
            #print "-"*10
            None
        elapsed = time.clock()-time0
        if elapsed > 1:
            print "\tA1:       \t%-8s %0.1e s" % (len(tuples), elapsed)
        return tuples
    
    #def getTupleLists(self):
    #    return (self.nonadjoinables,)
    #
    #def evaluate(self, dna, stage):
    #    if (stage == 1):
    #        self.initialize(dna)
    #    else:
    #        self.substitute(stage, dna)
    #        self.insert(stage, dna)
    #        
    #def initialize(self, dna):
    #    start = time.clock()
    #    for i in range(len(dna)):
    #        newTuple = Tuple((i,i+1,i+1,i+1), None, None, "A  ", Probabilities.A1_PROBABILITY)
    #        self.nonadjoinables.addToStage(newTuple)
    #    print "\tA1 initialize:\t%-8s %0.1e s" % (len(dna), time.clock()-start)
    #
    ## (i,j,j,j) => (i,j+1,j+1,j+1)
    #def substitute(self, stage, dna):
    #    start = time.clock()
    #    nSubstitutions = 0
    #    for tup in self.nonadjoinables.getTuplesOfYield(stage-1):
    #        (i,j,k,l) = tup.coordinates
    #        if j != k:  # Note: l may be > k because of insertions (and only for this reason)
    #            raise Exception()   # sanity check
    #        for terminalTuple in self.nonadjoinables.getTuplesOfIJKL((l,l+1,l+1,l+1)):
    #            newTuple = Tuple((i,l+1,l+1,l+1), tup, terminalTuple, "As ")
    #            newTuple.scaleProbability(Probabilities.Operations.A1_SUBSTITUTE)
    #            self.nonadjoinables.addToStage(newTuple)
    #            nSubstitutions = nSubstitutions + 1
    #    print "\tA1 substitute:\t%-8s %0.1e s" % (nSubstitutions, time.clock()-start)
    #
    ## Tree A2: Insert into "middle" sequence.
    #def insert(self, stage, dna):
    #    start = time.clock()
    #    nInsertions = 0
    #    for insertable in self.grammar.insertables.getAllTuples():
    #        (s,t,u,v) = insertable.coordinates
    #        if (t != u):
    #            raise Exception     # sanity check
    #        remainingYield = stage - insertable.getYield()
    #        for tup in self.nonadjoinables.getTuplesOfYieldAndL(remainingYield, s):
    #            (i,j,k,l) = tup.coordinates
    #            if (j != k):
    #                raise Exception     # sanity check
    #            newTuple = Tuple((i,j,k,v), insertable, tup, "Ai ")
    #            newTuple.scaleProbability(Probabilities.Operations.INSERT)
    #            self.nonadjoinables.addToStage(newTuple)
    #            #if self.grammar.isInsertable(newTuple, dna):
    #                #self.grammar.insertables.add(newTuple)                        # not insertable because not flanked by a TIR
    #    print "\tA1 insert:\t%-8s %0.1e s" % (nInsertions, time.clock()-start)



if __name__ == "__main__":
    import sys      # for command-line arguments
    import time     # for performance metrics
    import selfblastn01

    start = time.clock()
    if len(sys.argv) != 2:
        print "Usage: <program> DNA"
        exit(0)
    dna = sys.argv[1]
    print "\nInput:\n%s" % (dna)
    print "1234567890" * 3
    print "         1         2         3\n"
    print "Input length: %d" % len(dna)
    print "Probability of random sequence of this length: %s" % math.pow(0.25, len(dna))
    grammar = Grammar(selfblastn01.AllMatch(), dna)
    grammar.parse()
    end = time.clock()
    print; print "Elapsed time: %f s" % (end - start)
    #raw_input("Press <return> to exit: ")
    print
#
## End