#!/usr/bin/env python

"""
Simple linear tree adjoining grammar demo/test parser.
See Uemura et al. 1999
"""



import sys
import time



class BetaTree:
    """An auxiliary tree"""
    
    name = ""
    LU = 0
    LD = 0
    RD = 0
    RU = 0
    
    def L(self):
        return self.LU + self.LD
    
    def R(self):
        return self.RU + self.RD
    
    def __init__(self, aName, lu, ld, rd, ru):
        self.name = aName
        self.LU = lu
        self.LD = ld
        self.RD = rd
        self.RU = ru    



# 4D square matrix
def printMatrix(matrix):
    n = len(matrix[0][0][0])
    for i in range (0, n):
        for j in range (0, n):
            for k in range (0, n):
                for l in range (0, n):
                    if (matrix[i][j][k][l] != 0):
                        print "M(%d,%d,%d,%d) = %s" % (i,j,k,l, matrix[i][j][k][l])



def initialize(matrix, tree, dna):                           # TODO: multiple trees
    n = len(dna)
    for i in range(0, n+1):
        for k in range(i, n+1):
                                                        # TODO: for each beta tree
            if (initializeConditionsAreValid(tree, dna, i, k)):
                matrix[i][i+tree.L()][k][k+tree.R()] = tree.name



def initializeConditionsAreValid(tree, dna, i, k):
    # (1.1) Beta is mature                              # TODO
    # (1.2) i + |L(Beta)| <= k:
    if (i + tree.L() > k):   
        return False
    # (1.3) L(Beta) = a(i+1)...a(i+|L(Beta)|)
    # (1.4) R(Beta) = a(k+1)...a(k+|R(Beta)|)
    # combined 1.3 and 1.4 for treeB1 only:             # TODO: generalize, make part of tree classes
    if (dna[i:i+1] != dna[k:k+1]):
        return False
    return True



def construct(matrix, tree, dna):                       # TODO: multiple trees
    n = len(dna)
    for l in range(0, n+1):
        for i in range(l, -1, -1):
            for j in range(i, l+1):
                for k in range(l, j-1, -1):
                                                        # TODO: for each beta tree
                    if (constructConditionsAreValid(matrix, tree, dna, i, j, k, l)):
                        M[i][j][k][l] = tree.name



def constructConditionsAreValid(matrix, tree, dna, i, j, k, l):
    # (2.1) (See Urema et al. 1999)
    n = len(matrix[0][0][0])
    if (i+tree.LU < n and k+tree.RD < n and                           # !!! This part not specified in conditions... mistake in Uemura, or bug?
        matrix[i+tree.LU][j-tree.LD][k+tree.RD][l-tree.RU]) != "B1":
        return False
    #print "M(%d,%d,%d,%d) == B1" % (i+tree.LU, j-tree.LD, k+tree.RD, l-tree.RU)
    # (2.2) to (2.5):                                   # TODO: generalize, make part of tree classes
    if (dna[i:i+1] != dna[k:k+1]):
        return False
    return True



def recognize(matrix):
    n = len(matrix[0][0][0])
    for i in range(0, n):
        if (matrix[0][i][i][n-1] == "B1"):
            print "Accepted (M(%d,%d,%d,%d))" % (0,i,i,n)



# Main:

start = time.clock()
dna = sys.argv[1]
n = len(dna)
treeB1 = BetaTree("B1",1,0,1,0)  # direct repeat tree (LTR or TSD)
M = [[[[0 for x in xrange(n+1)] for x in xrange(n+1)] for x in xrange(n+1)] for x in xrange(n+1)]  # TODO: use NumPy
initialize(M, treeB1, dna)
#printMatrix(M); print; print
construct(M, treeB1, dna)
#printMatrix(M)
recognize(M)
end = time.clock()
print "dna length: ", n
print "elpsed time: ", end - start
raw_input("Press any key to exit: ")

# End