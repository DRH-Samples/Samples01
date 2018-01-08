#!/usr/bin/env python

'''
Adapted from http://en.wikipedia.org/wiki/Viterbi_algorithm
'''

import math
import sys
import argparse
from bedstats import BedStats
from weka_predictions_reader import WekaPredictionsReader

def viterbi(pred, states, startP, stateP, transP, adjustment):
    V = [{}]
    path = {}
 
    # Initialize base cases (t == 0)
    for s in states:
        V[0][s] =  startP[s] + pred[s][0]
        path[s] = [s]
 
    # Run Viterbi for t > 0
    n = len(pred[ pred.keys()[0] ])
    for t in range(1, n):
        V.append({})
        newpath = {}
 
        for s in states:
            (prob, state) = max( (V[t-1][s0] + transP[s0][s] + pred[s][t] - stateP[s] + adjustment , s0) for s0 in states )    # P(s|o)/P(s)
            #(prob, state) = max( (V[t-1][s0] + transP[s0][s] + pred[s0][t] - stateP[s0] + adjustment , s0) for s0 in states )    # P(s|o)/P(s)
            #(prob, state) = max( (V[t-1][s0] + transP[s0][s] + pred[s0][t] , s0) for s0 in states )    # P(s|o)
            V[t][s] = prob
            newpath[s] = path[state] + [s]
 
        # Don't need to remember the old paths
        path = newpath
        #if not t % 94:
        #    None
        #if not t % 10000:
        #    print(t)
    n = 0           # if only one element is observed max is sought in the initialization values
    if n != 1:
        n = t
    #print_dptable(V)
    (prob, state) = max((V[n][y], y) for y in states)
    return (prob, path[state])
 
# Print a table of the steps.
def print_dptable(V):
    s = "         " + " ".join(("%8d" % i) for i in range(len(V))) + "\n"
    for y in V[0]:
        s += "%9s: " % y
        s += " ".join("%.9s" % ("%.2e" % v[y]) for v in V)
        s += "\n"
    print(s)



def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('bedfile', type=argparse.FileType('r'),
        help='BED file from which to calculate transition probabilities')
    parser.add_argument('chromsizes', type=argparse.FileType('r'),
        help='chromosome sizes in a tab-separated file such as that generated by fetchChromSizes '
             '(see https://genome.ucsc.edu/goldenPath/help/bigBed.html)' )
    parser.add_argument('--background', type=str, default='None',
        help='Name of background state (i.e. what to call gaps) [default: \'None\']')
    parser.add_argument('predictionsFile', type=argparse.FileType('r'),
        help='BED file from which to calculate transition probabilities')
    parser.add_argument('--states', type=str, default='None,DNA,LINE,LTR,SINE,Unknown,Simple,Low',
        help='Name of background state (i.e. what to call gaps) '
            '[default: \'None,DNA,LINE,LTR,SINE,Unknown,Simple,Low\']')
    parser.add_argument('--fragX', type=float, default=1.0,
        help='Factor to adjust for fragmentation in calculating transition probabilities.'
            'The number of elements of each type will be divided by this number [default: 1.0].')
    parser.add_argument('--stateX', type=float, default=1.0,
        help='Factor to reduce influence of state probabilities. Range [0,1].'
            'State probabilities remain unadjusted at stateX=0 and approach 1 at stateX=1. [default: 0.0].')
    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_cl()
    bedStats = BedStats(args.bedfile, args.chromsizes, args.background)
    states = bedStats.states()
    
    totalStateP = 0.0
    for s in states:
        totalStateP += bedStats.statep(s, args.stateX)
    
    startPTable = {}
    statePTable = {}
    for s in states:
        startPTable[s] = math.log( 1.0/len(states) )
        statePTable[s] = math.log( bedStats.statep(s, args.stateX) / totalStateP )  # Normalized to sum = 1.
    
    stateP = []
    STATES = ("None", "DNA", "LINE", "LTR", "SINE", "Unknown", "Simple", "Low")
    for s in STATES:
        stateP.append(statePTable[s])
    print stateP
    
    transPTable = {}
    for i in states:
        transPTable[i] = {}
        for j in states:
            transPTable[i][j] = -1*sys.float_info.max
    
    for i in states:
        stateStats = bedStats.stateStats[i]
        transPTable[i][i] = math.log( stateStats.transp_self(args.fragX) )
        if i != bedStats.bgState:
            transPTable[i][bedStats.bgState] = math.log( stateStats.transp_to_bg(args.fragX) )
            transPTable[bedStats.bgState][i] = math.log( bedStats.transp_from_bg(i, args.fragX) )
            None
    #print(sum(map(math.exp, transPTable[bedStats.bgState].values())))
    
    print
    for i in STATES:
        transPi = []
        for j in STATES:
            transPi.append(transPTable[i][j])
        print transPi    
    None
    
    
    def printTable(table, name):
        print("\n" + name)
        for k in table:
            print(k + "\t" + str(math.exp(table[k])))
    
    def printTable2(table, name):
        print("\n" + name)
        for k in table:
            print
            for v in table[k]:
                print(k + "\t" + v + "\t" + str(math.exp(table[k][v])))
    
    #printTable(startPTable, 'Start')
    #printTable(statePTable, 'State')
    #printTable2(transPTable, 'Transition')
    #print
    
    adjustment = math.log( 1.0 / bedStats.scafSizes.genomeSize )
    print
    print adjustment
    None
    
    predReader = WekaPredictionsReader(args.predictionsFile, args.states)
    predictions = predReader.log_table()
    

    
    def example():
        return viterbi(predictions,
                       states,
                       startPTable,
                       statePTable,
                       transPTable,
                       adjustment)
    
    result = example()
    
    #print result[0]
    #i = 700001
    #for s in result[1]:
    #    print("{}: {}".format(i, s))
    #    i += 1
    
    def to_bed(result, scaf, offset, background):
        prevState = result[0]
        left = 0
        for i in range(0, len(result)-1):
            i += 1
            if result[i] != prevState:
                if prevState != background:
                    print( "\t".join((scaf, str(left+offset), str(i+offset), prevState)))
                prevState = result[i]
                left = i
        if prevState != background:
            print( "\t".join(scaf, left+offset, i+offset, prevState))
    to_bed(result[1], 'scaffold_1', 10700002, 'None')
    