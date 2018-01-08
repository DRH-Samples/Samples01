#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from subprocess import call
import argparse

VERSION = '02'

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedInFile', type=str,
                        default='/data/glenn.hickey/genomes/alyrata/tracks/manual/manualChaux.bed',
                        help='BED input file [default: manualChaux.bed]')
    parser.add_argument('--tracksInFile', type=str,
                        default='/data/douglas.hoen/piatea/models/hmm/experiments/eta02/tracks.xml',
                        help='tracks file [default: eta02/tracks.xml')
    parser.add_argument('--iter', type=int, default=100,
                        help='# EM iterations [default: 100]')
    parser.add_argument('--nTE', type=int, default=5,
                        help='# TE states (min. 2) [default: 5]')
    parser.add_argument('--nOther', type=int, default=5,
                        help='# non-TE states (min. 2) [default: 5]')
    parser.add_argument('--pXS', type=float, default=0.1,
                        help='P(transition) X-self [default: 0.1]')
    parser.add_argument('--pTS', type=float, default=0.99,
                        help='P(transition) TE-self [default: 0.99]')
    parser.add_argument('--pTX', type=float, default=0.001,
                        help='P(transition) TE-X [default: 0.001]')
    parser.add_argument('--pOS', type=float, default=0.99,
                        help='P(transition) Other-self [default: 0.99]')
    parser.add_argument('--pOX', type=float, default=0.001,
                        help='P(transition) Other-X [default: 0.001]')
    args = parser.parse_args()
    if args.nTE < 2:
        exit('--nTE must be at least 2')
    if args.nOther < 2:
        exit('--nOther must be at least 2') 
    return args

def run(args):
    name = init(args)
    #tFile = sys.stdout
    tFile = open('init_transitions', 'w+')
    print_transitions(tFile, args)
    tFile.close()
    #eFile = sys.stdout
    eFile = open('init_emissions', 'w+')
    print_emissions(eFile, args)
    eFile.close()
    train_and_eval(args, name, tFile, eFile)

def init(args):
    name='eta%s-%s-%s' % (VERSION, args.nTE, args.nOther) 
    if os.path.exists(name):
        exit(name + ' exists')
    os.mkdir(name)
    os.chdir(name)
    return name

def print_transitions(tFile, args):
    nTE = args.nTE;         assert nTE >= 2
    nOther = args.nOther;   assert nOther >= 2 
    pXS = args.pXS                  # self
    print('%s\t%s\t%s' % ('X', 'X', pXS), file=tFile)
    print('', file=tFile)
    
    pXA = (1-pXS) / (nTE+nOther)    # to (TE or Other) ('A')
    print_submodel_transitions('TE',    args.nTE,    pXA, args.pTX, args.pTS, tFile)
    print_submodel_transitions('Other', args.nOther, pXA, args.pOX, args.pOS, tFile)

def print_submodel_transitions(label, n, pXA, pAX, pAS, tFile):
    nAA = n-1               # number of Ai-Aj transitions
    pAA = (1-pAX-pAS)/nAA   # to other A, per transition
    for i in range(n):
        print('%s\t%s%d\t%s' % ('X', label, i, pXA), file=tFile)                # X-i
        print('%s%d\t%s\t%s' % (label, i, 'X', pAX), file=tFile)                # Ai-X
        print('%s%d\t%s%d\t%s' % (label, i, label, i, pAS), file=tFile)         # Ai-Ai
        for j in range(i+1, n):
            print('%s%d\t%s%d\t%s' % (label, i, label, j, pAA), file=tFile)     # Ai-Aj
            print('%s%d\t%s%d\t%s' % (label, j, label, i, pAA), file=tFile)     # Aj-Ai
        print('', file=tFile)

def print_emissions(eFile, args):
    print('%s\t%s\t%d\t%d' % ('X', 'copy_cv', 1, 1), file=eFile)
    for i in range(args.nTE):
        print('%s%d\t%s\t%d\t%d' % ('TE', i, 'copy', 0, 0), file=eFile)

def train_and_eval(args, name, tFile, eFile):
    out = open('run.out', 'w+')
    call(['date'], stdout=out)
    trainCommand = ['teHmmTrain.py', args.tracksInFile, args.bedInFile,         #teHmmTrain.py tracks.xml $bedFile model.mod --numStates $numStates --iter $iter
                    'model.mod', '--iter', str(args.iter),
                    '--initTransProbs', tFile.name, '--initEmProbs', eFile.name]
    call(['echo', ' '.join(trainCommand)], stdout=out)
    call(trainCommand, stdout=out) 
    call(['date'], stdout=out)
    evalCommand = ['teHmmEval.py', args.tracksInFile, 'model.mod',
                   args.bedInFile, '--bed', 'eval.bed']
    call(['echo', ' '.join(evalCommand)], stdout=out)
    call(evalCommand, stdout=out) 
    call(['date'], stdout=out)
    call(['addBedColours.py', 'eval.bed', name+'.bed'], stdout=out)
    call(['addTrackHeader.py', name+'.bed', name, name, '--rgb'], stdout=out)
    call(['ln', '-s', '/data/douglas.hoen/piatea/models/hmm/experiments/'+name+'/'+name+'.bed', '/data/douglas.hoen/sharebrowser/temp/'])
    call(['echo', 'http://MUGS:18transposon@mustang.biol.mcgill.ca:8190/sharebrowser/temp/'+name+'.bed'], stdout=out)

if __name__ == '__main__':
    args = parse_cl()
    run(args)