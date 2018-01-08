#!/usr/bin/env python

'''
Read weka predictions from a file
'''

import sys
import re
import math
import argparse
import datetime

class WekaPredictionsReader():
    
    def __init__(self, predFile, states):
        self.states = states.split(',')
        self.predictions = [[] for x in range( len(self.states) )]
        self._load(predFile)
    
    ''' Header lines:
        

    === Predictions on test data ===
    
     inst#     actual  predicted error distribution
    
    '''
    def _load(self, predFile):
        skipped = 0
        n = 0
        for line in predFile:
            if skipped < 5:                                                     # skip 5 header lines
                skipped += 1
                continue
            probs = self._parse(line.rstrip())
            for i in range(len(probs)):  
                self.predictions[i].append(probs[i])
            if n % 1000000 == 0:
                print str(n) + "\t" + str(datetime.datetime.now())
            n += 1
    
    # 700001     1:None     1:None       *0.986,0.002,0.002,0,0,0.008,0,0.002
    #    281     1:None   7:Simple   +   0.087,0,0,0,0,0.005,*0.716,0.192
    def _parse(self, line):
        m = re.search(r'[,\*.\d]*$', line, re.M)
        if not m:
            raise Exception('Unable to parse: ' + line)
        csv = re.sub(r'[*]', '', m.group())
        fields = csv.split(',')
        return map(float, fields)
    
    def log_table(self):
        d = {}
        i = 0
        for s in self.states:
            d[s] = map( math.log, map( lambda x: max(sys.float_info.min, x), self.predictions[i] ))
            i += 1
        return d
    
    def print_table(self, table):
        keys = table.keys()
        print(keys)
        for i in range(len(table[keys[0]])):
            row = []
            for state in keys:
                row.append(table[state][i])
            print ",".join(map(str, row))
            None


# From weka header:
# @attribute repeat_modeler {None, DNA, LINE, LTR, SINE, Unknown, Simple, Low}

def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('predictionsFile', type=argparse.FileType('r'),
        help='BED file from which to calculate transition probabilities')
    parser.add_argument('--states', type=str, default='None,DNA,LINE,LTR,SINE,Unknown,Simple,Low',
        help='Names of states'
            '[default: \'None,DNA,LINE,LTR,SINE,Unknown,Simple,Low\']')
    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_cl()
    reader = WekaPredictionsReader(args.predictionsFile, args.states)
    predictionsTable = reader.log_table()
    reader.print_table(predictionsTable)
    None
    