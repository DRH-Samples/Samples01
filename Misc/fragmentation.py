#!/usr/bin/env python

'''
Calculate and display various measures of the fragmentation of an annotation,
both inherent and relative to a reference annotation.

Intended to be used either as a stand-alone script, or as part of a larger
benchmarking system.
'''

from __future__ import print_function
from __future__ import division         # float division without casting
import sys
import os
import re
import argparse
import itertools
import pybedtools



class FragStats:

# TODO    
#- prec, rec, f1 at no merge and at 100 bp
#- plots
    
    def __init__(self, query, reference=None, referenceCoverage=0):
        self.merges = []
        self.md50 = -999
        self.precision0 = -999.0
        self.recall0 = -999.0
        self.f0 = -999.0
        self.precision100 = -999.0
        self.recall100 = -999.0
        self.f100 = -999.0
        self.deltaF = -999.0
        #for dist in itertools.chain( (-1,), range(10), range(10, 100, 10), range(100, 1001, 100) ):
        for dist in itertools.chain( (-1,), range(0, 100, 10) ):
            each = FragStatsElement(query.merge(d=dist), dist, reference, referenceCoverage)
            self.merges.append(each)
            if len(self.merges) > 1 and self.md50 == -999:
                first = self.merges[0]
                if (each.count / first.count <= 0.5):
                    self.md50 = each.mergeDist
            if reference is not None and dist == 0:
                self.precision0 = each.precision
                self.recall0 = each.recall
                self.f0 = each.f
            if reference is not None and dist == 100:
                self.precision100 = each.precision
                self.recall100 = each.recall
                self.f100 = each.f
                self.deltaF = self.f100-self.f0
    
    def __repr__(self):
        buf = ''
        buf = buf + ('{:>15s}'*8+'{}').format('MD50', 'Precision 0', 'Recall 0', 'F1 0',
                                              'Precision 100', 'Recall 100', 'F1 100',
                                              'Delta F1', '\n')
        buf = buf + ('{:>15d}'*1+'{:15.4f}'*7+'{}').format(
            self.md50, self.precision0, self.recall0, self.f0,
            self.precision100, self.recall100, self.f100, self.deltaF, '\n\n')
        buf = buf + FragStatsElement.header_str()
        prev = None
        for curr in self.merges:
            buf = buf + curr.to_str(prev) + '\n'
            prev = curr
        return buf



class FragStatsElement:
    def __init__(self, query, mergeDist, reference, referenceCoverage):
        self.query = query.sort()
        self.mergeDist = mergeDist
        self.count = query.count()
        self.coverage = query.total_coverage()
        self.meanLength = self.coverage / self.count
        self.precision = -999.0
        self.recall = -999.0
        self.f = -999.0
        if reference is not None:
            truePositives = query.intersect(reference).sort().total_coverage()
            self.precision = truePositives / self.coverage
            self.recall = truePositives / referenceCoverage
            self.f = 2 * self.precision * self.recall / (self.precision + self.recall)
    
    def delta(self, prev, attr):
        if prev is None:
            return 0
        else:
            return self.getattr(attr) - prev.getattr(attr)
    
    # Overload getattr()
    def getattr(self, attr):
        if hasattr(self, attr):
            return getattr(self, attr)
        else:
            exit('TO DO')
            #return getattr(self.bedtool, attr)()

    
    cols = ('coverage', 'count', 'meanLength', 'precision', 'recall', 'f')
    nCols = 2*len(cols)+1
    fmt = '{:15d}'*5+'{:15.0f}'*2+'{:15.4f}'*6
    
    @staticmethod
    def header_str():
        return ('{:>15s}'*FragStatsElement.nCols+'{}').format(
                'Merge Dist.', 'Coverage', '(Delta)', 'Fragments', '(Delta)',
                'Mean Length', '(Delta)', 'Precision', '(Delta)',
                'Recall', '(Delta)', 'F1', '(Delta)', '\n')
    
    def to_str(self, prev):
        vals = [self.mergeDist]
        #for attr in ('total_coverage', 'count', 'mean_length'):
        for attr in (self.cols):
            vals.append(self.getattr(attr))
            vals.append(self.delta(prev, attr))
        return self.fmt.format(*vals)
    
    
    
def parse_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=str, 
                            help='reference BED file')
    parser.add_argument('query', type=str, 
                            help='query BED file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_cl()
    
    tmpdir = './tmp'
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    pybedtools.set_tempdir(tmpdir)
    pybedtools.cleanup(remove_all=True)     # delete any files not deleted during previous runs
    
    query = pybedtools.BedTool(args.query).sort()
    reference = pybedtools.BedTool(args.reference).sort()
    
    print('\nReference: {}'.format(args.reference))
    print('Query: {}\n'.format(args.query))
    print('\nReference\n---------\n')
    print(FragStats(reference))
    print('\nQuery\n-----\n')
    print(FragStats(query, reference, reference.total_coverage()))
    

