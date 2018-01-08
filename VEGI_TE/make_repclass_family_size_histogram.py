#!/usr/bin/env python

import argparse
import re

argParser = argparse.ArgumentParser()
argParser.add_argument('--input', type = argparse.FileType('r')) # RepClass 'consensi.fa.classified' output file
args = argParser.parse_args()

rsRe = re.compile('^>rnd-1_family-\d+#(\w+)(/\S*)? \( RepeatScout Family Size = (\d+)')

# class -> [size] , e.g. 'LTR' => 130, 37, 204, ...
rsCounts = {}       # RepeatScout
reconCounts = {}    # RECON

for line in args.input:
    if line.startswith('>rnd-1'):   # round 1 is always RepeatScout, rounds 2+ are always RECON
        m = rsRe.match(line)
        if m:
            te = m.group(1)
            if not te in rsCounts:
                rsCounts[te] = []
            rsCounts[te].append(int(m.group(3)))
        else:
            raise Exception("Couldn't parse: " + line)
        
