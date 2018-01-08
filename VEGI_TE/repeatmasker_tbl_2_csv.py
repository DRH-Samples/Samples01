#!/usr/bin/env python

# Output headers:
# File Name, Library, Number of Elements(N) / Length Occupied (L) (bp) / Percentage of Sequence (P) (%),
# SINEs, LINEs, LTR elements, DNA elements, Unclassified, Small RNA, Satellites, Simple Repeats, Low Complexity

import argparse
import re

argParser = argparse.ArgumentParser()
argParser.add_argument('--input', type = argparse.FileType('r')) # 
args = argParser.parse_args()

stage = 0
parameters = []
counts = []
lengths = []
percents = []

reFile = re.compile('file name: (\w+)')                                                 # file name: AA  
reData = re.compile('[\w ]+:\s+(\d+)\s+(\d+) bp\s+([\d\.]+) %')                         # E.g.: SINEs:              100         9399 bp    0.01 %
reLib  = re.compile('The query was compared to classified sequences in "[\./]*(\S+)"') # The query was compared to classified sequences in ".../AL-RepClass.fa"

for line in args.input:               
        mFile = reFile.match(line)
        if mFile:
            parameters.append(mFile.group(1))
        mData = reData.match(line)
        if mData:
            counts.append(  mData.group(1))
            lengths.append( mData.group(2))
            percents.append(mData.group(3))
        mLib = reLib.match(line)
        if mLib:
            parameters.append(mLib.group(1))
        
print ', '.join(parameters + ['N'] + counts)
print ', '.join(parameters + ['L'] + lengths) 
print ', '.join(parameters + ['P'] + percents) 


# Example Input:
#
#==================================================
#file name: AA                       
#sequences:         59101
#total length:  199432295 bp  (165780879 bp excl N/X-runs)
#GC level:         33.97 %
#bases masked:   22999898 bp ( 13.87 %)
#==================================================
#               number of      length   percentage
#               elements*    occupied  of sequence
#--------------------------------------------------
#SINEs:              100         9399 bp    0.01 %
#      ALUs            0            0 bp    0.00 %
#      MIRs            0            0 bp    0.00 %
#
#LINEs:             2577       947780 bp    0.57 %
#      LINE1           0            0 bp    0.00 %
#      LINE2           0            0 bp    0.00 %
#      L3/CR1          0            0 bp    0.00 %
#
#LTR elements:     18334      6627154 bp    4.00 %
#      ERVL            0            0 bp    0.00 %
#      ERVL-MaLRs      0            0 bp    0.00 %
#      ERV_classI      0            0 bp    0.00 %
#      ERV_classII     0            0 bp    0.00 %
#
#DNA elements:      9283      2470810 bp    1.49 %
#     hAT-Charlie      0            0 bp    0.00 %
#     TcMar-Tigger     0            0 bp    0.00 %
#
#Unclassified:     11592      2026568 bp    1.22 %
#
#Total interspersed repeats: 12081711 bp    7.29 %
#
#
#Small RNA:            0            0 bp    0.00 %
#
#Satellites:           6          691 bp    0.00 %
#Simple repeats:   24858      1223543 bp    0.74 %
#Low complexity:  174707      9814803 bp    5.92 %
#==================================================
#
#* most repeats fragmented by insertions or deletions
#  have been counted as one element
#  Runs of >=20 X/Ns in query were excluded in % calcs
#
#
#The query species was assumed to be homo          
#RepeatMasker version open-3.3.0 , sensitive mode
#                                 
#run with rmblastn version : 2.2.23+
#The query was compared to classified sequences in ".../AL-RepClass.fa"
#RepBase Update 20120418, RM database version 20120418
