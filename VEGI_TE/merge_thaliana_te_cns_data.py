#!/usr/bin/env python

import re
import copy
from drhBio import Locatable
from drhBio import SimpleLocation
from drhBio import calculate_overlaps
from pprint import pprint

class Cns:
    name = None
    location = None
    trna = None
    lyrata_loc = None
    te_overlaps = None
    
    def __init__(self, name, loc):
        self.name = name
        self.location = loc
    
    @staticmethod
    def header_str():
        return '\t'.join(['ID', 'Type', 'tRNA', 'Thaliana Location', 'Lyrata Location',
                         'Thaliana Length', 'Lyrata Length', 'Max TE Overlap (bp)',
                         'Name', 'Family', 'Superfamily', '# TEs', 'TEs'])
    
    def __str__(self):
        cns_type = self.name.split('_')[0]
        te_str = ' ; '.join([str(o.b) for o in self.te_overlaps])
        return '\t'.join([self.name, cns_type, str(self.trna), str(self.location), 
                         str(self.lyrata_loc), str(self.location.length()),
                         str(self.lyrata_loc.length()), str(self.te_overlaps[0].length),
                         self.te_overlaps[0].b.name, self.te_overlaps[0].b.family,
                         self.te_overlaps[0].b.superfamily, str(len(self.te_overlaps)), te_str])
            
    def set_trna(self, trna):
        if (self.trna != None):
            raise Exception('trna for {0} is already set!'.format(self.name))
        self.trna = trna
        
    def set_lyrata_loc(self, loc):
        if (self.lyrata_loc != None):
            raise Exception('lyrata_loc for {0} is already set!'.format(self.name))
        self.lyrata_loc = lyrata_loc
        
    def set_te_overlaps(self, te_overlaps):
        if (self.te_overlaps != None):
            raise Exception('te_overlaps for {0} is already set!'.format(self.name))
        self.te_overlaps = sorted(te_overlaps, reverse=True, key=lambda o: o.length)

class Trna:
    kind = None
    location = None
    score = None
    
    def __init__(self, loc, kind, score):
        self.location = loc
        self.kind = kind
        self.score = score
    
    def __str__(self):
        return '{0} {1}'.format(self.location, self.kind)
        
class TeFragment:
    name = None
    family = None
    superfamily = None
    location = None
    
    def __init__(self, name, family, superfamily):
        self.name = name
        self.family = family
        self.superfamily = superfamily
    
    def __str__(self):
        return '{0} {1} {2} {3}'.format(self.name, self.location, self.family, self.superfamily)
        
    def set_location(self, loc):
        if (self.location != None):
            raise Exception('location for {0} is already set!'.format(self.name))
        self.location = loc

def main():
    print 'Hello. Here we go...'
    cnses = read_thaliana_cns()
    add_trnas(cnses)
    add_lyrata_cns(cnses)
    [te_index, te_table] = index_tes()
    cns_with_te = overlap_tes(cnses, te_index, te_table)
    print_te_table(te_table)
    print Cns.header_str()
    for cns in cns_with_te:
        print str(cns)
    print 'Thank-you for playing. Goodbye.'
    
def read_thaliana_cns(fname='FINAL_CNS4_AT2.bed'):
    print 'Processing ' + fname
    cnses = []
    with open(fname) as infile:
        for line in infile:
            line = line.rstrip()
            seq, start, end, name, score, strand, a,b,c = line.split('\t')      # ignore abc
            seq = seq.replace('AT_scaffold', 'Chr')                             # replace, e.g., AT_scaffold1 with Chr1, for consistency with other files
            if not seq.startswith('Chr'):
                continue                                                        # ignore "AT_mitochondria"
            loc = SimpleLocation(seq, start, end, strand)
            cns = Cns(name, loc)
            cnses.append(cns)
    #cnses.sort()
    #pprint([str(item) for item in cnses])
    return cnses

def add_trnas(cnses, fname='tRNA.bed'):
    print 'Processing ' + fname
    trnas = []
    with open(fname) as infile:
        for line in infile:
            line = line.rstrip()
            fields = line.split('\t')
            if fields[0].startswith('track'): continue                          # skip header line
            seq, start, end, kind, score, strand, a,b,c,d,e,f = fields          # ignore abcdef
            if not seq.startswith('Chr'):
                continue                                                        # ignore "mitochondria" & "chloroplast"
            loc = SimpleLocation(seq, start, end, strand)
            trna = Trna(loc, kind, score)
            trnas.append(trna)
    for cns, trna_overlaps in calculate_overlaps(cnses, trnas).items():
        if len(trna_overlaps)>0:
            cns.set_trna(trna_overlaps[0].b)
        if len(trna_overlaps)>1:
            print 'Warning: {0} has more than one overlapping tRNA. Ignoring all but the first.'.format(str(cns))

def add_lyrata_cns(cnses, fname='FINAL_CNS4_AL2.bed'):
    print 'Processing ' + fname
    cns_index = {}
    for cns in iter(cnses):
        cns_index[cns.name] = cns
    with open(fname) as infile:
        for line in infile:
            line = line.rstrip()
            seq, start, end, name, score, strand, a,b,c = line.split('\t')      # ignore abc
            loc = SimpleLocation(seq, start, end, strand)
            if cns_index.has_key(name):                                         # 2361 lyrata CNS not found in thaliana
                cns_index[name].lyrata_loc = loc
                
    raise Exception("TODO: fix lyrata location: from e.g. AL_scaffold4:12217751-12217815(+) to scaffold_4:12217751-12217815; Also get rid of strand on thaliana locations")

# note: see also http://www.arabidopsis.org/servlets/processor?type=transposonfamily&update_action=browse
# From random checks, the sums that I get match the values listed there.
def index_tes(fname='TAIR9_Transposable_Elements.txt'):
    te_table = {    'LTR/Copia':    [0, {}],                                    # superfamily -> sort_order, {family -> [te_count, overlap_count]}
                    'LTR/Gypsy':    [1, {}],  
                    'LINE/L1':      [2, {}],  
                    'LINE?':        [3, {}],  
                    'SINE':         [4, {}],  
                    'DNA/Tc1':      [5, {}],  
                    'DNA/Mariner':  [6, {}],  
                    'DNA/Pogo':     [7, {}],
                    'DNA/HAT':      [8, {}],
                    'DNA/MuDR':     [9, {}],
                    'DNA/Harbinger':[10, {}],
                    'DNA/En-Spm':   [11, {}],
                    'RC/Helitron':  [12, {}],
                    'DNA':          [13, {}],
                    'RathE1_cons':  [14, {}],
                    'RathE2_cons':  [15, {}],
                    'RathE3_cons':  [16, {}],
                    'Unassigned':   [17, {}] }
    
    print 'Processing ' + fname
    te_index = {}
    with open(fname) as infile:
        for line in infile:
            line = line.rstrip()
            fields = line.split('\t')
            if fields[0].startswith('Transposon_Name'): continue                # skip header line
            name, a,b,c, family, superfamily = fields                           # ignore strand, start, end: useless b/c no sequence ID in file!
            if te_index.has_key(name):
                print 'WARNING: {0} found twice in {1}'.format(name, fname)
            else:
                te_index[name] = TeFragment(name, family, superfamily)
            if not te_table.has_key(superfamily):
                raise Exception('Unexpected superfamily ' + superfamily)
            if not te_table[superfamily][1].has_key(family):
                te_table[superfamily][1][family] = [0,0]
            te_table[superfamily][1][family][0] += 1
            None
    return [te_index, te_table]

def overlap_tes(cnses, te_index, te_table, fname='TAIR9_TE-intersect-CNS4_liftover.bed'):
    print 'Processing ' + fname
    found = {}  # to see how many TEs found w/ >1 fragment, just for interest
    tes = []
    with open(fname) as infile:
        fragment_re = re.compile('transposon_fragment_Parent=(\w+)')
        for line in infile:
            line = line.rstrip()
            fields = line.split('\t')
            if fields[0].startswith('track'): continue                          # skip header line
            seq, start, end, details, score, strand = fields
            # only use fragments b/c 'parent' can erroneously join fragments of different TEs
            match = fragment_re.match(details)
            if match:
                name = match.group(1)
                if not te_index.has_key(name):
                    raise Exception('{0} not found in TE index'.format(name))
                if found.has_key(name):
                    print '{0} found more than once'.format(name)
                found[name] = ''
                te = copy.copy(te_index[name])                                  # b/c some TEs (4 of them) have more than 1 fragment 
                te.set_location(SimpleLocation(seq, start, end, strand))
                tes.append(te)
            else:
                None
    cns_with_te = []
    for cns, te_overlaps in calculate_overlaps(cnses, tes).items():
        if len(te_overlaps)>0:
            cns.set_te_overlaps(te_overlaps)
            cns_with_te.append(cns)
        for o in te_overlaps:
            te = o.b
            te_table[te.superfamily][1][te.family][1] += 1
    cns_with_te.sort(key=lambda x: x.location)
    return cns_with_te
    
def print_te_table(te_table):
    print 'Thaliana TEs'
    print '\t'.join(['Superfamily', 'Family', 'Total', 'CNS Overlaps'])
    for superfamily, [sort_id, families] in sorted(te_table.items(), key=lambda t: t[1][0]):
        print superfamily
        for family, [total, overlaps] in sorted(families.items(), reverse=True, key=lambda t: t[1][1]):
            print '\t'.join(['', family, str(total), str(overlaps)])
    
    
### Now, do it ###
main()


