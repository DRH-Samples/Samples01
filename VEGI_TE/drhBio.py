from pprint import pprint


#Input: two lists of Locatable objects
#Output: {a, [Overlap]}
def calculate_overlaps(locatablesA, locatablesB):
    omap = {}
    alist = sorted(locatablesA, key=lambda x: x.location)   # lambda creates an anonymous function
    blist = sorted(locatablesB, key=lambda x: x.location)
    i = 0
    for a in alist:
        aloc = a.location
        while (i>0 and blist[i-1].location.seq == aloc.seq and blist[i].location.end >= aloc.start):
            i=i-1 # back up in case a's are overlapping (allow one b to be assigned to >1 overlapping a, i.e. N-to-N)
        overs = []
        while i<len(blist):
            b = blist[i]
            bloc = b.location
            olength = aloc.overlap(bloc)
            if olength>0:
                overs.append(Overlap(a, b, olength))
                #print '{0} : {1}'.format(str(a), str(b))
            if bloc.start > aloc.end or bloc.seq != aloc.seq:
                break
            i=i+1
        if i==len(blist):
            i=i-1
        omap[a] = overs;
    return omap

#Abstract class: must set location (at least)
class Locatable:
    location = None
    def __cmp__(self, other):
        return cmp(self.location, other.location)
    #def __hash__(self):
    #    return id(self)
    
class Overlap:
    def __init__(self, a, b, length):
        self.a = a
        self.b = b
        self.length = length

class SimpleLocation:
    def __init__(self, seq, start, end, strand='.'):
        self.seq = seq
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        
    def __cmp__(self, other):
        return (cmp(self.seq, other.seq) or
                cmp(self.start, other.start) or
                cmp(self.end, other.end))
    
    def __str__(self):
        return '{0}:{1}-{2}({3})'.format(self.seq, self.start, self.end, self.strand)
        
    def has_overlap(self, loc):
        return self.overlap(self, loc) > 0
        #if (self.seq == loc.seq and
        #    self.start <= loc.end and
        #    self.end >= loc.start ) :
        #    return True
        #else:
        #    return False
        
    def overlap(self, loc):
        return min(self.end, loc.end) - max(self.start, loc.start) + 1
        
    def length(self):
        return self.end - self.start + 1