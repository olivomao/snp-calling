from itertools import izip
import time
import re

# precompile regular expression patterns
num_pattern = re.compile('[\d]+')
chr_pattern = re.compile('[^\d]+')
dna_pattern = re.compile(r'([acgtnACGTN]+)')
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
NUM_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

class Read:
    """Representation of an alignment line from a .SAM file.
    Spec: http://samtools.sourceforge.net/SAMv1.pdf
    """

    def __init__(self, line):
        # read out fields
        fields = line.split()
        self.line = fields
        self.id = fields[0]
        self.flag_num = int(fields[1])
        self.flag = '{0:08b}'.format(self.flag_num)
        self.pos = int(fields[3]) -1
        self.read = fields[9]

    def is_reversed(self):
        return (self.flag_num >> 4) & 1

    def is_first_segment(self):
        return (self.flag_num >> 6) & 1

class SNP:
    def __init__(self, pos):
        self.pos = pos # position on the target
        self.reads = []
        self.read_bases = []
        self.read_weights = []
        self.most_likely = ""

def get_sequence(address):
    """Loads the DNA sequence string at the address"""
    ref = ''
    with open(address,'rU') as ref_file:
        next(ref_file) # skip header lines
        for line in ref_file:
            ref = ref + line.strip()
    return ref

def write_sequence(address, sequence):
    with open(address, 'w+') as sequence_file:
        sequence_file.write(">" + address + "\n")
        sequence_file.write(sequence)

MAX_NOISE = 0.20 # threshold for random mutation of reads
def is_SNP(count):
    """Returns whether the count of reads at this locus looks like a SNP
    (41, 0, 0, 1) -> False
    (0, 80, 1, 23) -> True
    """
    counts = sum(count)
    return counts and float(counts - max(count)) / counts > MAX_NOISE

def distance(str1, str2):
    """Currently implements Hamming distance but anything can be used"""
    return sum(c1 != c2 for c1, c2 in izip(str1, str2))

def probability(distances):
    """Returns a list of the same length as distances where each element is
    probability that the read came from locus corresponding to that distance
    such that sum(list) = 1

    Currently adds 1 to each distance, takes the inverse then normalizes vector
    """
    v = [1.0/(d + 1) for d in distances]
    s = sum(v)
    return [i/s for i in v]

def distribution(bases, weights):
    """Return [P(A), P(C), P(G), P(T)] from the weighted reads"""
    dist = [0, 0, 0, 0]
    for base, weight in zip(bases, weights):
        dist[BASES[base]] += weight
    return dist